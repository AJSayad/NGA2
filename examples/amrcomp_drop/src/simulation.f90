!> AMR compressible drop test case
module simulation
   use precision,           only: WP
   use string,              only: str_medium
   use amrgrid_class,       only: amrgrid
   use amrmpcomp_class,     only: amrmpcomp
   use amrviz_class,        only: amrviz
   use amrdata_class,       only: amrdata
   use timetracker_class,   only: timetracker
   use event_class,         only: event
   use monitor_class,       only: monitor
   use amrio_class,         only: amrio
   use nasg_class,          only: nasg
   use ideal_gas_class,     only: ideal_gas
   use relax_ig_nasg_class, only: relax_ig_nasg, Prelax
   implicit none
   private

   public :: simulation_init,simulation_run,simulation_final

   !> AMR grid
   type(amrgrid), target :: amr

   !> Timetracker and compressible multiphase solver
   type(timetracker) :: time
   type(amrmpcomp), target :: fs
   type(amrdata) :: dQdt,drhoUdt,drhoVdt,drhoWdt,Umag,Mach

   !> Visualization
   type(event) :: viz_evt,plicviz_evt,tau_viz_evt
   type(amrviz) :: viz,plicviz,tau_viz

   ! Regrid parameters
   type(event) :: regrid_evt

   ! Restart parameters
   type(amrio) :: io
   type(event) :: save_evt
   character(len=str_medium) :: restart_dir
   logical :: restarted
   real(WP) :: restart_time

   !> Simulation monitoring
   type(monitor) :: mfile,consfile,cflfile,gridfile,tfile,rescfile

   !> Mixture total energy, the sum of the phasic total energies
   real(WP) :: Etot=0.0_WP

   !> Droplet QOI monitoring
   type(monitor) :: drop_QOIs
   real(WP) :: MLE_VF001,MLE_VF01,MLE_VF05                            !< Mist leading edge at VF>=0.01 and VF>=0.1 thresholds
   real(WP) :: drop_massVF001,drop_massVF01                           !< Droplet mass at VF>=0.01 and VF>=0.1 thresholds (smoothed)
   real(WP) :: MOI_xx,MOI_yy,MOI_zz,MOI_xy,MOI_xz,MOI_yz              !< Moment of inertia tensor about VF-weighted COM (VF>=0.1)
   real(WP) :: MOI_xx_001,MOI_yy_001,MOI_zz_001,MOI_xy_001,MOI_xz_001,MOI_yz_001 !< Moment of inertia tensor about VF-weighted COM (VF>=0.01)

   !> Materials
   type(nasg),      target :: water
   type(ideal_gas), target :: gas

   !> Relaxation model
   type(relax_ig_nasg), target :: relax_model

   !> Flow parameters
   real(WP) :: rhoG1,pG1,u1           !< Pre-shock gas state
   real(WP) :: rhoG2,pG2,u2           !< Post-shock gas state
   real(WP) :: rhoL1,pL1              !< Initial liquid state
   real(WP) :: M2,Xs                  !< Post-shock Mach and shock location
   real(WP) :: Ms                     !< Shock Mach number
   real(WP) :: density_ratio          !< rhoL1/rhoG1
   real(WP) :: ML                     !< Liquid Mach number
   real(WP) :: Reynolds,visc_ratio    !< Viscosity
   real(WP) :: Prandtl ,diff_ratio    !< Heat diffusivity
   real(WP) :: Weber                  !< Weber number

   !> Sutherland viscosity parameters: mu_g = (1+Suth_T)*T^Suth_n / (Re*(T+Suth_T))
   real(WP) :: Suth_n=1.5_WP          !< Sutherland exponent (1.0 for constant)
   real(WP) :: Suth_T=0.4042_WP       !< Sutherland temperature (0.0 for constant)

   !> Sponge parameters
   real(WP) :: R_spg=3.0_WP           !< Radius at which the sponge starts
   real(WP) :: L_spg=1.0_WP           !< Length of the sponge ramp
   real(WP) :: sig_spg=1.0_WP         !< Sponge relaxation rate [1/t] - set to 0 to disable
   real(WP) :: Ws                     !< Lab-frame shock speed - locates the sponge target state

   !> Spherical harmonics perturbation parameters
   integer :: nsh_modes=0
   integer , dimension(:), allocatable :: l_modes,m_modes
   real(WP), dimension(:), allocatable :: amp_modes,phase_modes

   !> Tagging parameter
   real(WP) :: Re_tag=huge(1.0_WP)
   real(WP) :: Rho_tag=huge(1.0_WP)
   real(WP) :: P_tag=huge(1.0_WP)
   real(WP) :: Ducros_tag=huge(1.0_WP)

contains

   !> Smooth Heaviside function
   real(WP) function Hshock(x,delta)
      real(WP), intent(in) :: x,delta
      Hshock=1.0_WP/(1.0_WP+exp(-x/delta))
   end function Hshock

   !> Levelset function for drop (centered at x=x_drop)
   function sphere_levelset(xyz,t) result(G)
      use mathtools, only: spherical_harmonic
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G,r,theta,phi,perturb
      integer :: i
      !G=0.5_WP-sqrt((xyz(1)-x_drop)**2+xyz(2)**2+xyz(3)**2)
      r=sqrt((xyz(1))**2+xyz(2)**2+xyz(3)**2)
      if (r.gt.1.0e-12_WP) then
         theta= acos(xyz(3)/r)
         phi  =atan2(xyz(2),xyz(1))
      else
         theta=0.0_WP
         phi  =0.0_WP
      end if
      ! Compute perturbation
      perturb=0.0_WP
      do i=1,nsh_modes
         perturb=perturb+amp_modes(i)*spherical_harmonic(l_modes(i),m_modes(i),theta,phi+phase_modes(i))
      end do
      G=0.5_WP+perturb-r
      if (amr%nz.eq.1) G=0.5_WP-sqrt((xyz(1))**2+xyz(2)**2) ! Enable quasi-2D runs
   end function sphere_levelset

   !> Compute viscosity: Sutherland for gas, VF-weighted blend with liquid
   subroutine get_viscosities()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pTG,pVF,pVisc,pBeta,pDiffL,pDiffG
      real(WP) :: mu_g,mu_l
      real(WP), parameter :: Tmax_visc=10.0_WP
      real(WP), parameter :: myeps=1.0e-15_WP
      ! Loop over levels
      do lvl=0,amr%clvl()
         ! Loop over domain
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pTG=>fs%TG%mf(lvl)%dataptr(mfi)
            pVF=>fs%VF%mf(lvl)%dataptr(mfi)
            pVisc=>fs%visc%mf(lvl)%dataptr(mfi)
            pBeta=>fs%beta%mf(lvl)%dataptr(mfi)
            pDiffL=>fs%diffL%mf(lvl)%dataptr(mfi)
            pDiffG=>fs%diffG%mf(lvl)%dataptr(mfi)
            ! Get tilebox grown by the viscosity ghost cells
            bx=mfi%growntilebox(min(fs%visc%ng,fs%beta%ng))
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                     ! Gas viscosity from Sutherland
                     mu_g=(1.0_WP+Suth_T)*min(pTG(i,j,k,1),Tmax_visc)**Suth_n/(Reynolds*(min(pTG(i,j,k,1),Tmax_visc)+Suth_T))
                     ! Liquid viscosity from ratio
                     mu_l=visc_ratio*Reynolds**(-1.0_WP)
                     ! Mixture viscosity
                     !pVisc(i,j,k,1)=pVF(i,j,k,1)*mu_l+(1.0_WP-pVF(i,j,k,1))*mu_g ! Arithmetic averaging
                     pVisc(i,j,k,1)=1.0_WP/(pVF(i,j,k,1)/max(mu_l,myeps)+(1.0_WP-pVF(i,j,k,1))/max(mu_g,myeps)) ! Harmonic averaging
                     ! Zero bulk viscosity
                     pBeta(i,j,k,1)=0.0_WP
                     ! Phasic heat diffusivities: gas k=cp*mu/Pr, liquid from ratio (no blending - solver uses phasic fields)
                     ! pDiffG(i,j,k,1)=gas%gamma*gas%cv*mu_g/Prandtl
                     ! pDiffL(i,j,k,1)=diff_ratio*gas%gamma*air%cv/(Reynolds*Prandtl)
                     ! Run with no heat transfer
                     pDiffG(i,j,k,1)=0.0_WP
                     pDiffL(i,j,k,1)=0.0_WP
                     ! No sponge viscosity here - the sponge is a relaxation term, see apply_sponge
                  end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine get_viscosities

   !> Undisturbed (drop-free) planar shock state at a given x and time
   subroutine freestream_state(x_cc,t,dx,rho,p,u)
      real(WP), intent(in) :: x_cc,t,dx
      real(WP), intent(out) :: rho,p,u
      real(WP) :: H
      ! Shock travels at Ws in the lab frame - once it exits, this is uniformly post-shock
      H=Hshock(x=Xs+Ws*t-x_cc,delta=0.5_WP*dx)
      rho=rhoG1+(rhoG2-rhoG1)*H
      p  =pG1  +(pG2  -pG1  )*H
      u  =u1   +(u2   -u1   )*H
   end subroutine freestream_state

   !> Sponge zone: relax the gas state toward the undisturbed solution outside r=R_spg.
   subroutine apply_sponge(dt,t)
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      integer :: lvl,i,j,k
      real(WP), intent(in) :: dt,t
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pRHOU,pRHOV,pRHOW
      real(WP) :: dx,dy,dz,fac,rho,p,u
      real(WP), parameter :: VF_spg=1.0e-6_WP ! Below this VF a cell is treated as pure gas
      ! Nothing to do if the sponge is off
      if (sig_spg.le.0.0_WP) return
      ! Loop over levels
      do lvl=0,amr%clvl()
         dx=amr%dx(lvl); dy=amr%dy(lvl); dz=amr%dz(lvl)
         ! Loop over domain
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pQ=>fs%Q%mf(lvl)%dataptr(mfi)
            pRHOU=>fs%rhoU%mf(lvl)%dataptr(mfi)
            pRHOV=>fs%rhoV%mf(lvl)%dataptr(mfi)
            pRHOW=>fs%rhoW%mf(lvl)%dataptr(mfi)
            ! Cells: relax the gas conserved variables toward the undisturbed state
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                     fac=sponge_factor(y=amr%ylo+(real(j,WP)+0.5_WP)*dy,z=amr%zlo+(real(k,WP)+0.5_WP)*dz)
                     if (fac.le.0.0_WP) cycle
                     ! Leave any cell carrying liquid alone
                     if (pQ(i,j,k,1).gt.VF_spg*rhoL1) cycle
                     ! Undisturbed state at the cell center
                     call freestream_state(x_cc=amr%xlo+(real(i,WP)+0.5_WP)*dx,t=t,dx=dx,rho=rho,p=p,u=u)
                     pQ(i,j,k,2)=pQ(i,j,k,2)+fac*(rho-pQ(i,j,k,2))
                     ! Q(4) is the gas total energy: relax it toward the undisturbed internal plus kinetic energy
                     pQ(i,j,k,4)=pQ(i,j,k,4)+fac*(rho*gas%get_e_from_p_rho(p=p,rho=rho,y=[1.0_WP])+0.5_WP*rho*u**2-pQ(i,j,k,4))
                  end do; end do; end do
            ! Faces: relax the mixture momentum toward the undisturbed state
            bx=mfi%nodaltilebox(1)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                     fac=sponge_factor(y=amr%ylo+(real(j,WP)+0.5_WP)*dy,z=amr%zlo+(real(k,WP)+0.5_WP)*dz)
                     if (fac.le.0.0_WP) cycle
                     if (max(pQ(i-1,j,k,1),pQ(i,j,k,1)).gt.VF_spg*rhoL1) cycle
                     ! Undisturbed state at the x-face
                     call freestream_state(x_cc=amr%xlo+real(i,WP)*dx,t=t,dx=dx,rho=rho,p=p,u=u)
                     pRHOU(i,j,k,1)=pRHOU(i,j,k,1)+fac*(rho*u-pRHOU(i,j,k,1))
                  end do; end do; end do
            bx=mfi%nodaltilebox(2)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                     fac=sponge_factor(y=amr%ylo+real(j,WP)*dy,z=amr%zlo+(real(k,WP)+0.5_WP)*dz)
                     if (fac.le.0.0_WP) cycle
                     if (max(pQ(i,j-1,k,1),pQ(i,j,k,1)).gt.VF_spg*rhoL1) cycle
                     pRHOV(i,j,k,1)=pRHOV(i,j,k,1)-fac*pRHOV(i,j,k,1)
                  end do; end do; end do
            bx=mfi%nodaltilebox(3)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                     fac=sponge_factor(y=amr%ylo+(real(j,WP)+0.5_WP)*dy,z=amr%zlo+real(k,WP)*dz)
                     if (fac.le.0.0_WP) cycle
                     if (max(pQ(i,j,k-1,1),pQ(i,j,k,1)).gt.VF_spg*rhoL1) cycle
                     pRHOW(i,j,k,1)=pRHOW(i,j,k,1)-fac*pRHOW(i,j,k,1)
                  end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   contains
      !> Relaxation factor over dt at (y,z): zero inside r=R_spg
      real(WP) function sponge_factor(y,z)
         real(WP), intent(in) :: y,z
         real(WP) :: r_cyl
         ! Distance from the drop axis
         r_cyl=sqrt(y**2+z**2)
         if (amr%nz.eq.1) r_cyl=abs(y) ! Enable quasi-2D runs
         ! Quadratic ramp - zero value and zero slope at the sponge entrance
         sponge_factor=1.0_WP-exp(-sig_spg*min(max(r_cyl-R_spg,0.0_WP)/L_spg,1.0_WP)**2*dt)
      end function sponge_factor
   end subroutine apply_sponge

   !> User init callback - set Q and VF/barycenters for a drop at rest with a shock
   subroutine shockdrop_init(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box
      use amrex_amr_module, only: amrex_mfiter_build,amrex_mfiter_destroy
      use mms_geom, only: initialize_volume_moments
      use amrmpcomp_class, only: VFlo
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF,pCL,pCG,pU
      real(WP), dimension(3) :: BL,BG
      real(WP) :: dx,dy,dz,myVF,IEL,x_cc,rhoG,pG,H
      integer :: i,j,k
      integer, parameter :: nref=3
      ! Get mesh size
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl); dz=solver%amr%dz(lvl)
      ! Get internal energy of liquid
      IEL=water%get_e_from_p_rho(p=pL1,rho=rhoL1,y=[1.0_WP])
      ! Use passed ba/dm since grid is being constructed
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         ! Get pointers to data
         pQ =>solver%Q%mf(lvl)%dataptr(mfi)
         pVF=>solver%VF%mf(lvl)%dataptr(mfi)
         if (lvl.eq.solver%amr%maxlvl) then
            pCL=>solver%CL%dataptr(mfi)
            pCG=>solver%CG%dataptr(mfi)
         end if
         ! Loop over grown tilebox
         bx=mfi%growntilebox(solver%nover)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Compute VF and barycenters from levelset
                  call initialize_volume_moments(lo=[solver%amr%xlo+real(i  ,WP)*dx,solver%amr%ylo+real(j  ,WP)*dy,solver%amr%zlo+real(k  ,WP)*dz], &
                  &                              hi=[solver%amr%xlo+real(i+1,WP)*dx,solver%amr%ylo+real(j+1,WP)*dy,solver%amr%zlo+real(k+1,WP)*dz], &
                  &                              levelset=sphere_levelset,time=time,level=nref,VFlo=VFlo,VF=myVF,BL=BL,BG=BG)
                  ! Store volume fraction
                  pVF(i,j,k,1)=myVF
                  ! Store barycenters
                  if (lvl.eq.solver%amr%maxlvl) then
                     pCL(i,j,k,:)=BL
                     pCG(i,j,k,:)=BG
                  end if
                  ! Compute local gas state from shock profile
                  x_cc=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
                  H=Hshock(x=Xs-x_cc,delta=0.5_WP*dx)
                  rhoG=rhoG1+(rhoG2-rhoG1)*H
                  pG  =pG1  +(pG2  -pG1  )*H
                  ! Set conserved variables: Q=(VF*rhoL, (1-VF)*rhoG, VF*rhoL*IL, (1-VF)*rhoG*IG); momentum is staggered.
                  ! Q(3:4) hold internal contents here: get_total_energy adds the kinetic energy after init_from_scratch
                  pQ(i,j,k,1)=(       myVF)*rhoL1
                  pQ(i,j,k,2)=(1.0_WP-myVF)*rhoG
                  pQ(i,j,k,3)=pQ(i,j,k,1)*IEL
                  pQ(i,j,k,4)=pQ(i,j,k,2)*gas%get_e_from_p_rho(p=pG,rho=rhoG,y=[1.0_WP])
               end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
      ! Set x-face velocity from the same shock profile, ghosts included. V and W are zero (set by on_init).
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         pU=>solver%U%mf(lvl)%dataptr(mfi)
         do k=lbound(pU,3),ubound(pU,3); do j=lbound(pU,2),ubound(pU,2); do i=lbound(pU,1),ubound(pU,1)
                  x_cc=solver%amr%xlo+real(i,WP)*dx
                  H=Hshock(x=Xs-x_cc,delta=0.5_WP*dx)
                  pU(i,j,k,1)=u1+(u2-u1)*H
               end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine shockdrop_init

   !> Apply inflow BC at low-x (face=1). Boundary momentum is not set here: apply_momentum_bc derives it from Q and U
   subroutine shock_dirichlet(solver,lvl,time,face,bx,comp,p)
      use amrex_amr_module, only: amrex_box
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      integer, intent(in) :: face
      type(amrex_box), intent(in) :: bx
      character(len=1), intent(in) :: comp
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      integer :: i,j,k
      select case (face)
       case (1)  ! X-LOW: Dirichlet inflow with post-shock (gas only, no liquid)
         select case (comp)
          case ('U')  ! Staggered U=u2
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                     p(i,j,k,1)=u2
                  end do; end do; end do
          case ('V','W')  ! Staggered V,W=0
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                     p(i,j,k,1)=0.0_WP
                  end do; end do; end do
          case ('Q')  ! Cell-centered Q=(0,rho2,0,E2); momentum is staggered and not part of Q
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                     p(i,j,k,1)=0.0_WP                  ! No liquid
                     p(i,j,k,2)=rhoG2                   ! Gas density
                     p(i,j,k,3)=0.0_WP                  ! No liquid energy
                     p(i,j,k,4)=rhoG2*gas%get_e_from_p_rho(p=pG2,rho=rhoG2,y=[1.0_WP])+0.5_WP*rhoG2*u2**2 ! Gas total energy (uniform u2)
                  end do; end do; end do
         end select
      end select
   end subroutine shock_dirichlet

   subroutine my_tagger(solver,lvl,time,tags_ptr)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      use amrtag,           only: lap_error,grd_error
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags_ptr
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx,gbx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pPL,pVF,pU,pV,pW,pC
      real(WP), dimension(:,:,:,:), allocatable :: pUVW !< Cell-centered velocity with one ghost layer, from the face velocities
      real(WP) :: dx,dy,dz,dxi,dyi,dzi,dxi2,dyi2,dzi2,delta,delta2
      real(WP) :: rho_cc,rho_xp,rho_xm,rho_yp,rho_ym,rho_zp,rho_zm
      real(WP) :: lapU,lapV,lapW,u_sgs,Re
      real(WP) :: divu,vortx,vorty,vortz,vort,Ducros,Deps
      real(WP) :: r_cyl
      logical  :: in_zone
      integer :: i,j,k
      real(WP), parameter :: Reps=1.0e-2_WP
      real(WP), parameter :: Peps=1.0e-2_WP
      real(WP), parameter :: Cduc=0.05_WP
      ! Get mesh size
      dx=solver%amr%dx(lvl); dxi=1.0_WP/dx; dxi2=1.0_WP/dx**2
      dy=solver%amr%dy(lvl); dyi=1.0_WP/dy; dyi2=1.0_WP/dy**2
      dz=solver%amr%dz(lvl); dzi=1.0_WP/dz; dzi2=1.0_WP/dz**2
      delta=solver%amr%min_meshsize(lvl); delta2=delta**2
      ! Recast tags
      tags=tags_ptr
      ! Compute tags
      call solver%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         ! Get pointers to data
         tagarr=>tags%dataPtr(mfi)
         pQ  =>solver%Q%mf(lvl)%dataptr(mfi)
         pPL =>solver%PL%mf(lvl)%dataptr(mfi)
         pVF =>solver%VF%mf(lvl)%dataptr(mfi)
         pU  =>solver%U%mf(lvl)%dataptr(mfi)
         pV  =>solver%V%mf(lvl)%dataptr(mfi)
         pW  =>solver%W%mf(lvl)%dataptr(mfi)
         pC  =>solver%C%mf(lvl)%dataptr(mfi)
         ! Solver UVW has no ghost cells: interpolate the face velocities (whose ghosts are filled) to cell centers
         ! over the tile and the one ghost layer read by the stencils below
         gbx=mfi%growntilebox(1)
         allocate(pUVW(gbx%lo(1):gbx%hi(1),gbx%lo(2):gbx%hi(2),gbx%lo(3):gbx%hi(3),3))
         do k=gbx%lo(3),gbx%hi(3); do j=gbx%lo(2),gbx%hi(2); do i=gbx%lo(1),gbx%hi(1)
                  pUVW(i,j,k,1)=0.5_WP*sum(pU(i:i+1,j,k,1))
                  pUVW(i,j,k,2)=0.5_WP*sum(pV(i,j:j+1,k,1))
                  pUVW(i,j,k,3)=0.5_WP*sum(pW(i,j,k:k+1,1))
               end do; end do; end do
         ! Loop over tile
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Refinement zone: away from sponge unless below maxlvl-1
                  r_cyl=sqrt((solver%amr%ylo+(real(j,WP)+0.5_WP)*dy)**2+(solver%amr%zlo+(real(k,WP)+0.5_WP)*dz)**2)
                  in_zone=(r_cyl.lt.R_spg+L_spg.or.lvl.lt.solver%amr%maxlvl-1)

                  ! Mixture density laplacian error
                  rho_cc=sum(pQ(i  ,j,  k,  1:2))
                  rho_xp=sum(pQ(i+1,j,  k,  1:2)); rho_xm=sum(pQ(i-1,j,  k,  1:2))
                  rho_yp=sum(pQ(i,  j+1,k,  1:2)); rho_ym=sum(pQ(i,  j-1,k,  1:2))
                  rho_zp=sum(pQ(i,  j,  k+1,1:2)); rho_zm=sum(pQ(i,  j,  k-1,1:2))
                  if (lap_error(rho_cc,rho_xm,rho_xp,rho_ym,rho_yp,rho_zm,rho_zp,Reps).gt.Rho_tag.and.in_zone) tagarr(i,j,k,1)=SETtag

                  ! Liquid pressure gradient
                  if (pVF(i,j,k,1).gt.0.0_WP) then
                     if (grd_error(pPL(i,j,k,1),pPL(i-1,j,k,1),pPL(i+1,j,k,1),pPL(i,j-1,k,1),pPL(i,j+1,k,1),pPL(i,j,k-1,1),pPL(i,j,k+1,1),Peps).gt.P_tag.and.in_zone) tagarr(i,j,k,1)=SETtag
                  end if

                  ! SGS cell Reynolds number
                  lapU=(pUVW(i+1,j,k,1)-2.0_WP*pUVW(i,j,k,1)+pUVW(i-1,j,k,1))*dxi2+(pUVW(i,j+1,k,1)-2.0_WP*pUVW(i,j,k,1)+pUVW(i,j-1,k,1))*dyi2+(pUVW(i,j,k+1,1)-2.0_WP*pUVW(i,j,k,1)+pUVW(i,j,k-1,1))*dzi2
                  lapV=(pUVW(i+1,j,k,2)-2.0_WP*pUVW(i,j,k,2)+pUVW(i-1,j,k,2))*dxi2+(pUVW(i,j+1,k,2)-2.0_WP*pUVW(i,j,k,2)+pUVW(i,j-1,k,2))*dyi2+(pUVW(i,j,k+1,2)-2.0_WP*pUVW(i,j,k,2)+pUVW(i,j,k-1,2))*dzi2
                  lapW=(pUVW(i+1,j,k,3)-2.0_WP*pUVW(i,j,k,3)+pUVW(i-1,j,k,3))*dxi2+(pUVW(i,j+1,k,3)-2.0_WP*pUVW(i,j,k,3)+pUVW(i,j-1,k,3))*dyi2+(pUVW(i,j,k+1,3)-2.0_WP*pUVW(i,j,k,3)+pUVW(i,j,k-1,3))*dzi2
                  u_sgs=0.2_WP*sqrt(lapU**2+lapV**2+lapW**2)*delta2
                  Re=Reynolds*u_sgs*delta
                  if (Re.gt.Re_tag.and.in_zone) tagarr(i,j,k,1)=SETtag

                  ! Ducros compression switch
                  divu =0.5_WP*dxi*(pUVW(i+1,j,k,1)-pUVW(i-1,j,k,1))+0.5_WP*dyi*(pUVW(i,j+1,k,2)-pUVW(i,j-1,k,2))+0.5_WP*dzi*(pUVW(i,j,k+1,3)-pUVW(i,j,k-1,3))
                  vortx=0.5_WP*dyi*(pUVW(i,j+1,k,3)-pUVW(i,j-1,k,3))-0.5_WP*dzi*(pUVW(i,j,k+1,2)-pUVW(i,j,k-1,2))
                  vorty=0.5_WP*dzi*(pUVW(i,j,k+1,1)-pUVW(i,j,k-1,1))-0.5_WP*dxi*(pUVW(i+1,j,k,3)-pUVW(i-1,j,k,3))
                  vortz=0.5_WP*dxi*(pUVW(i+1,j,k,2)-pUVW(i-1,j,k,2))-0.5_WP*dyi*(pUVW(i,j+1,k,1)-pUVW(i,j-1,k,1))
                  vort=sqrt(vortx**2+vorty**2+vortz**2)
                  Deps=(Cduc*pC(i,j,k,1)/delta)**2
                  Ducros=divu**2/max(divu**2+vort**2+Deps,tiny(1.0_WP))
                  if (divu.lt.0.0_WP.and.Ducros.gt.Ducros_tag.and.in_zone) tagarr(i,j,k,1)=SETtag
               end do; end do; end do
         deallocate(pUVW)
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine my_tagger

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      ! Initialize AMR grid
      create_amrgrid: block
         ! Set name
         amr%name='amrcomp_drop'
         ! Read in base grid size
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         ! Set domain
         amr%xlo=-05.0_WP; amr%xhi=+15.0_WP
         amr%ylo=-10.0_WP; amr%yhi=+10.0_WP
         amr%zlo=-10.0_WP; amr%zhi=+10.0_WP
         ! Set periodicity
         amr%xper=.false.; amr%yper=.true.; amr%zper=.true.
         ! Read in max level
         call param_read('Max level',amr%maxlvl)
         ! Enable quasi-2D
         if (amr%nz.eq.1) then
            amr%zlo=-0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny*2**amr%maxlvl,WP)
            amr%zhi=+0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny*2**amr%maxlvl,WP)
         end if
         ! Initialize
         call amr%initialize()
      end block create_amrgrid

      ! Read EoS and flow parameters
      init_eos_and_flow: block
         use messager, only: log,die
         use string,   only: str_long
         character(len=str_long) :: message
         real(WP) :: A,B,C
         real(WP) :: GammaL,PinfL,bL,CvL,qpL
         real(WP) :: GammaV,cvV,qV,qpV
         real(WP) :: GammaG,CvG
         real(WP) :: T_G
         ! Gas EoS parameters (ideal gas)
         call param_read('GammaG',GammaG)
         ! Liquid EoS: gamma only, PinfL is computed below
         call param_read('GammaL',GammaL)
         ! Shock parameters (gas phase, uses GammaG)
         call param_read('Gas Mach number',M2)
         call param_read('Shock location',Xs)
         ! Post-shock normalization: rhoG2=1, Deltau=1, T2=1
         rhoG2=1.0_WP
         pG2=1.0_WP/(GammaG*M2**2)
         ! Quadratic for rhoG1: A*rhoG1^2 - B*rhoG1 + C = 0
         A=2.0_WP*GammaG*pG2+(GammaG-1.0_WP)
         B=4.0_WP*GammaG*pG2+(GammaG+1.0_WP)
         C=2.0_WP*GammaG*pG2
         rhoG1=(B-sqrt(B**2-4.0_WP*A*C))/(2.0_WP*A)  ! smaller root for compression
         ! Shock-fixed frame velocities and pressure
         u1=1.0_WP/(1.0_WP-rhoG1)
         u2=u1-1.0_WP
         pG1=pG2-rhoG1/(1.0_WP-rhoG1)
         if (pG1.le.0.0_WP) call die('[simulation_init] Cannot achieve requested Mach number - negative pre-shock pressure')
         ! Shock Mach number
         Ms=u1/sqrt(GammaG*pG1/rhoG1)
         ! Lab-frame shock speed - pre-shock gas is at rest, so this is the shock-frame u1
         Ws=u1
         ! Shift to lab frame: pre-shock stationary
         u2=1.0_WP
         u1=0.0_WP
         ! CvG from T2=1
         CvG=pG2/(rhoG2*(GammaG-1.0_WP))
         ! Surface tension (set to 0 for this case)
         call param_read('Weber number',Weber)
         ! Liquid state from density ratio and liquid Mach number
         call param_read('Density ratio',density_ratio)
         call param_read('Liquid Mach number',ML)
         rhoL1=density_ratio
         pL1=pG1
         ! pL1=pG1+4.0_WP/Weber                   ! Force pressure equilibrium, accounting for 3D Laplace pressure
         ! if (amr%nz.eq.1) pL1=pG1+2.0_WP/Weber  ! Force pressure equilibrium, accounting for 2D Laplace pressure
         PinfL=rhoL1/(GammaL*ML**2)-pL1
         ! Pre-shock gas temperature (ideal gas, T = p/((gamma-1)*Cv*rho))
         T_G=pG1/(rhoG1*(GammaG-1.0_WP)*CvG)
         CvL=(pL1+PinfL)/(rhoL1*(GammaL-1.0_WP)*T_G) ! Force thermal equilibrium
         ! Build materials
         call gas%initialize  (gamma=GammaG,cv=CvG,q=0.0_WP,qp=0.0_WP,name='gas')
         call water%initialize(gamma=GammaL,pinf=PinfL,cv=CvL,q=0.0_WP,qp=0.0_WP,name='water')
         ! Viscous parameters
         call param_read('Reynolds number',Reynolds)
         call param_read('Prandtl number',Prandtl)
         call param_read('Viscosity ratio',visc_ratio)
         call param_read('Diffusivity ratio',diff_ratio)
         call param_read('Sutherland exponent',Suth_n)
         call param_read('Sutherland temperature',Suth_T)
         ! Log
         write(message,'("[Post-shock Mach] M2=",es12.5)') M2; call log(message)
         write(message,'("[Shock Mach]      Ms=",es12.5)') Ms; call log(message)
         write(message,'("[Pre-shock]  rhoG1=",es12.5," pG1=",es12.5)') rhoG1,pG1; call log(message)
         write(message,'("[Post-shock] rhoG2=",es12.5," pG2=",es12.5)') rhoG2,pG2; call log(message)
         write(message,'("[Liquid] rhoL1=",es12.5," pL1=",es12.5," ML=",es12.5)') rhoL1,pL1,ML; call log(message)
         write(message,'("[Temp]   TL=",es12.5," TG=",es12.5)') water%get_T_from_p_rho(p=pL1,rho=rhoL1,y=[1.0_WP]),T_G; call log(message)
         call water%print(); call gas%print()
         write(message,'("[Visc]   Re=",es12.5," mu*=",es12.5," Suth_n=",es12.5," Suth_T=",es12.5)') Reynolds,visc_ratio,Suth_n,Suth_T; call log(message)
         write(message,'("[Surface tension] We=",es12.5)') Weber; call log(message)
         ! Sponge parameters
         call param_read('Sponge radius',R_spg  ,default=3.0_WP)
         call param_read('Sponge length',L_spg  ,default=1.0_WP)
         call param_read('Sponge rate'  ,sig_spg,default=1.0_WP)
         write(message,'("[Sponge] R=",es12.5," L=",es12.5," sigma=",es12.5," Ws=",es12.5)') R_spg,L_spg,sig_spg,Ws; call log(message)
      end block init_eos_and_flow

      initialize_sph_modes: block
         use param,     only: param_read,param_exists
         use parallel,  only: amRoot,MPI_REAL_WP,comm
         use string,    only: str_long
         use messager,  only: log
         use mpi_f08,   only: MPI_BCAST,MPI_INTEGER
         use random,    only: random_initialize,random_uniform
         use mathtools, only: twoPi
         character(str_long) :: message
         integer  :: i,ierr
         real(WP) :: r
         ! Read number of modes (or set default)
         call param_read('SphHarm Nmode',nsh_modes,default=0)
         if (nsh_modes.gt.0) then
            allocate(l_modes(nsh_modes),m_modes(nsh_modes),amp_modes(nsh_modes),phase_modes(nsh_modes))
            if (param_exists('SphHarm l'))     call param_read('SphHarm l',    l_modes)
            if (param_exists('SphHarm m'))     call param_read('SphHarm m',    m_modes)
            if (param_exists('SphHarm amp'))   call param_read('SphHarm amp',  amp_modes)
            if (param_exists('SphHarm phase')) call param_read('SphHarm phase',phase_modes)
         else
            ! Hardcode no perturbation if no inputs are given
            nsh_modes=1
            allocate(l_modes(nsh_modes),m_modes(nsh_modes),amp_modes(nsh_modes),phase_modes(nsh_modes))
            l_modes(1)=1
            m_modes(1)=0
            amp_modes(1)=0.0_WP
            phase_modes(1)=0.0_WP
            ! ! Generate random modes if not provided
            ! nsh_modes=8
            ! allocate(l_modes(nsh_modes),m_modes(nsh_modes),amp_modes(nsh_modes),phase_modes(nsh_modes))
            ! if (amRoot) then
            !    call random_initialize()
            !    do i=1,nsh_modes
            !       l_modes    (i)=1+i
            !       m_modes    (i)=i-nsh_modes/2
            !       amp_modes  (i)=4.0e-3_WP
            !       phase_modes(i)=random_uniform(lo=0.0_WP,hi=twoPi)
            !    end do
            ! end if
            ! call MPI_BCAST(l_modes    ,nsh_modes,MPI_INTEGER,0,comm,ierr)
            ! call MPI_BCAST(m_modes    ,nsh_modes,MPI_INTEGER,0,comm,ierr)
            ! call MPI_BCAST(amp_modes  ,nsh_modes,MPI_REAL_WP,0,comm,ierr)
            ! call MPI_BCAST(phase_modes,nsh_modes,MPI_REAL_WP,0,comm,ierr)
         end if
         ! Log the selected modes
         if (amRoot) then
            write(message,'("[SphHarm] Nmode  =",i6)')         nsh_modes; call log(message)
            write(message,'("[SphHarm] l      =",1000(i4,x))') l_modes  ; call log(message)
            write(message,'("[SphHarm] m      =",1000(i4,x))') m_modes  ; call log(message)
            write(message,'("[SphHarm] amp    =",1000(es12.5,x))') amp_modes  ; call log(message)
            write(message,'("[SphHarm] phase  =",1000(es12.5,x))') phase_modes; call log(message)
         end if
      end block initialize_sph_modes

      ! Handle restart/saves here
      handle_restart: block
         integer :: restart_step
         ! Initialize IO object
         call io%initialize(amr=amr,nfiles=128)
         ! Check if restarting
         call param_read('Restart from',restart_dir,default='')
         restarted=(len_trim(restart_dir).gt.0)
         ! If restarting, read header
         if (restarted) call io%read_header(dirname=trim(restart_dir),time=restart_time,step=restart_step)
      end block handle_restart

      ! Initialize time tracker
      initialize_timetracker: block
         use messager, only: die
         time=timetracker(amRoot=amr%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max dt',time%dtmax)
         call param_read('Max CFL',time%cflmax)
         ! Iterative midpoint sub-iterations: 2 is explicit midpoint RK2, unstable on the imaginary axis;
         ! 3 or more is stable for |lambda*dt|<=2
         call param_read('Subiterations',time%itmax,default=2)
         if (time%itmax.lt.1) call die('Require at least one sub-iteration')
         call param_read('Max wall time',time%wtmax,default=huge(1.0_WP))
         time%dt=time%dtmax
         if (restarted) then
            call io%get_scalar('dt',time%dt)
            time%t=restart_time
         end if
      end block initialize_timetracker

      ! Initialize compressible multiphase solver
      create_solver: block
         use amrex_amr_module, only: amrex_bc_ext_dir,amrex_bc_foextrap
         use amrmpcomp_class,  only: BC_GAS
         ! Assign materials and create flow solver (initialize forces face-linear velocity interpolation)
         fs%liq=>water; fs%gas=>gas; call fs%initialize(amr=amr,name='drop')
         ! Set surface tension coefficient
         fs%sigma=1.0_WP/Weber; fs%sigma=0.0_WP
         ! Provide pressure relaxation model
         call relax_model%initialize(gas=gas,liq=water); fs%relax=>relax_model
         relax_model%model=Prelax
         relax_model%RHOGmin=0.0_WP
         ! Set initial conditions
         fs%user_init=>shockdrop_init
         ! Set BCs
         if (.not.amr%xper) then
            fs%lo_bc(1)=BC_GAS
            fs%Q%lo_bc(1,:)=amrex_bc_ext_dir; fs%Q%hi_bc(1,:)=amrex_bc_foextrap
            fs%U%lo_bc(1,:)=amrex_bc_ext_dir; fs%U%hi_bc(1,:)=amrex_bc_foextrap
            fs%V%lo_bc(1,:)=amrex_bc_ext_dir; fs%V%hi_bc(1,:)=amrex_bc_foextrap
            fs%W%lo_bc(1,:)=amrex_bc_ext_dir; fs%W%hi_bc(1,:)=amrex_bc_foextrap
            fs%user_bc=>shock_dirichlet
         end if
      end block create_solver

      ! Initialize workspaces
      create_workspace: block
         use amrdata_class, only: interp_none
         call dQdt%initialize(amr,name='dQdt',ncomp=fs%nQ,ng=0,interp=interp_none); call dQdt%register()
         call drhoUdt%initialize(amr,name='drhoUdt',ncomp=1,ng=0,nodal=[.true. ,.false.,.false.],interp=interp_none); call drhoUdt%register()
         call drhoVdt%initialize(amr,name='drhoVdt',ncomp=1,ng=0,nodal=[.false.,.true. ,.false.],interp=interp_none); call drhoVdt%register()
         call drhoWdt%initialize(amr,name='drhoWdt',ncomp=1,ng=0,nodal=[.false.,.false.,.true. ],interp=interp_none); call drhoWdt%register()
         call Umag%initialize(amr,name='Umag',ncomp=1    ,ng=0,interp=interp_none); call Umag%register()
         call Mach%initialize(amr,name='Mach',ncomp=1    ,ng=0,interp=interp_none); call Mach%register()
      end block create_workspace

      ! Initialize regridding
      init_regridding: block
         ! KnapSack load balancing
         amr%lb_strat=1
         ! Create regridding event
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         ! Set case-specific tagging
         fs%user_tagging=>my_tagger
         call param_read('Tag Reynolds value',Re_tag)
         call param_read('Tag density error' ,Rho_tag)
         call param_read('Tag pressure error',P_tag)
         call param_read('Tag Ducros value',Ducros_tag)
         ! Build the grid
         if (restarted) then
            ! Restore grid hierarchy from checkpoint
            call amr%init_from_checkpoint(dirname=trim(restart_dir),time=time%t)
            ! Restore solver state
            call fs%restore_checkpoint(io=io,dirname=trim(restart_dir),time=time%t)
         else
            ! Fresh start: post_regrid builds momentum from the face velocity and the primitive variables
            call amr%init_from_scratch(time=time%t)
            ! Phasic total energies from the internal contents written by user_init and the initial kinetic energy
            call fs%get_total_energy(time=time%t)
            ! Build PLIC
            call fs%build_plic(time%t)
         end if
         ! Compute viscosities
         call get_viscosities()
         ! Add SGS models
         call fs%add_viscartif(dt=time%dt,Cvisc=1.0e-2_WP)
         call fs%add_vreman(dt=time%dt)
         ! Compute cell-centered velocity, Umag, and Mach number
         call fs%get_cell_velocity()
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
      end block init_regridding

      ! Initialize checkpoint save event
      init_checkpoint: block
         ! Create checkpoint save event
         save_evt=event(time=time,name='Checkpoint')
         call param_read('Checkpoint period',save_evt%tper,default=-1.0_WP)
         ! Let solver self-register for checkpointing
         call fs%register_checkpoint(io)
         ! Add dt to checkpoint save
         call io%add_scalar(name='dt',value=time%dt)
      end block init_checkpoint

      ! Initialize visualization
      create_viz: block
         ! Create visualization object
         call viz%initialize(amr,'drop',use_hdf5=.false.)
         call viz%add_scalar(fs%VF,1,'VF')
         call viz%add_scalar(fs%RHOL,1,'RHOL')
         call viz%add_scalar(fs%RHOG,1,'RHOG')
         call viz%add_scalar(fs%PL,1,'PL')
         call viz%add_scalar(fs%PG,1,'PG')
         call viz%add_scalar(fs%TL,1,'TL')
         call viz%add_scalar(fs%TG,1,'TG')
         call viz%add_scalar(fs%UVW,1,'U')
         call viz%add_scalar(fs%UVW,2,'V')
         call viz%add_scalar(fs%UVW,3,'W')
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(Mach,1,'Mach')
         call viz%add_scalar(fs%visc,1,'visc')
         call viz%add_scalar(fs%beta,1,'beta')
         ! Create plic visualization object
         call plicviz%initialize(amr,'plic',use_hdf5=.false.)
         call plicviz%add_surfmesh(fs%smesh,'plic')
         ! Create total visualization at increments of tau of 0.25
         call tau_viz%initialize(amr,'tau_viz',use_hdf5=.false.)
         call tau_viz%add_scalar(fs%VF,1,'VF')
         call tau_viz%add_scalar(fs%RHOL,1,'RHOL')
         call tau_viz%add_scalar(fs%RHOG,1,'RHOG')
         call tau_viz%add_scalar(fs%PL,1,'PL')
         call tau_viz%add_scalar(fs%PG,1,'PG')
         call tau_viz%add_scalar(fs%TL,1,'TL')
         call tau_viz%add_scalar(fs%TG,1,'TG')
         call tau_viz%add_scalar(fs%UVW,1,'U')
         call tau_viz%add_scalar(fs%UVW,2,'V')
         call tau_viz%add_scalar(fs%UVW,3,'W')
         call tau_viz%add_scalar(Umag,1,'Umag')
         call tau_viz%add_scalar(Mach,1,'Mach')
         call tau_viz%add_surfmesh(fs%smesh,'plic')
         ! Create visualization output event
         viz_evt=event(time=time,name='Visualization output')
         plicviz_evt=event(time=time,name='PLIC Visualization output')
         tau_viz_evt=event(time=time,name='tau Visualization output')
         call param_read('Output period',viz_evt%tper)
         call param_read('PLIC Output period',plicviz_evt%tper)
         call param_read('tau Output period',tau_viz_evt%tper)
         ! Write initial state
         if (viz_evt%occurs()) call viz%write(time=time%t)
         if (plicviz_evt%occurs()) call plicviz%write(time=time%t)
         if (tau_viz_evt%occurs()) call tau_viz%write(time=time%t)
      end block create_viz

      ! Create monitors
      create_monitors: block
         ! Get solver info and cfl
         call fs%get_info(); Etot=fs%Qint(3)+fs%Qint(4)
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         ! Create simulation monitor
         mfile=monitor(amRoot=amr%amRoot,name='simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%RHOLmin,'rhoLmin')
         call mfile%add_column(fs%RHOLmax,'rhoLmax')
         call mfile%add_column(fs%PLmin,'PLmin')
         call mfile%add_column(fs%PLmax,'PLmax')
         call mfile%add_column(fs%TLmin,'TLmin')
         call mfile%add_column(fs%TLmax,'TLmax')
         call mfile%add_column(fs%RHOGmin,'rhoGmin')
         call mfile%add_column(fs%RHOGmax,'rhoGmax')
         call mfile%add_column(fs%PGmin,'PGmin')
         call mfile%add_column(fs%PGmax,'PGmax')
         call mfile%add_column(fs%TGmin,'TGmin')
         call mfile%add_column(fs%TGmax,'TGmax')
         call mfile%add_column(fs%VFmin,'VFmin')
         call mfile%add_column(fs%VFmax,'VFmax')
         call mfile%add_column(fs%VFint,'VFint')
         call mfile%add_column(fs%dPmax,'dPmax')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(amRoot=amr%amRoot,name='cfl')
         call cflfile%add_column(time%n,'Timestep')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(time%dt,'dt')
         call cflfile%add_column(fs%CFLc_x,'CFLc_x')
         call cflfile%add_column(fs%CFLc_y,'CFLc_y')
         call cflfile%add_column(fs%CFLc_z,'CFLc_z')
         call cflfile%add_column(fs%CFLa_x,'CFLa_x')
         call cflfile%add_column(fs%CFLa_y,'CFLa_y')
         call cflfile%add_column(fs%CFLa_z,'CFLa_z')
         call cflfile%add_column(fs%CFLv_x,'CFLv_x')
         call cflfile%add_column(fs%CFLv_y,'CFLv_y')
         call cflfile%add_column(fs%CFLv_z,'CFLv_z')
         call cflfile%add_column(fs%CFLst ,'CFLst' )
         call cflfile%write()
         ! Create conservation monitor
         consfile=monitor(amRoot=amr%amRoot,name='conservation')
         call consfile%add_column(time%n,'Timestep number')
         call consfile%add_column(time%t,'Time')
         call consfile%add_column(fs%VFint,'VFint')
         call consfile%add_column(fs%Qint(1),'Liquid Mass')
         call consfile%add_column(fs%Qint(2),'Gas Mass')
         call consfile%add_column(fs%Qint(3),'Liquid TotEnergy')
         call consfile%add_column(fs%Qint(4),'Gas TotEnergy')
         call consfile%add_column(Etot,'Total energy')
         call consfile%add_column(fs%rhoUint,'U Momentum')
         call consfile%add_column(fs%rhoVint,'V Momentum')
         call consfile%add_column(fs%rhoWint,'W Momentum')
         call consfile%add_column(fs%rhoKint,'Kinetic energy')
         call consfile%write()
         ! Create grid monitor
         gridfile=monitor(amRoot=amr%amRoot,name='grid')
         call gridfile%add_column(time%n,'Timestep')
         call gridfile%add_column(time%t,'Time')
         call gridfile%add_column(amr%nlevels,'Nlvl')
         call gridfile%add_column(amr%nboxes,'Nbox')
         call gridfile%add_column(amr%ncells,'Ncell')
         call gridfile%add_column(amr%compression,'Compression')
         call gridfile%add_column(amr%maxRSS,'Maximum RSS')
         call gridfile%add_column(amr%minRSS,'Minimum RSS')
         call gridfile%add_column(amr%avgRSS,'Average RSS')
         call gridfile%write()
         ! Create timing monitor
         tfile=monitor(amRoot=amr%amRoot,name='timing')
         call tfile%add_column(time%n,'Timestep')
         call tfile%add_column(time%t,'Time')
         ! Full routine times (max across ranks = wall-clock cost)
         call tfile%add_column(fs%wtmax_dQdt,'dQdt_max')
         call tfile%add_column(fs%wtmax_plic,'plic_max')
         call tfile%add_column(fs%wtmax_relax,'relax_max')
         call tfile%add_column(fs%wtmax_visc,'visc_max')
         ! Compute loop times (max = slowest rank, min = fastest rank)
         call tfile%add_column(fs%wtmax_prim,'prim_max')
         call tfile%add_column(fs%wtmin_prim,'prim_min')
         call tfile%add_column(fs%wtmax_slflux,'slflux_max')
         call tfile%add_column(fs%wtmin_slflux,'slflux_min')
         call tfile%add_column(fs%wtmax_fv,'fv_max')
         call tfile%add_column(fs%wtmin_fv,'fv_min')
         call tfile%add_column(fs%wtmax_div,'div_max')
         call tfile%add_column(fs%wtmin_div,'div_min')
         call tfile%add_column(fs%wtmax_plicnet,'plicnet_max')
         call tfile%add_column(fs%wtmin_plicnet,'plicnet_min')
         call tfile%add_column(fs%wtmax_polygon,'polygon_max')
         call tfile%add_column(fs%wtmin_polygon,'polygon_min')
         call tfile%add_column(fs%nmixed_max,'mixed_max')
         call tfile%add_column(fs%nmixed_min,'mixed_min')
         call tfile%write()
         ! Create rescue monitor (solver counters of the last step: relaxation failures, VF guard, slivers, positivity gate)
         rescfile=monitor(amRoot=amr%amRoot,name='rescue')
         call rescfile%add_column(time%n,'Timestep')
         call rescfile%add_column(time%t,'Time')
         call rescfile%add_column(fs%nrelax_fail,'Relax failures')
         call rescfile%add_column(fs%nvf_guard,'VF guard cells')
         call rescfile%add_column(fs%nsliver,'Sliver cells')
         call rescfile%add_column(fs%ngate,'Gate cells')
         call rescfile%add_column(fs%gate_mass,'Gate mass')
         call rescfile%add_column(fs%gate_energy,'Gate energy')
         call rescfile%add_column(fs%gate_rhoU,'Gate U momentum')
         call rescfile%add_column(fs%gate_rhoV,'Gate V momentum')
         call rescfile%add_column(fs%gate_rhoW,'Gate W momentum')
         call rescfile%write()
         ! Create droplet QOI monitor
         drop_QOIs=monitor(amRoot=amr%amRoot,name='drop_QOIs')
         call drop_QOIs%add_column(time%n,'Timestep')
         call drop_QOIs%add_column(time%t,'time')
         call drop_QOIs%add_column(MLE_VF001,'MLE_VF001')
         call drop_QOIs%add_column(MLE_VF01,'MLE_VF01')
         call drop_QOIs%add_column(MLE_VF05,'MLE_VF05')
         call drop_QOIs%add_column(drop_massVF001,'massVF001')
         call drop_QOIs%add_column(drop_massVF01,'massVF01')
         call drop_QOIs%add_column(MOI_xx,'Ixx')
         call drop_QOIs%add_column(MOI_yy,'Iyy')
         call drop_QOIs%add_column(MOI_zz,'Izz')
         call drop_QOIs%add_column(MOI_xy,'Ixy')
         call drop_QOIs%add_column(MOI_xz,'Ixz')
         call drop_QOIs%add_column(MOI_yz,'Iyz')
         call drop_QOIs%add_column(MOI_xx_001,'Ixx001')
         call drop_QOIs%add_column(MOI_yy_001,'Iyy001')
         call drop_QOIs%add_column(MOI_zz_001,'Izz001')
         call drop_QOIs%add_column(MOI_xy_001,'Ixy001')
         call drop_QOIs%add_column(MOI_xz_001,'Ixz001')
         call drop_QOIs%add_column(MOI_yz_001,'Iyz001')
         call drop_QOIs%write()
      end block create_monitors

   end subroutine simulation_init

   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Remember old state
         call fs%store_old()

         ! ======================= Iterative midpoint sub-iterations =======================
         do while (time%it.le.time%itmax)

            ! Midpoint state
            call fs%prepare_midpoint()

            ! Surface tension pressure jump at the midpoint
            call fs%add_surface_tension_jump()

            ! RHS at the midpoint: semi-Lagrangian band transport from the old state and bulk fluxes
            call fs%get_dQdt(dQdt=dQdt,drhoUdt=drhoUdt,drhoVdt=drhoVdt,drhoWdt=drhoWdt,dt=time%dt,time=time%t)

            ! Update conserved variables from the old state
            call fs%Q%lincomb   (a=1.0_WP,src1=fs%Qold   ,b=time%dt,src2=dQdt)
            call fs%rhoU%lincomb(a=1.0_WP,src1=fs%rhoUold,b=time%dt,src2=drhoUdt)
            call fs%rhoV%lincomb(a=1.0_WP,src1=fs%rhoVold,b=time%dt,src2=drhoVdt)
            call fs%rhoW%lincomb(a=1.0_WP,src1=fs%rhoWold,b=time%dt,src2=drhoWdt)
            ! Absorb outgoing waves in the sponge
            call apply_sponge(dt=time%dt,t=time%t)
            call fs%Q%average_down(); call fs%Q%fill(time=time%t)

            ! Phasic total energies of the new iterate: transported internal energy, phasic work, and kinetic energy
            call fs%get_energy(dt=time%dt,time=time%t)

            ! Rebuild PLIC of the new iterate (the band sweep only uses the old PLIC)
            call fs%build_plic(time=time%t)

            ! Get most up-to-date pressure: the last sub-iteration's relaxation is the end-of-step one
            call fs%apply_relax(dt=time%dt,time=time%t)

            ! Get primitive variables of the new iterate
            call fs%get_primitive(Q=fs%Q,rhoU=fs%rhoU,rhoV=fs%rhoV,rhoW=fs%rhoW,time=time%t)

            ! Increment sub-iteration counter
            time%it=time%it+1

         end do
         ! =================================================================================

         ! Regrid if event triggers
         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0,time=time%t)
            call gridfile%write()
         end if

         ! Compute viscosities
         call get_viscosities()

         ! Add SGS models
         call fs%add_viscartif(dt=time%dt,Cvisc=1.0e-2_WP)
         call fs%add_vreman(dt=time%dt)

         ! Compute cell-centered velocity, Umag, and Mach number
         call fs%get_cell_velocity()
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)

         ! Visualization output
         if (viz_evt%occurs()) call viz%write(time=time%t)
         ! Regrid clears the PLIC surface mesh when it remakes the finest level, and the relaxation moves the planes
         ! after build_plic extracted the polygons: rebuild them from the current PLIC before writing them
         if (plicviz_evt%occurs().or.tau_viz_evt%occurs()) call fs%build_polygons()
         if (plicviz_evt%occurs()) call plicviz%write(time=time%t)
         if (tau_viz_evt%occurs()) call tau_viz%write(time=time%t)

         ! Checkpoint save
         if (save_evt%occurs()) then
            save_checkpoint: block
               use string, only: rtoa
               call io%write(dirname='restart/drop_'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
            end block save_checkpoint
         end if

         ! Perform and output monitoring
         call fs%get_info(); Etot=fs%Qint(3)+fs%Qint(4)
         call mfile%write()
         call consfile%write()
         call cflfile%write()
         call tfile%write()
         call rescfile%write()

         ! Compute droplet metrics on a safe temporary copy of VF
         call compute_drop_metrics()
         call drop_QOIs%write()

      end do

      ! Force a final checkpoint on exit
      final_checkpoint: block
         use string, only: rtoa
         call io%write(dirname='restart/drop_'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
      end block final_checkpoint

   contains

      !> Compute droplet quantities of interest (MLE, mass, MOI)
      !> All operations are performed directly on fs%VF.
      subroutine compute_drop_metrics()
         use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
         use amrex_interface,  only: amrmask_make_fine
         use parallel,         only: MPI_REAL_WP,comm
         use mpi_f08,          only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM,MPI_MIN
         implicit none
         integer :: lvl,i,j,k,ierr
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         type(amrex_imultifab) :: mask
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pRHOL
         integer,  dimension(:,:,:,:), contiguous, pointer :: pMask
         real(WP) :: dx,dy,dz,vol,x_cc,y_cc,z_cc,vf_local,mi
         ! MOI accumulators (raw second moments about origin)
         real(WP) :: M_tot,Sx,Sy,Sz,Sxx,Syy,Szz,Sxy,Sxz,Syz,xc,yc,zc
         real(WP) :: M_tot_001,Sx_001,Sy_001,Sz_001,Sxx_001,Syy_001,Szz_001,Sxy_001,Sxz_001,Syz_001,xc_001,yc_001,zc_001
         ! Mass accumulators
         real(WP) :: mass_001,mass_01

         ! Initialize accumulators
         MLE_VF001=+huge(1.0_WP); MLE_VF01=+huge(1.0_WP); MLE_VF05=+huge(1.0_WP)
         mass_001=0.0_WP; mass_01=0.0_WP
         M_tot=0.0_WP; Sx=0.0_WP; Sy=0.0_WP; Sz=0.0_WP
         Sxx=0.0_WP; Syy=0.0_WP; Szz=0.0_WP
         Sxy=0.0_WP; Sxz=0.0_WP; Syz=0.0_WP
         M_tot_001=0.0_WP; Sx_001=0.0_WP; Sy_001=0.0_WP; Sz_001=0.0_WP
         Sxx_001=0.0_WP; Syy_001=0.0_WP; Szz_001=0.0_WP
         Sxy_001=0.0_WP; Sxz_001=0.0_WP; Syz_001=0.0_WP

         ! Loop over all levels
         do lvl=0,amr%clvl()
            dx=amr%dx(lvl); dy=amr%dy(lvl); dz=amr%dz(lvl)
            vol=dx*dy*dz
            ! Build fine mask for non-finest levels
            if (lvl.lt.amr%clvl()) then
               call amrex_imultifab_build(mask,amr%ba(lvl),amr%dm(lvl),1,0)
               call amrmask_make_fine(mask,amr%ba(lvl+1),[amr%rrefx(lvl),amr%rrefy(lvl),amr%rrefz(lvl)],0,1)
            end if
            call amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               pVF=>fs%VF%mf(lvl)%dataptr(mfi)
               pRHOL=>fs%RHOL%mf(lvl)%dataptr(mfi)
               if (lvl.lt.amr%clvl()) pMask=>mask%dataptr(mfi)
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                        ! Skip cells covered by finer levels
                        if (lvl.lt.amr%clvl()) then
                           if (pMask(i,j,k,1).eq.0) cycle
                        end if
                        vf_local=pVF(i,j,k,1)

                        ! Cell-center coordinates
                        x_cc=amr%xlo+(real(i,WP)+0.5_WP)*dx
                        y_cc=amr%ylo+(real(j,WP)+0.5_WP)*dy
                        z_cc=amr%zlo+(real(k,WP)+0.5_WP)*dz

                        ! MLE: track windward leading edge at VF thresholds
                        if (vf_local.ge.0.01_WP) MLE_VF001=min(MLE_VF001,x_cc-0.5_WP*dx)
                        if (vf_local.ge.0.1_WP)  MLE_VF01 =min(MLE_VF01, x_cc-0.5_WP*dx)
                        if (vf_local.ge.0.5_WP)  MLE_VF05 =min(MLE_VF05, x_cc-0.5_WP*dx)

                        ! Accumulate mass and moments
                        if (vf_local.ge.0.01_WP) then
                           mass_001=mass_001+vf_local*pRHOL(i,j,k,1)*vol
                           mi=vf_local*pRHOL(i,j,k,1)*vol
                           M_tot_001=M_tot_001+mi
                           Sx_001=Sx_001+mi*x_cc; Sy_001=Sy_001+mi*y_cc; Sz_001=Sz_001+mi*z_cc
                           Sxx_001=Sxx_001+mi*x_cc**2; Syy_001=Syy_001+mi*y_cc**2; Szz_001=Szz_001+mi*z_cc**2
                           Sxy_001=Sxy_001+mi*x_cc*y_cc; Sxz_001=Sxz_001+mi*x_cc*z_cc; Syz_001=Syz_001+mi*y_cc*z_cc
                        end if
                        if (vf_local.ge.0.1_WP) then
                           mass_01=mass_01+vf_local*pRHOL(i,j,k,1)*vol
                           mi=vf_local*pRHOL(i,j,k,1)*vol
                           M_tot=M_tot+mi
                           Sx=Sx+mi*x_cc; Sy=Sy+mi*y_cc; Sz=Sz+mi*z_cc
                           Sxx=Sxx+mi*x_cc**2; Syy=Syy+mi*y_cc**2; Szz=Szz+mi*z_cc**2
                           Sxy=Sxy+mi*x_cc*y_cc; Sxz=Sxz+mi*x_cc*z_cc; Syz=Syz+mi*y_cc*z_cc
                        end if
                     end do; end do; end do
            end do
            call amr%mfiter_destroy(mfi)
            if (lvl.lt.amr%clvl()) call amrex_imultifab_destroy(mask)
         end do

         ! MPI reductions
         call MPI_ALLREDUCE(MPI_IN_PLACE,MLE_VF001,1,MPI_REAL_WP,MPI_MIN,comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,MLE_VF01, 1,MPI_REAL_WP,MPI_MIN,comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,MLE_VF05, 1,MPI_REAL_WP,MPI_MIN,comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,mass_001,1,MPI_REAL_WP,MPI_SUM,comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,mass_01, 1,MPI_REAL_WP,MPI_SUM,comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,M_tot,1,MPI_REAL_WP,MPI_SUM,comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Sx,1,MPI_REAL_WP,MPI_SUM,comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Sy,1,MPI_REAL_WP,MPI_SUM,comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Sz,1,MPI_REAL_WP,MPI_SUM,comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Sxx,1,MPI_REAL_WP,MPI_SUM,comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Syy,1,MPI_REAL_WP,MPI_SUM,comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Szz,1,MPI_REAL_WP,MPI_SUM,comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Sxy,1,MPI_REAL_WP,MPI_SUM,comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Sxz,1,MPI_REAL_WP,MPI_SUM,comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Syz,1,MPI_REAL_WP,MPI_SUM,comm,ierr)

         call MPI_ALLREDUCE(MPI_IN_PLACE,M_tot_001,1,MPI_REAL_WP,MPI_SUM,comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Sx_001,1,MPI_REAL_WP,MPI_SUM,comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Sy_001,1,MPI_REAL_WP,MPI_SUM,comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Sz_001,1,MPI_REAL_WP,MPI_SUM,comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Sxx_001,1,MPI_REAL_WP,MPI_SUM,comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Syy_001,1,MPI_REAL_WP,MPI_SUM,comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Szz_001,1,MPI_REAL_WP,MPI_SUM,comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Sxy_001,1,MPI_REAL_WP,MPI_SUM,comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Sxz_001,1,MPI_REAL_WP,MPI_SUM,comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Syz_001,1,MPI_REAL_WP,MPI_SUM,comm,ierr)

         ! Store mass
         drop_massVF001=mass_001
         drop_massVF01 =mass_01

         ! Compute inertia tensor about VF-weighted COM via parallel axis theorem
         if (M_tot.gt.0.0_WP) then
            xc=Sx/M_tot; yc=Sy/M_tot; zc=Sz/M_tot
            ! Diagonal: Ixx = sum(m*(y^2+z^2)) - M*(yc^2+zc^2)
            MOI_xx=(Syy+Szz)-M_tot*(yc**2+zc**2)
            MOI_yy=(Sxx+Szz)-M_tot*(xc**2+zc**2)
            MOI_zz=(Sxx+Syy)-M_tot*(xc**2+yc**2)
            ! Off-diagonal: Ixy = -(sum(m*x*y) - M*xc*yc)
            MOI_xy=-(Sxy-M_tot*xc*yc)
            MOI_xz=-(Sxz-M_tot*xc*zc)
            MOI_yz=-(Syz-M_tot*yc*zc)
         else
            MOI_xx=0.0_WP; MOI_yy=0.0_WP; MOI_zz=0.0_WP
            MOI_xy=0.0_WP; MOI_xz=0.0_WP; MOI_yz=0.0_WP
         end if

         if (M_tot_001.gt.0.0_WP) then
            xc_001=Sx_001/M_tot_001; yc_001=Sy_001/M_tot_001; zc_001=Sz_001/M_tot_001
            ! Diagonal: Ixx = sum(m*(y^2+z^2)) - M*(yc^2+zc^2)
            MOI_xx_001=(Syy_001+Szz_001)-M_tot_001*(yc_001**2+zc_001**2)
            MOI_yy_001=(Sxx_001+Szz_001)-M_tot_001*(xc_001**2+zc_001**2)
            MOI_zz_001=(Sxx_001+Syy_001)-M_tot_001*(xc_001**2+yc_001**2)
            ! Off-diagonal: Ixy = -(sum(m*x*y) - M*xc*yc)
            MOI_xy_001=-(Sxy_001-M_tot_001*xc_001*yc_001)
            MOI_xz_001=-(Sxz_001-M_tot_001*xc_001*zc_001)
            MOI_yz_001=-(Syz_001-M_tot_001*yc_001*zc_001)
         else
            MOI_xx_001=0.0_WP; MOI_yy_001=0.0_WP; MOI_zz_001=0.0_WP
            MOI_xy_001=0.0_WP; MOI_xz_001=0.0_WP; MOI_yz_001=0.0_WP
         end if

      end subroutine compute_drop_metrics

   end subroutine simulation_run

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Finalize time
      call time%finalize()
      ! Finalize grid
      call amr%finalize()
      call regrid_evt%finalize()
      ! Finalize solver
      call fs%finalize()
      call dQdt%finalize()
      call drhoUdt%finalize()
      call drhoVdt%finalize()
      call drhoWdt%finalize()
      call Umag%finalize()
      call Mach%finalize()
      ! Finalize materials
      call water%finalize()
      call gas%finalize()
      ! Finalize visualization
      call viz%finalize()
      call viz_evt%finalize()
      call plicviz%finalize()
      call plicviz_evt%finalize()
      call tau_viz%finalize()
      call tau_viz_evt%finalize()
      ! Finalize checkpoint
      call save_evt%finalize()
      call io%finalize()
      ! Finalize monitoring
      call mfile%finalize()
      call cflfile%finalize()
      call consfile%finalize()
      call gridfile%finalize()
      call tfile%finalize()
      call rescfile%finalize()
   end subroutine simulation_final

end module simulation
