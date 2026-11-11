!> Sembian blastwave case: a high-pressure gas slab releases a blast wave onto a water cylinder (Sembian et al. 2016)
!> Dimensional (SI) setup on an AMR grid in a closed slip-wall channel
module simulation
   use precision,           only: WP
   use string,              only: str_medium
   use amrgrid_class,       only: amrgrid
   use amrmpcomp_class, only: amrmpcomp
   use amrviz_class,        only: amrviz
   use amrdata_class,       only: amrdata
   use timetracker_class,   only: timetracker
   use event_class,         only: event
   use monitor_class,       only: monitor
   use amrio_class,         only: amrio
   use nasg_class,          only: nasg
   use ideal_gas_class,     only: ideal_gas
   use relax_ig_nasg_class, only: relax_ig_nasg,Prelax
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
   type(event) :: viz_evt
   type(amrviz) :: viz

   ! Regrid parameters
   type(event) :: regrid_evt

   ! Restart parameters
   type(amrio) :: io
   type(event) :: save_evt
   character(len=str_medium) :: restart_dir
   logical :: restarted
   real(WP) :: restart_time

   !> Simulation monitoring
   type(monitor) :: mfile,consfile,cflfile,gridfile,tfile

   !> Droplet QOI monitoring
   type(monitor) :: drop_QOIs
   real(WP) :: MLE_VF001,MLE_VF01,MLE_VF05                            !< Mist leading edge at VF>=0.01 and VF>=0.1 thresholds
   real(WP) :: drop_massVF001,drop_massVF01                           !< Droplet mass at VF>=0.01 and VF>=0.1 thresholds (smoothed)
   real(WP) :: MOI_xx,MOI_yy,MOI_zz,MOI_xy,MOI_xz,MOI_yz              !< Moment of inertia tensor about VF-weighted COM (VF>=0.1)
   real(WP) :: MOI_xx_001,MOI_yy_001,MOI_zz_001,MOI_xy_001,MOI_xz_001,MOI_yz_001 !< Moment of inertia tensor about VF-weighted COM (VF>=0.01)

   !> Mixture total energy, the sum of the phasic total energies
   real(WP) :: Etot=0.0_WP

   !> Materials
   type(nasg),      target :: water
   type(ideal_gas), target :: gas

   !> Relaxation model
   type(relax_ig_nasg), target :: relax_model

   !> Flow parameters - all dimensional (SI)
   real(WP) :: rhoG1,pG1              !< Ambient (pre-shock) gas state at rest [kg/m^3,Pa]
   real(WP) :: rhoHP,pHP              !< High-pressure gas state at rest [kg/m^3,Pa]
   real(WP) :: xHP,wHP                !< High-pressure region center and thickness [m]
   real(WP) :: rhoL1,pL1              !< Initial liquid state [kg/m^3,Pa]
   real(WP) :: d_drop,x_drop          !< Drop diameter [m] and streamwise center location [m]
   real(WP) :: visc_l,visc_g          !< Liquid viscosity and gas Sutherland reference viscosity [Pa.s]
   real(WP) :: Suth_T0,Suth_S         !< Sutherland reference and Sutherland temperatures [K]
   real(WP) :: cond_g,cond_l          !< Gas and liquid thermal conductivities [W/(m K)]
   real(WP) :: sigma                  !< Surface tension coefficient [N/m]
   real(WP) :: T1                     !< Ambient (pre-shock) gas temperature [K]

   !> Tagging parameters
   real(WP) :: Re_tag=huge(1.0_WP)
   real(WP) :: Rho_tag=huge(1.0_WP)
   real(WP) :: P_tag=huge(1.0_WP)
   real(WP) :: Ducros_tag=huge(1.0_WP)

   !> Wall-clock timers of this rank [s], written to the timing monitor (solver-internal timers are in amrmpcomp)
   real(WP) :: wt_total=0.0_WP        !< Whole time step
   real(WP) :: wt_cfl=0.0_WP          !< CFL evaluation and time step adjustment
   real(WP) :: wt_iter=0.0_WP         !< Sub-iteration loop
   real(WP) :: wt_regrid=0.0_WP       !< Regridding
   real(WP) :: wt_post=0.0_WP         !< Viscosities and SGS models after the sub-iterations
   real(WP) :: wt_out=0.0_WP          !< Visualization and checkpoint output
   real(WP) :: wt_info=0.0_WP         !< Solver info and monitor writes
   real(WP) :: wt_metrics=0.0_WP      !< Droplet metrics

   !> Artificial bulk viscosity coefficient and viscous CFL clipping
   real(WP) :: Cartif,max_visc_cfl

contains

   !> Levelset function for the drop (cylinder in quasi-2D, sphere in 3D)
   function sphere_levelset(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=0.5_WP*d_drop-sqrt((xyz(1)-x_drop)**2+xyz(2)**2+xyz(3)**2)
      if (amr%nz.eq.1) G=0.5_WP*d_drop-sqrt((xyz(1)-x_drop)**2+xyz(2)**2) ! Enable quasi-2D runs
   end function sphere_levelset

   !> Compute viscosity: Sutherland for gas, harmonic blend with liquid; constant phasic conductivities
   subroutine get_viscosities()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pTG,pVF,pVisc,pBeta,pDiffL,pDiffG
      real(WP) :: Ts,mu_g
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
            bx=mfi%growntilebox(min(fs%visc%ng,fs%beta%ng,fs%diffL%ng,fs%diffG%ng))
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                     ! Gas viscosity from Sutherland, with the temperature capped at 10 times the reference as in amrcomp_drop
                     Ts=min(pTG(i,j,k,1)/Suth_T0,10.0_WP)
                     mu_g=visc_g*(1.0_WP+Suth_S/Suth_T0)*Ts**1.5_WP/(Ts+Suth_S/Suth_T0)
                     ! Mixture viscosity
                     !pVisc(i,j,k,1)=pVF(i,j,k,1)*visc_l+(1.0_WP-pVF(i,j,k,1))*mu_g ! Arithmetic averaging
                     pVisc(i,j,k,1)=1.0_WP/(pVF(i,j,k,1)/max(visc_l,myeps)+(1.0_WP-pVF(i,j,k,1))/max(mu_g,myeps)) ! Harmonic averaging
                     ! Zero bulk viscosity
                     pBeta(i,j,k,1)=0.0_WP
                     ! Phasic thermal conductivities
                     pDiffG(i,j,k,1)=cond_g
                     pDiffL(i,j,k,1)=cond_l
                  end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine get_viscosities

   !> User init callback - set Q and VF/barycenters for a drop and gas at rest with a high-pressure gas slab
   subroutine blastwave_init(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box
      use amrex_amr_module, only: amrex_mfiter_build,amrex_mfiter_destroy
      use mms_geom, only: initialize_volume_moments
      use amrmpcomp_class, only: VFlo
      use mathtools, only: twoPi
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF,pCL,pCG,pC
      real(WP), dimension(3) :: BL,BG
      real(WP) :: dx,dy,dz,myVF,IEL,x_cc,rhoG,pG,cL,cG
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
                  ! Local gas state: high-pressure slab if the cell center lies inside it, ambient otherwise
                  x_cc=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
                  if (abs(x_cc-xHP).lt.0.5_WP*wHP) then
                     rhoG=rhoHP; pG=pHP
                  else
                     rhoG=rhoG1; pG=pG1
                  end if
                  ! Set conserved variables: Q=(VF*rhoL, (1-VF)*rhoG, VF*rhoL*IL, (1-VF)*rhoG*IG); momentum is staggered
                  pQ(i,j,k,1)=(       myVF)*rhoL1
                  pQ(i,j,k,2)=(1.0_WP-myVF)*rhoG
                  pQ(i,j,k,3)=pQ(i,j,k,1)*IEL
                  pQ(i,j,k,4)=pQ(i,j,k,2)*gas%get_e_from_p_rho(p=pG,rho=rhoG,y=[1.0_WP])
               end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
      ! Face velocities stay zero (set by on_init); face momentum is built from them and Q in post_regrid
      ! Set mixture sound speed: the tagger's Ducros sensor reads it before get_primitive first sets it
      cL=water%get_c_from_p_rho(p=pL1,rho=rhoL1,y=[1.0_WP])
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         pQ=>solver%Q%mf(lvl)%dataptr(mfi)
         pC=>solver%C%mf(lvl)%dataptr(mfi)
         do k=lbound(pC,3),ubound(pC,3); do j=lbound(pC,2),ubound(pC,2); do i=lbound(pC,1),ubound(pC,1)
                  x_cc=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
                  if (abs(x_cc-xHP).lt.0.5_WP*wHP) then
                     rhoG=rhoHP; pG=pHP
                  else
                     rhoG=rhoG1; pG=pG1
                  end if
                  cG=gas%get_c_from_p_rho(p=pG,rho=rhoG,y=[1.0_WP])
                  pC(i,j,k,1)=sqrt((pQ(i,j,k,1)*cL**2+pQ(i,j,k,2)*cG**2)/max(pQ(i,j,k,1)+pQ(i,j,k,2),solver%rho_floor))
               end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine blastwave_init

   !> Tagger from amrcomp_drop: mixture density curvature, liquid pressure gradient, SGS cell Reynolds number
   !> and Ducros compression switch, with velocity derivatives taken directly from the face velocities
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
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pPL,pVF,pU,pV,pW,pC
      real(WP) :: dx,dy,dz,dxi,dyi,dzi,dxi2,dyi2,dzi2,delta,delta2
      real(WP) :: rho_cc,rho_xp,rho_xm,rho_yp,rho_ym,rho_zp,rho_zm
      real(WP) :: lapU,lapV,lapW,u_sgs,mu,Re
      real(WP) :: divu,vort2,Ducros,Deps
      integer :: i,j,k,ii,jj,kk
      real(WP), parameter :: Reps=1.0e-2_WP
      real(WP), parameter :: Cduc=0.05_WP
      real(WP), parameter :: myeps=1.0e-15_WP
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
         ! Loop over tile
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Mixture density laplacian error
                  rho_cc=sum(pQ(i  ,j,  k,  1:2))
                  rho_xp=sum(pQ(i+1,j,  k,  1:2)); rho_xm=sum(pQ(i-1,j,  k,  1:2))
                  rho_yp=sum(pQ(i,  j+1,k,  1:2)); rho_ym=sum(pQ(i,  j-1,k,  1:2))
                  rho_zp=sum(pQ(i,  j,  k+1,1:2)); rho_zm=sum(pQ(i,  j,  k-1,1:2))
                  if (lap_error(rho_cc,rho_xm,rho_xp,rho_ym,rho_yp,rho_zm,rho_zp,Reps).gt.Rho_tag) tagarr(i,j,k,1)=SETtag

                  ! Liquid pressure gradient, with a noise floor of 1% of the ambient pressure
                  if (pVF(i,j,k,1).gt.0.0_WP) then
                     if (grd_error(pPL(i,j,k,1),pPL(i-1,j,k,1),pPL(i+1,j,k,1),pPL(i,j-1,k,1),pPL(i,j+1,k,1),pPL(i,j,k-1,1),pPL(i,j,k+1,1),1.0e-2_WP*pG1).gt.P_tag) tagarr(i,j,k,1)=SETtag
                  end if

                  ! SGS cell Reynolds number from the local mixture density and molecular viscosity, with each velocity
                  ! Laplacian taken on the cell's two faces normal to that component and averaged to the cell center
                  lapU=0.5_WP*(face_lap(pU,i,j,k)+face_lap(pU,i+1,j,k))
                  lapV=0.5_WP*(face_lap(pV,i,j,k)+face_lap(pV,i,j+1,k))
                  lapW=0.5_WP*(face_lap(pW,i,j,k)+face_lap(pW,i,j,k+1))
                  u_sgs=0.2_WP*sqrt(lapU**2+lapV**2+lapW**2)*delta2
                  mu=1.0_WP/(pVF(i,j,k,1)/max(visc_l,myeps)+(1.0_WP-pVF(i,j,k,1))/max(visc_g,myeps))
                  Re=rho_cc*u_sgs*delta/mu
                  if (Re.gt.Re_tag) tagarr(i,j,k,1)=SETtag

                  ! Ducros compression switch: compact divergence from the face velocities, and each curl component
                  ! squared on the four surrounding parallel edges then averaged, as in the artificial bulk viscosity
                  divu=dxi*(pU(i+1,j,k,1)-pU(i,j,k,1))+dyi*(pV(i,j+1,k,1)-pV(i,j,k,1))+dzi*(pW(i,j,k+1,1)-pW(i,j,k,1))
                  vort2=0.0_WP
                  do kk=0,1; do jj=0,1
                     vort2=vort2+(dyi*(pW(i,j+jj,k+kk,1)-pW(i,j+jj-1,k+kk,1))-dzi*(pV(i,j+jj,k+kk,1)-pV(i,j+jj,k+kk-1,1)))**2
                  end do; end do
                  do kk=0,1; do ii=0,1
                     vort2=vort2+(dzi*(pU(i+ii,j,k+kk,1)-pU(i+ii,j,k+kk-1,1))-dxi*(pW(i+ii,j,k+kk,1)-pW(i+ii-1,j,k+kk,1)))**2
                  end do; end do
                  do jj=0,1; do ii=0,1
                     vort2=vort2+(dxi*(pV(i+ii,j+jj,k,1)-pV(i+ii-1,j+jj,k,1))-dyi*(pU(i+ii,j+jj,k,1)-pU(i+ii,j+jj-1,k,1)))**2
                  end do; end do
                  vort2=0.25_WP*vort2
                  Deps=(Cduc*pC(i,j,k,1)/delta)**2
                  Ducros=divu**2/max(divu**2+vort2+Deps,tiny(1.0_WP))
                  if (divu.lt.0.0_WP.and.Ducros.gt.Ducros_tag) tagarr(i,j,k,1)=SETtag
               end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   contains
      !> 7-point Laplacian of a face velocity component at face (fi,fj,fk)
      real(WP) function face_lap(p,fi,fj,fk)
         real(WP), dimension(:,:,:,:), contiguous, pointer, intent(in) :: p
         integer, intent(in) :: fi,fj,fk
         face_lap=(p(fi+1,fj,fk,1)-2.0_WP*p(fi,fj,fk,1)+p(fi-1,fj,fk,1))*dxi2 &
         &       +(p(fi,fj+1,fk,1)-2.0_WP*p(fi,fj,fk,1)+p(fi,fj-1,fk,1))*dyi2 &
         &       +(p(fi,fj,fk+1,1)-2.0_WP*p(fi,fj,fk,1)+p(fi,fj,fk-1,1))*dzi2
      end function face_lap
   end subroutine my_tagger

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      ! Read EoS and flow parameters - everything below is dimensional (SI)
      init_eos_and_flow: block
         use messager, only: log,die
         use string,   only: str_long
         character(len=str_long) :: message
         real(WP) :: GammaL,PinfL,CvL
         real(WP) :: GammaG,CvG
         real(WP) :: Ms,Us,rhoG2,pG2,u2
         ! Gas EoS parameters (stiffened gas with Pinf=0, i.e. ideal gas)
         call param_read('Gas gamma',GammaG)
         call param_read('Gas Cv'   ,CvG)
         ! Liquid EoS parameters (stiffened gas)
         call param_read('Liquid gamma',GammaL)
         call param_read('Liquid Pinf' ,PinfL)
         call param_read('Liquid Cv'   ,CvL)
         ! Pre-shock (ambient) gas state - stationary
         call param_read('Pre-shock density' ,rhoG1)
         call param_read('Pre-shock pressure',pG1)
         ! High-pressure gas slab - stationary
         call param_read('HP density'         ,rhoHP)
         call param_read('HP pressure'        ,pHP)
         call param_read('HP center location' ,xHP)
         call param_read('HP region thickness',wHP)
         ! Drop geometry
         call param_read('Drop diameter',d_drop)
         call param_read('Drop location',x_drop)
         ! Initial liquid state - at rest and in mechanical equilibrium with the ambient gas
         call param_read('Liquid density',rhoL1)
         pL1=pG1
         if (pL1+PinfL.le.0.0_WP) call die('[simulation_init] Non-physical liquid state - pL1+PinfL must be positive')
         ! Pre-shock gas temperature (ideal gas, T=p/((gamma-1)*Cv*rho))
         T1=pG1/(rhoG1*(GammaG-1.0_WP)*CvG)
         ! Liquid Cv forced so that the drop starts in thermal equilibrium with the ambient gas
         ! CvL=(pL1+PinfL)/(rhoL1*(GammaL-1.0_WP)*T1)
         ! Build materials
         call gas%initialize  (gamma=GammaG,cv=CvG,q=0.0_WP,qp=0.0_WP,name='air')
         call water%initialize(gamma=GammaL,pinf=PinfL,cv=CvL,q=0.0_WP,qp=0.0_WP,name='water')
         ! Diagnostics: Rankine-Hugoniot state behind the target shock Mach number (logged only, not imposed)
         call param_read('Mach number of shock',Ms,default=0.0_WP)
         pG2=pG1*(2.0_WP*GammaG*Ms**2-(GammaG-1.0_WP))/(GammaG+1.0_WP)
         rhoG2=rhoG1*(GammaG+1.0_WP)*Ms**2/((GammaG-1.0_WP)*Ms**2+2.0_WP)
         Us=Ms*sqrt(GammaG*pG1/rhoG1)
         u2=Us*(1.0_WP-rhoG1/max(rhoG2,epsilon(1.0_WP)))
         ! Liquid viscosity, and the gas Sutherland reference viscosity and temperatures
         call param_read('Liquid viscosity',visc_l)
         call param_read('Gas viscosity'   ,visc_g)
         call param_read('Sutherland reference temperature',Suth_T0)
         call param_read('Sutherland temperature'          ,Suth_S)
         ! Thermal conductivities
         call param_read('Gas conductivity'   ,cond_g)
         call param_read('Liquid conductivity',cond_l)
         ! Surface tension
         call param_read('Surface tension coefficient',sigma,default=0.0_WP)
         ! Log
         write(message,'("[Pre-shock]  rhoG1=",es12.5," pG1=",es12.5)') rhoG1,pG1; call log(message)
         write(message,'("[HP region]  rhoHP=",es12.5," pHP=",es12.5," center=",es12.5," m  thickness=",es12.5," m")') rhoHP,pHP,xHP,wHP; call log(message)
         if (Ms.gt.0.0_WP) then
            write(message,'("[Target shock] Ms=",es12.5," Us=",es12.5," m/s")') Ms,Us; call log(message)
            write(message,'("[Post-shock] rhoG2=",es12.5," pG2=",es12.5," u2=",es12.5)') rhoG2,pG2,u2; call log(message)
         end if
         write(message,'("[Liquid] rhoL1=",es12.5," pL1=",es12.5," CvL=",es12.5)') rhoL1,pL1,CvL; call log(message)
         write(message,'("[Temp]   TL=",es12.5," TG=",es12.5)') water%get_T_from_p_rho(p=pL1,rho=rhoL1,y=[1.0_WP]),T1; call log(message)
         call water%print(); call gas%print()
         write(message,'("[Visc]   muL=",es12.5," muG=",es12.5," Pa.s")') visc_l,visc_g; call log(message)
         write(message,'("[Surface tension] sigma=",es12.5," N/m")') sigma; call log(message)
         if (sigma.gt.0.0_WP.and.Ms.gt.0.0_WP) then
            write(message,'("[Scales] We=",es12.5)') rhoG2*u2**2*d_drop/sigma; call log(message)
         end if
      end block init_eos_and_flow

      ! Initialize the AMR grid
      create_amrgrid: block
         use messager, only: log
         use string,   only: str_long
         character(len=str_long) :: message
         real(WP) :: Lx,Ly
         ! Set name
         amr%name='amr_sembian_blastwave'
         ! Read in max level and base grid size
         call param_read('Max level',amr%maxlvl)
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         ! Set domain: x in [0,Lx], y in [-Ly/2,+Ly/2]
         call param_read('Lx',Lx)
         call param_read('Ly',Ly)
         amr%xlo=0.0_WP;        amr%xhi=+Lx
         amr%ylo=-0.5_WP*Ly;    amr%yhi=+0.5_WP*Ly
         amr%zlo=-0.5_WP*Ly;    amr%zhi=+0.5_WP*Ly
         ! Set periodicity: all directions are slip walls, except z in quasi-2D runs
         amr%xper=.false.; amr%yper=.false.; amr%zper=.false.
         ! Enable quasi-2D: one finest cell deep, dz=dy on the finest level, periodic
         if (amr%nz.eq.1) then
            amr%zlo=-0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny*2**amr%maxlvl,WP)
            amr%zhi=+0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny*2**amr%maxlvl,WP)
            amr%zper=.true.
         end if
         ! Initialize
         call amr%initialize()
         ! Log drop resolution on the finest level (needs the initialized grid)
         write(message,'("[Drop]   d0=",es12.5," m  x0=",es12.5," m  CPD=",es12.5)') d_drop,x_drop,d_drop/amr%dy(amr%maxlvl); call log(message)
      end block create_amrgrid

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
         use amrex_amr_module, only: amrex_bc_reflect_even,amrex_bc_reflect_odd
         use amrmpcomp_class, only: BC_REFLECT
         ! Surface tension coefficient
         fs%sigma=sigma
         ! Band momentum flux: CTU below the staggered phase fraction CTU alpha, MC elsewhere. Set before initialize, which validates them
         call param_read('CTU alpha',fs%alpha_ctu,default=0.5_WP)
         ! Assign materials and create flow solver (initialize forces face-linear velocity interpolation)
         fs%liq=>water; fs%gas=>gas; call fs%initialize(amr=amr,name='drop')
         ! Provide pressure relaxation model (mechanical only, no thermal relaxation)
         call relax_model%initialize(gas=gas,liq=water); fs%relax=>relax_model
         relax_model%model=Prelax
         ! Refuse relaxations that leave the gas below 1e-2; the model's own RHOGmin keeps its 1e-2 default
         fs%rhoG_min=1.0e-2_WP
         ! Positivity gate temperature floors (0: undeclared, a pure cell outside its law's domain stops the run)
         call param_read('Liquid Tmin',fs%TminL,default=0.0_WP)
         call param_read('Gas Tmin'   ,fs%TminG,default=0.0_WP)
         ! Set initial conditions
         fs%user_init=>blastwave_init
         ! Set BCs: slip walls on every non-periodic boundary. Staggered solver conventions:
         !  - reflect_odd on the normal velocity zeroes the wall face and mirrors the ghosts; apply_momentum_bc then zeroes the wall face momentum
         !  - reflect_even on the tangential velocities and on Q (masses and energies only, momentum is staggered) mirrors the ghosts
         ! VOF: mirrored across the walls
         if (.not.amr%xper) then
            fs%lo_bc(1)=BC_REFLECT; fs%hi_bc(1)=BC_REFLECT
            fs%Q%lo_bc(1,:)=amrex_bc_reflect_even; fs%Q%hi_bc(1,:)=amrex_bc_reflect_even
            fs%U%lo_bc(1,:)=amrex_bc_reflect_odd ; fs%U%hi_bc(1,:)=amrex_bc_reflect_odd
            fs%V%lo_bc(1,:)=amrex_bc_reflect_even; fs%V%hi_bc(1,:)=amrex_bc_reflect_even
            fs%W%lo_bc(1,:)=amrex_bc_reflect_even; fs%W%hi_bc(1,:)=amrex_bc_reflect_even
         end if
         if (.not.amr%yper) then
            fs%lo_bc(2)=BC_REFLECT; fs%hi_bc(2)=BC_REFLECT
            fs%Q%lo_bc(2,:)=amrex_bc_reflect_even; fs%Q%hi_bc(2,:)=amrex_bc_reflect_even
            fs%U%lo_bc(2,:)=amrex_bc_reflect_even; fs%U%hi_bc(2,:)=amrex_bc_reflect_even
            fs%V%lo_bc(2,:)=amrex_bc_reflect_odd ; fs%V%hi_bc(2,:)=amrex_bc_reflect_odd
            fs%W%lo_bc(2,:)=amrex_bc_reflect_even; fs%W%hi_bc(2,:)=amrex_bc_reflect_even
         end if
         if (.not.amr%zper) then
            fs%lo_bc(3)=BC_REFLECT; fs%hi_bc(3)=BC_REFLECT
            fs%Q%lo_bc(3,:)=amrex_bc_reflect_even; fs%Q%hi_bc(3,:)=amrex_bc_reflect_even
            fs%U%lo_bc(3,:)=amrex_bc_reflect_even; fs%U%hi_bc(3,:)=amrex_bc_reflect_even
            fs%V%lo_bc(3,:)=amrex_bc_reflect_even; fs%V%hi_bc(3,:)=amrex_bc_reflect_even
            fs%W%lo_bc(3,:)=amrex_bc_reflect_odd ; fs%W%hi_bc(3,:)=amrex_bc_reflect_odd
         end if
      end block create_solver

      ! Initialize workspaces
      create_workspace: block
         use amrdata_class, only: interp_none
         call dQdt%initialize(amr,name='dQdt',ncomp=fs%nQ,ng=0,interp=interp_none); call dQdt%register()
         call Umag%initialize(amr,name='Umag',ncomp=1    ,ng=0,interp=interp_none); call Umag%register()
         call drhoUdt%initialize(amr,name='drhoUdt',ncomp=1,ng=0,nodal=[.true. ,.false.,.false.],interp=interp_none); call drhoUdt%register()
         call drhoVdt%initialize(amr,name='drhoVdt',ncomp=1,ng=0,nodal=[.false.,.true. ,.false.],interp=interp_none); call drhoVdt%register()
         call drhoWdt%initialize(amr,name='drhoWdt',ncomp=1,ng=0,nodal=[.false.,.false.,.true. ],interp=interp_none); call drhoWdt%register()
         call Mach%initialize(amr,name='Mach',ncomp=1    ,ng=0,interp=interp_none); call Mach%register()
      end block create_workspace

      ! Initialize regridding and build the initial state
      init_regridding: block
         ! Artificial bulk viscosity for the blast wave: off by default
         call param_read('Artificial viscosity coefficient',Cartif,default=0.0_WP)
         call param_read('Max viscous CFL',max_visc_cfl,default=0.5_WP)
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
         call param_read('Tag Ducros value'  ,Ducros_tag)
         if (restarted) then
            ! Restore grid hierarchy from checkpoint
            call amr%init_from_checkpoint(dirname=trim(restart_dir),time=time%t)
            ! Restore solver state
            call fs%restore_checkpoint(io=io,dirname=trim(restart_dir),time=time%t)
         else
            ! Fresh start
            call amr%init_from_scratch(time=time%t)
            ! Phasic total energies from the internal contents written by user_init and the initial kinetic energy
            call fs%get_total_energy(time=time%t)
            ! Build PLIC
            call fs%build_plic(time%t)
         end if
         ! Compute viscosities
         call get_viscosities()
         ! Add SGS models
         call fs%add_viscartif(dt=time%dt,Cartif=Cartif,max_cfl=max_visc_cfl)
         call fs%add_vreman(dt=time%dt)
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
         call viz%add_surfmesh(fs%smesh,'plic')
         ! Create visualization output event
         viz_evt=event(time=time,name='Visualization output')
         call param_read('Output period',viz_evt%tper)
         ! Write initial state
         if (viz_evt%occurs()) then
            ! Cell-centered velocity, Umag, and Mach number are only computed for visualization
            call fs%get_cell_velocity()
            call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
            call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
            call viz%write(time=time%t)
         end if
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
         call consfile%add_column(fs%nrelax_fail,'Relax failures')
         call consfile%add_column(fs%nvf_guard,'VF guard cells')
         call consfile%add_column(fs%nsliver,'Sliver cells')
         call consfile%add_column(fs%ngate,'Gate cells')
         call consfile%add_column(fs%gate_mass,'Gate mass')
         call consfile%add_column(fs%gate_energy,'Gate energy')
         call consfile%add_column(fs%gate_rhoU,'Gate U momentum')
         call consfile%add_column(fs%gate_rhoV,'Gate V momentum')
         call consfile%add_column(fs%ngate_floor,'Gate floors')
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
         ! Step-level timers of this rank
         call tfile%add_column(wt_total  ,'step_total')
         call tfile%add_column(wt_cfl    ,'cfl')
         call tfile%add_column(wt_iter   ,'iterations')
         call tfile%add_column(wt_regrid ,'regrid')
         call tfile%add_column(wt_post   ,'visc_post')
         call tfile%add_column(wt_out    ,'output')
         call tfile%add_column(wt_info   ,'info')
         call tfile%add_column(wt_metrics,'drop_metrics')
         ! Full routine times (max across ranks = wall-clock cost)
         call tfile%add_column(fs%wtmax_dQdt,'dQdt_max')
         call tfile%add_column(fs%wtmax_plic,'plic_max')
         call tfile%add_column(fs%wtmax_relax,'relax_max')
         call tfile%add_column(fs%wtmax_visc,'visc_max')
         call tfile%add_column(fs%wtmax_mid,'midpoint_max')
         ! Compute loop times (max = slowest rank, min = fastest rank)
         call tfile%add_column(fs%wtmax_prim,'prim_max')
         call tfile%add_column(fs%wtmin_prim,'prim_min')
         call tfile%add_column(fs%wtmax_fvel,'fvel_max')
         call tfile%add_column(fs%wtmin_fvel,'fvel_min')
         ! get_dQdt breakdown
         call tfile%add_column(fs%wtmax_alloc,'alloc_max')
         call tfile%add_column(fs%wtmin_alloc,'alloc_min')
         call tfile%add_column(fs%wtmax_slflux,'slflux_max')
         call tfile%add_column(fs%wtmin_slflux,'slflux_min')
         call tfile%add_column(fs%wtmax_vmom,'vmom_max')
         call tfile%add_column(fs%wtmin_vmom,'vmom_min')
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
      use mpi_f08, only: MPI_Wtime
      implicit none
      real(WP) :: tw0,tw1

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         tw0=MPI_Wtime(); tw1=tw0
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         call time%adjust_dt()
         call time%increment()
         wt_cfl=MPI_Wtime()-tw1; tw1=MPI_Wtime()

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
         wt_iter=MPI_Wtime()-tw1; tw1=MPI_Wtime()

         ! Regrid if event triggers
         wt_regrid=0.0_WP
         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0,time=time%t)
            call gridfile%write()
            wt_regrid=MPI_Wtime()-tw1; tw1=MPI_Wtime()
         end if

         ! Refresh accepted-state coefficients for the next CFL check and output.
         call get_viscosities()

         ! Add SGS models
         call fs%add_viscartif(dt=time%dt,Cartif=Cartif,max_cfl=max_visc_cfl)
         call fs%add_vreman(dt=time%dt)

         wt_post=MPI_Wtime()-tw1; tw1=MPI_Wtime()

         ! Visualization output: cell-centered velocity, Umag, and Mach number are only computed when written
         if (viz_evt%occurs()) then
            ! Regrid clears the PLIC surface mesh when it remakes the finest level: rebuild it from the current PLIC
            call fs%build_polygons()
            call fs%get_cell_velocity()
            call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
            call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
            call viz%write(time=time%t)
         end if

         ! Checkpoint save
         if (save_evt%occurs()) then
            save_checkpoint: block
               use string, only: rtoa
               call io%write(dirname='restart/drop_'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
            end block save_checkpoint
         end if

         wt_out=MPI_Wtime()-tw1; tw1=MPI_Wtime()

         ! Perform and output monitoring
         call fs%get_info(); Etot=fs%Qint(3)+fs%Qint(4)
         call mfile%write()
         call consfile%write()
         call cflfile%write()
         wt_info=MPI_Wtime()-tw1; tw1=MPI_Wtime()

         ! Compute droplet metrics on a safe temporary copy of VF
         call compute_drop_metrics()
         call drop_QOIs%write()
         wt_metrics=MPI_Wtime()-tw1
         wt_total=MPI_Wtime()-tw0

         ! Timing monitor is written last so that it covers the whole step
         call tfile%write()

      end do

      ! Force a final checkpoint on exit
      final_checkpoint: block
         use string, only: rtoa
         call io%write(dirname='restart/drop_'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
      end block final_checkpoint

   contains

      !> Compute droplet QOIs
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
      ! Finalize checkpoint
      call save_evt%finalize()
      call io%finalize()
      ! Finalize monitoring
      call mfile%finalize()
      call cflfile%finalize()
      call consfile%finalize()
      call gridfile%finalize()
      call tfile%finalize()
   end subroutine simulation_final

end module simulation
