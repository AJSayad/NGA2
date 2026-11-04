!> Simplified Sembian shock-drop case: Ms=2.4 planar shock over a water cylinder
!> Dimensional (SI) setup on a static uniform grid - no AMR, no sponge, no IB wall
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
   use relax_isentropic_class, only: relax_isentropic
   ! use safe_relax_class,    only: safe_relax
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

   ! Restart parameters
   type(amrio) :: io
   type(event) :: save_evt
   character(len=str_medium) :: restart_dir
   logical :: restarted
   real(WP) :: restart_time

   !> Simulation monitoring
   type(monitor) :: mfile,consfile,cflfile,gridfile,tfile,rescfile
   !> Relaxation-model census (relax_model%acc reduced across ranks for the rescue monitor)
   real(WP) :: diss_n=0.0_WP,diss_m=0.0_WP
   real(WP) :: quad_n=0.0_WP,swap_n=0.0_WP,flr_n=0.0_WP,flr_e=0.0_WP,stuck_n=0.0_WP

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
   type(relax_isentropic), target :: relax_model

   !> Flow parameters - all dimensional (SI)
   real(WP) :: rhoG1,pG1,u1           !< Pre-shock gas state [kg/m^3,Pa,m/s]
   real(WP) :: rhoG2,pG2,u2           !< Post-shock gas state [kg/m^3,Pa,m/s]
   real(WP) :: rhoL1,pL1              !< Initial liquid state [kg/m^3,Pa]
   real(WP) :: Ms,Us,Xs               !< Shock Mach number, shock speed [m/s] and location [m]
   real(WP) :: d_drop,x_drop          !< Drop diameter [m] and streamwise center location [m]
   real(WP) :: visc_l,visc_g          !< Constant dynamic viscosities [Pa.s]
   real(WP) :: sigma                  !< Surface tension coefficient [N/m]
   real(WP) :: T1                     !< Ambient (pre-shock) gas temperature [K]
   real(WP) :: relax_start            !< Time at which relaxation starts [s]

   !> Wall-clock timers of this rank [s], written to the timing monitor (solver-internal timers are in amrmpcomp)
   real(WP) :: wt_total=0.0_WP        !< Whole time step
   real(WP) :: wt_cfl=0.0_WP          !< CFL evaluation and time step adjustment
   real(WP) :: wt_iter=0.0_WP         !< Sub-iteration loop
   real(WP) :: wt_post=0.0_WP         !< Viscosities and SGS models after the sub-iterations
   real(WP) :: wt_out=0.0_WP          !< Visualization and checkpoint output
   real(WP) :: wt_info=0.0_WP         !< Solver info and monitor writes
   real(WP) :: wt_metrics=0.0_WP      !< Droplet metrics

   !> Artificial bulk viscosity coefficient and viscous CFL clipping
   real(WP) :: Cartif,max_visc_cfl

contains

   !> Smooth Heaviside function
   real(WP) function Hshock(x,delta)
      real(WP), intent(in) :: x,delta
      Hshock=1.0_WP/(1.0_WP+exp(-x/delta))
   end function Hshock

   !> Levelset function for the drop (cylinder in quasi-2D, sphere in 3D)
   function sphere_levelset(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=0.5_WP*d_drop-sqrt((xyz(1)-x_drop)**2+xyz(2)**2+xyz(3)**2)
      if (amr%nz.eq.1) G=0.5_WP*d_drop-sqrt((xyz(1)-x_drop)**2+xyz(2)**2) ! Enable quasi-2D runs
   end function sphere_levelset

   !> Compute viscosity: constant per phase, VF-weighted blend in mixed cells
   subroutine get_viscosities()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pVisc,pBeta
      ! real(WP), dimension(:,:,:,:), contiguous, pointer :: pDiffL,pDiffG (not available in the staggered solver)
      real(WP), parameter :: myeps=1.0e-15_WP
      ! Loop over levels
      do lvl=0,amr%clvl()
         ! Loop over domain
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pVF=>fs%VF%mf(lvl)%dataptr(mfi)
            pVisc=>fs%visc%mf(lvl)%dataptr(mfi)
            pBeta=>fs%beta%mf(lvl)%dataptr(mfi)
            ! pDiffL=>fs%diffL%mf(lvl)%dataptr(mfi)
            ! pDiffG=>fs%diffG%mf(lvl)%dataptr(mfi)
            ! Get tilebox grown by the viscosity ghost cells
            bx=mfi%growntilebox(min(fs%visc%ng,fs%beta%ng))
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                     ! Mixture viscosity from constant phasic values
                     !pVisc(i,j,k,1)=pVF(i,j,k,1)*visc_l+(1.0_WP-pVF(i,j,k,1))*visc_g ! Arithmetic averaging
                     pVisc(i,j,k,1)=1.0_WP/(pVF(i,j,k,1)/max(visc_l,myeps)+(1.0_WP-pVF(i,j,k,1))/max(visc_g,myeps)) ! Harmonic averaging
                     ! Zero bulk viscosity
                     pBeta(i,j,k,1)=0.0_WP
                     ! Run with no heat transfer
                     ! pDiffG(i,j,k,1)=0.0_WP
                     ! pDiffL(i,j,k,1)=0.0_WP
                  end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine get_viscosities

   !> Refresh only LAD after prepare_midpoint; molecular and Vreman shear viscosity
   !> remain at the accepted start-of-step values to isolate bulk-viscosity timing.
   subroutine update_midpoint_lad()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      implicit none
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pRHOL,pRHOG,pPL,pPG,pC
      real(WP) :: cL,cG,mL,mG

      ! prepare_midpoint updates the phasic densities and pressures, but not C.
      ! Use those midpoint primitives for the LAD sound-speed floor. Use the latest
      ! iterate's VF, matching add_viscartif's density conversion; do not modify VF
      ! or the start-of-step geometry used by the subsequent geometric flux sweep.
      do lvl=0,amr%clvl()
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pVF=>fs%VF%mf(lvl)%dataptr(mfi)
            pRHOL=>fs%RHOL%mf(lvl)%dataptr(mfi)
            pRHOG=>fs%RHOG%mf(lvl)%dataptr(mfi)
            pPL=>fs%PL%mf(lvl)%dataptr(mfi)
            pPG=>fs%PG%mf(lvl)%dataptr(mfi)
            pC=>fs%C%mf(lvl)%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               cL=0.0_WP; cG=0.0_WP
               if (pRHOL(i,j,k,1).gt.0.0_WP) cL=fs%liq%get_c_from_p_rho(p=pPL(i,j,k,1),rho=pRHOL(i,j,k,1),y=[1.0_WP])
               if (pRHOG(i,j,k,1).gt.0.0_WP) cG=fs%gas%get_c_from_p_rho(p=pPG(i,j,k,1),rho=pRHOG(i,j,k,1),y=[1.0_WP])
               mL=pVF(i,j,k,1)*pRHOL(i,j,k,1)
               mG=(1.0_WP-pVF(i,j,k,1))*pRHOG(i,j,k,1)
               pC(i,j,k,1)=sqrt((mL*cL**2+mG*cG**2)/max(mL+mG,fs%rho_floor))
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do

      ! add_viscartif accumulates: clear beta, including ghosts, on every iterate.
      ! Cvisc=0 leaves the existing molecular+Vreman shear viscosity unchanged.
      call fs%beta%setval(val=0.0_WP)
      call fs%add_viscartif(dt=time%dt,Cartif=Cartif,Cvisc=0.0_WP,max_cfl=max_visc_cfl)
   end subroutine update_midpoint_lad

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

   !> User init callback - set Q and VF/barycenters for a drop at rest with a shock
   subroutine shockdrop_init(solver,lvl,time,ba,dm)
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
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF,pCL,pCG
      real(WP), dimension(3) :: BL,BG
      real(WP) :: dx,dy,dz,myVF,IEL,x_cc,rhoG,pG,H
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU
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
                  ! Compute local gas state from shock profile (post-shock state sits at low x)
                  x_cc=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
                  H=Hshock(x=Xs-x_cc,delta=0.5_WP*dx)
                  rhoG=rhoG1+(rhoG2-rhoG1)*H
                  pG  =pG1  +(pG2  -pG1  )*H
                  ! Set conserved variables: Q=(VF*rhoL, (1-VF)*rhoG, VF*rhoL*IL, (1-VF)*rhoG*IG); momentum is staggered
                  pQ(i,j,k,1)=(       myVF)*rhoL1
                  pQ(i,j,k,2)=(1.0_WP-myVF)*rhoG
                  pQ(i,j,k,3)=pQ(i,j,k,1)*IEL
                  pQ(i,j,k,4)=pQ(i,j,k,2)*gas%get_e_from_p_rho(p=pG,rho=rhoG,y=[1.0_WP])
               end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
      ! Set x-face velocity from the same shock profile, ghosts included. V and W are zero (set by on_init).
      ! Face momentum is built from this velocity and Q in post_regrid.
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
         u1=0.0_WP
         ! Post-shock gas state
         call param_read('Post-shock density' ,rhoG2)
         call param_read('Post-shock pressure',pG2)
         call param_read('Post-shock velocity',u2)
         ! Shock location
         call param_read('Shock location',Xs)
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
         ! Diagnostics: shock speed and Mach number recovered from the imposed jump
         Us=u2/(1.0_WP-rhoG1/rhoG2)
         Ms=Us/sqrt(GammaG*pG1/rhoG1)
         ! Constant dynamic viscosities
         call param_read('Liquid viscosity',visc_l)
         call param_read('Gas viscosity'   ,visc_g)
         ! Surface tension
         call param_read('Surface tension coefficient',sigma,default=0.0_WP)
         ! Log
         write(message,'("[Shock] Ms=",es12.5," Us=",es12.5," m/s  Xs=",es12.5," m")') Ms,Us,Xs; call log(message)
         write(message,'("[Pre-shock]  rhoG1=",es12.5," pG1=",es12.5," u1=",es12.5)') rhoG1,pG1,u1; call log(message)
         write(message,'("[Post-shock] rhoG2=",es12.5," pG2=",es12.5," u2=",es12.5)') rhoG2,pG2,u2; call log(message)
         write(message,'("[Liquid] rhoL1=",es12.5," pL1=",es12.5," CvL=",es12.5)') rhoL1,pL1,CvL; call log(message)
         write(message,'("[Temp]   TL=",es12.5," TG=",es12.5)') water%get_T_from_p_rho(p=pL1,rho=rhoL1,y=[1.0_WP]),T1; call log(message)
         call water%print(); call gas%print()
         write(message,'("[Visc]   muL=",es12.5," muG=",es12.5," Pa.s")') visc_l,visc_g; call log(message)
         write(message,'("[Surface tension] sigma=",es12.5," N/m")') sigma; call log(message)
         if (sigma.gt.0.0_WP) then
            write(message,'("[Scales] We=",es12.5)') rhoG2*u2**2*d_drop/sigma; call log(message)
         end if
      end block init_eos_and_flow

      ! Initialize the grid - single static level, no refinement
      create_amrgrid: block
         use messager, only: log
         use string,   only: str_long
         character(len=str_long) :: message
         real(WP) :: Lx,Ly
         ! Set name
         amr%name='amr_sembian_simplified'
         ! Static uniform grid: level 0 only
         amr%maxlvl=0
         ! Read in grid size
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         ! Set domain: x in [0,Lx], y in [-Ly/2,+Ly/2]
         call param_read('Lx',Lx)
         call param_read('Ly',Ly)
         amr%xlo=0.0_WP;        amr%xhi=+Lx
         amr%ylo=-0.5_WP*Ly;    amr%yhi=+0.5_WP*Ly
         amr%zlo=-0.5_WP*Ly;    amr%zhi=+0.5_WP*Ly
         ! Set periodicity: x and y are physical (constant extrapolation), z is the quasi-2D direction
         amr%xper=.false.; amr%yper=.false.; amr%zper=.true.
         ! Enable quasi-2D: one cell deep, dz=dy
         if (amr%nz.eq.1) then
            amr%zlo=-0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny,WP)
            amr%zhi=+0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny,WP)
         end if
         ! Initialize
         call amr%initialize()
         ! Log drop resolution (needs the initialized grid)
         write(message,'("[Drop]   d0=",es12.5," m  x0=",es12.5," m  CPD=",es12.5)') d_drop,x_drop,d_drop/amr%dy(0); call log(message)
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
         use amrex_amr_module, only: amrex_bc_ext_dir,amrex_bc_foextrap
         use amrmpcomp_class, only: BC_GAS
         ! Surface tension coefficient
         fs%sigma=sigma
         ! Band momentum flux: CTU below the staggered phase fraction CTU alpha, MC elsewhere. Set before initialize, which validates them
         call param_read('CTU alpha',fs%alpha_ctu,default=0.5_WP)
         ! Relaxation is applied from t=0 unless a later start time is requested
         call param_read('Relaxation start time',relax_start,default=0.0_WP)
         ! Assign materials and create flow solver (initialize forces face-linear velocity interpolation)
         fs%liq=>water; fs%gas=>gas; call fs%initialize(amr=amr,name='drop')
         ! Provide pressure relaxation model (gas-isentropic)
         call relax_model%initialize(gas=gas,liq=water); fs%relax=>relax_model
         ! Set initial conditions
         fs%user_init=>shockdrop_init
         ! Set BCs. Staggered solver conventions:
         !  - ext_dir on a velocity component calls shock_dirichlet (normal: boundary face + ghosts; tangential: ghosts)
         !  - ext_dir on Q calls shock_dirichlet to fill the ghosts; boundary face momentum is then derived by apply_momentum_bc
         !  - foextrap on a normal velocity leaves the boundary face as an evolved unknown; ghosts copy it
         ! x-low: post-shock inflow. x-high and y: constant extrapolation (outflow). z stays periodic (quasi-2D).
         ! VOF: gas ghosts at the inflow; at the outflows the ghosts are also gas so that liquid leaves the domain freely
         if (.not.amr%xper) then
            fs%lo_bc(1)=BC_GAS; fs%hi_bc(1)=BC_GAS
            fs%Q%lo_bc(1,:)=amrex_bc_ext_dir; fs%Q%hi_bc(1,:)=amrex_bc_foextrap
            fs%U%lo_bc(1,:)=amrex_bc_ext_dir; fs%U%hi_bc(1,:)=amrex_bc_foextrap
            fs%V%lo_bc(1,:)=amrex_bc_ext_dir; fs%V%hi_bc(1,:)=amrex_bc_foextrap
            fs%W%lo_bc(1,:)=amrex_bc_ext_dir; fs%W%hi_bc(1,:)=amrex_bc_foextrap
            fs%user_bc=>shock_dirichlet
         end if
         if (.not.amr%yper) then
            fs%lo_bc(2)=BC_GAS; fs%hi_bc(2)=BC_GAS
            fs%Q%lo_bc(2,:)=amrex_bc_foextrap; fs%Q%hi_bc(2,:)=amrex_bc_foextrap
            fs%U%lo_bc(2,:)=amrex_bc_foextrap; fs%U%hi_bc(2,:)=amrex_bc_foextrap
            fs%V%lo_bc(2,:)=amrex_bc_foextrap; fs%V%hi_bc(2,:)=amrex_bc_foextrap
            fs%W%lo_bc(2,:)=amrex_bc_foextrap; fs%W%hi_bc(2,:)=amrex_bc_foextrap
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

      ! Build the static grid and the initial state
      init_grid: block
         ! Artificial bulk viscosity: off by default, there are no shocks
         call param_read('Artificial viscosity coefficient',Cartif,default=0.0_WP)
         call param_read('Max viscous CFL',max_visc_cfl,default=0.5_WP)
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
      end block init_grid

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
         ! Create rescue-census monitor (cumulative counters/amounts per mechanism)
         rescfile=monitor(amRoot=amr%amRoot,name='rescue')
         call rescfile%add_column(time%n,'Timestep')
         call rescfile%add_column(time%t,'Time')
         ! call rescfile%add_column(fs%resc_nl,'LiqResc n')
         ! call rescfile%add_column(fs%resc_ml,'LiqResc dm')
         ! call rescfile%add_column(fs%resc_el,'LiqResc dE')
         ! call rescfile%add_column(fs%resc_ng,'GasResc n')
         ! call rescfile%add_column(fs%resc_mg,'GasResc dm')
         ! call rescfile%add_column(fs%resc_eg,'GasResc dE')
         ! call rescfile%add_column(diss_n,'Diss n')
         ! call rescfile%add_column(diss_m,'Diss dm')
         ! call rescfile%add_column(quad_n,'Quad n')
         ! call rescfile%add_column(swap_n,'Swap n')
         ! call rescfile%add_column(flr_n,'Floor n')
         ! call rescfile%add_column(flr_e,'Floor dE')
         ! call rescfile%add_column(stuck_n,'Stuck n')
         ! call rescfile%add_column(fs%pool_n,'Pool n')
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
            call update_midpoint_lad()

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
            if (time%t.ge.relax_start) call fs%apply_relax(dt=time%dt,time=time%t)

            ! Get primitive variables of the new iterate
            call fs%get_primitive(Q=fs%Q,rhoU=fs%rhoU,rhoV=fs%rhoV,rhoW=fs%rhoW,time=time%t)

            ! Increment sub-iteration counter
            time%it=time%it+1

         end do
         ! =================================================================================
         wt_iter=MPI_Wtime()-tw1; tw1=MPI_Wtime()

         ! Refresh accepted-state coefficients for the next CFL check and output.
         call get_viscosities()

         ! Add SGS models
         call fs%add_viscartif(dt=time%dt,Cartif=Cartif,max_cfl=max_visc_cfl)
         call fs%add_vreman(dt=time%dt)

         wt_post=MPI_Wtime()-tw1; tw1=MPI_Wtime()

         ! Visualization output: cell-centered velocity, Umag, and Mach number are only computed when written
         if (viz_evt%occurs()) then
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
         ! relax_census: block
         !    use mpi_f08,  only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM
         !    use parallel, only: MPI_REAL_WP
         !    real(WP), dimension(7) :: tmp
         !    integer :: ierr
         !    tmp=relax_model%acc
         !    call MPI_ALLREDUCE(MPI_IN_PLACE,tmp,7,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
         !    diss_n=tmp(1); diss_m=tmp(2); quad_n=tmp(3); swap_n=tmp(4)
         !    flr_n=tmp(5); flr_e=tmp(6); stuck_n=tmp(7)
         ! end block relax_census
         call mfile%write()
         call consfile%write()
         call cflfile%write()
         call rescfile%write()
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
      call rescfile%finalize()
   end subroutine simulation_final

end module simulation
