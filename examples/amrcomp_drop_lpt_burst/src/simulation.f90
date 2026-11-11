!> AMR compressible drop case with Lagrangian spray transfer
module simulation
   use precision,           only: WP,I8
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
   use relax_ig_nasg_class, only: relax_ig_nasg, Prelax !, PThybrid
   ! use safe_relax_class,    only: safe_relax
   use amrccl_class,        only: amrccl
   use amrlpt_class,        only: amrlpt,part,particle_finite,PART_MOVES,PART_COLLIDES,PART_EXCHANGES,AMRLPT_OPEN
   use amrlptviz_class,     only: amrlptviz
   use timer_class,         only: timer
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

   !> Characteristic time tau=t/tc_sim, tc_sim=sqrt(rhoL1/rhoG2)*d_drop/u2, with tau=0 when the shock reaches the drop leading edge (t=t_hit)
   real(WP) :: tc_sim=1.0_WP,t_hit=0.0_WP

   !> Burst output (camera exposure replication, see burst_implementation.md): PLIC+spray every step inside [burst_t0,burst_t1]
   integer :: nburst=0
   real(WP), dimension(:), allocatable :: burst_tau,burst_t0,burst_t1   !< Window centre (tau) and bounds (simulation time)
   logical,  dimension(:), allocatable :: burst_open                    !< Writers of this burst have been initialized
   type(amrviz),    dimension(:), allocatable :: burst_viz             !< PLIC writer per burst
   type(amrlptviz), dimension(:), allocatable :: burst_lpviz           !< Spray writer per burst

   ! Regrid parameters
   type(event) :: regrid_evt

   ! Restart parameters
   type(amrio) :: io
   type(event) :: save_evt
   character(len=str_medium) :: restart_dir
   logical :: restarted
   integer :: restart_step=0
   real(WP) :: restart_time

   !> Simulation monitoring
   type(monitor) :: mfile,consfile,cflfile,gridfile,tfile,rescfile

   !> Mixture total energy, the sum of the phasic total energies
   real(WP) :: Etot=0.0_WP

   !> Droplet QOI monitoring, evaluated on metrics_evt rather than every step
   type(monitor) :: drop_QOIs
   type(event)   :: metrics_evt
   !> Isocontour-based quantities, over the whole domain
   real(WP) :: MLE_VF01=0.0_WP,MLE_VF05=0.0_WP                        !< Mist leading edge at VF>=0.1 and VF>=0.5
   real(WP) :: drop_massVF01=0.0_WP                                   !< Droplet mass at VF>=0.1
   real(WP) :: MOI_xx=0.0_WP,MOI_yy=0.0_WP,MOI_zz=0.0_WP              !< Inertia tensor about the VF-weighted COM (VF>=0.1)
   real(WP) :: MOI_xy=0.0_WP,MOI_xz=0.0_WP,MOI_yz=0.0_WP
   !> Primary-drop quantities, over every cell in the selected core
   real(WP) :: MOI_xx_core=0.0_WP,MOI_yy_core=0.0_WP,MOI_zz_core=0.0_WP
   real(WP) :: MOI_xy_core=0.0_WP,MOI_xz_core=0.0_WP,MOI_yz_core=0.0_WP
   real(WP) :: core_mass=0.0_WP,core_vol=0.0_WP                       !< Primary-drop mass and liquid volume
   real(WP) :: MLE_core=0.0_WP                                        !< Windward edge of the primary drop (any liquid)
   real(WP) :: MLE_core001=0.0_WP,MLE_core01=0.0_WP,MLE_core05=0.0_WP !< ... restricted to VF>=0.01/0.1/0.5
   !> Monitors for primary droplet
   real(WP) :: core_valid=0.0_WP,core_open=0.0_WP,core_wrap=0.0_WP

   !> Sub-grid drop transfer: connected-component labelling hands unresolved liquid structures to a Lagrangian spray
   type(amrccl), target :: ccl
   type(amrdata)   :: spray_energy,spray_force(3),spray_cells !< Frozen cell energy / face momentum rates
   type(amrdata)   :: particle_stress    !< Resolved effective stress divergence, force/volume
   type(amrdata)   :: particle_thermal  !< Gas temperature and molecular conductivity
   type(amrdata)   :: particle_gas      !< Molecular gas viscosity and gas sound speed
   type(amrdata)   :: particle_gas_fraction !< Valid gas fraction for normalized interpolation
   type(amrlpt)    :: lp
   type(amrlptviz) :: lpviz
   type(event)     :: lpviz_evt,xfer_evt
   type(monitor)   :: sprayfile
   logical  :: do_transfer=.false.        !< Master switch, off by default
   real(WP) :: xfer_dmax_ratio=1.0_WP     !< Maximum equivalent diameter / coupling-cap spacing
   real(WP) :: xfer_diam_threshold=1.0e-1_WP !< Minimum particle diameter / finest-grid spacing
   real(WP) :: xfer_margin=-1.0_WP        !< Core-box padding for wake tagging only; <0 means 0.5*d_drop
   real(WP) :: xfer_dmax=0.0_WP           !< Equal-volume diameter of the largest transferable drop
   real(WP) :: lp_dx=0.0_WP               !< Cell size on the finest allowed particle level
   real(WP) :: lp_dres=0.0_WP             !< xfer_dmax/lp_dx - drop diameter in coupling-mesh cells.
   real(WP) :: lp_dres_limit=1.0_WP       !< Largest permitted particle diameter / particle-cell size
   real(WP), dimension(3) :: corelo=0.0_WP,corehi=0.0_WP  !< Core bounding box, reported to the monitor
   integer  :: lp_lvl_off=3               !< Levels below amr%maxlvl that particles are capped to
   logical  :: lp_collide=.false.         !< Soft-sphere collisions between transferred drops
   !> Number of per-structure moments used by label_and_moments
   integer, parameter :: NMOM=24 !< Last three: donor gas volume, mass, internal energy
   integer :: core_id=0                  !< Largest liquid-volume label
   real(WP), dimension(NMOM) :: core_mom=0.0_WP
   !> Windward face x of the primary drop: any liquid, then VF>=0.01/0.1/0.5
   real(WP), dimension(4)    :: core_mle=0.0_WP
   logical  :: core_mom_valid=.false.     !< core_mom/core_mle hold a labelling from THIS step
   !> Transfer bookkeeping - per-step counts and running totals for the spray monitor
   integer     :: ccl_nstruct=0,xfer_n_step=0,flot_n_step=0
   integer(I8) :: lp_np=0
   real(WP)    :: lp_tprev=0.0_WP   !< Flow time at the last particle advance
   real(WP)    :: dt_src=1.0_WP     !< Interval lp%src was accumulated over (measured, not assumed)
   !> Candidate census of the last labelling (see label_and_moments)
   integer  :: xfer_ncand=0
   real(WP), dimension(2) :: xfer_rej=0.0_WP
   real(WP), dimension(4) :: xfer_fin=0.0_WP !< Size, coarse, no-donor, and invalid-birth rejections
   real(WP) :: core_fallbacks=0.0_WP      !< Labellings whose core could not be certified from the top-K sets
   !> Per-step coupling tracker
   real(WP)    :: src_momx=0.0_WP,src_ener=0.0_WP,src_emin=0.0_WP
   !> Largest velocity change the coupling imposes on the gas in one step, max|src(1:3)|/rho.
   real(WP)    :: src_dumax=0.0_WP
   !> Per-phase wall-clock profiling.
   type(timer)   :: tmr_step_all,tmr_prep,tmr_lpt,tmr_lptinfo,tmr_divstress,tmr_loop,tmr_rk
   type(timer)   :: tmr_xfer,tmr_ccl,tmr_msweep,tmr_mcomm,tmr_xpart,tmr_surg,tmr_resync
   type(timer)   :: tmr_diag,tmr_regrid,tmr_props,tmr_sgs,tmr_viz,tmr_ckpt,tmr_mon
   type(monitor) :: perffile
   real(WP)    :: t_step=0.0_WP,t_prep=0.0_WP,t_lpt=0.0_WP,t_lptinfo=0.0_WP,t_divstress=0.0_WP
   real(WP)    :: t_lptstep=0.0_WP,t_lptsrc=0.0_WP,t_lptvf=0.0_WP
   real(WP)    :: t_loop=0.0_WP,t_rk=0.0_WP
   real(WP)    :: t_xfer=0.0_WP,t_ccl=0.0_WP,t_msweep=0.0_WP,t_mcomm=0.0_WP,t_xpart=0.0_WP,t_surg=0.0_WP,t_resync=0.0_WP
   real(WP)    :: t_diag=0.0_WP,t_regrid=0.0_WP,t_props=0.0_WP,t_sgs=0.0_WP,t_viz=0.0_WP,t_ckpt=0.0_WP,t_mon=0.0_WP
   real(WP)    :: t_other=0.0_WP
   real(WP)    :: xfer_vol=0.0_WP,xfer_mass=0.0_WP,xfer_momx=0.0_WP,xfer_momy=0.0_WP,xfer_momz=0.0_WP,xfer_energy=0.0_WP
   !> Cumulative numerical liquid discarded by the minimum particle-diameter threshold
   real(WP)    :: flot_n_total=0.0_WP,flot_vol=0.0_WP,flot_mass=0.0_WP
   real(WP) :: flot_mom(3)=0.0_WP,flot_energy=0.0_WP
   !> Gas mass the surgery adds when filling vacated volume
   real(WP)    :: xfer_gasadd=0.0_WP
   !> Numerical refill residual: mesh change plus removed liquid momentum/energy (including flotsam)
   real(WP), dimension(4) :: xfer_refill=0.0_WP
   real(WP) :: xfer_refused=0.0_WP,xfer_badthermal=0.0_WP
   !> Cumulative change made by regrids: liquid mass, gas mass, momentum, total energy
   real(WP), dimension(6) :: regrid_book=0.0_WP
   !> Per-FAB pointers rebound by ccl_prepare and read by ccl_make / ccl_same
   real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF=>null(),pPLIC=>null()
   real(WP) :: cdx,cdy,cdz                !< Mesh size at the labeled level

   !> Materials
   type(nasg),      target :: water
   type(ideal_gas), target :: gas

   !> Relaxation model
   type(relax_ig_nasg), target :: relax_model

   !> Flow parameters
   real(WP) :: rhoG1,pG1,u1             !< Pre-shock gas state
   real(WP) :: rhoG2,pG2,u2             !< Post-shock gas state
   real(WP) :: rhoL1,pL1                !< Initial liquid state
   real(WP) :: M2,Xs                    !< Post-shock Mach and shock location
   real(WP) :: Ms                       !< Shock Mach number
   real(WP) :: density_ratio            !< rhoL1/rhoG2
   real(WP), parameter :: d_drop=1.0_WP !> Initial drop diameter
   real(WP) :: ML                       !< Liquid Mach number
   real(WP) :: Reynolds,visc_ratio      !< Viscosity
   real(WP) :: Prandtl ,diff_ratio      !< Heat diffusivity
   real(WP) :: Weber                    !< Weber number

   !> Sutherland viscosity parameters: mu_g = (1+Suth_T)*T^Suth_n / (Re*(T+Suth_T))
   real(WP) :: Suth_n=1.5_WP            !< Sutherland exponent (1.0 for constant)
   real(WP) :: Suth_T=0.4042_WP         !< Sutherland temperature (0.0 for constant)

   !> Sponge parameters
   real(WP) :: R_spg=3.0_WP             !< Radius at which the sponge starts
   real(WP) :: L_spg=2.0_WP             !< Length of the sponge ramp
   real(WP) :: sig_spg=1.0_WP           !< Sponge relaxation rate [1/t] - set to 0 to disable
   real(WP) :: Ws                       !< Lab-frame shock speed - locates the sponge target state

   !> Spherical harmonics perturbation parameters
   integer :: nsh_modes=0
   integer , dimension(:), allocatable :: l_modes,m_modes
   real(WP), dimension(:), allocatable :: amp_modes,phase_modes

   !> Tagging parameter
   real(WP) :: Re_tag=huge(1.0_WP)
   real(WP) :: Rho_tag=huge(1.0_WP)
   real(WP) :: P_tag=huge(1.0_WP)
   real(WP) :: Ducros_tag=huge(1.0_WP)
   !> Far-field coarsening, beyond the sponge ramp (r_cyl > R_spg+L_spg)
   integer, parameter :: spg_bandcells=16 !< Width of each band, in cells of that band's own level
   integer :: spg_minlvl=2                !< Coarsest level the graded cap will descend to
   !> Sponge edge
   real(WP), dimension(:), allocatable :: spg_redge
   !> Global refinement FLOOR. 
   integer :: min_lvl=2
   !> Wake coarsening
   integer  :: wake_lvl=huge(1)           !< SGS-Re threshold starts ramping above this level
   real(WP) :: wake_ramp=1.0_WP           !< Threshold multiplier per level above wake_lvl
   integer  :: wake_maxlvl_off=0           !< Hard wake cap: maxlvl-offset (0 preserves old behavior)
   !> Downstream station where the distance ramp starts
   real(WP) :: wake_x0=huge(1.0_WP)
   real(WP) :: wake_len=1.0_WP            !< Distance over which the threshold doubles
   !> Also ramp the density-Laplacian criterion
   logical  :: wake_rho=.false.

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
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pTG,pVF,pVisc,pBeta
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
            ! Get tilebox (TG has one ghost layer, so visc and beta ghosts are filled below)
            bx=mfi%tilebox()
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
                     ! Resolved-fluid conduction remains off: diffL/diffG stay zero; particle heat is separate
                  end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
      ! Fill the viscosity ghost cells read by the momentum fluxes
      call fs%visc%fill(time=time%t); call fs%beta%fill(time=time%t)
   end subroutine get_viscosities

   !> Rebuild molecular gas properties from the beginning-of-step state, over the particle ghost width
   subroutine get_particle_gas_properties()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      use messager, only: die
      implicit none
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pRHOG,pPG,pGas,pFraction,pTherm
      real(WP) :: Tgas,Tmu
      integer :: lvl,i,j,k
      do lvl=0,min(lp%maxlvl,amr%clvl())
         call amr%mfiter_build(lvl,mfi,tiling=.false.)
         do while (mfi%next())
            pVF=>fs%VF%mf(lvl)%dataptr(mfi)
            pRHOG=>fs%RHOG%mf(lvl)%dataptr(mfi)
            pPG=>fs%PG%mf(lvl)%dataptr(mfi)
            pGas=>particle_gas%mf(lvl)%dataptr(mfi)
            pTherm=>particle_thermal%mf(lvl)%dataptr(mfi)
            pFraction=>particle_gas_fraction%mf(lvl)%dataptr(mfi)
            bx=mfi%growntilebox(lp%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pGas(i,j,k,:)=0.0_WP; pFraction(i,j,k,1)=0.0_WP; pTherm(i,j,k,:)=0.0_WP
               if (.not.particle_finite(pVF(i,j,k,1))) call die('[particle gas properties] non-finite liquid fraction')
               if (pVF(i,j,k,1).ge.1.0_WP) cycle
               if (.not.particle_finite(pRHOG(i,j,k,1)).or..not.particle_finite(pPG(i,j,k,1))) &
               & call die('[particle gas properties] non-finite gas density or pressure')
               ! The solver treats the gas as absent where its mass or internal energy is not positive
               if (pRHOG(i,j,k,1).le.0.0_WP.or.pPG(i,j,k,1).le.0.0_WP) cycle
               Tgas=gas%get_T_from_p_rho(p=pPG(i,j,k,1),rho=pRHOG(i,j,k,1),y=[1.0_WP])
               if (.not.particle_finite(Tgas).or.Tgas.le.0.0_WP) &
               & call die('[particle gas properties] invalid gas temperature')
               Tmu=min(Tgas,10.0_WP) ! Same molecular-viscosity cap as get_viscosities
               pGas(i,j,k,1)=(1.0_WP+Suth_T)*Tmu**Suth_n/(Reynolds*(Tmu+Suth_T))
               pGas(i,j,k,2)=sqrt(gas%gamma*(gas%gamma-1.0_WP)*gas%cv*Tgas)
               pTherm(i,j,k,1)=Tgas
               ! Nondimensional k: mu*cp/Pr, with cp in the same energy/temperature units as the EOS.
               pTherm(i,j,k,2)=pGas(i,j,k,1)*gas%cp/Prandtl
               pFraction(i,j,k,1)=max(0.0_WP,1.0_WP-pVF(i,j,k,1))
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine get_particle_gas_properties

   !> Resolved effective stress divergence for particle coupling, in force/volume, on levels up to maxlvl
   subroutine get_div_stress(stress_div,time,maxlvl)
      use amrex_amr_module, only: amrex_multifab,amrex_multifab_destroy,amrex_mfiter,amrex_box
      use amrex_interface,  only: amrmfab_average_down_face
      use messager,         only: die
      implicit none
      type(amrdata), intent(inout) :: stress_div
      real(WP), intent(in) :: time
      integer, intent(in) :: maxlvl
      type(amrex_multifab), dimension(0:amr%maxlvl) :: Fx,Fy,Fz
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx,gbx,fbx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pFx,pFy,pFz,pU,pV,pW,pP,pVisc,pBeta,pRHS
      real(WP), dimension(:,:,:,:), allocatable :: uc
      real(WP), dimension(1:3,1:3) :: dUdx
      real(WP) :: dxi,dyi,dzi,div,visc_f,beta_f
      integer :: top,lvl,i,j,k,n

      if (any(stress_div%nodal).or.stress_div%ncomp.ne.3.or.stress_div%ng.lt.1) &
      & call die('[get_div_stress] output needs 3 cell-centered components and ghosts')
      top=min(amr%clvl(),maxlvl)
      if (top.lt.0) call die('[get_div_stress] negative level ceiling')
      call stress_div%setval(0.0_WP)

      ! Allocate face-centered stress fluxes (x,y,z momentum components)
      do lvl=0,top
         call amr%mfab_build(lvl=lvl,mfab=Fx(lvl),ncomp=3,nover=0,atface=[.true. ,.false.,.false.]); call Fx(lvl)%setval(0.0_WP)
         call amr%mfab_build(lvl=lvl,mfab=Fy(lvl),ncomp=3,nover=0,atface=[.false.,.true. ,.false.]); call Fy(lvl)%setval(0.0_WP)
         call amr%mfab_build(lvl=lvl,mfab=Fz(lvl),ncomp=3,nover=0,atface=[.false.,.false.,.true. ]); call Fz(lvl)%setval(0.0_WP)
      end do

      ! Pressure and viscous stresses on faces
      do lvl=0,top
         dxi=1.0_WP/amr%dx(lvl); if (amr%nx.eq.1) dxi=0.0_WP
         dyi=1.0_WP/amr%dy(lvl); if (amr%ny.eq.1) dyi=0.0_WP
         dzi=1.0_WP/amr%dz(lvl); if (amr%nz.eq.1) dzi=0.0_WP
         call amr%mfiter_build(lvl,mfi,tiling=.false.)
         do while (mfi%next())
            pFx=>Fx(lvl)%dataptr(mfi); pFy=>Fy(lvl)%dataptr(mfi); pFz=>Fz(lvl)%dataptr(mfi)
            pU=>fs%U%mf(lvl)%dataptr(mfi); pV=>fs%V%mf(lvl)%dataptr(mfi); pW=>fs%W%mf(lvl)%dataptr(mfi)
            pP=>fs%P%mf(lvl)%dataptr(mfi)
            pVisc=>fs%visc%mf(lvl)%dataptr(mfi); pBeta=>fs%beta%mf(lvl)%dataptr(mfi)
            ! Cell velocity over one ghost layer, from the face velocities (two ghost layers)
            gbx=mfi%growntilebox(1)
            allocate(uc(gbx%lo(1):gbx%hi(1),gbx%lo(2):gbx%hi(2),gbx%lo(3):gbx%hi(3),3))
            do k=gbx%lo(3),gbx%hi(3); do j=gbx%lo(2),gbx%hi(2); do i=gbx%lo(1),gbx%hi(1)
               uc(i,j,k,1)=0.5_WP*sum(pU(i:i+1,j,k,1))
               uc(i,j,k,2)=0.5_WP*sum(pV(i,j:j+1,k,1))
               uc(i,j,k,3)=0.5_WP*sum(pW(i,j,k:k+1,1))
            end do; end do; end do
            ! X-fluxes
            if (amr%nx.gt.1) then
               fbx=mfi%nodaltilebox(1)
               do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
                  ! Add pressure stress
                  pFx(i,j,k,1)=pFx(i,j,k,1)-0.5_WP*sum(pP(i-1:i,j,k,1))
                  ! Velocity gradients at x-face
                  dUdx(1,1)=dxi*(uc(i,j,k,1)-uc(i-1,j,k,1))
                  dUdx(2,1)=0.25_WP*dyi*(uc(i-1,j+1,k,1)-uc(i-1,j-1,k,1)+uc(i,j+1,k,1)-uc(i,j-1,k,1))
                  dUdx(3,1)=0.25_WP*dzi*(uc(i-1,j,k+1,1)-uc(i-1,j,k-1,1)+uc(i,j,k+1,1)-uc(i,j,k-1,1))
                  dUdx(1,2)=dxi*(uc(i,j,k,2)-uc(i-1,j,k,2))
                  dUdx(2,2)=0.25_WP*dyi*(uc(i-1,j+1,k,2)-uc(i-1,j-1,k,2)+uc(i,j+1,k,2)-uc(i,j-1,k,2))
                  dUdx(3,2)=0.25_WP*dzi*(uc(i-1,j,k+1,2)-uc(i-1,j,k-1,2)+uc(i,j,k+1,2)-uc(i,j,k-1,2))
                  dUdx(1,3)=dxi*(uc(i,j,k,3)-uc(i-1,j,k,3))
                  dUdx(2,3)=0.25_WP*dyi*(uc(i-1,j+1,k,3)-uc(i-1,j-1,k,3)+uc(i,j+1,k,3)-uc(i,j-1,k,3))
                  dUdx(3,3)=0.25_WP*dzi*(uc(i-1,j,k+1,3)-uc(i-1,j,k-1,3)+uc(i,j,k+1,3)-uc(i,j,k-1,3))
                  div=dUdx(1,1)+dUdx(2,2)+dUdx(3,3)
                  ! Viscosities at x-face
                  visc_f=2.0_WP*product(pVisc(i-1:i,j,k,1))/(sum(pVisc(i-1:i,j,k,1))+tiny(1.0_WP))
                  beta_f=2.0_WP*product(pBeta(i-1:i,j,k,1))/(sum(pBeta(i-1:i,j,k,1))+tiny(1.0_WP))
                  ! Viscous stress at x-face
                  pFx(i,j,k,1)=pFx(i,j,k,1)+visc_f*(dUdx(1,1)+dUdx(1,1))+(beta_f-2.0_WP/3.0_WP*visc_f)*div
                  pFx(i,j,k,2)=pFx(i,j,k,2)+visc_f*(dUdx(2,1)+dUdx(1,2))
                  pFx(i,j,k,3)=pFx(i,j,k,3)+visc_f*(dUdx(3,1)+dUdx(1,3))
               end do; end do; end do
            end if
            ! Y-fluxes
            if (amr%ny.gt.1) then
               fbx=mfi%nodaltilebox(2)
               do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
                  ! Add pressure stress
                  pFy(i,j,k,2)=pFy(i,j,k,2)-0.5_WP*sum(pP(i,j-1:j,k,1))
                  ! Velocity gradients at y-face
                  dUdx(1,1)=0.25_WP*dxi*(uc(i+1,j-1,k,1)-uc(i-1,j-1,k,1)+uc(i+1,j,k,1)-uc(i-1,j,k,1))
                  dUdx(2,1)=dyi*(uc(i,j,k,1)-uc(i,j-1,k,1))
                  dUdx(3,1)=0.25_WP*dzi*(uc(i,j-1,k+1,1)-uc(i,j-1,k-1,1)+uc(i,j,k+1,1)-uc(i,j,k-1,1))
                  dUdx(1,2)=0.25_WP*dxi*(uc(i+1,j-1,k,2)-uc(i-1,j-1,k,2)+uc(i+1,j,k,2)-uc(i-1,j,k,2))
                  dUdx(2,2)=dyi*(uc(i,j,k,2)-uc(i,j-1,k,2))
                  dUdx(3,2)=0.25_WP*dzi*(uc(i,j-1,k+1,2)-uc(i,j-1,k-1,2)+uc(i,j,k+1,2)-uc(i,j,k-1,2))
                  dUdx(1,3)=0.25_WP*dxi*(uc(i+1,j-1,k,3)-uc(i-1,j-1,k,3)+uc(i+1,j,k,3)-uc(i-1,j,k,3))
                  dUdx(2,3)=dyi*(uc(i,j,k,3)-uc(i,j-1,k,3))
                  dUdx(3,3)=0.25_WP*dzi*(uc(i,j-1,k+1,3)-uc(i,j-1,k-1,3)+uc(i,j,k+1,3)-uc(i,j,k-1,3))
                  div=dUdx(1,1)+dUdx(2,2)+dUdx(3,3)
                  ! Viscosities at y-face
                  visc_f=2.0_WP*product(pVisc(i,j-1:j,k,1))/(sum(pVisc(i,j-1:j,k,1))+tiny(1.0_WP))
                  beta_f=2.0_WP*product(pBeta(i,j-1:j,k,1))/(sum(pBeta(i,j-1:j,k,1))+tiny(1.0_WP))
                  ! Viscous stress at y-face
                  pFy(i,j,k,1)=pFy(i,j,k,1)+visc_f*(dUdx(1,2)+dUdx(2,1))
                  pFy(i,j,k,2)=pFy(i,j,k,2)+visc_f*(dUdx(2,2)+dUdx(2,2))+(beta_f-2.0_WP/3.0_WP*visc_f)*div
                  pFy(i,j,k,3)=pFy(i,j,k,3)+visc_f*(dUdx(3,2)+dUdx(2,3))
               end do; end do; end do
            end if
            ! Z-fluxes
            if (amr%nz.gt.1) then
               fbx=mfi%nodaltilebox(3)
               do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
                  ! Add pressure stress
                  pFz(i,j,k,3)=pFz(i,j,k,3)-0.5_WP*sum(pP(i,j,k-1:k,1))
                  ! Velocity gradients at z-face
                  dUdx(1,1)=0.25_WP*dxi*(uc(i+1,j,k-1,1)-uc(i-1,j,k-1,1)+uc(i+1,j,k,1)-uc(i-1,j,k,1))
                  dUdx(2,1)=0.25_WP*dyi*(uc(i,j+1,k-1,1)-uc(i,j-1,k-1,1)+uc(i,j+1,k,1)-uc(i,j-1,k,1))
                  dUdx(3,1)=dzi*(uc(i,j,k,1)-uc(i,j,k-1,1))
                  dUdx(1,2)=0.25_WP*dxi*(uc(i+1,j,k-1,2)-uc(i-1,j,k-1,2)+uc(i+1,j,k,2)-uc(i-1,j,k,2))
                  dUdx(2,2)=0.25_WP*dyi*(uc(i,j+1,k-1,2)-uc(i,j-1,k-1,2)+uc(i,j+1,k,2)-uc(i,j-1,k,2))
                  dUdx(3,2)=dzi*(uc(i,j,k,2)-uc(i,j,k-1,2))
                  dUdx(1,3)=0.25_WP*dxi*(uc(i+1,j,k-1,3)-uc(i-1,j,k-1,3)+uc(i+1,j,k,3)-uc(i-1,j,k,3))
                  dUdx(2,3)=0.25_WP*dyi*(uc(i,j+1,k-1,3)-uc(i,j-1,k-1,3)+uc(i,j+1,k,3)-uc(i,j-1,k,3))
                  dUdx(3,3)=dzi*(uc(i,j,k,3)-uc(i,j,k-1,3))
                  div=dUdx(1,1)+dUdx(2,2)+dUdx(3,3)
                  ! Viscosities at z-face
                  visc_f=2.0_WP*product(pVisc(i,j,k-1:k,1))/(sum(pVisc(i,j,k-1:k,1))+tiny(1.0_WP))
                  beta_f=2.0_WP*product(pBeta(i,j,k-1:k,1))/(sum(pBeta(i,j,k-1:k,1))+tiny(1.0_WP))
                  ! Viscous stress at z-face
                  pFz(i,j,k,1)=pFz(i,j,k,1)+visc_f*(dUdx(1,3)+dUdx(3,1))
                  pFz(i,j,k,2)=pFz(i,j,k,2)+visc_f*(dUdx(2,3)+dUdx(3,2))
                  pFz(i,j,k,3)=pFz(i,j,k,3)+visc_f*(dUdx(3,3)+dUdx(3,3))+(beta_f-2.0_WP/3.0_WP*visc_f)*div
               end do; end do; end do
            end if
            deallocate(uc)
         end do
         call amr%mfiter_destroy(mfi)
      end do

      ! Average down face stresses for C/F conservation
      do lvl=top,1,-1
         call amrmfab_average_down_face(fmf=Fx(lvl),cmf=Fx(lvl-1),rr=[amr%rrefx(lvl-1),amr%rrefy(lvl-1),amr%rrefz(lvl-1)],cgeom=amr%geom(lvl-1))
         call amrmfab_average_down_face(fmf=Fy(lvl),cmf=Fy(lvl-1),rr=[amr%rrefx(lvl-1),amr%rrefy(lvl-1),amr%rrefz(lvl-1)],cgeom=amr%geom(lvl-1))
         call amrmfab_average_down_face(fmf=Fz(lvl),cmf=Fz(lvl-1),rr=[amr%rrefx(lvl-1),amr%rrefy(lvl-1),amr%rrefz(lvl-1)],cgeom=amr%geom(lvl-1))
      end do

      ! Stress divergence on the coupling hierarchy
      do lvl=0,top
         dxi=1.0_WP/amr%dx(lvl); if (amr%nx.eq.1) dxi=0.0_WP
         dyi=1.0_WP/amr%dy(lvl); if (amr%ny.eq.1) dyi=0.0_WP
         dzi=1.0_WP/amr%dz(lvl); if (amr%nz.eq.1) dzi=0.0_WP
         call amr%mfiter_build(lvl,mfi,tiling=.false.)
         do while (mfi%next())
            pRHS=>stress_div%mf(lvl)%dataptr(mfi)
            pFx=>Fx(lvl)%dataptr(mfi); pFy=>Fy(lvl)%dataptr(mfi); pFz=>Fz(lvl)%dataptr(mfi)
            bx=mfi%tilebox()
            do n=1,3; do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pRHS(i,j,k,n)=dxi*(pFx(i+1,j,k,n)-pFx(i,j,k,n)) &
               &          +dyi*(pFy(i,j+1,k,n)-pFy(i,j,k,n)) &
               &          +dzi*(pFz(i,j,k+1,n)-pFz(i,j,k,n))
            end do; end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do

      ! Synchronize the output and release flux storage
      do lvl=top-1,0,-1
         call stress_div%average_downto(lvl)
      end do
      do lvl=0,top
         call stress_div%fill_lvl(lvl,time)
         call amrex_multifab_destroy(Fx(lvl))
         call amrex_multifab_destroy(Fy(lvl))
         call amrex_multifab_destroy(Fz(lvl))
      end do
   end subroutine get_div_stress

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
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pRhoU,pRhoV,pRhoW
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
            pRhoU=>fs%rhoU%mf(lvl)%dataptr(mfi)
            pRhoV=>fs%rhoV%mf(lvl)%dataptr(mfi)
            pRhoW=>fs%rhoW%mf(lvl)%dataptr(mfi)
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
                     pRhoU(i,j,k,1)=pRhoU(i,j,k,1)+fac*(rho*u-pRhoU(i,j,k,1))
                  end do; end do; end do
            bx=mfi%nodaltilebox(2)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                     fac=sponge_factor(y=amr%ylo+real(j,WP)*dy,z=amr%zlo+(real(k,WP)+0.5_WP)*dz)
                     if (fac.le.0.0_WP) cycle
                     if (max(pQ(i,j-1,k,1),pQ(i,j,k,1)).gt.VF_spg*rhoL1) cycle
                     pRhoV(i,j,k,1)=pRhoV(i,j,k,1)-fac*pRhoV(i,j,k,1)
                  end do; end do; end do
            bx=mfi%nodaltilebox(3)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                     fac=sponge_factor(y=amr%ylo+(real(j,WP)+0.5_WP)*dy,z=amr%zlo+real(k,WP)*dz)
                     if (fac.le.0.0_WP) cycle
                     if (max(pQ(i,j,k-1,1),pQ(i,j,k,1)).gt.VF_spg*rhoL1) cycle
                     pRhoW(i,j,k,1)=pRhoW(i,j,k,1)-fac*pRhoW(i,j,k,1)
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
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF,pCL,pCG,pU,pC
      real(WP), dimension(3) :: BL,BG
      real(WP) :: dx,dy,dz,myVF,IEL,x_cc,rhoG,pG,H,cL,cG
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
                  ! Set conserved variables
                  pQ(i,j,k,1)=(       myVF)*rhoL1
                  pQ(i,j,k,2)=(1.0_WP-myVF)*rhoG
                  pQ(i,j,k,3)=pQ(i,j,k,1)*IEL
                  pQ(i,j,k,4)=pQ(i,j,k,2)*gas%get_e_from_p_rho(p=pG,rho=rhoG,y=[1.0_WP])
               end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
      ! Set x-face velocity from the same shock profile
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
      ! Set mixture sound speed
      cL=water%get_c_from_p_rho(p=pL1,rho=rhoL1,y=[1.0_WP])
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         pQ=>solver%Q%mf(lvl)%dataptr(mfi)
         pC=>solver%C%mf(lvl)%dataptr(mfi)
         do k=lbound(pC,3),ubound(pC,3); do j=lbound(pC,2),ubound(pC,2); do i=lbound(pC,1),ubound(pC,1)
                  x_cc=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
                  H=Hshock(x=Xs-x_cc,delta=0.5_WP*dx)
                  rhoG=rhoG1+(rhoG2-rhoG1)*H
                  pG  =pG1  +(pG2  -pG1  )*H
                  cG=gas%get_c_from_p_rho(p=pG,rho=rhoG,y=[1.0_WP])
                  pC(i,j,k,1)=sqrt((pQ(i,j,k,1)*cL**2+pQ(i,j,k,2)*cG**2)/max(pQ(i,j,k,1)+pQ(i,j,k,2),solver%rho_floor))
               end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine shockdrop_init

   !> Apply the inflow BC at low-x (boundary momentum comes from apply_momentum_bc)
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

   !> Build the coarsening ramp outside sponge zone
   subroutine build_sponge_bands()
      use messager, only: log,die
      use string,   only: str_long
      implicit none
      character(len=str_long) :: message
      real(WP) :: r,rmax
      integer  :: lvl
      allocate(spg_redge(0:amr%maxlvl)); spg_redge=huge(1.0_WP)
      if (spg_minlvl.lt.0.or.spg_minlvl.gt.amr%maxlvl) call die('[amrcomp_drop_lpt_prod] Minimum sponge level must lie in [0,Max level]')
      ! Largest r_cyl any cell can reach (z collapses in quasi-2D, as in apply_sponge)
      rmax=max(abs(amr%ylo),abs(amr%yhi))
      if (amr%nz.gt.1) rmax=sqrt(rmax**2+max(abs(amr%zlo),abs(amr%zhi))**2)
      ! Walk outward from the end of the ramp, one level coarser per band
      r=R_spg+L_spg
      do lvl=amr%maxlvl,spg_minlvl+1,-1
         r=r+real(spg_bandcells,WP)*amr%dy(lvl)
         spg_redge(lvl)=r
      end do
      ! The floor owns everything beyond the last band, however far that reaches
      spg_redge(0:spg_minlvl)=huge(1.0_WP)
      ! A band the blocking factor will swallow is a band that does not exist
      if (spg_bandcells.lt.amr%nbloc) then
         write(message,'("[Sponge bands] WARNING: ",i0," cells per band is below the blocking factor ",i0,"; bands may be absorbed by grid clustering")') spg_bandcells,amr%nbloc
         call log(message)
      end if
      if (r.gt.rmax) then
         write(message,'("[Sponge bands] WARNING: bands reach r=",es12.5," beyond the domain edge r=",es12.5,"; level ",i0," is never reached")') r,rmax,spg_minlvl
         call log(message)
      end if
      do lvl=amr%maxlvl,spg_minlvl+1,-1
         write(message,'("[Sponge bands] level ",i0," capped out to r=",es12.5)') lvl,spg_redge(lvl)
         call log(message)
      end do
      write(message,'("[Sponge bands] level ",i0," floor from r=",es12.5," to the domain edge r=",es12.5)') spg_minlvl,min(r,rmax),rmax
      call log(message)
   end subroutine build_sponge_bands

   !> Refinement cap at this radius
   integer function sponge_lvlcap(r_cyl,maxlvl)
      implicit none
      real(WP), intent(in) :: r_cyl
      integer,  intent(in) :: maxlvl
      integer :: lvl
      if (r_cyl.le.R_spg+L_spg) then
         sponge_lvlcap=maxlvl; return
      end if
      ! MUST descend
      sponge_lvlcap=spg_minlvl
      do lvl=maxlvl,spg_minlvl+1,-1
         if (r_cyl.le.spg_redge(lvl)) then; sponge_lvlcap=lvl; return; end if
      end do
   end function sponge_lvlcap

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
      real(WP) :: rho_cc,rho_xp,rho_xm,rho_yp,rho_ym,rho_zp,rho_zm,rho_err
      real(WP) :: lapU,lapV,lapW,u_sgs,Re,wramp,x_cc
      real(WP) :: divu,vortx,vorty,vortz,vort,Ducros,Deps
      real(WP) :: r_cyl
      logical  :: in_zone,in_wake,wake_can_refine
      integer :: i,j,k,wake_maxlvl
      real(WP), parameter :: Reps=1.0e-2_WP
      real(WP), parameter :: Peps=1.0e-2_WP
      real(WP), parameter :: Cduc=0.05_WP
      real(WP) :: wake_xstart
      ! Get mesh size
      dx=solver%amr%dx(lvl); dxi=1.0_WP/dx; dxi2=1.0_WP/dx**2
      dy=solver%amr%dy(lvl); dyi=1.0_WP/dy; dyi2=1.0_WP/dy**2
      dz=solver%amr%dz(lvl); dzi=1.0_WP/dz; dzi2=1.0_WP/dz**2
      delta=solver%amr%min_meshsize(lvl); delta2=delta**2
      ! Begin wake coarsening downstream of the padded core box
      wake_xstart=wake_x0
      if (do_transfer.and.corehi(1).gt.corelo(1)) then
         if (wake_x0.ge.0.5_WP*huge(1.0_WP)) then
            wake_xstart=corehi(1)
         else
            wake_xstart=max(wake_x0,corehi(1))
         end if
      end if
      ! AMReX tags level L to create level L+1, so the wake criteria stop tagging at level maxlvl-offset
      wake_maxlvl=max(0,solver%amr%maxlvl-wake_maxlvl_off)
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
         ! Cell-centered velocity over the tile and one ghost layer, from the face velocities
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
                  ! Refinement floor: force tags below min_lvl
                  if (lvl.lt.min_lvl) then
                     tagarr(i,j,k,1)=SETtag
                  end if

                  ! Far-field cap: unrestricted inside the sponge ramp, then one level coarser per band beyond it
                  r_cyl=sqrt((solver%amr%ylo+(real(j,WP)+0.5_WP)*dy)**2+(solver%amr%zlo+(real(k,WP)+0.5_WP)*dz)**2)
                  if (solver%amr%nz.eq.1) r_cyl=abs(solver%amr%ylo+(real(j,WP)+0.5_WP)*dy) ! match apply_sponge in quasi-2D
                  in_zone=(lvl.lt.sponge_lvlcap(r_cyl,solver%amr%maxlvl))

                  ! Mixture density laplacian error
                  rho_cc=sum(pQ(i  ,j,  k,  1:2))
                  rho_xp=sum(pQ(i+1,j,  k,  1:2)); rho_xm=sum(pQ(i-1,j,  k,  1:2))
                  rho_yp=sum(pQ(i,  j+1,k,  1:2)); rho_ym=sum(pQ(i,  j-1,k,  1:2))
                  rho_zp=sum(pQ(i,  j,  k+1,1:2)); rho_zm=sum(pQ(i,  j,  k-1,1:2))
                  rho_err=lap_error(rho_cc,rho_xm,rho_xp,rho_ym,rho_yp,rho_zm,rho_zp,Reps)
                  ! Wake ramp factor, shared by the two turbulence-sensitive criteria
                  wramp=1.0_WP
                  x_cc=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
                  in_wake=x_cc.gt.wake_xstart
                  wake_can_refine=(.not.in_wake.or.lvl.lt.wake_maxlvl)
                  if (in_wake) then
                     if (lvl.gt.wake_lvl) wramp=wramp*wake_ramp**(lvl-wake_lvl)
                     wramp=wramp*2.0_WP**((x_cc-wake_xstart)/wake_len)
                  end if
                  ! Density criterion, capped beyond the sponge and wake-ramped only when wake_rho is set
                  if (rho_err.gt.Rho_tag*merge(wramp,1.0_WP,wake_rho).and.in_zone) then
                     if (.not.(wake_rho.and.in_wake).or.wake_can_refine) tagarr(i,j,k,1)=SETtag
                  end if

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
                  ! SGS criterion, ramped with level and downstream distance in the wake
                  if (Re.gt.Re_tag*wramp.and.in_zone.and.wake_can_refine) tagarr(i,j,k,1)=SETtag

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

   !> ================= Sub-grid drop -> Lagrangian particle transfer =================

   !> Label the current Eulerian liquid, which invalidates the cached core moments
   subroutine build_current_ccl()
      implicit none
      core_mom_valid=.false.
      core_id=0
      call tmr_ccl%start()
      call ccl%build_composite(prepare_fab=ccl_prepare,make_label=ccl_make,same_label=ccl_same)
      call tmr_ccl%stop()
   end subroutine build_current_ccl

   !> Rebind our array pointers for the FAB amrccl is about to scan
   subroutine ccl_prepare(lvl,mfi)
      use amrex_amr_module, only: amrex_mfiter
      implicit none
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(inout) :: mfi
      pVF  =>fs%VF%mf(lvl)%dataptr(mfi)
      cdx=amr%dx(lvl); cdy=amr%dy(lvl); cdz=amr%dz(lvl)
      nullify(pPLIC)
      if (lvl.eq.amr%maxlvl) pPLIC=>fs%PLIC%dataptr(mfi)
   end subroutine ccl_prepare

   !> A cell belongs to a structure if it holds any liquid
   logical function ccl_make(i1,j1,k1)
      use amrmpcomp_class, only: VFlo
      implicit none
      integer, intent(in) :: i1,j1,k1
      ccl_make=(pVF(i1,j1,k1,1).gt.VFlo)
   end function ccl_make

   !> Thin-film split test after shock_droplet: opposed interfaces that face each other get different labels
   logical function ccl_same(i1,j1,k1,i2,j2,k2)
      implicit none
      integer, intent(in) :: i1,j1,k1,i2,j2,k2
      real(WP), dimension(3) :: N1,N2,O1,O2,X1,X2
      ! Big default, assume same label
      ccl_same=.true.
      ! PLIC exists only at maxlvl, so coarser liquid uses face connectivity
      if (.not.associated(pPLIC)) return
      ! Both cells must carry a real interface
      if (abs(pPLIC(i1,j1,k1,4)).gt.1.0e9_WP) return
      if (abs(pPLIC(i2,j2,k2,4)).gt.1.0e9_WP) return
      N1=pPLIC(i1,j1,k1,1:3)
      N2=pPLIC(i2,j2,k2,1:3)
      ! If not at least ~72 degrees apart, return
      if (dot_product(N1,N2).ge.0.3_WP) return
      ! Foot of the perpendicular from each cell center onto its own plane
      X1=[amr%xlo+(real(i1,WP)+0.5_WP)*cdx,amr%ylo+(real(j1,WP)+0.5_WP)*cdy,amr%zlo+(real(k1,WP)+0.5_WP)*cdz]
      X2=[amr%xlo+(real(i2,WP)+0.5_WP)*cdx,amr%ylo+(real(j2,WP)+0.5_WP)*cdy,amr%zlo+(real(k2,WP)+0.5_WP)*cdz]
      O1=X1-(dot_product(N1,X1)-pPLIC(i1,j1,k1,4))*N1
      O2=X2-(dot_product(N2,X2)-pPLIC(i2,j2,k2,4))*N2
      ! If pointing towards one another, use different labels
      if (dot_product(O1-O2,N1).lt.0.0_WP.and.dot_product(O2-O1,N2).lt.0.0_WP) ccl_same=.false.
   end function ccl_same


   !> Mixture kinetic energy of a cell, a copy of amrmpcomp_class's private cell_ke (keep identical)
   pure real(WP) function cell_ke(pQ,pU,pV,pW,i,j,k,rho_floor) result(ke)
      implicit none
      real(WP), dimension(:,:,:,:), contiguous, pointer, intent(in) :: pQ,pU,pV,pW
      integer, intent(in) :: i,j,k
      real(WP), intent(in) :: rho_floor
      ke=0.25_WP*(max(0.5_WP*sum(pQ(i-1:i,j,k,1:2)),rho_floor)*pU(i,j,k,1)**2+max(0.5_WP*sum(pQ(i:i+1,j,k,1:2)),rho_floor)*pU(i+1,j,k,1)**2 &
      &          +max(0.5_WP*sum(pQ(i,j-1:j,k,1:2)),rho_floor)*pV(i,j,k,1)**2+max(0.5_WP*sum(pQ(i,j:j+1,k,1:2)),rho_floor)*pV(i,j+1,k,1)**2 &
      &          +max(0.5_WP*sum(pQ(i,j,k-1:k,1:2)),rho_floor)*pW(i,j,k,1)**2+max(0.5_WP*sum(pQ(i,j,k:k+1,1:2)),rho_floor)*pW(i,j,k+1,1)**2)
   end function cell_ke

   !> Maximum liquid volume admitted by the equivalent-sphere diameter limit.
   real(WP) function transfer_volume_limit() result(vmax)
      use mathtools, only: Pi
      vmax=Pi/6.0_WP*xfer_dmax**3
   end function transfer_volume_limit

   real(WP) function transfer_min_diameter() result(dmin)
      dmin=xfer_diam_threshold*amr%min_meshsize(amr%maxlvl)
   end function transfer_min_diameter

   !> Label the liquid and reduce the moments of the transfer candidates and the core; call before any regrid
   subroutine label_and_moments(cid,mom,cmap,ns)
      use mpi_f08,          only: MPI_ALLREDUCE,MPI_ALLGATHER,MPI_ALLGATHERV,MPI_IN_PLACE,MPI_SUM,MPI_MIN,MPI_MAX,MPI_INTEGER
      use parallel,         only: MPI_REAL_WP
      use mathtools,        only: twoPi
      use amrmpcomp_class,  only: VFhi
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      implicit none
      integer,  dimension(:),   allocatable, intent(out) :: cid
      real(WP), dimension(:,:), allocatable, intent(out) :: mom
      integer,  dimension(:),   allocatable, intent(out) :: cmap
      integer, intent(out) :: ns
      integer, parameter :: KTOP=4
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pId,pVF,pQ,pUVW,pCL,pU,pV,pW
      ! Rank-local table: slot(label) -> local slot, sid(slot) -> label
      integer,  dimension(:),   allocatable :: slot,sid
      real(WP), dimension(:,:), allocatable :: lmom,lbox,lle
      logical,  dimension(:),   allocatable :: lcoarse
      integer :: nloc,cap
      ! Selection
      integer,  dimension(KTOP) :: topid
      real(WP), dimension(KTOP) :: topvol
      integer,  dimension(:), allocatable :: alltop,tid,nom,allnom,counts,displs
      real(WP), dimension(:), allocatable :: tbuf,vall
      real(WP), dimension(7) :: minbuf
      real(WP), dimension(2) :: rejbuf
      real(WP), dimension(3) :: L,Olo,X,cell_second
      real(WP) :: vol,theta,mi,xface,boxmargin,volmax,v,ag,ig
      real(WP), dimension(2) :: masses
      integer :: lvl,i,j,k,n,d,s,m,p,nT,nnom,ntot,nU,ierr

      ! Coarse liquid is included even when the hierarchy has not yet reached maxlvl.
      core_mom=0.0_WP; core_mle=+huge(1.0_WP); core_mom_valid=.false.; core_id=0
      L=[amr%xhi-amr%xlo,amr%yhi-amr%ylo,amr%zhi-amr%zlo]
      Olo=[amr%xlo,amr%ylo,amr%zlo]

      ! Label the liquid
      call build_current_ccl()
      ns=ccl%nstruct; ccl_nstruct=ns
      if (ns.eq.0) then
         allocate(cid(0),mom(NMOM,0),cmap(0))
         xfer_ncand=0
         return
      end if

      ! Accumulate rank-local moments per label (rows 9:14 are circular sums, kept only for wrapped directions)
      call tmr_msweep%start()
      allocate(slot(ns)); slot=0
      cap=1024; nloc=0
      allocate(sid(cap),lmom(NMOM,cap),lbox(6,cap),lle(4,cap),lcoarse(cap))
      ! Composite IDs are zero on covered coarse cells, so each physical volume contributes once.
      do lvl=0,amr%clvl()
         vol=amr%cell_vol(lvl)
         cdx=amr%dx(lvl); cdy=amr%dy(lvl); cdz=amr%dz(lvl)
         call amr%mfiter_build(lvl,mfi,tiling=.false.)
         do while (mfi%next())
            pId =>ccl%id%mf(lvl)%dataptr(mfi)
            pVF =>fs%VF%mf(lvl)%dataptr(mfi)
            pQ  =>fs%Q%mf(lvl)%dataptr(mfi)
            pUVW=>fs%UVW%mf(lvl)%dataptr(mfi)
            pU=>fs%U%mf(lvl)%dataptr(mfi); pV=>fs%V%mf(lvl)%dataptr(mfi); pW=>fs%W%mf(lvl)%dataptr(mfi)
            if (lvl.eq.amr%maxlvl) pCL =>fs%CL%dataptr(mfi)
            bx=mfi%validbox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                     n=nint(pId(i,j,k,1))
                     if (n.le.0) cycle
                     s=slot(n); if (s.eq.0) s=new_slot(n)
                     X=[amr%xlo+(real(i,WP)+0.5_WP)*cdx, &
                        amr%ylo+(real(j,WP)+0.5_WP)*cdy,amr%zlo+(real(k,WP)+0.5_WP)*cdz]
                     if (lvl.eq.amr%maxlvl) X=pCL(i,j,k,1:3)
                     ! Add the exact within-cell second moment for pure-liquid cells
                     cell_second=0.0_WP
                     if (pVF(i,j,k,1).ge.VFhi) cell_second=[cdx**2,cdy**2,cdz**2]/12.0_WP
                     mi=vol*pQ(i,j,k,1)
                     lmom( 1,s)=lmom( 1,s)+vol*pVF(i,j,k,1)
                     lmom( 2,s)=lmom( 2,s)+mi
                     lmom( 3:5 ,s)=lmom( 3:5 ,s)+mi*X
                     lmom( 6:8 ,s)=lmom( 6:8 ,s)+mi*pUVW(i,j,k,1:3)
                     do d=1,3
                        if (ccl%struct(n)%per(d).ne.1) cycle
                        theta=twoPi*(X(d)-Olo(d))/L(d)
                        lmom( 8+d,s)=lmom( 8+d,s)+mi*cos(theta)
                        lmom(11+d,s)=lmom(11+d,s)+mi*sin(theta)
                     end do
                     lmom(15,s)=lmom(15,s)+vol*pQ(i,j,k,3)
                     lmom(16,s)=lmom(16,s)+mi*(X(1)**2+cell_second(1))
                     lmom(17,s)=lmom(17,s)+mi*(X(2)**2+cell_second(2))
                     lmom(18,s)=lmom(18,s)+mi*(X(3)**2+cell_second(3))
                     lmom(19,s)=lmom(19,s)+mi*X(1)*X(2)
                     lmom(20,s)=lmom(20,s)+mi*X(1)*X(3)
                     lmom(21,s)=lmom(21,s)+mi*X(2)*X(3)
                     ! Donor mean: gas-volume-weighted density and internal-energy density.
                     ! Use exactly the same staggered kinetic-energy subtraction
                     ag=1.0_WP-pVF(i,j,k,1)
                     if (ag.ge.0.5_WP.and.pQ(i,j,k,2).gt.0.0_WP) then
                        masses=max(pQ(i,j,k,1:2),0.0_WP)
                        ig=pQ(i,j,k,4)-masses(2)/max(sum(masses),fs%rho_floor)*cell_ke(pQ,pU,pV,pW,i,j,k,fs%rho_floor)
                        if (ig.gt.0.0_WP) lmom(22:24,s)=lmom(22:24,s)+vol*[ag,pQ(i,j,k,2),ig]
                     end if
                     lbox(1:3,s)=min(lbox(1:3,s),[amr%xlo+real(i  ,WP)*cdx,amr%ylo+real(j  ,WP)*cdy,amr%zlo+real(k  ,WP)*cdz])
                     lbox(4:6,s)=max(lbox(4:6,s),[amr%xlo+real(i+1,WP)*cdx,amr%ylo+real(j+1,WP)*cdy,amr%zlo+real(k+1,WP)*cdz])
                     ! The leading edge is the windward cell face
                     xface=amr%xlo+real(i,WP)*cdx
                     lle(1,s)=min(lle(1,s),xface)
                     if (pVF(i,j,k,1).ge.0.01_WP) lle(2,s)=min(lle(2,s),xface)
                     if (pVF(i,j,k,1).ge.0.10_WP) lle(3,s)=min(lle(3,s),xface)
                     if (pVF(i,j,k,1).ge.0.50_WP) lle(4,s)=min(lle(4,s),xface)
                     if (lvl.lt.amr%maxlvl) lcoarse(s)=.true.
                  end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
      call tmr_msweep%stop()

      call tmr_mcomm%start()

      ! 1) Core: largest liquid volume, certified from the per-rank top-KTOP sets
      topid=0; topvol=-1.0_WP
      do s=1,nloc
         v=lmom(1,s)
         if (v.le.topvol(KTOP)) cycle
         p=KTOP
         do while (p.gt.1)
            if (topvol(p-1).ge.v) exit
            topvol(p)=topvol(p-1); topid(p)=topid(p-1); p=p-1
         end do
         topvol(p)=v; topid(p)=sid(s)
      end do
      allocate(alltop(KTOP*amr%nproc))
      call MPI_ALLGATHER(topid,KTOP,MPI_INTEGER,alltop,KTOP,MPI_INTEGER,amr%comm,ierr)
      ! cmap doubles as a mark array here; it is cleared again before the union is built
      allocate(cmap(ns)); cmap=0
      nT=0
      do m=1,size(alltop)
         n=alltop(m)
         if (n.le.0) cycle
         if (cmap(n).gt.0) cycle
         nT=nT+1; cmap(n)=nT
      end do
      allocate(tid(nT),tbuf(nT+1))
      do m=1,size(alltop)
         n=alltop(m)
         if (n.gt.0) tid(cmap(n))=n
      end do
      tbuf=0.0_WP
      do m=1,nT
         s=slot(tid(m))
         if (s.gt.0) tbuf(m)=lmom(1,s)
      end do
      ! Upper bound contribution for labels outside the set: this rank's KTOP-th largest local volume
      tbuf(nT+1)=max(topvol(KTOP),0.0_WP)
      call MPI_ALLREDUCE(MPI_IN_PLACE,tbuf,nT+1,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
      v=-1.0_WP
      do m=1,nT
         if (tbuf(m).gt.v.or.(tbuf(m).eq.v.and.tid(m).lt.core_id)) then
            v=tbuf(m); core_id=tid(m)
         end if
      end do
      do m=1,nT
         cmap(tid(m))=0
      end do
      if (.not.(v.gt.tbuf(nT+1))) then
         ! Not certified: reduce every label's volume, as the full reduction did
         allocate(vall(ns)); vall=0.0_WP
         do s=1,nloc
            vall(sid(s))=lmom(1,s)
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,vall,ns,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
         core_id=maxloc(vall,dim=1)
         deallocate(vall)
         core_fallbacks=core_fallbacks+1.0_WP
      end if

      ! 2) Core box and leading edges
      s=slot(core_id)
      if (s.gt.0) then
         minbuf(1:3)=lbox(1:3,s); minbuf(4:7)=lle(:,s); corehi=lbox(4:6,s)
      else
         minbuf=+huge(1.0_WP); corehi=-huge(1.0_WP)
      end if
      call MPI_ALLREDUCE(MPI_IN_PLACE,minbuf,7,MPI_REAL_WP,MPI_MIN,amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,corehi,3,MPI_REAL_WP,MPI_MAX,amr%comm,ierr)
      corelo=minbuf(1:3); core_mle=minbuf(4:7)
      ! Preserve the padded core box for wake tagging; it no longer restricts transfer
      boxmargin=xfer_margin; if (boxmargin.lt.0.0_WP) boxmargin=0.5_WP*d_drop
      corelo=corelo-boxmargin; corehi=corehi+boxmargin

      ! 3) Local nomination: a local volume is a lower bound on the component volume
      volmax=transfer_volume_limit()
      rejbuf=0.0_WP
      allocate(nom(max(nloc,1))); nnom=0
      do s=1,nloc
         n=sid(s)
         if (n.eq.core_id) cycle
         if (lcoarse(s)) then
            rejbuf(1)=rejbuf(1)+1.0_WP; cycle
         end if
         if (lmom(1,s).gt.volmax) then
            rejbuf(2)=rejbuf(2)+1.0_WP; cycle
         end if
         nnom=nnom+1; nom(nnom)=n
      end do

      ! 4) Union of the nominations, plus the core, in ascending label order
      allocate(counts(amr%nproc),displs(amr%nproc))
      call MPI_ALLGATHER(nnom,1,MPI_INTEGER,counts,1,MPI_INTEGER,amr%comm,ierr)
      displs(1)=0
      do m=2,amr%nproc
         displs(m)=displs(m-1)+counts(m-1)
      end do
      ntot=sum(counts)
      allocate(allnom(max(ntot,1)))
      call MPI_ALLGATHERV(nom,nnom,MPI_INTEGER,allnom,counts,displs,MPI_INTEGER,amr%comm,ierr)
      do m=1,ntot
         cmap(allnom(m))=1
      end do
      cmap(core_id)=1
      nU=0
      do n=1,ns
         if (cmap(n).eq.0) cycle
         nU=nU+1; cmap(n)=nU
      end do
      allocate(cid(nU))
      do n=1,ns
         if (cmap(n).gt.0) cid(cmap(n))=n
      end do

      ! 5) Packed reduction of the union's moments, and the prefilter census
      allocate(mom(NMOM,nU)); mom=0.0_WP
      do m=1,nU
         s=slot(cid(m))
         if (s.gt.0) mom(:,m)=lmom(:,s)
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,mom,NMOM*nU,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,rejbuf,2,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
      call tmr_mcomm%stop()
      xfer_ncand=nU; xfer_rej=rejbuf

      ! Cache the core's row for the metrics, since transfer_drops deallocates mom
      core_mom=mom(:,cmap(core_id))
      core_mom_valid=.true.

      deallocate(slot,sid,lmom,lbox,lle,lcoarse,alltop,tid,tbuf,nom,allnom,counts,displs)

   contains

      !> Open a local slot for label lab
      integer function new_slot(lab) result(sn)
         integer, intent(in) :: lab
         if (nloc.eq.cap) call grow()
         nloc=nloc+1; sn=nloc
         slot(lab)=sn; sid(sn)=lab
         lmom(:,sn)=0.0_WP
         lbox(1:3,sn)=+huge(1.0_WP); lbox(4:6,sn)=-huge(1.0_WP)
         lle(:,sn)=+huge(1.0_WP)
         lcoarse(sn)=.false.
      end function new_slot

      !> Double the local table
      subroutine grow()
         integer,  dimension(:),   allocatable :: itmp
         real(WP), dimension(:,:), allocatable :: rtmp
         logical,  dimension(:),   allocatable :: ltmp
         cap=2*cap
         allocate(itmp(cap)); itmp(1:nloc)=sid(1:nloc); call move_alloc(itmp,sid)
         allocate(rtmp(NMOM,cap)); rtmp(:,1:nloc)=lmom(:,1:nloc); call move_alloc(rtmp,lmom)
         allocate(rtmp(6,cap)); rtmp(:,1:nloc)=lbox(:,1:nloc); call move_alloc(rtmp,lbox)
         allocate(rtmp(4,cap)); rtmp(:,1:nloc)=lle(:,1:nloc); call move_alloc(rtmp,lle)
         allocate(ltmp(cap)); ltmp(1:nloc)=lcoarse(1:nloc); call move_alloc(ltmp,lcoarse)
      end subroutine grow

   end subroutine label_and_moments

   !> Validate the new particle against the liquid EOS before changing mesh or particle state.
   logical function particle_birth_state(mass,volume,momentum,energy,eint) result(valid)
      real(WP), intent(in) :: mass,volume,momentum(3),energy
      real(WP), intent(out) :: eint
      real(WP) :: rho,Tp,velocity(3)
      valid=.false.; eint=0.0_WP
      if (.not.all(particle_finite([mass,volume,momentum,energy]))) return
      if (mass.le.0.0_WP.or.volume.le.0.0_WP) return
      velocity=momentum/mass; rho=mass/volume
      eint=energy/mass-0.5_WP*dot_product(velocity,velocity)
      if (.not.all(particle_finite([rho,eint]))) return
      Tp=water%get_T_from_rho_e(rho=rho,e=eint,y=[1.0_WP])
      valid=particle_finite(Tp).and.Tp.gt.0.0_WP
   end function particle_birth_state

   !> Remove unresolved liquid structures, turning them into particles or flotsam, and rebuild the gas state
   subroutine transfer_drops()
      use mpi_f08,          only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM
      use parallel,         only: MPI_REAL_WP
      use mathtools,        only: Pi
      use messager,         only: log
      use string,           only: str_long
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_multifab
      implicit none
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx,fbx
      type(amrex_multifab) :: work
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pId,pVF,pQ,pCL,pCG,pPLIC,pWk
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW,pRhoU,pRhoV,pRhoW
      integer,  dimension(:),   allocatable :: cid,cmap  !< Reduced labels and label -> index map, from label_and_moments
      real(WP), dimension(:,:), allocatable :: mom       !< (1:NMOM,1:size(cid)), filled by label_and_moments
      logical , dimension(:  ), allocatable :: remove,xfer
      real(WP), dimension(:), allocatable :: birth_e
      type(part), dimension(:), allocatable, target :: pnew
      real(WP), dimension(3) :: L,Olo,X,Bd,Ud
      real(WP), dimension(2) :: mk
      real(WP) :: vol,volmax,Vd,Md,gasadd,ag
      real(WP), dimension(4) :: before,after,removed
      integer :: lvl,ilvl,i,j,k,n,m,nc,ns,np_new,nremove,ncreate,ierr
      character(len=str_long) :: msg
      xfer_n_step=0; flot_n_step=0
      ! amrvof only carries PLIC/CL/CG at maxlvl, and only once the hierarchy has grown into it
      if (amr%clvl().lt.amr%maxlvl) return
      lvl=amr%maxlvl
      vol=amr%cell_vol(lvl)
      L=[amr%xhi-amr%xlo,amr%yhi-amr%ylo,amr%zhi-amr%zlo]
      Olo=[amr%xlo,amr%ylo,amr%zlo]

      ! Label the liquid and measure the candidates (also sets the core box and core_mom)
      call label_and_moments(cid,mom,cmap,ns)
      if (ns.eq.0) then
         deallocate(cid,mom,cmap)
         return
      end if
      nc=size(cid)
      cdx=amr%dx(lvl); cdy=amr%dy(lvl); cdz=amr%dz(lvl)

      ! Decide what goes; every rank holds the same global moments, so all ranks agree
      allocate(remove(nc),xfer(nc),birth_e(nc)); remove=.false.; xfer=.false.; birth_e=0.0_WP
      volmax=transfer_volume_limit()
      do m=1,nc
         n=cid(m)
         if (n.eq.core_id) cycle                    ! never transfer the core
         if (mom(1,m).gt.volmax) then               ! still resolved
            xfer_fin(1)=xfer_fin(1)+1.0_WP; cycle
         end if
         if (mom(2,m).le.0.0_WP) cycle              ! nothing to transfer
         if (ccl%struct(n)%open) cycle
         if (ccl%struct(n)%has_coarse) then         ! removal below operates only at maxlvl
            xfer_fin(2)=xfer_fin(2)+1.0_WP; cycle
         end if
         ! Refuse atomically, before creating a particle or booking a flotsam deletion.
         if (mom(22,m).le.0.0_WP) then
            xfer_fin(3)=xfer_fin(3)+1.0_WP
            xfer_refused=xfer_refused+1.0_WP
            cycle
         end if
         remove(m)=.true.
         ! Create a particle at or above the minimum diameter; strictly smaller structures are flotsam
         Vd=mom(1,m)
         if ((6.0_WP*Vd/Pi)**(1.0_WP/3.0_WP).ge.transfer_min_diameter()) xfer(m)=.true.
         if (xfer(m)) then
            if (.not.particle_birth_state(mom(2,m),Vd,mom(6:8,m),mom(15,m),birth_e(m))) then
               remove(m)=.false.; xfer(m)=.false.
               xfer_fin(4)=xfer_fin(4)+1.0_WP; xfer_badthermal=xfer_badthermal+1.0_WP
            end if
         end if
      end do
      nremove=count(remove)
      ncreate=count(xfer)
      xfer_n_step=ncreate
      flot_n_step=nremove-ncreate
      flot_n_total=flot_n_total+real(flot_n_step,WP)

      ! No births, redistribution, or mesh writes when every component is refused.
      if (nremove.eq.0) then
         deallocate(cid,mom,cmap,remove,xfer,birth_e)
         return
      end if
      call transfer_mesh_integrals(before)
      removed=0.0_WP
      do m=1,nc
         if (remove(m)) removed=removed+[mom(6:8,m),mom(15,m)]
      end do

      ! Build the new particles on root, then append collectively
      call tmr_xpart%start()
      allocate(pnew(max(ncreate,1)))
      np_new=0
      do m=1,nc
         if (.not.remove(m)) cycle
         Vd=mom(1,m); Md=mom(2,m)
         if (.not.xfer(m)) then
            flot_vol =flot_vol +Vd
            flot_mass=flot_mass+Md
            flot_mom=flot_mom+mom(6:8,m)
            flot_energy=flot_energy+mom(15,m)
            cycle
         end if
         Ud=mom(6:8,m)/Md
         Bd=struct_barycenter(m,cid(m))
         ! Running totals for the conservation monitor: what leaves the Eulerian field
         xfer_vol   =xfer_vol   +Vd
         xfer_mass  =xfer_mass  +Md
         xfer_momx  =xfer_momx  +Md*Ud(1)
         xfer_momy  =xfer_momy  +Md*Ud(2)
         xfer_momz  =xfer_momz  +Md*Ud(3)
         xfer_energy=xfer_energy+mom(15,m)   ! Q(3) is the liquid total energy, its kinetic part included
         if (.not.amr%amRoot) cycle
         np_new=np_new+1
         pnew(np_new)%pos   =Bd
         pnew(np_new)%d     =(6.0_WP*Vd/Pi)**(1.0_WP/3.0_WP)   ! Preserve volume independently of mass
         pnew(np_new)%mass  =Md
         pnew(np_new)%eint  =birth_e(m)
         pnew(np_new)%vel   =Ud
         pnew(np_new)%angVel=0.0_WP
         pnew(np_new)%Acol  =0.0_WP
         pnew(np_new)%Tcol  =0.0_WP
         pnew(np_new)%dt    =0.0_WP                            ! let it find its own sub-step
         pnew(np_new)%flag  =PART_MOVES+PART_EXCHANGES
         if (lp_collide) pnew(np_new)%flag=pnew(np_new)%flag+PART_COLLIDES
      end do
      call lp%append(pnew,int(np_new,I8))
      ! The default maxlvl keeps new particles at the particle level cap
      call lp%redistribute()
      deallocate(pnew)
      call tmr_xpart%stop()

      call tmr_surg%start()

      ! Remove each cell's mass share of the kinetic energy, so the surgery works on internal energies
      do ilvl=0,amr%clvl()
         call amr%mfiter_build(ilvl,mfi)
         do while (mfi%next())
            pQ=>fs%Q%mf(ilvl)%dataptr(mfi)
            pU=>fs%U%mf(ilvl)%dataptr(mfi); pV=>fs%V%mf(ilvl)%dataptr(mfi); pW=>fs%W%mf(ilvl)%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                     mk=max(pQ(i,j,k,1:2),0.0_WP)
                     pQ(i,j,k,3:4)=pQ(i,j,k,3:4)-mk/max(sum(mk),fs%rho_floor)*cell_ke(pQ,pU,pV,pW,i,j,k,fs%rho_floor)
                  end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do

      ! One mask marks cells emptied by surgery for the face-momentum rebuild.
      call amr%mfab_build(lvl=lvl,mfab=work,ncomp=1,nover=1)
      call work%setval(0.0_WP)
      gasadd=0.0_WP

      ! Pass A: delete the eligible liquid from the Eulerian fields
      call amr%mfiter_build(lvl,mfi,tiling=.false.)
      do while (mfi%next())
         pId   =>ccl%id%mf(lvl)%dataptr(mfi)
         pVF   =>fs%VF%mf(lvl)%dataptr(mfi)
         pQ    =>fs%Q%mf(lvl)%dataptr(mfi)
         pCL   =>fs%CL%dataptr(mfi)
         pCG   =>fs%CG%dataptr(mfi)
         pPLIC =>fs%PLIC%dataptr(mfi)
         pWk   =>work%dataptr(mfi)
         bx=mfi%validbox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  n=nint(pId(i,j,k,1))
                  if (n.le.0) cycle
                  m=cmap(n)
                  if (m.eq.0) cycle
                  if (.not.remove(m)) cycle
                  pWk(i,j,k,1)=1.0_WP
                  ag=1.0_WP-pVF(i,j,k,1)
                  gasadd=gasadd-vol*pQ(i,j,k,2)
                  if (ag.ge.0.5_WP.and.pQ(i,j,k,2).gt.0.0_WP.and.pQ(i,j,k,4).gt.0.0_WP) then
                     ! Keep the cell's own gas density and specific internal energy.
                     pQ(i,j,k,2)=pQ(i,j,k,2)/ag
                     pQ(i,j,k,4)=pQ(i,j,k,4)/ag
                  else
                     ! All ranks have the same component donor mean; no neighbor flood is needed.
                     pQ(i,j,k,2)=mom(23,m)/mom(22,m)
                     pQ(i,j,k,4)=mom(24,m)/mom(22,m)
                  end if
                  gasadd=gasadd+vol*pQ(i,j,k,2)
                  pQ(i,j,k,1)=0.0_WP
                  pQ(i,j,k,3)=0.0_WP
                  ! Clear the liquid volume moments and plane, as apply_relax does for a cell that turns pure gas
                  pVF(i,j,k,1)=0.0_WP
                  X=[amr%xlo+(real(i,WP)+0.5_WP)*cdx,amr%ylo+(real(j,WP)+0.5_WP)*cdy,amr%zlo+(real(k,WP)+0.5_WP)*cdz]
                  pCL(i,j,k,1:3)=X
                  pCG(i,j,k,1:3)=X
                  pPLIC(i,j,k,:)=[0.0_WP,0.0_WP,0.0_WP,-1.0e10_WP]
               end do; end do; end do
      end do
      call amr%mfiter_destroy(mfi)
      ! Roll up the gas the surgery had to invent
      call MPI_ALLREDUCE(MPI_IN_PLACE,gasadd,1,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
      xfer_gasadd=xfer_gasadd+gasadd

      call tmr_surg%stop()
      call tmr_resync%start()

      ! Average down and fill VF and Q
      call fs%VF%average_down(); call fs%Q%average_down()
      call fs%fill(time=time%t); call fs%Q%fill(time=time%t)

      ! Rebuild: on faces next to emptied cells, keep the face velocity and rebuild momentum from the new face density
      call work%fill_boundary(amr%geom(lvl))
      call amr%mfiter_build(lvl,mfi,tiling=.false.)
      do while (mfi%next())
         pQ=>fs%Q%mf(lvl)%dataptr(mfi)
         pWk=>work%dataptr(mfi)
         pU=>fs%U%mf(lvl)%dataptr(mfi); pRhoU=>fs%rhoU%mf(lvl)%dataptr(mfi)
         pV=>fs%V%mf(lvl)%dataptr(mfi); pRhoV=>fs%rhoV%mf(lvl)%dataptr(mfi)
         pW=>fs%W%mf(lvl)%dataptr(mfi); pRhoW=>fs%rhoW%mf(lvl)%dataptr(mfi)
         fbx=mfi%nodaltilebox(1)
         do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
                  if (max(pWk(i-1,j,k,1),pWk(i,j,k,1)).lt.0.5_WP) cycle
                  pRhoU(i,j,k,1)=max(0.5_WP*sum(pQ(i-1:i,j,k,1:2)),fs%rho_floor)*pU(i,j,k,1)
               end do; end do; end do
         fbx=mfi%nodaltilebox(2)
         do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
                  if (max(pWk(i,j-1,k,1),pWk(i,j,k,1)).lt.0.5_WP) cycle
                  pRhoV(i,j,k,1)=max(0.5_WP*sum(pQ(i,j-1:j,k,1:2)),fs%rho_floor)*pV(i,j,k,1)
               end do; end do; end do
         fbx=mfi%nodaltilebox(3)
         do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
                  if (max(pWk(i,j,k-1,1),pWk(i,j,k,1)).lt.0.5_WP) cycle
                  pRhoW(i,j,k,1)=max(0.5_WP*sum(pQ(i,j,k-1:k,1:2)),fs%rho_floor)*pW(i,j,k,1)
               end do; end do; end do
      end do
      call amr%mfiter_destroy(mfi)
      call amr%mfab_destroy(work)

      ! Back to total energies (this also recovers face velocities and primitives), then the cell velocity
      call fs%get_total_energy(time=time%t)
      call fs%get_cell_velocity()
      call transfer_mesh_integrals(after)
      xfer_refill=xfer_refill+after-before+removed

      ! Rebuild polygons and particle VF here, unless the regrid below does it
      if (.not.regrid_evt%occurs()) then
         call fs%build_polygons()
         call lp%update_VF()
      end if

      call tmr_resync%stop()

      if (amr%amRoot) then
         write(msg,'(">> Removed ",i0," Eulerian structure(s): transferred ",i0," to spray, discarded ",i0," flotsam(s) at t=",es12.5,"; cumulative mass ",es12.5," kg")') &
         &   nremove,ncreate,flot_n_step,time%t,xfer_mass
         call log(trim(msg))
      end if

      deallocate(cid,mom,cmap,remove,xfer,birth_e)

   contains

      !> Mass-weighted barycenter of candidate m, using the circular mean in wrapped directions
      function struct_barycenter(m,n) result(X0)
         use mathtools, only: twoPi
         implicit none
         integer, intent(in) :: m,n
         real(WP), dimension(3) :: X0
         real(WP) :: th
         integer :: dd
         do dd=1,3
            if (ccl%struct(n)%per(dd).eq.1) then
               th=atan2(mom(11+dd,m),mom(8+dd,m))
               X0(dd)=Olo(dd)+L(dd)*modulo(th/twoPi,1.0_WP)
            else
               X0(dd)=mom(2+dd,m)/mom(2,m)
            end if
         end do
      end function struct_barycenter

   end subroutine transfer_drops


   !> Composite cell-centered momentum quadrature and total energy accounting.
   !> Includes gas neighbors whose face masses/kinetic energies change during refill.
   subroutine transfer_mesh_integrals(totals,mass)
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
      use amrex_interface, only: amrmask_make_fine
      use mpi_f08, only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM
      use parallel, only: MPI_REAL_WP
      real(WP), intent(out) :: totals(4)
      real(WP), intent(out), optional :: mass(2)
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_imultifab) :: mask
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pRhoU,pRhoV,pRhoW
      integer, dimension(:,:,:,:), contiguous, pointer :: pMask
      real(WP) :: vol,buf(6)
      integer :: lvl,i,j,k,ierr
      buf=0.0_WP
      do lvl=0,amr%clvl()
         vol=amr%cell_vol(lvl)
         if (lvl.lt.amr%clvl()) then
            call amrex_imultifab_build(mask,amr%ba(lvl),amr%dm(lvl),1,0)
            call amrmask_make_fine(mask,amr%ba(lvl+1),[amr%rrefx(lvl),amr%rrefy(lvl),amr%rrefz(lvl)],0,1)
         end if
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pQ=>fs%Q%mf(lvl)%dataptr(mfi)
            pRhoU=>fs%rhoU%mf(lvl)%dataptr(mfi); pRhoV=>fs%rhoV%mf(lvl)%dataptr(mfi); pRhoW=>fs%rhoW%mf(lvl)%dataptr(mfi)
            if (lvl.lt.amr%clvl()) pMask=>mask%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               if (lvl.lt.amr%clvl()) then
                  if (pMask(i,j,k,1).eq.0) cycle
               end if
               buf=buf+vol*[0.5_WP*sum(pRhoU(i:i+1,j,k,1)),0.5_WP*sum(pRhoV(i,j:j+1,k,1)), &
                  0.5_WP*sum(pRhoW(i,j,k:k+1,1)),sum(pQ(i,j,k,3:4)),pQ(i,j,k,1),pQ(i,j,k,2)]
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
         if (lvl.lt.amr%clvl()) call amrex_imultifab_destroy(mask)
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,buf,6,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
      totals=buf(1:4)
      if (present(mass)) mass=buf(5:6)
   end subroutine transfer_mesh_integrals

   !> ================= Burst output (see burst_implementation.md) =================

   !> Output folder of burst ib, inside amrviz/
   function burst_name(ib) result(name)
      use string, only: rtoa
      implicit none
      integer, intent(in) :: ib
      character(len=str_medium) :: name
      name='burst_tau_'//trim(adjustl(rtoa(burst_tau(ib))))
   end function burst_name

   !> True if the current time lies inside the window of burst ib
   logical function burst_in(ib)
      implicit none
      integer, intent(in) :: ib
      burst_in=(time%t.ge.burst_t0(ib).and.time%t.le.burst_t1(ib))
   end function burst_in

   !> True if the current time lies inside any burst window
   logical function burst_any()
      implicit none
      integer :: ib
      burst_any=.false.
      do ib=1,nburst
         if (burst_in(ib)) burst_any=.true.
      end do
   end function burst_any

   !> Write one burst frame (PLIC with ccl_id, spray, frames.txt) for every window containing the current time.
   !> Polygons, labels and surface variables must already be current (done in the viz section).
   subroutine write_bursts()
      use messager, only: log
      use string,   only: str_long
      implicit none
      character(len=str_long) :: message
      integer :: ib,n,iunit
      real(WP) :: dtf
      do ib=1,nburst
         if (.not.burst_in(ib)) cycle
         ! Open the writers when the window is first reached (on restart this resumes the existing series)
         if (.not.burst_open(ib)) then
            call burst_viz(ib)%initialize(amr,trim(burst_name(ib)),use_hdf5=.false.)
            call burst_viz(ib)%add_surfmesh(fs%smesh,'plic')
            if (do_transfer) then
               call burst_lpviz(ib)%initialize(lp,trim(burst_name(ib)))
               call burst_lpviz(ib)%select_comp('d'   ,.true.)
               call burst_lpviz(ib)%select_comp('vel' ,.true.)
               call burst_lpviz(ib)%select_comp('flag',.true.)
               call burst_lpviz(ib)%select_comp('mass',.true.)
               call burst_lpviz(ib)%select_comp('eint',.true.)
            end if
            burst_open(ib)=.true.
            write(message,'("[burst] opened amrviz/",a," at t=",es12.5)') trim(burst_name(ib)),time%t; call log(message)
         end if
         ! Write the frame
         call burst_viz(ib)%write(time=time%t)
         if (do_transfer) call burst_lpviz(ib)%write(time=time%t)
         ! Rewrite the frame log from the writer's own time series (consistent with plic.pvd, also after a restart)
         if (amr%amRoot) then
            open(newunit=iunit,file='amrviz/'//trim(burst_name(ib))//'/frames.txt',status='replace',action='write')
            write(iunit,'("# Burst tau=",es24.16,"  window t=[",es24.16,",",es24.16,"]  tc_sim=",es24.16,"  t_hit=",es24.16)') &
            &  burst_tau(ib),burst_t0(ib),burst_t1(ib),tc_sim,t_hit
            write(iunit,'("#",a7,3(1x,a24))') 'frame','t','tau','dt_prev'
            do n=1,burst_viz(ib)%ntime
               dtf=0.0_WP; if (n.gt.1) dtf=burst_viz(ib)%time(n)-burst_viz(ib)%time(n-1)
               write(iunit,'(i8,3(1x,es24.16))') n,burst_viz(ib)%time(n),(burst_viz(ib)%time(n)-t_hit)/tc_sim,dtf
            end do
            close(iunit)
         end if
      end do
   end subroutine write_bursts

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read,param_exists,param_getsize
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
         ! Read in max level for fluid
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
         call param_read('Sponge length',L_spg  ,default=2.0_WP)
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
         ! Midpoint sub-iterations (2 is explicit midpoint RK2, unstable on the imaginary axis; use 3 or more)
         call param_read('Subiterations',time%itmax,default=2)
         if (time%itmax.lt.1) call die('Require at least one sub-iteration')
         call param_read('Max wall time',time%wtmax,default=huge(1.0_WP))
         time%dt=time%dtmax
         if (restarted) then
            call io%get_scalar('dt',time%dt)
            time%t=restart_time
            time%n=restart_step
         end if
      end block initialize_timetracker

      ! Initialize compressible multiphase solver
      create_solver: block
         use amrex_amr_module, only: amrex_bc_ext_dir,amrex_bc_foextrap
         use amrmpcomp_class,  only: BC_GAS
         ! Assign materials and create flow solver (initialize forces face-linear velocity interpolation)
         fs%liq=>water; fs%gas=>gas; call fs%initialize(amr=amr,name='drop')
         ! Surface tension is off (fs%sigma defaults to 0)
         ! Provide pressure relaxation model
         call relax_model%initialize(gas=gas,liq=water); fs%relax=>relax_model
         relax_model%model=Prelax
         ! Refuse relaxations that leave the gas below 1e-2 (torrens); the model's own RHOGmin keeps its 1e-2 default
         fs%rhoG_min=1.0e-2_WP
         ! Positivity gate temperature floors (0: undeclared, a pure cell outside its law's domain stops the run)
         call param_read('Liquid Tmin',fs%TminL,default=0.0_WP)
         call param_read('Gas Tmin'   ,fs%TminG,default=0.0_WP)
         ! relax_model%vol=amr%cell_vol(amr%maxlvl)
         ! Small-cell agglomeration and clean_Q rescue limits are not available in the staggered solver
         ! fs%merge_sick=100.0_WP
         ! relax_model%diss_P=200.0_WP
         ! fs%Pmin_liq=-0.98_WP*water%pinf
         ! fs%Tmin_liq=0.1_WP
         ! fs%Pmin_gas=1.0e-4_WP
         ! fs%Tmin_gas=0.1_WP
         ! relax_model%Pmin_liq=fs%Pmin_liq; relax_model%Tmin_liq=fs%Tmin_liq
         ! relax_model%Pmin_gas=fs%Pmin_gas; relax_model%Tmin_gas=fs%Tmin_gas
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

      ! Drop transfer objects register regrid callbacks, so they must exist before the grid is built
      create_transfer: block
         use mathtools, only: Pi
         use amrdata_class, only: interp_none,interp_const
         use messager, only: die,log
         call param_read('Drop transfer',do_transfer,default=.false.)
         ! CCL output is mandatory, including while transfer is disabled for debugging.
         call ccl%initialize(amr=amr,name='ccl')
         call fs%add_surfvar(ccl%id,1,'ccl_id')
         ! Particle setup only when transfer is on, so amrlpt's cost callback is not installed otherwise
         if (do_transfer) then
            ! Set density and drag model (Osnes et al. 2023) before initialize(), which prints them
            lp%rho=rhoL1
            call param_read('Particle drag model',lp%drag_model,default='Osnes')
            ! This case uses isolated-sphere drag; retain deposited VF for diagnostics and coupling.
            lp%drag_vf_override=0.0_WP
            lp%total_energy_coupling=.true.
            lp%cull_gas_free=.true.
            lp%heat_material=>water
            lp%heat_prandtl=Prandtl
            if (.not.particle_finite(Prandtl).or.Prandtl.le.0.0_WP) call die('[particle heat] Prandtl must be positive')
            ! Lagrangian spray
            call lp%initialize(amr,name='spray',gas_gamma=gas%gamma)
            if (lp%nover.gt.min(fs%RHOG%ng,fs%PG%ng)) call die('[particle gas properties] insufficient primitive ghost width')
            ! Source and normal cell widths for conservative averaging across refinement boundaries.
            call spray_cells%initialize(amr,name='spray_cell_projection',ncomp=6,ng=1,interp=interp_const)
            call spray_cells%register()
            spray_cells%lo_bc(:,1:3)=lp%src%lo_bc(:,1:3); spray_cells%hi_bc(:,1:3)=lp%src%hi_bc(:,1:3)
            spray_cells%lo_bc(:,4:6)=lp%src%lo_bc(:,1:3); spray_cells%hi_bc(:,4:6)=lp%src%hi_bc(:,1:3)
            call spray_energy%initialize(amr,name='spray_energy_rate',ncomp=1,ng=0,interp=interp_none)
            call spray_energy%register()
            call spray_force(1)%initialize(amr,name='spray_xforce_rate',ncomp=1,ng=0,nodal=[.true.,.false.,.false.],interp=interp_none)
            call spray_force(2)%initialize(amr,name='spray_yforce_rate',ncomp=1,ng=0,nodal=[.false.,.true.,.false.],interp=interp_none)
            call spray_force(3)%initialize(amr,name='spray_zforce_rate',ncomp=1,ng=0,nodal=[.false.,.false.,.true.],interp=interp_none)
            call spray_force(1)%register(); call spray_force(2)%register(); call spray_force(3)%register()
            call particle_gas%initialize(amr,name='particle_gas',ncomp=2,ng=lp%nover)
            call particle_gas%register()
            call particle_thermal%initialize(amr,name='particle_thermal',ncomp=2,ng=lp%nover)
            call particle_thermal%register()
            call particle_gas_fraction%initialize(amr,name='particle_gas_fraction',ncomp=1,ng=lp%nover)
            call particle_gas_fraction%register()
            call particle_stress%initialize(amr,name='particle_stress',ncomp=3,ng=lp%nover)
            call particle_stress%register()
            particle_stress%lo_bc=lp%src%lo_bc(:,1:3)
            particle_stress%hi_bc=lp%src%hi_bc(:,1:3)
            ! Cap particle sampling and deposition at maxlvl-offset
            call param_read('Particle level offset',lp_lvl_off,default=3)
            lp_lvl_off=max(0,min(lp_lvl_off,amr%maxlvl))
            ! Particle level ceiling, independent of the fluid maximum
            lp%maxlvl=amr%maxlvl-lp_lvl_off
            call param_read('Particle max CFL'   ,lp%cflmax    ,default=0.5_WP)
            ! Optional trace of the particles with maximum Re, Ma and Kn (diagnostic restarts only)
            call param_read('Particle drag diagnostics',lp%drag_diagnostics,default=.false.)
            ! Diameter limits: maximum uses the coupling cap; minimum always uses the finest grid.
            if (param_exists('Drop transfer volume coeff')) &
            &   call die('[transfer] replace Drop transfer volume coeff with Drop transfer maximum diameter ratio (coupling-cell units)')
            call param_read('Drop transfer maximum diameter ratio',xfer_dmax_ratio,default=1.0_WP)
            call param_read('Drop transfer diameter threshold',xfer_diam_threshold,default=0.1_WP)
            if (xfer_dmax_ratio.le.0.0_WP.or.xfer_diam_threshold.lt.0.0_WP) &
            &   call die('[transfer] maximum diameter ratio must be positive; minimum threshold must be non-negative')
            lp_dx=amr%min_meshsize(lp%maxlvl)
            xfer_dmax=xfer_dmax_ratio*lp_dx
            if (transfer_min_diameter().gt.xfer_dmax) call die('[transfer] minimum diameter exceeds maximum diameter')
            lp_dres=xfer_dmax/lp_dx
            call log('[transfer] minimum diameter uses finest-grid spacing; maximum uses coupling-cap spacing')
            transfer_limits: block
               character(len=160) :: limits
               write(limits,'("[transfer] dmin=",es14.6," dmax=",es14.6)') transfer_min_diameter(),xfer_dmax
               call log(trim(limits))
               write(limits,'("[transfer] particle maxlvl=",i0," coupling dx=",es14.6," maximum d/dx=",es14.6)') lp%maxlvl,lp_dx,lp_dres
               call log(trim(limits))
            end block transfer_limits
            call param_read('Particle maximum d/dx',lp_dres_limit,default=1.0_WP)
            if (lp_dres_limit.le.0.0_WP) call die('[amrcomp_drop_lpt] Particle maximum d/dx must be positive')
            if (lp_dres.gt.lp_dres_limit*(1.0_WP+100.0_WP*epsilon(1.0_WP))) &
            &   call die('[amrcomp_drop_lpt] transfer threshold admits particles larger than the configured particle-cell limit; increase Particle level offset or Particle maximum d/dx')
            ! Source filter width, off by default
            call param_read('Particle filter width',lp%filter_width,default=0.0_WP)
            ! filter_ncell=3 puts a 3.5*dx filter on the particle level
            lp%filter_ncell=3
            call param_read('Particle max dt'    ,lp%dtmax     ,default=time%dtmax)
            call param_read('Particle collisions',lp_collide   ,default=.false.)
            if (lp_collide) call param_read('Particle collision time',lp%tau_col,default=10.0_WP*time%dtmax)
            lp%dt=lp%dtmax
            call param_read('Particle max substeps',lp%nsub_max,default=1000)

            ! Core-box padding for wake tagging only; <0 selects 0.5*d_drop
            call param_read('Drop transfer margin'      ,xfer_margin,default=-1.0_WP)
            ! Particles leave through the open x/y boundaries; z is the periodic quasi-2D direction
            lp%lo_bc=AMRLPT_OPEN; lp%hi_bc=AMRLPT_OPEN
            ! Run the transfer every N steps
            xfer_evt=event(time=time,name='Drop transfer')
            call param_read('Drop transfer nsteps',xfer_evt%nper,default=1)
         end if
      end block create_transfer

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
         ! Graded far-field coarsening beyond the sponge ramp
         call param_read('Minimum level',min_lvl,default=2)
         min_lvl=max(0,min(min_lvl,amr%maxlvl))
         call param_read('Sponge min level' ,spg_minlvl   ,default=2)
         spg_minlvl=max(min_lvl,min(spg_minlvl,amr%maxlvl))
         call build_sponge_bands()
         ! Wake coarsening (see the declarations); defaults leave the tagger unchanged
         call param_read('Wake tag level' ,wake_lvl ,default=huge(1))
         call param_read('Wake tag ramp'  ,wake_ramp,default=1.0_WP)
         call param_read('Wake max level offset',wake_maxlvl_off,default=0)
         wake_maxlvl_off=max(0,min(wake_maxlvl_off,amr%maxlvl))
         call param_read('Wake tag x0'    ,wake_x0  ,default=huge(1.0_WP))
         call param_read('Wake tag length',wake_len ,default=1.0_WP)
         call param_read('Wake tag density',wake_rho ,default=.false.)
         ! Build the grid
         if (restarted) then
            ! Restore grid hierarchy from checkpoint
            call amr%init_from_checkpoint(dirname=trim(restart_dir),time=time%t)
            ! Restore solver state (rebuilds primitives and face velocities from the restored momentum)
            call fs%restore_checkpoint(io=io,dirname=trim(restart_dir),time=time%t)
            ! Fire the post-regrid callbacks that init_from_checkpoint skips, so amrlpt sees the restored grids
            call amr%fire_postregrid(time=time%t)
         else
            ! Fresh start: post_regrid builds momentum from the face velocity and the primitive variables
            call amr%init_from_scratch(time=time%t)
            ! Phasic total energies from the internal contents written by user_init and the initial kinetic energy
            call fs%get_total_energy(time=time%t)
            ! Build PLIC (no build_subVF: the staggered solver does not read the sub-cell volume fractions)
            call fs%build_plic(time%t)
         end if
         ! Compute viscosities
         call get_viscosities()
         ! Add SGS models
         call fs%add_viscartif(dt=time%dt,Cvisc=2.5e-2_WP)
         call fs%add_vreman(dt=time%dt)
         ! Compute cell-centered velocity, Umag and Mach number
         call fs%get_cell_velocity()
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
         if (do_transfer) then
            ! Anchor the particle clock to the flow clock (matters on restart)
            lp%t=time%t; lp_tprev=time%t
         end if
         ! Restore the spray and the transfer ledger from the checkpoint
         if (restarted.and.do_transfer) then
            ! read() enforces the particle level cap and rebuilds lp%VF
            call lp%read(dirname=trim(restart_dir)//'_spray')
            call lp%get_info(); lp_np=lp%np
            restore_ledger: block
               real(WP) :: val
               logical  :: found
               call io%get_scalar('xfer_vol'   ,val,found); if (found) xfer_vol   =val
               call io%get_scalar('xfer_mass'  ,val,found); if (found) xfer_mass  =val
               call io%get_scalar('xfer_momx'  ,val,found); if (found) xfer_momx  =val
               call io%get_scalar('xfer_energy',val,found); if (found) xfer_energy=val
               call io%get_scalar('xfer_gasadd',val,found); if (found) xfer_gasadd=val
               call io%get_scalar('refill_energy',val,found); if (found) xfer_refill(4)=val
               call io%get_scalar('refill_momz',val,found); if (found) xfer_refill(3)=val
               call io%get_scalar('refill_momy',val,found); if (found) xfer_refill(2)=val
               call io%get_scalar('refill_momx',val,found); if (found) xfer_refill(1)=val
               call io%get_scalar('xfer_refused',val,found); if (found) xfer_refused=val
               call io%get_scalar('xfer_badthermal',val,found); if (found) xfer_badthermal=val
               call io%get_scalar('cull_count',val,found); if (found) lp%np_culled_tot=val
               call io%get_scalar('cull_mass',val,found); if (found) lp%culled(1)=val
               call io%get_scalar('cull_momx',val,found); if (found) lp%culled(2)=val
               call io%get_scalar('cull_momy',val,found); if (found) lp%culled(3)=val
               call io%get_scalar('cull_momz',val,found); if (found) lp%culled(4)=val
               call io%get_scalar('cull_energy',val,found); if (found) lp%culled(5)=val
               call io%get_scalar('exchange_heat',val,found); if (found) lp%heat_exchange=val
               call io%get_scalar('exchange_px',val,found); if (found) lp%exchange(1)=val
               call io%get_scalar('exchange_py',val,found); if (found) lp%exchange(2)=val
               call io%get_scalar('exchange_pz',val,found); if (found) lp%exchange(3)=val
               call io%get_scalar('exchange_mech',val,found); if (found) lp%exchange(4)=val
               call io%get_scalar('drag_dissipation',val,found); if (found) lp%exchange(5)=val
               call io%get_scalar('out_mass',val,found); if (found) lp%outflow(1)=val
               call io%get_scalar('out_momx',val,found); if (found) lp%outflow(2)=val
               call io%get_scalar('out_momy',val,found); if (found) lp%outflow(3)=val
               call io%get_scalar('out_momz',val,found); if (found) lp%outflow(4)=val
               call io%get_scalar('out_energy',val,found); if (found) lp%outflow(5)=val
               call io%get_scalar('flot_n_total',val,found); if (found) flot_n_total=val
               call io%get_scalar('flot_vol'    ,val,found); if (found) flot_vol    =val
               call io%get_scalar('xfer_momy',val,found); if (found) xfer_momy=val
               call io%get_scalar('xfer_momz',val,found); if (found) xfer_momz=val
               call io%get_scalar('flot_momx',val,found); if (found) flot_mom(1)=val
               call io%get_scalar('flot_momy',val,found); if (found) flot_mom(2)=val
               call io%get_scalar('flot_momz',val,found); if (found) flot_mom(3)=val
               call io%get_scalar('flot_energy',val,found); if (found) flot_energy=val
               call io%get_scalar('flot_mass'   ,val,found); if (found) flot_mass   =val
            end block restore_ledger
         end if
         ! Restore the regrid ledger from the checkpoint
         if (restarted) then
            restore_regrid: block
               real(WP) :: val
               logical  :: found
               call io%get_scalar('regrid_massL' ,val,found); if (found) regrid_book(1)=val
               call io%get_scalar('regrid_massG' ,val,found); if (found) regrid_book(2)=val
               call io%get_scalar('regrid_momx'  ,val,found); if (found) regrid_book(3)=val
               call io%get_scalar('regrid_momy'  ,val,found); if (found) regrid_book(4)=val
               call io%get_scalar('regrid_momz'  ,val,found); if (found) regrid_book(5)=val
               call io%get_scalar('regrid_energy',val,found); if (found) regrid_book(6)=val
            end block restore_regrid
         end if
      end block init_regridding

      ! Initialize checkpoint save event
      init_checkpoint: block
         ! Create checkpoint save event
         save_evt=event(time=time,name='Checkpoint')
         call param_read('Checkpoint period',save_evt%tper,default=-1.0_WP)
         ! Let solver self-register for checkpointing
         call fs%register_checkpoint(io)
         ! Checkpoint the transfer ledger (particles go to a sibling directory)
         if (do_transfer) then
            call io%add_scalar(name='xfer_vol'   ,value=xfer_vol)
            call io%add_scalar(name='xfer_mass'  ,value=xfer_mass)
            call io%add_scalar(name='xfer_momx'  ,value=xfer_momx)
            call io%add_scalar(name='xfer_energy',value=xfer_energy)
            call io%add_scalar(name='xfer_gasadd',value=xfer_gasadd)
            call io%add_scalar(name='refill_energy',value=xfer_refill(4))
            call io%add_scalar(name='refill_momz',value=xfer_refill(3))
            call io%add_scalar(name='refill_momy',value=xfer_refill(2))
            call io%add_scalar(name='refill_momx',value=xfer_refill(1))
            call io%add_scalar(name='xfer_refused',value=xfer_refused)
            call io%add_scalar(name='xfer_badthermal',value=xfer_badthermal)
            call io%add_scalar(name='cull_count',value=lp%np_culled_tot)
            call io%add_scalar(name='cull_mass',value=lp%culled(1))
            call io%add_scalar(name='cull_momx',value=lp%culled(2))
            call io%add_scalar(name='cull_momy',value=lp%culled(3))
            call io%add_scalar(name='cull_momz',value=lp%culled(4))
            call io%add_scalar(name='cull_energy',value=lp%culled(5))
            call io%add_scalar(name='exchange_heat',value=lp%heat_exchange)
            call io%add_scalar(name='exchange_px',value=lp%exchange(1))
            call io%add_scalar(name='exchange_py',value=lp%exchange(2))
            call io%add_scalar(name='exchange_pz',value=lp%exchange(3))
            call io%add_scalar(name='exchange_mech',value=lp%exchange(4))
            call io%add_scalar(name='drag_dissipation',value=lp%exchange(5))
            call io%add_scalar(name='out_mass',value=lp%outflow(1))
            call io%add_scalar(name='out_momx',value=lp%outflow(2))
            call io%add_scalar(name='out_momy',value=lp%outflow(3))
            call io%add_scalar(name='out_momz',value=lp%outflow(4))
            call io%add_scalar(name='out_energy',value=lp%outflow(5))
            call io%add_scalar(name='flot_n_total',value=flot_n_total)
            call io%add_scalar(name='flot_vol'    ,value=flot_vol)
            call io%add_scalar(name='xfer_momy',value=xfer_momy)
            call io%add_scalar(name='xfer_momz',value=xfer_momz)
            call io%add_scalar(name='flot_momx',value=flot_mom(1))
            call io%add_scalar(name='flot_momy',value=flot_mom(2))
            call io%add_scalar(name='flot_momz',value=flot_mom(3))
            call io%add_scalar(name='flot_energy',value=flot_energy)
            call io%add_scalar(name='flot_mass'   ,value=flot_mass)
         end if
         ! Add dt to checkpoint save
         call io%add_scalar(name='dt',value=time%dt)
         ! Checkpoint the regrid ledger
         call io%add_scalar(name='regrid_massL' ,value=regrid_book(1))
         call io%add_scalar(name='regrid_massG' ,value=regrid_book(2))
         call io%add_scalar(name='regrid_momx'  ,value=regrid_book(3))
         call io%add_scalar(name='regrid_momy'  ,value=regrid_book(4))
         call io%add_scalar(name='regrid_momz'  ,value=regrid_book(5))
         call io%add_scalar(name='regrid_energy',value=regrid_book(6))
      end block init_checkpoint

      ! Per-phase timers, created before the initial labelling that uses some of them
      create_timers: block
         tmr_step_all =timer(comm=amr%comm,name='step'     ,sync=.false.)
         tmr_prep     =timer(comm=amr%comm,name='prep'     ,sync=.false.)
         tmr_lpt      =timer(comm=amr%comm,name='lpt'      ,sync=.false.)
         tmr_lptinfo  =timer(comm=amr%comm,name='lptinfo'  ,sync=.false.)
         tmr_divstress=timer(comm=amr%comm,name='divstress',sync=.false.)
         tmr_loop     =timer(comm=amr%comm,name='loop'     ,sync=.false.)
         tmr_rk       =timer(comm=amr%comm,name='rk'       ,sync=.false.)
         tmr_xfer     =timer(comm=amr%comm,name='xfer'     ,sync=.false.)
         tmr_ccl      =timer(comm=amr%comm,name='ccl'      ,sync=.false.)
         tmr_msweep   =timer(comm=amr%comm,name='msweep'   ,sync=.false.)
         tmr_mcomm    =timer(comm=amr%comm,name='mcomm'    ,sync=.false.)
         tmr_xpart    =timer(comm=amr%comm,name='xpart'    ,sync=.false.)
         tmr_surg     =timer(comm=amr%comm,name='surgery'  ,sync=.false.)
         tmr_resync   =timer(comm=amr%comm,name='resync'   ,sync=.false.)
         tmr_diag     =timer(comm=amr%comm,name='diag'     ,sync=.false.)
         tmr_regrid   =timer(comm=amr%comm,name='regrid'   ,sync=.false.)
         tmr_props    =timer(comm=amr%comm,name='props'    ,sync=.false.)
         tmr_sgs      =timer(comm=amr%comm,name='sgs'      ,sync=.false.)
         tmr_viz      =timer(comm=amr%comm,name='viz'      ,sync=.false.)
         tmr_ckpt     =timer(comm=amr%comm,name='ckpt'     ,sync=.false.)
         tmr_mon      =timer(comm=amr%comm,name='mon'      ,sync=.false.)
      end block create_timers

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
         ! Characteristic time and tau origin, from the flow parameters (see burst_implementation.md)
         tc_sim=sqrt(rhoL1/rhoG2)*d_drop/u2
         t_hit=(-0.5_WP*d_drop-Xs)/Ws
         tau_period: block
            real(WP) :: tau_per
            call param_read('tau Output period [tau]',tau_per)
            tau_viz_evt%tper=tau_per*tc_sim
         end block tau_period
         ! Burst windows: tau centres and camera exposure converted to simulation time
         init_bursts: block
            use messager, only: log,die
            use string,   only: str_long,rtoa
            character(len=str_long) :: message
            real(WP) :: exposure,tc_exp,L_tau,L_sim
            integer :: ib
            write(message,'("[tau] tc_sim=",es12.5," t_hit=",es12.5," (tau=0: shock at drop leading edge)")') tc_sim,t_hit; call log(message)
            nburst=0
            if (param_exists('Burst tau')) nburst=param_getsize('Burst tau')
            allocate(burst_tau(nburst),burst_t0(nburst),burst_t1(nburst),burst_open(nburst),burst_viz(nburst),burst_lpviz(nburst))
            burst_open=.false.
            if (nburst.gt.0) then
               call param_read('Burst tau',burst_tau)
               call param_read('Burst exposure [s]',exposure)
               call param_read('Experiment t_c [s]',tc_exp)
               if (exposure.le.0.0_WP) call die('[burst] Burst exposure [s] must be positive')
               if (tc_exp.le.0.0_WP) call die('[burst] Experiment t_c [s] must be positive')
               L_tau=exposure/tc_exp
               L_sim=L_tau*tc_sim
               burst_t0=t_hit+burst_tau*tc_sim-0.5_WP*L_sim
               burst_t1=t_hit+burst_tau*tc_sim+0.5_WP*L_sim
               write(message,'("[burst] exposure=",es12.5," s, tc_exp=",es12.5," s -> L_tau=",es12.5," L_sim=",es12.5)') exposure,tc_exp,L_tau,L_sim; call log(message)
               if (L_sim.lt.3.0_WP*time%dtmax) then
                  write(message,'("[burst] WARNING: window L_sim=",es12.5," is shorter than 3*Max dt=",es12.5,"; few frames per burst")') L_sim,3.0_WP*time%dtmax; call log(message)
               end if
               do ib=1,nburst
                  if (burst_t0(ib).lt.0.0_WP) call die('[burst] window of Burst tau='//trim(adjustl(rtoa(burst_tau(ib))))//' starts before t=0')
                  write(message,'("[burst] tau=",es12.5," t=[",es12.5,",",es12.5,"] tau=[",es12.5,",",es12.5,"] -> amrviz/",a)') &
                  &  burst_tau(ib),burst_t0(ib),burst_t1(ib),(burst_t0(ib)-t_hit)/tc_sim,(burst_t1(ib)-t_hit)/tc_sim,trim(burst_name(ib)); call log(message)
                  if (burst_t1(ib).gt.time%tmax) then
                     write(message,'("[burst] WARNING: window of tau=",es12.5," ends after Max time=",es12.5)') burst_tau(ib),time%tmax; call log(message)
                  end if
               end do
            end if
         end block init_bursts
         ! Structure labels in every series
         call viz%add_scalar(ccl%id,1,'ccl_id')
         call plicviz%add_scalar(ccl%id,1,'ccl_id')
         call tau_viz%add_scalar(ccl%id,1,'ccl_id')
         if (do_transfer) then
            ! Coupling source: src(1:3) is the fluid momentum gain per volume over a step, src(4) the drag dissipation
            call viz%add_scalar(lp%src,1,'src_momx')
            call viz%add_scalar(lp%src,2,'src_momy')
            call viz%add_scalar(lp%src,4,'src_energy')
            call viz%add_scalar(lp%VF,1,'part_VF')
            ! Spray particles get their own plotfile series
            call lpviz%initialize(lp,'spray')
            call lpviz%select_comp('d'   ,.true.)
            call lpviz%select_comp('vel' ,.true.)
            call lpviz%select_comp('flag',.true.)
            call lpviz%select_comp('mass',.true.)
            call lpviz%select_comp('eint',.true.)
            lpviz_evt=event(time=time,name='Spray visualization output')
            call param_read('Spray output period',lpviz_evt%tper,default=viz_evt%tper)
         end if
         ! Label the liquid and compute the QOIs (no transfer), then write the initial state
         initial_core: block
            integer,  dimension(:),   allocatable :: cid,cmap
            real(WP), dimension(:,:), allocatable :: mom
            integer :: ns
            call label_and_moments(cid,mom,cmap,ns)
            deallocate(cid,mom,cmap)
            call compute_drop_metrics()
         end block initial_core
         if (plicviz_evt%occurs().or.tau_viz_evt%occurs().or.burst_any()) call fs%build_polygons()
         if (viz_evt%occurs().or.plicviz_evt%occurs().or.tau_viz_evt%occurs().or.burst_any()) then
            call fs%update_surfvars()
         end if
         if (viz_evt%occurs()) call viz%write(time=time%t)
         if (plicviz_evt%occurs()) call plicviz%write(time=time%t)
         if (tau_viz_evt%occurs()) call tau_viz%write(time=time%t)
         call write_bursts()
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
         call mfile%add_column(u2,'u2_exact')
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
         call consfile%add_column(regrid_book(1),'Regrid liquid mass')
         call consfile%add_column(regrid_book(2),'Regrid gas mass')
         call consfile%add_column(regrid_book(3),'Regrid U momentum')
         call consfile%add_column(regrid_book(4),'Regrid V momentum')
         call consfile%add_column(regrid_book(5),'Regrid W momentum')
         call consfile%add_column(regrid_book(6),'Regrid energy')
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
         ! Create timing monitor: solver routine times, slowest rank (max) and fastest rank (min)
         tfile=monitor(amRoot=amr%amRoot,name='timing')
         call tfile%add_column(time%n,'Timestep')
         call tfile%add_column(time%t,'Time')
         ! Full routine times (max across ranks = wall-clock cost)
         call tfile%add_column(fs%wtmax_dQdt,'dQdt_max')
         call tfile%add_column(fs%wtmax_plic,'plic_max')
         call tfile%add_column(fs%wtmax_relax,'relax_max')
         call tfile%add_column(fs%wtmax_visc,'visc_max')
         call tfile%add_column(fs%wtmax_mid,'mid_max')
         call tfile%add_column(fs%wtmax_fvel,'fvel_max')
         call tfile%add_column(fs%wtmin_fvel,'fvel_min')
         ! Compute loop times (max = slowest rank, min = fastest rank)
         call tfile%add_column(fs%wtmax_prim,'prim_max')
         call tfile%add_column(fs%wtmin_prim,'prim_min')
         call tfile%add_column(fs%wtmax_alloc,'alloc_max')
         call tfile%add_column(fs%wtmax_slflux,'slflux_max')
         call tfile%add_column(fs%wtmin_slflux,'slflux_min')
         call tfile%add_column(fs%wtmax_vmom,'vmom_max')
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
         ! Create rescue monitor (solver counters of the last step)
         rescfile=monitor(amRoot=amr%amRoot,name='rescue')
         call rescfile%add_column(time%n,'Timestep')
         call rescfile%add_column(time%t,'Time')
         ! The clean_Q rescue and small-cell pooling of the collocated solver do not exist here
         call rescfile%add_column(fs%nrelax_fail,'Relax fail')
         call rescfile%add_column(fs%nvf_guard,'VF guard cells')
         call rescfile%add_column(fs%nsliver,'Sliver cells')
         call rescfile%add_column(fs%ngate,'Gate cells')
         call rescfile%add_column(fs%ngate_floor,'Gate floors')
         call rescfile%add_column(fs%gate_mass,'Gate mass')
         call rescfile%add_column(fs%gate_energy,'Gate energy')
         call rescfile%add_column(fs%gate_rhoU,'Gate U momentum')
         call rescfile%add_column(fs%gate_rhoV,'Gate V momentum')
         call rescfile%add_column(fs%gate_rhoW,'Gate W momentum')
         call rescfile%write()
         ! Droplet QOI monitor, evaluated on the transfer interval by default so it reuses the transfer's labels
         metrics_evt=event(time=time,name='Droplet QOIs')
         if (do_transfer) then
            call param_read('Drop QOI nsteps',metrics_evt%nper,default=xfer_evt%nper)
         else
            call param_read('Drop QOI nsteps',metrics_evt%nper,default=10)
         end if
         drop_QOIs=monitor(amRoot=amr%amRoot,name='drop_QOIs')
         call drop_QOIs%add_column(time%n,'Timestep')
         call drop_QOIs%add_column(time%t,'time')
         ! Isocontour family: every cell in the domain above the threshold, fragments included
         call drop_QOIs%add_column(MLE_VF01,'MLE_VF01')
         call drop_QOIs%add_column(MLE_VF05,'MLE_VF05')
         call drop_QOIs%add_column(drop_massVF01,'massVF01')
         call drop_QOIs%add_column(MOI_xx,'Ixx')
         call drop_QOIs%add_column(MOI_yy,'Iyy')
         call drop_QOIs%add_column(MOI_zz,'Izz')
         call drop_QOIs%add_column(MOI_xy,'Ixy')
         call drop_QOIs%add_column(MOI_xz,'Ixz')
         call drop_QOIs%add_column(MOI_yz,'Iyz')
         ! Core family: largest liquid-volume structure, every labelled cell, no VF cut
         call drop_QOIs%add_column(MLE_core,'MLE_core')
         call drop_QOIs%add_column(MLE_core001,'MLE_core001')
         call drop_QOIs%add_column(MLE_core01,'MLE_core01')
         call drop_QOIs%add_column(MLE_core05,'MLE_core05')
         call drop_QOIs%add_column(core_mass,'core_mass')
         call drop_QOIs%add_column(core_vol,'core_vol')
         call drop_QOIs%add_column(MOI_xx_core,'Ixx_core')
         call drop_QOIs%add_column(MOI_yy_core,'Iyy_core')
         call drop_QOIs%add_column(MOI_zz_core,'Izz_core')
         call drop_QOIs%add_column(MOI_xy_core,'Ixy_core')
         call drop_QOIs%add_column(MOI_xz_core,'Ixz_core')
         call drop_QOIs%add_column(MOI_yz_core,'Iyz_core')
         ! Core health flags: measured, incomplete label, straddles a periodic boundary
         call drop_QOIs%add_column(core_valid,'core_valid')
         call drop_QOIs%add_column(core_open,'core_open')
         call drop_QOIs%add_column(core_wrap,'core_wrap')
         call drop_QOIs%add_column(ccl_nstruct,'core_nstruct')
         call drop_QOIs%write()
         ! Spray monitor: transfer ledger, flotsam, particle statistics and coupling diagnostics
         if (do_transfer) then
            sprayfile=monitor(amRoot=amr%amRoot,name='spray')
            call sprayfile%add_column(time%n,'Timestep number')
            call sprayfile%add_column(time%t,'Time')
            call sprayfile%add_column(ccl_nstruct,'Nstruct')
            call sprayfile%add_column(xfer_n_step,'Nxfer_step')
            call sprayfile%add_column(flot_n_step,'Nflot_step')
            call sprayfile%add_column(flot_n_total,'Nflot total')
            call sprayfile%add_column(lp_np,'Np')
            call sprayfile%add_column(xfer_vol,'Xfer volume')
            call sprayfile%add_column(xfer_mass,'Xfer mass')
            call sprayfile%add_column(flot_vol,'Flot volume')
            call sprayfile%add_column(xfer_momy,'Xfer momy')
            call sprayfile%add_column(xfer_momz,'Xfer momz')
            call sprayfile%add_column(flot_mom(1),'Flot momx')
            call sprayfile%add_column(flot_mom(2),'Flot momy')
            call sprayfile%add_column(flot_mom(3),'Flot momz')
            call sprayfile%add_column(flot_energy,'Flot energy')
            call sprayfile%add_column(flot_mass,'Flot mass')
            call sprayfile%add_column(xfer_momx,'Xfer momx')
            call sprayfile%add_column(xfer_energy,'Xfer energy')
            call sprayfile%add_column(xfer_gasadd,'Gas added')
            call sprayfile%add_column(xfer_refill(1),'Refill momx')
            call sprayfile%add_column(xfer_refill(2),'Refill momy')
            call sprayfile%add_column(xfer_refill(3),'Refill momz')
            call sprayfile%add_column(xfer_refill(4),'Refill energy')
            call sprayfile%add_column(xfer_refused,'Refused total')
            call sprayfile%add_column(lp%occupied_minlvl,'lp lvl min')
            call sprayfile%add_column(lp%occupied_maxlvl,'lp lvl max')
            call sprayfile%add_column(lp%np_above_maxlvl,'Np above cap')
            call sprayfile%add_column(lp%np_demoted,'Np demoted')
            call sprayfile%add_column(lp%dmin,'dmin')
            call sprayfile%add_column(lp%dmax,'dmax')
            call sprayfile%add_column(lp%dmean,'dmean')
            call sprayfile%add_column(lp%Umean,'Umean')
            call sprayfile%add_column(lp%Vp_tot,'Vp total')
            call sprayfile%add_column(lp%mass_tot,'Particle mass')
            call sprayfile%add_column(lp%momentum_tot(1),'Particle momx')
            call sprayfile%add_column(lp%momentum_tot(2),'Particle momy')
            call sprayfile%add_column(lp%momentum_tot(3),'Particle momz')
            call sprayfile%add_column(lp%kinetic_tot,'Particle KE')
            call sprayfile%add_column(lp%internal_tot,'Particle IE')
            call sprayfile%add_column(lp%outflow(1),'Exit mass')
            call sprayfile%add_column(lp%outflow(2),'Exit momx')
            call sprayfile%add_column(lp%outflow(3),'Exit momy')
            call sprayfile%add_column(lp%outflow(4),'Exit momz')
            call sprayfile%add_column(lp%outflow(5),'Exit energy')
            ! Coupling ledger and drag-regime diagnostics
            call sprayfile%add_column(src_momx,'Src momx')
            call sprayfile%add_column(lp%np_culled,'Cull step')
            call sprayfile%add_column(lp%np_culled_tot,'Cull total')
            call sprayfile%add_column(lp%culled(1),'Cull mass')
            call sprayfile%add_column(lp%culled(2),'Cull momx')
            call sprayfile%add_column(lp%culled(3),'Cull momy')
            call sprayfile%add_column(lp%culled(4),'Cull momz')
            call sprayfile%add_column(lp%culled(5),'Cull energy')
            call sprayfile%add_column(lp%heat_exchange,'Fluid heat E')
            call sprayfile%add_column(lp%exchange(1),'Fluid Jx')
            call sprayfile%add_column(lp%exchange(2),'Fluid Jy')
            call sprayfile%add_column(lp%exchange(3),'Fluid Jz')
            call sprayfile%add_column(lp%exchange(4),'Fluid mech E')
            call sprayfile%add_column(lp%exchange(5),'Drag dissip')
            call sprayfile%add_column(src_ener,'Src energy')
            call sprayfile%add_column(src_emin,'Src emin')
            call sprayfile%add_column(lp%Remax,'Re max')
            call sprayfile%add_column(lp%Mamax,'Ma max')
            call sprayfile%add_column(lp%Knmax,'Kn max')
            call sprayfile%add_column(corelo(1),'Core xlo')
            call sprayfile%add_column(corehi(1),'Core xhi')
            call sprayfile%add_column(src_dumax,'Src du max')
            call sprayfile%add_column(lp%np_substep_capped,'Np capped')
            ! Candidate census (see label_and_moments and transfer_drops)
            call sprayfile%add_column(xfer_ncand ,'Ncand')
            call sprayfile%add_column(xfer_rej(1),'Rej coarse')
            call sprayfile%add_column(xfer_rej(2),'Rej volume')
            call sprayfile%add_column(xfer_fin(1),'Fin volume')
            call sprayfile%add_column(xfer_fin(2),'Fin coarse')
            call sprayfile%add_column(xfer_fin(3),'Fin no donor')
            call sprayfile%add_column(xfer_fin(4),'Fin bad thermal')
            call sprayfile%add_column(xfer_badthermal,'Bad thermal total')
            call sprayfile%add_column(core_fallbacks,'Core fallback')
            call sprayfile%write()
         end if
         ! Per-phase profiling: slowest-rank seconds per phase, with 't other' the untimed remainder
         perffile=monitor(amRoot=amr%amRoot,name='performance')
         call perffile%add_column(time%n,'Timestep')
         call perffile%add_column(time%t,'Time')
         call perffile%add_column(t_step,'t step')
         call perffile%add_column(t_prep,'t prep')
         call perffile%add_column(t_lpt,'t lpt')
         call perffile%add_column(t_lptstep,'t lpt step')
         call perffile%add_column(t_lptsrc,'t lpt src')
         call perffile%add_column(t_lptvf,'t lpt VF')
         call perffile%add_column(t_lptinfo,'t lpt info')
         call perffile%add_column(t_divstress,'t divstress')
         call perffile%add_column(t_loop,'t loop')
         call perffile%add_column(t_rk,'t rk')
         call perffile%add_column(t_xfer,'t transfer')
         call perffile%add_column(t_ccl,'t ccl')
         call perffile%add_column(t_msweep,'t mom sweep')
         call perffile%add_column(t_mcomm,'t mom comm')
         call perffile%add_column(t_xpart,'t xfer part')
         call perffile%add_column(t_surg,'t surgery')
         call perffile%add_column(t_resync,'t resync')
         call perffile%add_column(t_diag,'t diag')
         call perffile%add_column(t_regrid,'t regrid')
         call perffile%add_column(t_props,'t props')
         call perffile%add_column(t_sgs,'t sgs')
         call perffile%add_column(t_viz,'t viz')
         call perffile%add_column(t_ckpt,'t ckpt')
         call perffile%add_column(t_mon,'t mon')
         call perffile%add_column(t_other,'t other')
         call perffile%write()
      end block create_monitors

   end subroutine simulation_init

   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      logical :: spray_source_ready,regridded

      ! Perform time integration
      do while (.not.time%done())

         ! Whole-step wall clock, and every phase timer zeroed for this step
         call reset_timers()
         call tmr_step_all%start()
         ! No label survives into this step
         core_mom_valid=.false.

         ! Increment time and remember old state
         call tmr_prep%start()
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         call time%adjust_dt()
         call time%increment()
         call fs%store_old()
         call tmr_prep%stop()

         ! Advance the spray first; its source is frozen over the step and applied in every sub-iteration
         spray_source_ready=.false.
         t_lptstep=0.0_WP; t_lptsrc=0.0_WP; t_lptvf=0.0_WP
         if (do_transfer) then
            call tmr_lpt%start()
            dt_src=max(time%t-lp_tprev,epsilon(1.0_WP))
            if (lp%np.gt.0_I8) then
               call get_particle_gas_properties()
               ! As in shock_droplet, sample effective resolved stress on the coupling mesh.
               call tmr_divstress%start()
               call get_div_stress(particle_stress,time=lp_tprev,maxlvl=lp%maxlvl)
               call tmr_divstress%stop()
               ! Sample the face velocities directly (amrlpt interpolates face-centered fields)
               call lp%advance_to(time=time%t,do_collide=lp_collide, &
               &                  U=fs%U,Ucomp=1,V=fs%V,Vcomp=1,W=fs%W,Wcomp=1, &
               &                  rho=fs%RHOG,rhocomp=1,visc=particle_gas,visccomp=1,snd=particle_gas,sndcomp=2, &
               &                  stress_div=particle_stress,gas_fraction=particle_gas_fraction,thermal=particle_thermal)
               ! Keep this step's deposited source even if the last particle has just exited.
               spray_source_ready=.true.
               call prepare_spray_rates()
               ! src(4) now includes mechanical work and minus the accepted particle heat increment.
               ! get_info also resets amrlpt's internal timers, so harvest them here
               call tmr_lptinfo%start(); call lp%get_info(); call tmr_lptinfo%stop(); lp_np=lp%np
               t_lptstep=lp%tmr_step%tmax; t_lptsrc=lp%tmr_src%tmax; t_lptvf=lp%tmr_vf%tmax
            else
               ! Empty spray: skip advance_to and its deposit pipeline, which is costly on the whole hierarchy
               call lp%src%setval(0.0_WP)
               lp%t=time%t
            end if
            lp_tprev=time%t
            call tmr_lpt%stop()
         end if

         ! ======================= Iterative midpoint sub-iterations =======================
         call tmr_loop%start()
         do while (time%it.le.time%itmax)

            ! Midpoint state
            call fs%prepare_midpoint()

            ! Surface tension pressure jump at the midpoint
            call fs%add_surface_tension_jump()

            ! RHS at the midpoint: semi-Lagrangian band transport from the old state and bulk fluxes
            call fs%get_dQdt(dQdt=dQdt,drhoUdt=drhoUdt,drhoVdt=drhoVdt,drhoWdt=drhoWdt,dt=time%dt,time=time%t)

            ! Spray reaction, conserved-variable update from the old state, sponge, and Q ghosts
            call tmr_rk%start()
            call add_spray_source()
            call fs%Q%lincomb   (a=1.0_WP,src1=fs%Qold   ,b=time%dt,src2=dQdt)
            call fs%rhoU%lincomb(a=1.0_WP,src1=fs%rhoUold,b=time%dt,src2=drhoUdt)
            call fs%rhoV%lincomb(a=1.0_WP,src1=fs%rhoVold,b=time%dt,src2=drhoVdt)
            call fs%rhoW%lincomb(a=1.0_WP,src1=fs%rhoWold,b=time%dt,src2=drhoWdt)
            ! Absorb outgoing waves in the sponge
            call apply_sponge(dt=time%dt,t=time%t)
            call fs%Q%average_down(); call fs%Q%fill(time=time%t)
            call tmr_rk%stop()

            ! Phasic total energies of the new iterate: transported internal energy, phasic work, and kinetic energy
            if (spray_source_ready) then
               call fs%get_energy(dt=time%dt,time=time%t,energy_source=spray_energy,force_source=spray_force)
            else
               call fs%get_energy(dt=time%dt,time=time%t)
            end if

            ! Rebuild PLIC of the new iterate (the band sweep only uses the old PLIC)
            call fs%build_plic(time=time%t)

            ! Get most up-to-date pressure: the last sub-iteration's relaxation is the end-of-step one
            call fs%apply_relax(dt=time%dt,time=time%t)

            ! Get primitive variables of the new iterate
            call fs%get_primitive(Q=fs%Q,rhoU=fs%rhoU,rhoV=fs%rhoV,rhoW=fs%rhoW,time=time%t)

            ! Increment sub-iteration counter
            time%it=time%it+1

         end do
         call tmr_loop%stop()
         ! =================================================================================

         ! Cell velocity for the transfer moments and the metrics (get_primitive does not update it)
         call tmr_props%start(); call fs%get_cell_velocity(); call tmr_props%stop()

         ! Transfer after the full step and before the regrid, which invalidates the labels
         if (do_transfer) then
            ! Reset the per-step counts, since the transfer runs only on event steps
            xfer_n_step=0; flot_n_step=0; xfer_fin=0.0_WP
            if (xfer_evt%occurs()) then
               call tmr_xfer%start()
               call transfer_drops()
               call tmr_xpart%start(); call lp%get_info(); lp_np=lp%np; call tmr_xpart%stop()
               call tmr_xfer%stop()
            end if
         end if

         ! Droplet QOIs before the regrid, reusing the transfer's labels when available
         if (metrics_evt%occurs()) then
            call tmr_diag%start()
            ! No fresh label from a transfer this step, so build one
            if (.not.core_mom_valid) then
               build_label_for_metrics: block
                  integer,  dimension(:),   allocatable :: cid,cmap
                  real(WP), dimension(:,:), allocatable :: mom
                  integer :: ns
                  call label_and_moments(cid,mom,cmap,ns)
                  deallocate(cid,mom,cmap)
               end block build_label_for_metrics
            end if
            call compute_drop_metrics()
            call drop_QOIs%write()
            call tmr_diag%stop()
         end if

         ! Regrid if event triggers
         regridded=regrid_evt%occurs()
         if (regridded) then
            call tmr_regrid%start()
            ! Book the change the regrid makes to the uncovered-cell mass, momentum and total energy
            book_regrid: block
               real(WP) :: mass0(2),mass1(2),tot0(4),tot1(4)
               call transfer_mesh_integrals(tot0,mass0)
               call amr%regrid(baselvl=0,time=time%t)
               call transfer_mesh_integrals(tot1,mass1)
               regrid_book=regrid_book+[mass1-mass0,tot1-tot0]
            end block book_regrid
            call fs%build_polygons()
            call gridfile%write()
            call tmr_regrid%stop()
         end if

         ! Compute viscosities
         call tmr_props%start()
         call get_viscosities()
         call tmr_props%stop()

         ! Add SGS models
         call tmr_sgs%start()
         call fs%add_viscartif(dt=time%dt,Cvisc=2.5e-2_WP)
         call fs%add_vreman(dt=time%dt)
         call tmr_sgs%stop()

         ! Compute Umag and Mach number (the regrid resets the cell velocity on remade levels)
         call tmr_props%start()
         if (regridded) call fs%get_cell_velocity()
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
         call tmr_props%stop()

         ! Visualization output: relabel the liquid and rebuild polygons from the current PLIC first
         call tmr_viz%start()
         if (plicviz_evt%occurs().or.tau_viz_evt%occurs().or.burst_any()) call fs%build_polygons()
         if (viz_evt%occurs().or.plicviz_evt%occurs().or.tau_viz_evt%occurs().or.burst_any()) then
            call build_current_ccl()
            call fs%update_surfvars()
         end if
         if (viz_evt%occurs()) call viz%write(time=time%t)
         if (plicviz_evt%occurs()) call plicviz%write(time=time%t)
         if (tau_viz_evt%occurs()) call tau_viz%write(time=time%t)
         if (do_transfer) then
            if (lpviz_evt%occurs()) call lpviz%write(time=time%t)
         end if
         ! Burst frames (PLIC with ccl_id, spray) inside any exposure window
         call write_bursts()
         call tmr_viz%stop()

         ! Checkpoint save
         call tmr_ckpt%start()
         if (save_evt%occurs()) then
            save_checkpoint: block
               use string, only: rtoa
               character(len=str_medium) :: ckdir
               ckdir='restart/drop_'//trim(adjustl(rtoa(time%t)))
               call io%write(dirname=trim(ckdir),time=time%t,step=time%n)
               ! The spray goes to a sibling directory, since AMReX's particle checkpoint cleans its parent
               if (do_transfer) call lp%write(dirname=trim(ckdir)//'_spray')
            end block save_checkpoint
         end if
         call tmr_ckpt%stop()

         ! Perform and output monitoring
         call tmr_mon%start()
         call fs%get_info(); Etot=fs%Qint(3)+fs%Qint(4)
         call mfile%write()
         call consfile%write()
         call cflfile%write()
         call tfile%write()
         call rescfile%write()

         ! Droplet QOIs are written above, before the regrid, where the amrccl labels are live

         ! Spray monitoring
         if (do_transfer) then
            call spray_ledger()
            call sprayfile%write()
         end if
         call tmr_mon%stop()

         ! Per-phase profiling: slowest rank for each phase
         call tmr_step_all%stop()
         call perf_stats()
         call perffile%write()

      end do

      ! Force a final checkpoint on exit
      final_checkpoint: block
         use string, only: rtoa
         character(len=str_medium) :: ckdir
         ckdir='restart/drop_'//trim(adjustl(rtoa(time%t)))
         call io%write(dirname=trim(ckdir),time=time%t,step=time%n)
         if (do_transfer) call lp%write(dirname=trim(ckdir)//'_spray')
      end block final_checkpoint

   contains

      !> Zero every phase timer at the start of a step
      subroutine reset_timers()
         implicit none
         call tmr_step_all%reset(); call tmr_prep%reset(); call tmr_lpt%reset(); call tmr_lptinfo%reset()
         call tmr_divstress%reset(); call tmr_loop%reset(); call tmr_rk%reset(); call tmr_xfer%reset()
         call tmr_ccl%reset(); call tmr_msweep%reset(); call tmr_mcomm%reset(); call tmr_xpart%reset()
         call tmr_surg%reset(); call tmr_resync%reset(); call tmr_diag%reset(); call tmr_regrid%reset()
         call tmr_props%reset(); call tmr_sgs%reset(); call tmr_viz%reset(); call tmr_ckpt%reset()
         call tmr_mon%reset()
      end subroutine reset_timers

      !> Slowest-rank time of every phase, in a single reduction
      subroutine perf_stats()
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_MAX
         use parallel, only: MPI_REAL_WP
         implicit none
         integer, parameter :: ntmr=21
         real(WP), dimension(ntmr) :: tt
         integer :: ierr
         tt=[tmr_step_all%time,tmr_prep%time,tmr_lpt%time,tmr_lptinfo%time,tmr_divstress%time, &
         &   tmr_loop%time,tmr_rk%time,tmr_xfer%time,tmr_ccl%time,tmr_msweep%time,tmr_mcomm%time, &
         &   tmr_xpart%time,tmr_surg%time,tmr_resync%time,tmr_diag%time,tmr_regrid%time, &
         &   tmr_props%time,tmr_sgs%time,tmr_viz%time,tmr_ckpt%time,tmr_mon%time]
         call MPI_ALLREDUCE(MPI_IN_PLACE,tt,ntmr,MPI_REAL_WP,MPI_MAX,amr%comm,ierr)
         t_step  =tt( 1); t_prep  =tt( 2); t_lpt   =tt( 3); t_lptinfo=tt( 4); t_divstress=tt( 5)
         t_loop  =tt( 6); t_rk    =tt( 7); t_xfer  =tt( 8); t_ccl    =tt( 9); t_msweep   =tt(10)
         t_mcomm =tt(11); t_xpart =tt(12); t_surg  =tt(13); t_resync =tt(14); t_diag     =tt(15)
         t_regrid=tt(16); t_props =tt(17); t_sgs   =tt(18); t_viz    =tt(19); t_ckpt     =tt(20)
         t_mon   =tt(21)
         t_other=t_step-(t_prep+t_lpt+t_loop+t_xfer+t_diag+t_regrid+t_props+t_sgs+t_viz+t_ckpt+t_mon)
      end subroutine perf_stats

      !> Integrate lp%src over the domain so the coupling can be audited from the monitor
      subroutine spray_ledger()
         use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
         use amrex_interface,  only: amrmask_make_fine
         use mpi_f08,          only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM,MPI_MIN,MPI_MAX
         use parallel,         only: MPI_REAL_WP
         implicit none
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         type(amrex_imultifab) :: mask
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pS,pQ
         integer,  dimension(:,:,:,:), contiguous, pointer :: pMask
         real(WP) :: vol,rho
         integer :: lvl,i,j,k,ierr
         src_momx=0.0_WP; src_ener=0.0_WP; src_emin=huge(1.0_WP); src_dumax=0.0_WP
         do lvl=0,amr%clvl()
            vol=amr%cell_vol(lvl)
            if (lvl.lt.amr%clvl()) then
               call amrex_imultifab_build(mask,amr%ba(lvl),amr%dm(lvl),1,0)
               call amrmask_make_fine(mask,amr%ba(lvl+1),[amr%rrefx(lvl),amr%rrefy(lvl),amr%rrefz(lvl)],0,1)
            end if
            call amr%mfiter_build(lvl,mfi,tiling=.false.)
            do while (mfi%next())
               pS=>lp%src%mf(lvl)%dataptr(mfi)
               pQ=>fs%Q%mf(lvl)%dataptr(mfi)
               if (lvl.lt.amr%clvl()) pMask=>mask%dataptr(mfi)
               bx=mfi%validbox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                        ! Count only cells not covered by a finer level
                        if (lvl.lt.amr%clvl()) then; if (pMask(i,j,k,1).eq.0) cycle; end if
                        src_momx=src_momx+vol*pS(i,j,k,1)
                        src_ener=src_ener+vol*pS(i,j,k,4)
                        src_emin=min(src_emin,pS(i,j,k,4))
                        ! Velocity change the coupling imposes in one step
                        rho=max(pQ(i,j,k,1)+pQ(i,j,k,2),tiny(1.0_WP))
                        src_dumax=max(src_dumax,norm2(pS(i,j,k,1:3))/rho)
                     end do; end do; end do
            end do
            call amr%mfiter_destroy(mfi)
            if (lvl.lt.amr%clvl()) call amrex_imultifab_destroy(mask)
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,src_dumax,1,MPI_REAL_WP,MPI_MAX,amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,src_momx,1,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,src_ener,1,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,src_emin,1,MPI_REAL_WP,MPI_MIN,amr%comm,ierr)
      end subroutine spray_ledger

      !> Apply the frozen rates to each new RHS, which updates from the old flow state.
      subroutine add_spray_source()
         if (.not.spray_source_ready) return
         call drhoUdt%add(spray_force(1))
         call drhoVdt%add(spray_force(2))
         call drhoWdt%add(spray_force(3))
      end subroutine add_spray_source

   end subroutine simulation_run

   !> Convert accumulated particle increments to rates once; reuse these exact faces in get_energy.
   subroutine prepare_spray_rates()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pS,pF
      integer :: lvl,i,j,k
      call spray_cells%copy(lp%src,srccomp=1,dstcomp=1,ncomp=3)
      do lvl=0,amr%clvl()
         call spray_cells%mf(lvl)%setval(amr%dx(lvl),icomp=4,ncomp=1)
         call spray_cells%mf(lvl)%setval(amr%dy(lvl),icomp=5,ncomp=1)
         call spray_cells%mf(lvl)%setval(amr%dz(lvl),icomp=6,ncomp=1)
      end do
      ! Piecewise-constant coarse ghost values retain the adjacent coarse cell's source and width.
      ! A face owns half of each adjacent cell volume: F=(hL*SL+hR*SR)/(hL+hR).
      ! This reduces to the original arithmetic average away from a coarse/fine interface.
      call spray_cells%fill(time=time%t)
      call spray_energy%copy(lp%src,srccomp=4,dstcomp=1,ncomp=1)
      call spray_energy%mult(1.0_WP/dt_src)
      do lvl=0,amr%clvl()
         call amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            pS=>spray_cells%mf(lvl)%dataptr(mfi)
            pF=>spray_force(1)%mf(lvl)%dataptr(mfi); bx=mfi%nodaltilebox(1)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pF(i,j,k,1)=(pS(i-1,j,k,4)*pS(i-1,j,k,1)+pS(i,j,k,4)*pS(i,j,k,1))/(pS(i-1,j,k,4)+pS(i,j,k,4))/dt_src
            end do; end do; end do
            pF=>spray_force(2)%mf(lvl)%dataptr(mfi); bx=mfi%nodaltilebox(2)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pF(i,j,k,1)=(pS(i,j-1,k,5)*pS(i,j-1,k,2)+pS(i,j,k,5)*pS(i,j,k,2))/(pS(i,j-1,k,5)+pS(i,j,k,5))/dt_src
            end do; end do; end do
            pF=>spray_force(3)%mf(lvl)%dataptr(mfi); bx=mfi%nodaltilebox(3)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pF(i,j,k,1)=(pS(i,j,k-1,6)*pS(i,j,k-1,3)+pS(i,j,k,6)*pS(i,j,k,3))/(pS(i,j,k-1,6)+pS(i,j,k,6))/dt_src
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
      ! Momentum recovery restricts/overrides faces too; use the identical synchronized force
      ! in the momentum RHS and energy work, including the coarse side of refinement interfaces.
      call fs%sync_momentum(spray_force(1),spray_force(2),spray_force(3))
   end subroutine prepare_spray_rates


   !> Droplet QOIs for all liquid (isocontour family) and for the largest structure (core family)
   subroutine compute_drop_metrics()
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
      use amrex_interface,  only: amrmask_make_fine
      use parallel,         only: MPI_REAL_WP,comm
      use mpi_f08,          only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM,MPI_MIN
      implicit none
      integer :: lvl,i,j,k,ierr,d
      integer, dimension(3) :: ncell
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_imultifab) :: mask
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pRHOL
      integer,  dimension(:,:,:,:), contiguous, pointer :: pMask
      real(WP) :: dx,dy,dz,vol,x_cc,y_cc,z_cc,vf_local,mi
      real(WP) :: M_tot,Sx,Sy,Sz,Sxx,Syy,Szz,Sxy,Sxz,Syz,xc,yc,zc
      real(WP) :: mass_01,Mc

      ! Core family: cached moments of the selected core, no mesh sweep
      core_valid=0.0_WP; core_open=0.0_WP; core_wrap=0.0_WP
      MOI_xx_core=0.0_WP; MOI_yy_core=0.0_WP; MOI_zz_core=0.0_WP
      MOI_xy_core=0.0_WP; MOI_xz_core=0.0_WP; MOI_yz_core=0.0_WP
      core_mass=0.0_WP; core_vol=0.0_WP
      MLE_core=0.0_WP; MLE_core001=0.0_WP; MLE_core01=0.0_WP; MLE_core05=0.0_WP
      if (core_mom_valid) then
         Mc=core_mom(2)
         if (Mc.gt.0.0_WP) then
            core_valid=1.0_WP
            core_vol =core_mom(1)
            core_mass=Mc
            if (ccl%struct(core_id)%open) core_open=1.0_WP
            ! Wrap flag only for directions the mesh resolves
            ncell=[amr%nx,amr%ny,amr%nz]
            do d=1,3
               if (ncell(d).le.1) cycle
               if (ccl%struct(core_id)%per(d).eq.1) core_wrap=1.0_WP
            end do
            ! Arithmetic barycenter, consistent with the raw second moments
            xc=core_mom(3)/Mc; yc=core_mom(4)/Mc; zc=core_mom(5)/Mc
            ! Parallel axis: I_xx = sum(m*(y^2+z^2)) - M*(yc^2+zc^2), and so on
            MOI_xx_core=(core_mom(17)+core_mom(18))-Mc*(yc**2+zc**2)
            MOI_yy_core=(core_mom(16)+core_mom(18))-Mc*(xc**2+zc**2)
            MOI_zz_core=(core_mom(16)+core_mom(17))-Mc*(xc**2+yc**2)
            MOI_xy_core=-(core_mom(19)-Mc*xc*yc)
            MOI_xz_core=-(core_mom(20)-Mc*xc*zc)
            MOI_yz_core=-(core_mom(21)-Mc*yc*zc)
            ! Unreached VF thresholds fall back to the drop's own edge
            MLE_core   =core_mle(1)
            MLE_core001=merge(core_mle(1),core_mle(2),core_mle(2).ge.0.5_WP*huge(1.0_WP))
            MLE_core01 =merge(core_mle(1),core_mle(3),core_mle(3).ge.0.5_WP*huge(1.0_WP))
            MLE_core05 =merge(core_mle(1),core_mle(4),core_mle(4).ge.0.5_WP*huge(1.0_WP))
         end if
      end if

      ! Isocontour family: hierarchy sweep masked against finer levels
      MLE_VF01=+huge(1.0_WP); MLE_VF05=+huge(1.0_WP)
      mass_01=0.0_WP
      M_tot=0.0_WP; Sx=0.0_WP; Sy=0.0_WP; Sz=0.0_WP
      Sxx=0.0_WP; Syy=0.0_WP; Szz=0.0_WP
      Sxy=0.0_WP; Sxz=0.0_WP; Syz=0.0_WP

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
                     if (vf_local.lt.0.1_WP) cycle

                     ! Cell-center coordinates
                     x_cc=amr%xlo+(real(i,WP)+0.5_WP)*dx
                     y_cc=amr%ylo+(real(j,WP)+0.5_WP)*dy
                     z_cc=amr%zlo+(real(k,WP)+0.5_WP)*dz

                     ! MLE: windward leading face at each VF threshold
                     MLE_VF01=min(MLE_VF01,x_cc-0.5_WP*dx)
                     if (vf_local.ge.0.5_WP) MLE_VF05=min(MLE_VF05,x_cc-0.5_WP*dx)

                     ! Mass and raw moments about the origin
                     mi=vf_local*pRHOL(i,j,k,1)*vol
                     mass_01=mass_01+mi
                     M_tot=M_tot+mi
                     Sx=Sx+mi*x_cc; Sy=Sy+mi*y_cc; Sz=Sz+mi*z_cc
                     Sxx=Sxx+mi*x_cc**2; Syy=Syy+mi*y_cc**2; Szz=Szz+mi*z_cc**2
                     Sxy=Sxy+mi*x_cc*y_cc; Sxz=Sxz+mi*x_cc*z_cc; Syz=Syz+mi*y_cc*z_cc
                  end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
         if (lvl.lt.amr%clvl()) call amrex_imultifab_destroy(mask)
      end do

      ! MPI reductions
      call MPI_ALLREDUCE(MPI_IN_PLACE,MLE_VF01,1,MPI_REAL_WP,MPI_MIN,comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,MLE_VF05,1,MPI_REAL_WP,MPI_MIN,comm,ierr)
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

      drop_massVF01=mass_01

      ! Inertia tensor about the VF-weighted COM, via the parallel axis theorem
      if (M_tot.gt.0.0_WP) then
         xc=Sx/M_tot; yc=Sy/M_tot; zc=Sz/M_tot
         MOI_xx=(Syy+Szz)-M_tot*(yc**2+zc**2)
         MOI_yy=(Sxx+Szz)-M_tot*(xc**2+zc**2)
         MOI_zz=(Sxx+Syy)-M_tot*(xc**2+yc**2)
         MOI_xy=-(Sxy-M_tot*xc*yc)
         MOI_xz=-(Sxz-M_tot*xc*zc)
         MOI_yz=-(Syz-M_tot*yc*zc)
      else
         MOI_xx=0.0_WP; MOI_yy=0.0_WP; MOI_zz=0.0_WP
         MOI_xy=0.0_WP; MOI_xz=0.0_WP; MOI_yz=0.0_WP
         MLE_VF01=0.0_WP; MLE_VF05=0.0_WP
      end if
      ! An unreached threshold must not leave huge() in the monitor
      if (MLE_VF05.ge.0.5_WP*huge(1.0_WP)) MLE_VF05=MLE_VF01

   end subroutine compute_drop_metrics

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
      ! Finalize burst writers
      finalize_bursts: block
         integer :: ib
         do ib=1,nburst
            if (.not.burst_open(ib)) cycle
            call burst_viz(ib)%finalize()
            if (do_transfer) call burst_lpviz(ib)%finalize()
         end do
         if (allocated(burst_tau)) deallocate(burst_tau,burst_t0,burst_t1,burst_open,burst_viz,burst_lpviz)
      end block finalize_bursts
      ! Finalize drop transfer
      call ccl%finalize()
      if (do_transfer) then
         call spray_cells%finalize()
         call spray_energy%finalize()
         call spray_force(1)%finalize(); call spray_force(2)%finalize(); call spray_force(3)%finalize()
         call particle_stress%finalize()
         call particle_thermal%finalize()
         call particle_gas%finalize()
         call particle_gas_fraction%finalize()
         call lp%finalize()
         call lpviz%finalize()
         call lpviz_evt%finalize()
         call sprayfile%finalize()
         call xfer_evt%finalize()
      end if
      call perffile%finalize()
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
      call drop_QOIs%finalize()
      call metrics_evt%finalize()
      if (allocated(spg_redge)) deallocate(spg_redge)
   end subroutine simulation_final

end module simulation
