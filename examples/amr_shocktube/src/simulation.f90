!> 1D multiphase shock tube simulation with amr code, using the staggered compressible multiphase solver
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

   !> Previous sub-iterate, for the convergence monitor
   type(amrdata) :: Qprev,rhoUprev

   !> Visualization
   type(event) :: viz_evt
   type(amrviz) :: viz

   !> Simulation monitoring
   type(monitor) :: mfile,consfile,cflfile,gridfile,cvgfile

   !> Sub-iteration convergence: max over components of |X^{k+1}-X^k|/|X^{k+1}| (L-infinity)
   real(WP) :: cvgQ,cvgU

   !> Mixture total energy, the sum of the phasic total energies
   real(WP) :: Etot=0.0_WP,Etot0=0.0_WP,dEtot=0.0_WP   !< Mixture total energy, initial value, relative drift

   !> Materials
   type(nasg),      target :: water   !< Stiffened gas: NASG with zero co-volume
   type(ideal_gas), target :: air

   !> Relaxation model
   type(relax_ig_nasg), target :: relax_model

   !> Flow parameters
   real(WP) :: rhoG,pG_init            !< Gas state (right side)
   real(WP) :: rhoL,pL_init            !< Liquid state (left side)
   real(WP) :: x_int                   !< Interface position

   !> Artificial bulk viscosity coefficient and viscous CFL clipping
   real(WP) :: Cartif,max_visc_cfl

contains

   !> Levelset function for planar interface at x = x_int
   !> Returns positive for liquid (x < x_int), negative for gas (x >= x_int)
   function planar_levelset(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=x_int-xyz(1)
   end function planar_levelset

   !> Inviscid: zero out transport coefficients (artificial viscosity is added separately)
   subroutine get_viscosities()
      implicit none
      call fs%visc%setval(val=0.0_WP)
      call fs%beta%setval(val=0.0_WP)
   end subroutine get_viscosities

   !> User init callback - set Q, VF, and barycenters for 1D shocktube; face velocities are zero (set by on_init)
   !> Liquid on left (x < x_int), gas on right (x >= x_int)
   subroutine shocktube_init(solver,lvl,time,ba,dm)
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
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF,pCL,pCG
      real(WP), dimension(3) :: BL,BG
      real(WP) :: dx,dy,dz,myVF,IEL,IEG
      integer :: i,j,k
      integer, parameter :: nref=3
      ! Get mesh size
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl); dz=solver%amr%dz(lvl)
      ! Get internal energies
      IEL=water%get_e_from_p_rho(p=pL_init,rho=rhoL,y=[1.0_WP])
      IEG=air%get_e_from_p_rho  (p=pG_init,rho=rhoG,y=[1.0_WP])
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
                  ! Compute VF from planar levelset
                  call initialize_volume_moments(lo=[solver%amr%xlo+real(i  ,WP)*dx,solver%amr%ylo+real(j  ,WP)*dy,solver%amr%zlo+real(k  ,WP)*dz], &
                  &                              hi=[solver%amr%xlo+real(i+1,WP)*dx,solver%amr%ylo+real(j+1,WP)*dy,solver%amr%zlo+real(k+1,WP)*dz], &
                  &                              levelset=planar_levelset,time=time,level=nref,VFlo=VFlo,VF=myVF,BL=BL,BG=BG)
                  ! Store volume fraction
                  pVF(i,j,k,1)=myVF
                  ! Store barycenters
                  if (lvl.eq.solver%amr%maxlvl) then
                     pCL(i,j,k,:)=BL
                     pCG(i,j,k,:)=BG
                  end if
                  ! Set conserved variables: Q=(VF*rhoL,(1-VF)*rhoG,VF*rhoL*IL,(1-VF)*rhoG*IG); momentum lives on faces
                  pQ(i,j,k,1)=(       myVF)*rhoL
                  pQ(i,j,k,2)=(1.0_WP-myVF)*rhoG
                  pQ(i,j,k,3)=pQ(i,j,k,1)*IEL
                  pQ(i,j,k,4)=pQ(i,j,k,2)*IEG
               end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine shocktube_init

   !> Initialization of problem solver
   subroutine simulation_init
      use param,    only: param_read
      use messager, only: die
      implicit none

      ! Read EoS and flow parameters
      init_eos_and_flow: block
         use messager, only: log
         use string,   only: str_long
         character(len=str_long) :: message
         real(WP) :: GammaL,PinfL,CvL,GammaG,CvG
         ! Gas EoS parameters (ideal gas)
         call param_read('GammaG',GammaG)
         ! Liquid EoS parameters (stiffened gas)
         call param_read('GammaL',GammaL)
         call param_read('PinfL',PinfL)
         ! Specific heats
         call param_read('CvG',CvG)
         call param_read('CvL',CvL)
         ! Build materials
         call air%initialize  (gamma=GammaG,cv=CvG,q=0.0_WP,qp=0.0_WP,name='air')
         call water%initialize(gamma=GammaL,pinf=PinfL,cv=CvL,q=0.0_WP,qp=0.0_WP,name='water')
         ! Gas state (right side)
         call param_read('Gas density',rhoG)
         call param_read('Gas pressure',pG_init)
         ! Liquid state (left side)
         call param_read('Liquid density',rhoL)
         call param_read('Liquid pressure',pL_init)
         ! Interface location
         call param_read('Interface location',x_int)
         ! Log
         call water%print(); call air%print()
         write(message,'("[Gas]    rhoG=",es12.5," pG=",es12.5)') rhoG,pG_init; call log(message)
         write(message,'("[Liquid] rhoL=",es12.5," pL=",es12.5)') rhoL,pL_init; call log(message)
         write(message,'("[Interface] x_int=",es12.5)') x_int; call log(message)
      end block init_eos_and_flow

      ! Initialize AMR grid (1D domain: [-1.5, 1.5])
      create_amrgrid: block
         amr%name='shocktube'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         amr%xlo=-1.5_WP; amr%xhi=+1.5_WP
         amr%ylo=-1.0_WP; amr%yhi=+1.0_WP
         amr%zlo=-1.0_WP; amr%zhi=+1.0_WP
         amr%xper=.false.; amr%yper=.true.; amr%zper=.true.
         call param_read('Max level',amr%maxlvl)
         ! Enable quasi-1D: set y extent to match one cell at finest level
         if (amr%ny.eq.1) then
            amr%ylo=-0.5_WP*(amr%xhi-amr%xlo)/real(amr%nx*2**amr%maxlvl,WP)
            amr%yhi=+0.5_WP*(amr%xhi-amr%xlo)/real(amr%nx*2**amr%maxlvl,WP)
         end if
         ! Enable quasi-2D/1D: set z extent to match one cell at finest level
         if (amr%nz.eq.1) then
            amr%zlo=-0.5_WP*(amr%xhi-amr%xlo)/real(amr%nx*2**amr%maxlvl,WP)
            amr%zhi=+0.5_WP*(amr%xhi-amr%xlo)/real(amr%nx*2**amr%maxlvl,WP)
         end if
         amr%default_tiling=.false.
         call amr%initialize()
      end block create_amrgrid

      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=amr%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max dt',time%dtmax)
         call param_read('Max CFL',time%cflmax)
         ! Iterative midpoint sub-iterations
         call param_read('Subiterations',time%itmax,default=3)
         if (time%itmax.lt.1) call die('Require at least one sub-iteration')
         time%dt=time%dtmax
      end block initialize_timetracker

      ! Initialize compressible multiphase solver
      create_solver: block
         use amrex_amr_module, only: amrex_bc_foextrap
         use amrmpcomp_class,  only: BC_GAS,BC_LIQ
         ! Assign materials
         fs%liq=>water; fs%gas=>air
         ! Provide pressure relaxation model (mechanical only, no thermal relaxation)
         call relax_model%initialize(gas=air,liq=water); fs%relax=>relax_model
         relax_model%model=Prelax
         ! Refuse relaxations that leave the gas below 1e-2; the model's own RHOGmin keeps its 1e-2 default
         fs%rhoG_min=1.0e-2_WP
         ! Positivity gate temperature floors (0: undeclared, a pure cell outside its law's domain stops the run)
         call param_read('Liquid Tmin',fs%TminL,default=0.0_WP)
         call param_read('Gas Tmin'   ,fs%TminG,default=0.0_WP)
         ! Band momentum flux: CTU below the staggered phase fraction CTU alpha, MC elsewhere
         call param_read('CTU alpha',fs%alpha_ctu,default=0.5_WP)
         ! No VF treatment: build_plic and relaxation keep the transported VF
         fs%reset_vf=.false.
         ! Create flow solver
         call fs%initialize(amr=amr,name='shocktube')
         ! Set initial conditions
         fs%user_init=>shocktube_init
         ! Set BCs: outflow (foextrap) on both x boundaries
         if (.not.amr%xper) then
            fs%lo_bc(1)=BC_LIQ
            fs%hi_bc(1)=BC_GAS
            fs%Q%lo_bc(1,:)=amrex_bc_foextrap
            fs%Q%hi_bc(1,:)=amrex_bc_foextrap
            fs%U%lo_bc(1,:)=amrex_bc_foextrap; fs%U%hi_bc(1,:)=amrex_bc_foextrap
            fs%V%lo_bc(1,:)=amrex_bc_foextrap; fs%V%hi_bc(1,:)=amrex_bc_foextrap
            fs%W%lo_bc(1,:)=amrex_bc_foextrap; fs%W%hi_bc(1,:)=amrex_bc_foextrap
         end if
      end block create_solver

      ! Initialize workspaces
      create_workspace: block
         use amrdata_class, only: interp_none
         call dQdt%initialize   (amr,name='dQdt'   ,ncomp=fs%nQ,ng=0,interp=interp_none); call dQdt%register()
         call drhoUdt%initialize(amr,name='drhoUdt',ncomp=1,ng=0,nodal=[.true. ,.false.,.false.],interp=interp_none); call drhoUdt%register()
         call drhoVdt%initialize(amr,name='drhoVdt',ncomp=1,ng=0,nodal=[.false.,.true. ,.false.],interp=interp_none); call drhoVdt%register()
         call drhoWdt%initialize(amr,name='drhoWdt',ncomp=1,ng=0,nodal=[.false.,.false.,.true. ],interp=interp_none); call drhoWdt%register()
         call Qprev%initialize   (amr,name='Qprev'   ,ncomp=fs%nQ,ng=0,interp=interp_none); call Qprev%register()
         call rhoUprev%initialize(amr,name='rhoUprev',ncomp=1,ng=0,nodal=[.true. ,.false.,.false.],interp=interp_none); call rhoUprev%register()
         call Umag%initialize(amr,name='Umag',ncomp=1,ng=0,interp=interp_none); call Umag%register()
         call Mach%initialize(amr,name='Mach',ncomp=1,ng=0,interp=interp_none); call Mach%register()
      end block create_workspace

      ! Initialize grid (no regridding for uniform mesh)
      init_grid: block
         ! Artificial bulk viscosity (shock capturing)
         call param_read('Artificial viscosity coefficient',Cartif,default=5.0_WP)
         call param_read('Max viscous CFL',max_visc_cfl,default=0.5_WP)
         if (Cartif.lt.0.0_WP.or.max_visc_cfl.le.0.0_WP) call die('Invalid artificial viscosity parameters')
         ! Fresh start: post_regrid builds momentum from the (zero) face velocity and the primitive variables
         call amr%init_from_scratch(time=time%t)
         ! Phasic total energies from the internal contents written by user_init and the initial kinetic energy
         call fs%get_total_energy(time=time%t)
         ! Build PLIC
         call fs%build_plic(time%t)
         ! Set viscosities to zero
         call get_viscosities()
         ! Add artificial viscosity only
         call fs%add_viscartif(dt=time%dt,Cartif=Cartif,max_cfl=max_visc_cfl)
      end block init_grid

      ! Initialize visualization
      create_viz: block
         ! Create visualization object
         call viz%initialize(amr,'shocktube',use_hdf5=.false.)
         call viz%add_scalar(fs%VF,1,'VF')
         call viz%add_scalar(fs%RHOL,1,'RHOL')
         call viz%add_scalar(fs%RHOG,1,'RHOG')
         call viz%add_scalar(fs%PL,1,'PL')
         call viz%add_scalar(fs%PG,1,'PG')
         call viz%add_scalar(fs%P,1,'P')
         call viz%add_scalar(fs%UVW,1,'U')
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(Mach,1,'Mach')
         call viz%add_scalar(fs%beta,1,'beta')
         call viz%add_scalar(fs%C,1,'C')
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
         call fs%get_info(); Etot=fs%Qint(3)+fs%Qint(4); if (time%n.eq.0) Etot0=Etot; dEtot=(Etot-Etot0)/Etot0
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
         call mfile%add_column(fs%RHOGmin,'rhoGmin')
         call mfile%add_column(fs%RHOGmax,'rhoGmax')
         call mfile%add_column(fs%PGmin,'PGmin')
         call mfile%add_column(fs%PGmax,'PGmax')
         call mfile%add_column(fs%VFmin,'VFmin')
         call mfile%add_column(fs%VFmax,'VFmax')
         call mfile%add_column(fs%VFint,'VFint')
         call mfile%add_column(fs%dPmax,'dPmax')
         call mfile%add_column(fs%TLmin,'TLmin')
         call mfile%add_column(fs%TLmax,'TLmax')
         call mfile%add_column(fs%TGmin,'TGmin')
         call mfile%add_column(fs%TGmax,'TGmax')
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
         call consfile%add_column(fs%Qint(1),'Liquid Mass')
         call consfile%add_column(fs%Qint(2),'Gas Mass')
         call consfile%add_column(fs%Qint(3),'Liquid TotEnergy')
         call consfile%add_column(fs%Qint(4),'Gas TotEnergy')
         call consfile%add_column(Etot,'Total energy')
         call consfile%add_column(dEtot,'Energy drift')
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
         ! Create sub-iteration convergence monitor
         cvgQ=0.0_WP; cvgU=0.0_WP
         cvgfile=monitor(amRoot=amr%amRoot,name='convergence')
         call cvgfile%add_column(time%n,'Timestep number')
         call cvgfile%add_column(time%it,'Iteration')
         call cvgfile%add_column(time%t,'Time')
         call cvgfile%add_column(cvgQ,'dQ/Q')
         call cvgfile%add_column(cvgU,'drhoU/rhoU')
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

         ! Transport coefficients from the start-of-step state: inviscid, with artificial bulk viscosity for shock capturing
         call get_viscosities()
         call fs%add_viscartif(dt=time%dt,Cartif=Cartif,max_cfl=max_visc_cfl)

         ! ======================= Iterative midpoint sub-iterations =======================
         do while (time%it.le.time%itmax)

            ! Midpoint state
            call fs%prepare_midpoint()

            ! RHS at the midpoint: semi-Lagrangian band transport from the old state and bulk fluxes
            call fs%get_dQdt(dQdt=dQdt,drhoUdt=drhoUdt,drhoVdt=drhoVdt,drhoWdt=drhoWdt,dt=time%dt,time=time%t)

            ! Update conserved variables from the old state
            call Qprev%copy(src=fs%Q); call rhoUprev%copy(src=fs%rhoU)
            call fs%Q%lincomb   (a=1.0_WP,src1=fs%Qold   ,b=time%dt,src2=dQdt)
            call fs%rhoU%lincomb(a=1.0_WP,src1=fs%rhoUold,b=time%dt,src2=drhoUdt)
            call fs%rhoV%lincomb(a=1.0_WP,src1=fs%rhoVold,b=time%dt,src2=drhoVdt)
            call fs%rhoW%lincomb(a=1.0_WP,src1=fs%rhoWold,b=time%dt,src2=drhoWdt)
            call fs%Q%average_down(); call fs%Q%fill(time=time%t)

            ! Phasic total energies of the new iterate: transported internal energy, phasic work, and kinetic energy
            call fs%get_energy(dt=time%dt,time=time%t)

            ! Rebuild PLIC of the new iterate (the band sweep only uses the old PLIC)
            call fs%build_plic(time=time%t)

            ! Add surface tension term
            ! call fs%add_surface_tension(scale=time%dt)

            ! Get most up-to-date pressure: the last sub-iteration's relaxation is the end-of-step one
            call fs%apply_relax(dt=time%dt,time=time%t)

            ! Get primitive variables of the new iterate
            call fs%get_primitive(Q=fs%Q,rhoU=fs%rhoU,rhoV=fs%rhoV,rhoW=fs%rhoW,time=time%t)

            ! Monitor convergence of the sub-iterations
            call get_convergence()
            call cvgfile%write()

            ! Increment sub-iteration counter
            time%it=time%it+1

         end do
         ! =================================================================================

         ! Stop on nonfinite or negative state
         call check_state()

         ! Visualization output
         if (viz_evt%occurs()) then
            ! Cell-centered velocity, Umag, and Mach number are only computed for visualization
            call fs%get_cell_velocity()
            call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
            call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
            call viz%write(time=time%t)
         end if

         ! Perform and output monitoring
         call fs%get_info(); Etot=fs%Qint(3)+fs%Qint(4); if (time%n.eq.0) Etot0=Etot; dEtot=(Etot-Etot0)/Etot0
         call mfile%write()
         call consfile%write()
         call cflfile%write()

      end do

      ! Extract mixture data at final time for x in [0, 1.5]
      mixture_data: block
         use mpi_f08
         use parallel,         only: MPI_REAL_WP
         use string,           only: str_medium
         use filesys,          only: makedir,isdir
         use amrex_amr_module, only: amrex_mfiter,amrex_box
         integer :: i,j,k,lvl,iunit,ierr,nlocal,ntotal
         integer, dimension(:), allocatable :: recvcount,displs
         real(WP), dimension(:), allocatable :: local_x,local_VF,local_RHO,local_P,local_UVW
         real(WP), dimension(:), allocatable :: global_x,global_VF,global_RHO,global_P,global_UVW
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pQ,pP,pUVW
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP) :: x_cc,vf_val,dx
         character(len=str_medium) :: filename,timestamp
         lvl=0 ! Single level, no AMR
         dx=amr%dx(lvl)
         ! Cell-centered velocity
         call fs%get_cell_velocity()
         ! First pass: count local cells in range [0, 1.5]
         nlocal=0
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            bx=mfi%tilebox()
            do i=bx%lo(1),bx%hi(1)
               x_cc=amr%xlo+(real(i,WP)+0.5_WP)*dx
               if (x_cc.ge.0.0_WP.and.x_cc.le.1.5_WP) nlocal=nlocal+1
            end do
         end do
         call amr%mfiter_destroy(mfi)
         ! Allocate local arrays
         allocate(local_x(nlocal),local_VF(nlocal),local_RHO(nlocal),local_P(nlocal),local_UVW(nlocal))
         ! Second pass: fill local arrays - mixture density is Q(1)+Q(2), mixture pressure is P=VF*PL+(1-VF)*PG
         nlocal=0
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pVF  =>fs%VF%mf(lvl)%dataptr(mfi)
            pQ   =>fs%Q%mf(lvl)%dataptr(mfi)
            pP   =>fs%P%mf(lvl)%dataptr(mfi)
            pUVW =>fs%UVW%mf(lvl)%dataptr(mfi)
            bx=mfi%tilebox()
            j=bx%lo(2); k=bx%lo(3)
            do i=bx%lo(1),bx%hi(1)
               x_cc=amr%xlo+(real(i,WP)+0.5_WP)*dx
               if (x_cc.ge.0.0_WP.and.x_cc.le.1.5_WP) then
                  nlocal=nlocal+1
                  local_x(nlocal)  =x_cc
                  local_VF(nlocal) =pVF(i,j,k,1)
                  local_RHO(nlocal)=pQ(i,j,k,1)+pQ(i,j,k,2)
                  local_P(nlocal)  =pP(i,j,k,1)
                  local_UVW(nlocal)=pUVW(i,j,k,1)
               end if
            end do
         end do
         call amr%mfiter_destroy(mfi)
         ! Gather counts from all ranks
         allocate(recvcount(amr%nproc),displs(amr%nproc))
         call MPI_GATHER(nlocal,1,MPI_INTEGER,recvcount,1,MPI_INTEGER,0,amr%comm,ierr)
         ! Compute displacements and total count on root
         ntotal=0
         if (amr%amRoot) then
            displs(1)=0
            do i=2,amr%nproc
               displs(i)=displs(i-1)+recvcount(i-1)
            end do
            ntotal=displs(amr%nproc)+recvcount(amr%nproc)
         end if
         ! Allocate global arrays on root
         allocate(global_x(ntotal),global_VF(ntotal),global_RHO(ntotal),global_P(ntotal),global_UVW(ntotal))
         ! Gather data to root
         call MPI_GATHERV(local_x,  nlocal,MPI_REAL_WP,global_x,  recvcount,displs,MPI_REAL_WP,0,amr%comm,ierr)
         call MPI_GATHERV(local_VF, nlocal,MPI_REAL_WP,global_VF, recvcount,displs,MPI_REAL_WP,0,amr%comm,ierr)
         call MPI_GATHERV(local_RHO,nlocal,MPI_REAL_WP,global_RHO,recvcount,displs,MPI_REAL_WP,0,amr%comm,ierr)
         call MPI_GATHERV(local_P,  nlocal,MPI_REAL_WP,global_P,  recvcount,displs,MPI_REAL_WP,0,amr%comm,ierr)
         call MPI_GATHERV(local_UVW,nlocal,MPI_REAL_WP,global_UVW,recvcount,displs,MPI_REAL_WP,0,amr%comm,ierr)
         ! Sort by position on root (insertion sort — small array)
         if (amr%amRoot) then
            do i=2,ntotal
               if (global_x(i).lt.global_x(i-1)) then
                  ! Find insertion point
                  j=i-1
                  do while (j.ge.1)
                     if (global_x(j).le.global_x(i)) exit ! Separate test: .and. does not short-circuit, global_x(0) is out of bounds
                     j=j-1
                  end do
                  j=j+1
                  ! Rotate element i into position j
                  vf_val=global_x(i);   global_x(j+1:i)  =global_x(j:i-1);   global_x(j)  =vf_val
                  vf_val=global_VF(i);  global_VF(j+1:i) =global_VF(j:i-1);  global_VF(j) =vf_val
                  vf_val=global_RHO(i); global_RHO(j+1:i)=global_RHO(j:i-1); global_RHO(j)=vf_val
                  vf_val=global_P(i);   global_P(j+1:i)  =global_P(j:i-1);   global_P(j)  =vf_val
                  vf_val=global_UVW(i); global_UVW(j+1:i)=global_UVW(j:i-1); global_UVW(j)=vf_val
               end if
            end do
         end if
         ! Only root outputs to file
         if (amr%amRoot) then
            if (.not.isdir('data')) call makedir('data')
            filename='data_'; write(timestamp,'(es12.5)') time%t
            open(newunit=iunit,file='data/'//trim(adjustl(filename))//trim(adjustl(timestamp))//'.csv', &
            &    form='formatted',status='replace',iostat=ierr)
            write(iunit,'(a)') 'Position,VF,RHO,P,Velocity'
            do i=1,ntotal
               write(iunit,'(g0.17,",",g0.17,",",g0.17,",",g0.17,",",g0.17)') &
               &    global_x(i),global_VF(i),global_RHO(i),global_P(i),global_UVW(i)
            end do
            close(iunit)
         end if
         ! Clean up
         deallocate(local_x,local_VF,local_RHO,local_P,local_UVW)
         deallocate(global_x,global_VF,global_RHO,global_P,global_UVW)
         deallocate(recvcount,displs)
      end block mixture_data

   end subroutine simulation_run

   !> Change between two sub-iterates (level 0, L-infinity), relative to the new iterate
   subroutine get_convergence()
      implicit none
      integer :: n
      real(WP) :: ref
      call Qprev%lincomb(a=1.0_WP,src1=fs%Q,b=-1.0_WP,src2=Qprev)
      call rhoUprev%lincomb(a=1.0_WP,src1=fs%rhoU,b=-1.0_WP,src2=rhoUprev)
      cvgQ=0.0_WP
      do n=1,fs%nQ
         ref=fs%Q%norm0(lvl=0,comp=n)
         if (ref.gt.0.0_WP) cvgQ=max(cvgQ,Qprev%norm0(lvl=0,comp=n)/ref)
      end do
      ref=fs%rhoU%norm0(lvl=0)
      cvgU=0.0_WP; if (ref.gt.0.0_WP) cvgU=rhoUprev%norm0(lvl=0)/ref
   end subroutine get_convergence

   !> Stop if any cell has a nonfinite state or negative phasic mass
   subroutine check_state()
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_INTEGER,MPI_MAX
      use messager, only: die
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), contiguous, pointer :: pQ(:,:,:,:),pP(:,:,:,:)
      integer :: bad,ierr,i,j,k,lvl
      bad=0
      do lvl=0,amr%clvl()
         call amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            bx=mfi%tilebox(); pQ=>fs%Q%mf(lvl)%dataptr(mfi); pP=>fs%P%mf(lvl)%dataptr(mfi)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               if (.not.all(ieee_is_finite(pQ(i,j,k,:))).or..not.ieee_is_finite(pP(i,j,k,1))) bad=1
               if (min(pQ(i,j,k,1),pQ(i,j,k,2)).lt.0.0_WP) bad=1
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,bad,1,MPI_INTEGER,MPI_MAX,amr%comm,ierr)
      if (bad.ne.0) call die('Nonfinite state or negative phasic mass: reduce dt/CFL')
   end subroutine check_state

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Finalize time
      call time%finalize()
      ! Finalize solver
      call fs%finalize()
      call dQdt%finalize(); call drhoUdt%finalize(); call drhoVdt%finalize(); call drhoWdt%finalize()
      call Qprev%finalize(); call rhoUprev%finalize()
      call Umag%finalize()
      call Mach%finalize()
      ! Finalize materials
      call water%finalize()
      call air%finalize()
      ! Finalize visualization
      call viz%finalize()
      call viz_evt%finalize()
      ! Finalize monitoring
      call mfile%finalize()
      call cflfile%finalize()
      call consfile%finalize()
      call gridfile%finalize()
      call cvgfile%finalize()
      ! Finalize grid
      call amr%finalize()
   end subroutine simulation_final

end module simulation
