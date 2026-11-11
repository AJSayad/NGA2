!> Inviscid compressible Taylor-Green vortex with a slab of a second fluid, as Jain & Moin (J. Comput. Phys. 2022)
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
   use stiffened_gas_class, only: stiffened_gas
   use nasg_class,          only: nasg
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

   ! Restart parameters
   type(amrio) :: io
   type(event) :: save_evt
   character(len=str_medium) :: restart_dir
   logical :: restarted
   integer :: restart_step=0
   real(WP) :: restart_time

   !> Simulation monitoring
   type(monitor) :: mfile,consfile,cflfile,gridfile

   !> Taylor-Green diagnostics: kinetic energy and its initial value, phase entropies and their initial values, relative
   !> changes K/K0 and (S_l-S_l0)/|S_l0|
   real(WP) :: Kin,Kin0,Kratio
   real(WP), dimension(2) :: Sph,Sph0,dSph

   !> Mixture total energy, the sum of the phasic total energies
   real(WP) :: Etot=0.0_WP,Etot0=0.0_WP,dEtot=0.0_WP   !< Mixture total energy, initial value, relative drift

   !> Materials: fluid 1 (slab, the solver's liquid) and fluid 2 (the solver's gas)
   type(nasg),          target :: water   !< Stiffened gas: NASG 
   type(stiffened_gas), target :: air

   !> Relaxation model
   type(relax_ig_nasg), target :: relax_model

   !> Flow parameters
   real(WP) :: rho1,rho2               !< Densities of fluid 1 (slab) and fluid 2
   real(WP) :: p0                      !< Base pressure
   real(WP) :: mach0,U0                !< Mach number and velocity amplitude U0=mach0*c2(p0,rho2)
   real(WP), dimension(3) :: center    !< Centre of the box (the slab is centred on x=center(1))
   real(WP) :: radius                  !< Slab half-width (0: no slab)

   !> Artificial bulk viscosity coefficient and viscous CFL clipping
   real(WP) :: Cartif,max_visc_cfl

contains

   !> Levelset function for the slab of half-width radius centred on x=center(1): positive in fluid 1
   function slab_levelset(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=radius-abs(xyz(1)-center(1))
   end function slab_levelset

   !> Inviscid: zero out transport coefficients
   subroutine get_viscosities()
      implicit none
      call fs%visc%setval(val=0.0_WP)
      call fs%beta%setval(val=0.0_WP)
   end subroutine get_viscosities

   !> User init callback - set Q, VF, barycenters, and face velocities of the Taylor-Green vortex with the slab
   subroutine tgv_init(solver,lvl,time,ba,dm)
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
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF,pCL,pCG,pU,pV
      real(WP), dimension(3) :: BL,BG,xc
      real(WP) :: dx,dy,dz,myVF,rho,p,x,y,z
      integer :: i,j,k
      integer, parameter :: nref=3
      ! Get mesh size
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl); dz=solver%amr%dz(lvl)
      ! Use passed ba/dm since grid is being constructed
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         ! Get pointers to data
         pQ =>solver%Q%mf(lvl)%dataptr(mfi)
         pVF=>solver%VF%mf(lvl)%dataptr(mfi)
         pU =>solver%U%mf(lvl)%dataptr(mfi)
         pV =>solver%V%mf(lvl)%dataptr(mfi)
         if (lvl.eq.solver%amr%maxlvl) then
            pCL=>solver%CL%dataptr(mfi)
            pCG=>solver%CG%dataptr(mfi)
         end if
         ! Cells, ghosts included
         bx=mfi%growntilebox(solver%nover)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            xc=[solver%amr%xlo+(real(i,WP)+0.5_WP)*dx,solver%amr%ylo+(real(j,WP)+0.5_WP)*dy,solver%amr%zlo+(real(k,WP)+0.5_WP)*dz]
            ! Volume fraction and barycenters of the sharp slab
            if (radius.gt.0.0_WP) then
               call initialize_volume_moments(lo=[solver%amr%xlo+real(i  ,WP)*dx,solver%amr%ylo+real(j  ,WP)*dy,solver%amr%zlo+real(k  ,WP)*dz], &
               &                              hi=[solver%amr%xlo+real(i+1,WP)*dx,solver%amr%ylo+real(j+1,WP)*dy,solver%amr%zlo+real(k+1,WP)*dz], &
               &                              levelset=slab_levelset,time=time,level=nref,VFlo=VFlo,VF=myVF,BL=BL,BG=BG)
            else
               myVF=0.0_WP; BL=xc; BG=xc
            end if
            pVF(i,j,k,1)=myVF
            if (lvl.eq.solver%amr%maxlvl) then
               pCL(i,j,k,:)=BL
               pCG(i,j,k,:)=BG
            end if
            ! Pressure from the mixture density, shared by both phases; Q=(VF*rho1,(1-VF)*rho2,VF*rho1*e1,(1-VF)*rho2*e2)
            rho=myVF*rho1+(1.0_WP-myVF)*rho2
            p=p0+rho*U0**2/16.0_WP*(cos(2.0_WP*xc(1))+cos(2.0_WP*xc(2)))*(cos(2.0_WP*xc(3))+2.0_WP)
            pQ(i,j,k,1)=(       myVF)*rho1
            pQ(i,j,k,2)=(1.0_WP-myVF)*rho2
            pQ(i,j,k,3)=pQ(i,j,k,1)*water%get_e_from_p_rho(p=p,rho=rho1,y=[1.0_WP])
            pQ(i,j,k,4)=pQ(i,j,k,2)*air%get_e_from_p_rho  (p=p,rho=rho2,y=[1.0_WP])
         end do; end do; end do
         ! Face velocities at the face centres, ghosts included (w=0 is set by on_init)
         do k=lbound(pU,3),ubound(pU,3); do j=lbound(pU,2),ubound(pU,2); do i=lbound(pU,1),ubound(pU,1)
            x=solver%amr%xlo+real(i,WP)*dx; y=solver%amr%ylo+(real(j,WP)+0.5_WP)*dy; z=solver%amr%zlo+(real(k,WP)+0.5_WP)*dz
            pU(i,j,k,1)=+U0*sin(x)*cos(y)*cos(z)
         end do; end do; end do
         do k=lbound(pV,3),ubound(pV,3); do j=lbound(pV,2),ubound(pV,2); do i=lbound(pV,1),ubound(pV,1)
            x=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx; y=solver%amr%ylo+real(j,WP)*dy; z=solver%amr%zlo+(real(k,WP)+0.5_WP)*dz
            pV(i,j,k,1)=-U0*cos(x)*sin(y)*cos(z)
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine tgv_init

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
         real(WP) :: Gamma1,Pinf1,Cv1,Gamma2,Pinf2,Cv2
         ! Fluid 1 (slab) and fluid 2 (stiffened gases, ideal for Pinf=0)
         call param_read('Gamma1',Gamma1); call param_read('Pinf1',Pinf1,default=0.0_WP); call param_read('Cv1',Cv1)
         call param_read('Gamma2',Gamma2); call param_read('Pinf2',Pinf2,default=0.0_WP); call param_read('Cv2',Cv2)
         call water%initialize(gamma=Gamma1,pinf=Pinf1,cv=Cv1,q=0.0_WP,qp=0.0_WP,name='fluid1')
         call air%initialize  (gamma=Gamma2,pinf=Pinf2,cv=Cv2,q=0.0_WP,qp=0.0_WP,name='fluid2')
         ! Densities, base pressure, Mach number, and slab half-width
         call param_read('Density 1',rho1)
         call param_read('Density 2',rho2)
         call param_read('Base pressure',p0,default=1.0_WP/Gamma2)
         call param_read('Mach number',mach0)
         call param_read('Slab half-width',radius)
         if (min(rho1,rho2,p0,mach0).le.0.0_WP.or.radius.lt.0.0_WP) call die('Invalid flow parameters')
         ! Velocity amplitude from the sound speed of fluid 2
         U0=mach0*air%get_c_from_p_rho(p=p0,rho=rho2,y=[1.0_WP])
         ! Log
         call water%print(); call air%print()
         write(message,'("[Fluids] rho1=",es12.5," rho2=",es12.5," p0=",es12.5)') rho1,rho2,p0; call log(message)
         write(message,'("[Vortex] Mach=",es12.5," U0=",es12.5," slab half-width=",es12.5)') mach0,U0,radius; call log(message)
      end block init_eos_and_flow

      ! Initialize AMR grid (triply periodic box [0,2pi]^3)
      create_amrgrid: block
         use mathtools, only: twoPi
         amr%name='taylor_green'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         amr%xlo=0.0_WP; amr%xhi=twoPi
         amr%ylo=0.0_WP; amr%yhi=twoPi
         amr%zlo=0.0_WP; amr%zhi=twoPi
         amr%xper=.true.; amr%yper=.true.; amr%zper=.true.
         call param_read('Max level',amr%maxlvl)
         ! Tiling must stay off: the fused flux loops of the staggered solver write outside their tile
         amr%default_tiling=.false.
         call amr%initialize()
         ! Centre of the box: the slab is centred on x=pi
         center=0.5_WP*[amr%xlo+amr%xhi,amr%ylo+amr%yhi,amr%zlo+amr%zhi]
      end block create_amrgrid

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
         time=timetracker(amRoot=amr%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max dt',time%dtmax)
         call param_read('Max CFL',time%cflmax)
         ! Iterative midpoint sub-iterations
         call param_read('Subiterations',time%itmax,default=3)
         if (time%itmax.lt.1) call die('Require at least one sub-iteration')
         time%dt=time%dtmax
         if (restarted) then
            call io%get_scalar('dt',time%dt)
            time%t=restart_time
            time%n=restart_step
         end if
      end block initialize_timetracker

      ! Initialize compressible multiphase solver
      create_solver: block
         ! Assign materials
         fs%liq=>water; fs%gas=>air
         ! Provide pressure relaxation model
         if (air%pinf.ne.0.0_WP) call die('[simulation_init] Prelax relaxation requires Pinf2=0')
         call relax_model%initialize(gas=air,liq=water); fs%relax=>relax_model
         relax_model%model=Prelax
         ! Refuse relaxations that leave the gas below 1e-2
         fs%rhoG_min=1.0e-2_WP
         ! Positivity gate temperature floors
         call param_read('Liquid Tmin',fs%TminL,default=0.0_WP)
         call param_read('Gas Tmin'   ,fs%TminG,default=0.0_WP)
         ! Band momentum flux: CTU below the staggered phase fraction CTU alpha, MC elsewhere
         call param_read('CTU alpha',fs%alpha_ctu,default=0.5_WP)
         ! No VF treatment: build_plic and relaxation keep the transported VF
         fs%reset_vf=.false.
         ! Create flow solver
         call fs%initialize(amr=amr,name='taylor_green')
         ! Set initial conditions
         fs%user_init=>tgv_init
         ! All directions are periodic - no BCs to set
      end block create_solver

      ! Initialize workspaces
      create_workspace: block
         use amrdata_class, only: interp_none
         call dQdt%initialize   (amr,name='dQdt'   ,ncomp=fs%nQ,ng=0,interp=interp_none); call dQdt%register()
         call drhoUdt%initialize(amr,name='drhoUdt',ncomp=1,ng=0,nodal=[.true. ,.false.,.false.],interp=interp_none); call drhoUdt%register()
         call drhoVdt%initialize(amr,name='drhoVdt',ncomp=1,ng=0,nodal=[.false.,.true. ,.false.],interp=interp_none); call drhoVdt%register()
         call drhoWdt%initialize(amr,name='drhoWdt',ncomp=1,ng=0,nodal=[.false.,.false.,.true. ],interp=interp_none); call drhoWdt%register()
         call Umag%initialize(amr,name='Umag',ncomp=1,ng=0,interp=interp_none); call Umag%register()
         call Mach%initialize(amr,name='Mach',ncomp=1,ng=0,interp=interp_none); call Mach%register()
      end block create_workspace

      ! Initialize grid (no regridding for uniform mesh)
      init_grid: block
         ! Artificial bulk viscosity: off by default
         call param_read('Artificial viscosity coefficient',Cartif,default=0.0_WP)
         call param_read('Max viscous CFL',max_visc_cfl,default=0.5_WP)
         if (Cartif.lt.0.0_WP.or.max_visc_cfl.le.0.0_WP) call die('Invalid artificial viscosity parameters')
         if (restarted) then
            ! Restore grid hierarchy from checkpoint
            call amr%init_from_checkpoint(dirname=trim(restart_dir),time=time%t)
            ! Restore solver state
            call fs%restore_checkpoint(io=io,dirname=trim(restart_dir),time=time%t)
            ! Restore the monitors' reference values, which are only set at the first step
            call io%get_scalar('Etot0',Etot0); call io%get_scalar('Kin0',Kin0)
            call io%get_scalar('Sph0_1',Sph0(1)); call io%get_scalar('Sph0_2',Sph0(2))
         else
            ! Fresh start: post_regrid builds momentum from the face velocity and the primitive variables
            call amr%init_from_scratch(time=time%t)
            ! Phasic total energies from the internal contents written by user_init and the initial kinetic energy
            call fs%get_total_energy(time=time%t)
            ! Build PLIC
            call fs%build_plic(time%t)
         end if
         ! Set viscosities to zero
         call get_viscosities()
         ! Add artificial viscosity only
         call fs%add_viscartif(dt=time%dt,Cartif=Cartif,max_cfl=max_visc_cfl)
      end block init_grid

      ! Initialize checkpoint save event
      init_checkpoint: block
         ! Create checkpoint save event
         save_evt=event(time=time,name='Checkpoint')
         call param_read('Checkpoint period',save_evt%tper,default=-1.0_WP)
         ! Let solver self-register for checkpointing
         call fs%register_checkpoint(io)
         ! Add dt and the monitors' reference values to checkpoint save
         call io%add_scalar(name='dt',value=time%dt)
         call io%add_scalar(name='Etot0',value=Etot0); call io%add_scalar(name='Kin0',value=Kin0)
         call io%add_scalar(name='Sph0_1',value=Sph0(1)); call io%add_scalar(name='Sph0_2',value=Sph0(2))
      end block init_checkpoint

      ! Initialize visualization
      create_viz: block
         ! Create visualization object
         call viz%initialize(amr,'taylor_green',use_hdf5=.false.)
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
         call fs%get_info(); Etot=fs%Qint(3)+fs%Qint(4); if (time%n.eq.0) Etot0=Etot; dEtot=(Etot-Etot0)/Etot0; call get_tgv()
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         ! Create simulation monitor, with the Taylor-Green diagnostics
         mfile=monitor(amRoot=amr%amRoot,name='simulation',restart=restarted)
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(Kin,'Kinetic energy')
         call mfile%add_column(Kratio,'K/K0')
         call mfile%add_column(Sph(1),'Entropy 1')
         call mfile%add_column(Sph(2),'Entropy 2')
         call mfile%add_column(dSph(1),'dS1/|S1_0|')
         call mfile%add_column(dSph(2),'dS2/|S2_0|')
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
         call mfile%add_column(fs%TLmin,'TLmin')
         call mfile%add_column(fs%TLmax,'TLmax')
         call mfile%add_column(fs%TGmin,'TGmin')
         call mfile%add_column(fs%TGmax,'TGmax')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(amRoot=amr%amRoot,name='cfl',restart=restarted)
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
         consfile=monitor(amRoot=amr%amRoot,name='conservation',restart=restarted)
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
         gridfile=monitor(amRoot=amr%amRoot,name='grid',restart=restarted)
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
      end block create_monitors

   end subroutine simulation_init

   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time - constant timestep, so no dt adjustment
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         call time%increment()

         ! Remember old state
         call fs%store_old()

         ! Transport coefficients from the start-of-step state: inviscid, with optional artificial bulk viscosity
         call get_viscosities()
         call fs%add_viscartif(dt=time%dt,Cartif=Cartif,max_cfl=max_visc_cfl)

         ! ======================= Iterative midpoint sub-iterations =======================
         do while (time%it.le.time%itmax)

            ! Midpoint state
            call fs%prepare_midpoint()

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

         ! Visualization output
         if (viz_evt%occurs()) then
            ! Cell-centered velocity, Umag, and Mach number are only computed for visualization
            call fs%get_cell_velocity()
            call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
            call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
            call viz%write(time=time%t)
         end if

         ! Checkpoint save
         if (save_evt%occurs()) then
            save_checkpoint: block
               use string, only: rtoa
               call io%write(dirname='restart/tgv_'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
            end block save_checkpoint
         end if

         ! Perform and output monitoring
         call fs%get_info(); Etot=fs%Qint(3)+fs%Qint(4); if (time%n.eq.0) Etot0=Etot; dEtot=(Etot-Etot0)/Etot0; call get_tgv()
         call mfile%write()
         call consfile%write()
         call cflfile%write()

      end do

      ! Force a final checkpoint on exit
      final_checkpoint: block
         use string, only: rtoa
         call io%write(dirname='restart/tgv_'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
      end block final_checkpoint

   end subroutine simulation_run

   !> Taylor-Green diagnostics
   subroutine get_tgv()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM
      use parallel, only: MPI_REAL_WP
      implicit none
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), contiguous, pointer :: pQ(:,:,:,:),pRHOL(:,:,:,:),pRHOG(:,:,:,:),pPL(:,:,:,:),pPG(:,:,:,:)
      integer :: i,j,k,ierr
      Sph=0.0_WP
      call amr%mfiter_build(0,mfi)
      do while(mfi%next())
         pQ=>fs%Q%mf(0)%dataptr(mfi)
         pRHOL=>fs%RHOL%mf(0)%dataptr(mfi); pPL=>fs%PL%mf(0)%dataptr(mfi)
         pRHOG=>fs%RHOG%mf(0)%dataptr(mfi); pPG=>fs%PG%mf(0)%dataptr(mfi)
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            if (pQ(i,j,k,1).gt.0.0_WP) Sph(1)=Sph(1)+pQ(i,j,k,1)*water%cv*log((pPL(i,j,k,1)+water%pinf)/pRHOL(i,j,k,1)**water%gamma)
            if (pQ(i,j,k,2).gt.0.0_WP) Sph(2)=Sph(2)+pQ(i,j,k,2)*air%cv  *log((pPG(i,j,k,1)+air%pinf  )/pRHOG(i,j,k,1)**air%gamma  )
         end do; end do; end do
      end do
      call amr%mfiter_destroy(mfi)
      call MPI_ALLREDUCE(MPI_IN_PLACE,Sph,2,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
      Sph=Sph*amr%cell_vol(0)
      Kin=fs%rhoKint
      if (time%n.eq.0) then
         Kin0=Kin; Sph0=Sph
      end if
      Kratio=Kin/Kin0
      dSph=0.0_WP
      where (abs(Sph0).gt.0.0_WP) dSph=(Sph-Sph0)/abs(Sph0)
   end subroutine get_tgv

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Finalize time
      call time%finalize()
      ! Finalize solver
      call fs%finalize()
      call dQdt%finalize(); call drhoUdt%finalize(); call drhoVdt%finalize(); call drhoWdt%finalize()
      call Umag%finalize()
      call Mach%finalize()
      ! Finalize materials
      call water%finalize()
      call air%finalize()
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
      ! Finalize grid
      call amr%finalize()
   end subroutine simulation_final

end module simulation
