!> Supersonic compressible Taylor-Green vortex from Chapelier et al. (Phys. Fluids 2024)
module simulation
   use precision,           only: WP
   use amrgrid_class,       only: amrgrid
   use amrmpcomp_class, only: amrmpcomp
   use amrviz_class,        only: amrviz
   use amrdata_class,       only: amrdata
   use timetracker_class,   only: timetracker
   use event_class,         only: event
   use monitor_class,       only: monitor
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

   !> Molecular shear viscosity (fs%visc also holds the Vreman part), kept for the monitors
   type(amrdata) :: viscm

   !> Visualization
   type(event) :: viz_evt
   type(amrviz) :: viz

   !> Simulation monitoring
   type(monitor) :: mfile,consfile,cflfile

   !> Monitors of the reference: integrals of volume, mass, internal and kinetic energy, their rates, the resolved and
   !> subgrid solenoidal and dilatational dissipations, and the pressure dilatation
   real(WP) :: vol_total,mass,IE,KE,IEold,KEold,dkdt,dedt
   real(WP) :: sdiss,ddiss,ddiss_sgsb,sdiss_sgss,ddiss_sgss,pdil

   !> Total energy (sum of the phasic total energies) and relative drifts of total energy and mass
   real(WP) :: Etot,Etot0,dEtot,mass0,dmass

   !> Materials: the gas, and the unused liquid
   type(nasg),          target :: liq   !< Stiffened gas: NASG with zero co-volume
   type(stiffened_gas), target :: gas

   !> Relaxation model (inactive without mixed cells)
   type(relax_ig_nasg), target :: relax_model

   !> Equation of state and flow conditions
   real(WP) :: Gamma,Mach0,Reynolds,Prandtl,bulk2shear

   !> Subgrid-scale modeling: artificial bulk viscosity (bSGS) and Vreman eddy viscosity (sSGS)
   logical  :: use_sgsb,use_sgss
   real(WP) :: Cartif,max_visc_cfl

contains

   !> Transport coefficients from the gas temperature, as the reference: Sutherland's law for the shear viscosity, bulk
   !> viscosity from the bulk-to-shear ratio, the optional subgrid models, then the gas heat conductivity cp*visc/Pr (visc
   !> with its Vreman part, as the reference).
   subroutine get_viscosities()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      implicit none
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pTG,pVisc,pViscm,pBeta
      integer :: lvl,i,j,k
      do lvl=0,amr%clvl()
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pTG   =>fs%TG%mf(lvl)%dataptr(mfi)
            pVisc =>fs%visc%mf(lvl)%dataptr(mfi)
            pViscm=>viscm%mf(lvl)%dataptr(mfi)
            pBeta =>fs%beta%mf(lvl)%dataptr(mfi)
            ! The gas temperature holds the viscosity ghost layer, so no ghost fill is needed
            bx=mfi%growntilebox(fs%visc%ng)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pViscm(i,j,k,1)=Reynolds**(-1.0_WP)*(1.4042_WP*pTG(i,j,k,1)**1.5_WP)/(pTG(i,j,k,1)+0.4042_WP)
               pVisc(i,j,k,1)=pViscm(i,j,k,1)
               pBeta(i,j,k,1)=bulk2shear*pViscm(i,j,k,1)
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
      ! Subgrid-scale models
      if (use_sgsb) call fs%add_viscartif(dt=time%dt,Cartif=Cartif,max_cfl=max_visc_cfl)
      if (use_sgss) call fs%add_vreman(dt=time%dt)
      ! Heat conductivities
      call fs%diffG%copy(src=fs%visc); call fs%diffG%mult(val=gas%cv*gas%gamma/Prandtl)
      call fs%diffL%setval(val=0.0_WP)
   end subroutine get_viscosities

   !> Monitor conservation properties
   subroutine analyze_conservation()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM
      use parallel, only: MPI_REAL_WP
      implicit none
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW,pP,pVisc,pViscm,pBeta
      real(WP), dimension(6) :: buf
      real(WP) :: dxi,dyi,dzi,div,wz2,wx2,wy2
      integer :: i,j,k,ierr
      ! Integrals and their rates
      mass=fs%Qint(1)+fs%Qint(2)
      KE=fs%rhoKint
      Etot=fs%Qint(3)+fs%Qint(4)
      IE=Etot-KE
      if (time%n.eq.0) then
         IEold=IE; KEold=KE; Etot0=Etot; mass0=mass
      end if
      dkdt=(KE-KEold)/time%dt
      dedt=(IE-IEold)/time%dt
      dEtot=(Etot-Etot0)/Etot0
      dmass=(mass-mass0)/mass0
      ! Dissipations and pressure dilatation
      dxi=1.0_WP/amr%dx(0); dyi=1.0_WP/amr%dy(0); dzi=1.0_WP/amr%dz(0)
      buf=0.0_WP
      call amr%mfiter_build(0,mfi)
      do while (mfi%next())
         pU=>fs%U%mf(0)%dataptr(mfi); pV=>fs%V%mf(0)%dataptr(mfi); pW=>fs%W%mf(0)%dataptr(mfi)
         pP=>fs%P%mf(0)%dataptr(mfi); pBeta=>fs%beta%mf(0)%dataptr(mfi)
         pVisc=>fs%visc%mf(0)%dataptr(mfi); pViscm=>viscm%mf(0)%dataptr(mfi)
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Cell-centered divergence: total, artificial bulk, and Vreman dilatational dissipation, and pressure dilatation
            div=dxi*(pU(i+1,j,k,1)-pU(i,j,k,1))+dyi*(pV(i,j+1,k,1)-pV(i,j,k,1))+dzi*(pW(i,j,k+1,1)-pW(i,j,k,1))
            buf(2)=buf(2)+(pBeta(i,j,k,1)+4.0_WP/3.0_WP*pVisc(i,j,k,1))*div**2
            buf(3)=buf(3)+(pBeta(i,j,k,1)-bulk2shear*pViscm(i,j,k,1))*div**2
            buf(5)=buf(5)+4.0_WP/3.0_WP*(pVisc(i,j,k,1)-pViscm(i,j,k,1))*div**2
            buf(6)=buf(6)+pP(i,j,k,1)*div
            ! Squared vorticity on the three edges whose lower corner is (i,j,k), each edge counted once: total and Vreman
            ! solenoidal dissipation
            wz2=(dxi*(pV(i,j,k,1)-pV(i-1,j,k,1))-dyi*(pU(i,j,k,1)-pU(i,j-1,k,1)))**2
            wx2=(dyi*(pW(i,j,k,1)-pW(i,j-1,k,1))-dzi*(pV(i,j,k,1)-pV(i,j,k-1,1)))**2
            wy2=(dzi*(pU(i,j,k,1)-pU(i,j,k-1,1))-dxi*(pW(i,j,k,1)-pW(i-1,j,k,1)))**2
            buf(1)=buf(1)+0.25_WP*sum(pVisc(i-1:i,j-1:j,k,1))*wz2+0.25_WP*sum(pVisc(i,j-1:j,k-1:k,1))*wx2+0.25_WP*sum(pVisc(i-1:i,j,k-1:k,1))*wy2
            buf(4)=buf(4)+0.25_WP*sum(pVisc(i-1:i,j-1:j,k,1)-pViscm(i-1:i,j-1:j,k,1))*wz2 &
            &            +0.25_WP*sum(pVisc(i,j-1:j,k-1:k,1)-pViscm(i,j-1:j,k-1:k,1))*wx2 &
            &            +0.25_WP*sum(pVisc(i-1:i,j,k-1:k,1)-pViscm(i-1:i,j,k-1:k,1))*wy2
         end do; end do; end do
      end do
      call amr%mfiter_destroy(mfi)
      buf=buf*amr%cell_vol(0)
      call MPI_ALLREDUCE(MPI_IN_PLACE,buf,6,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
      sdiss_sgss=buf(4); ddiss_sgsb=buf(3); ddiss_sgss=buf(5); pdil=buf(6)
      ! Resolved dissipations
      sdiss=buf(1)-sdiss_sgss
      ddiss=buf(2)-ddiss_sgsb-ddiss_sgss
   end subroutine analyze_conservation

   !> User init callback - set Q, VF, barycenters, and face velocities of the Taylor-Green vortex
   subroutine tgv_init(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box
      use amrex_amr_module, only: amrex_mfiter_build,amrex_mfiter_destroy
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF,pCL,pCG,pU,pV
      real(WP), dimension(3) :: xc
      real(WP) :: dx,dy,dz,E0,rho,x,y,z
      integer :: i,j,k
      ! Get mesh size and the uniform internal energy
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl); dz=solver%amr%dz(lvl)
      E0=1.0_WP/((Gamma-1.0_WP)*Gamma*Mach0**2)
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
         ! Cells, ghosts included: gas only, with Q=(0,rho,0,rho*E0)
         bx=mfi%growntilebox(solver%nover)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            xc=[solver%amr%xlo+(real(i,WP)+0.5_WP)*dx,solver%amr%ylo+(real(j,WP)+0.5_WP)*dy,solver%amr%zlo+(real(k,WP)+0.5_WP)*dz]
            pVF(i,j,k,1)=0.0_WP
            if (lvl.eq.solver%amr%maxlvl) then
               pCL(i,j,k,:)=xc
               pCG(i,j,k,:)=xc
            end if
            ! Density from the pressure and the uniform internal energy
            rho=(1.0_WP/(Gamma*Mach0**2)+(cos(2.0_WP*xc(1))+cos(2.0_WP*xc(2)))*(cos(2.0_WP*xc(3))+2.0_WP)/16.0_WP)/(E0*(Gamma-1.0_WP))
            pQ(i,j,k,1)=0.0_WP
            pQ(i,j,k,2)=rho
            pQ(i,j,k,3)=0.0_WP
            pQ(i,j,k,4)=rho*E0
         end do; end do; end do
         ! Face velocities at the face centres, ghosts included (w=0 is set by on_init)
         do k=lbound(pU,3),ubound(pU,3); do j=lbound(pU,2),ubound(pU,2); do i=lbound(pU,1),ubound(pU,1)
            x=solver%amr%xlo+real(i,WP)*dx; y=solver%amr%ylo+(real(j,WP)+0.5_WP)*dy; z=solver%amr%zlo+(real(k,WP)+0.5_WP)*dz
            pU(i,j,k,1)=+sin(x)*cos(y)*cos(z)
         end do; end do; end do
         do k=lbound(pV,3),ubound(pV,3); do j=lbound(pV,2),ubound(pV,2); do i=lbound(pV,1),ubound(pV,1)
            x=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx; y=solver%amr%ylo+real(j,WP)*dy; z=solver%amr%zlo+(real(k,WP)+0.5_WP)*dz
            pV(i,j,k,1)=-cos(x)*sin(y)*cos(z)
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine tgv_init

   !> Initialization of problem solver
   subroutine simulation_init
      use param,    only: param_read
      use messager, only: die
      implicit none

      ! Read EoS and flow conditions
      init_eos_and_flow: block
         call param_read('Gamma'   ,Gamma   )
         call param_read('Mach'    ,Mach0   )
         call param_read('Reynolds',Reynolds)
         call param_read('Prandtl' ,Prandtl )
         call param_read('Bulk to shear ratio',bulk2shear,default=0.0_WP)
         ! Ideal gas with cv such that T=1, and the same for the unused liquid
         call gas%initialize(gamma=Gamma,pinf=0.0_WP,cv=1.0_WP/(Gamma*(Gamma-1.0_WP)*Mach0**2),q=0.0_WP,qp=0.0_WP,name='gas')
         call liq%initialize(gamma=Gamma,pinf=0.0_WP,cv=1.0_WP/(Gamma*(Gamma-1.0_WP)*Mach0**2),q=0.0_WP,qp=0.0_WP,name='liquid (unused)')
         call gas%print()
         ! Subgrid-scale modeling: both off by default, as the reference
         call param_read('Use bSGS model',use_sgsb,default=.false.)
         call param_read('Use sSGS model',use_sgss,default=.false.)
         call param_read('Artificial viscosity coefficient',Cartif,default=5.0_WP)
         call param_read('Max viscous CFL',max_visc_cfl,default=0.5_WP)
         if (Cartif.lt.0.0_WP.or.max_visc_cfl.le.0.0_WP) call die('Invalid artificial viscosity parameters')
      end block init_eos_and_flow

      ! Initialize AMR grid (triply periodic box [-pi,pi]^3, uniform)
      create_amrgrid: block
         use mathtools, only: twoPi
         amr%name='supersonic_taylor_green'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         call param_read('Max level',amr%maxlvl,default=0)
         call param_read('Max grid size'  ,amr%nmax ,default=32)
         call param_read('Blocking factor',amr%nbloc,default=8 )
         if (amr%maxlvl.ne.0) call die('This case requires Max level=0')
         amr%xlo=-0.5_WP*twoPi; amr%xhi=+0.5_WP*twoPi
         amr%ylo=-0.5_WP*twoPi; amr%yhi=+0.5_WP*twoPi
         amr%zlo=-0.5_WP*twoPi; amr%zhi=+0.5_WP*twoPi
         amr%xper=.true.; amr%yper=.true.; amr%zper=.true.
         ! Tiling must stay off: the fused flux loops of the staggered solver write outside their tile
         amr%default_tiling=.false.
         call amr%initialize()
         vol_total=(amr%xhi-amr%xlo)*(amr%yhi-amr%ylo)*(amr%zhi-amr%zlo)
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

      ! Initialize compressible multiphase solver: without an interface, the band setting (CTU alpha) keeps the solver default
      create_solver: block
         ! Assign materials
         fs%liq=>liq; fs%gas=>gas
         ! Provide pressure relaxation model (not used for this case)
         call relax_model%initialize(gas=gas,liq=liq); fs%relax=>relax_model
         relax_model%model=Prelax
         ! Refuse relaxations that leave the gas below 1e-2; the model's own RHOGmin keeps its 1e-2 default
         fs%rhoG_min=1.0e-2_WP
         ! Positivity gate temperature floors (0: undeclared, a pure cell outside its law's domain stops the run)
         call param_read('Liquid Tmin',fs%TminL,default=0.0_WP)
         call param_read('Gas Tmin'   ,fs%TminG,default=0.0_WP)
         ! No VF treatment: build_plic and relaxation keep the transported VF
         fs%reset_vf=.false.
         ! Create flow solver
         call fs%initialize(amr=amr,name='supersonic_taylor_green')
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
         call viscm%initialize(amr,name='viscm',ncomp=1,ng=fs%visc%ng,interp=interp_none); call viscm%register()
      end block create_workspace

      ! Initialize grid (no regridding for uniform mesh)
      init_grid: block
         ! Fresh start: post_regrid builds momentum from the face velocity and the primitive variables
         call amr%init_from_scratch(time=time%t)
         ! Phasic total energies from the internal contents written by user_init and the initial kinetic energy
         call fs%get_total_energy(time=time%t)
         ! Build PLIC
         call fs%build_plic(time%t)
         ! Transport coefficients
         call get_viscosities()
      end block init_grid

      ! Initialize visualization
      create_viz: block
         ! Create visualization object
         call viz%initialize(amr,'supersonic_taylor_green',use_hdf5=.false.)
         call viz%add_scalar(fs%RHOG,1,'density')
         call viz%add_scalar(fs%P,1,'pressure')
         call viz%add_scalar(fs%UVW,1,'U')
         call viz%add_scalar(fs%UVW,2,'V')
         call viz%add_scalar(fs%UVW,3,'W')
         call viz%add_scalar(fs%IG,1,'energy')
         call viz%add_scalar(fs%TG,1,'temperature')
         call viz%add_scalar(Mach,1,'Mach')
         call viz%add_scalar(fs%visc,1,'visc')
         call viz%add_scalar(fs%beta,1,'beta')
         ! Create visualization output event
         viz_evt=event(time=time,name='Visualization output')
         call param_read('Output period',viz_evt%tper)
         ! Write initial state
         if (viz_evt%occurs()) call write_viz()
      end block create_viz

      ! Create monitors
      create_monitors: block
         ! Get solver info, cfl, and the reference monitors
         call fs%get_info()
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         call analyze_conservation()
         ! Create simulation monitor
         mfile=monitor(amRoot=amr%amRoot,name='simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%PGmax,'Pmax')
         call mfile%add_column(fs%PGmin,'Pmin')
         call mfile%add_column(fs%IGmax,'Emax')
         call mfile%add_column(fs%IGmin,'Emin')
         call mfile%add_column(fs%RHOGmax,'RHOmax')
         call mfile%add_column(fs%RHOGmin,'RHOmin')
         call mfile%add_column(fs%TGmax,'Tmax')
         call mfile%add_column(fs%TGmin,'Tmin')
         call mfile%add_column(fs%VFmax,'VFmax')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(amRoot=amr%amRoot,name='cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLa_x,'Acoustic xCFL')
         call cflfile%add_column(fs%CFLa_y,'Acoustic yCFL')
         call cflfile%add_column(fs%CFLa_z,'Acoustic zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
         ! Create conservation monitor: the reference's columns, then the conservation checks
         consfile=monitor(amRoot=amr%amRoot,name='conservation')
         call consfile%add_column(time%n,'Timestep number')
         call consfile%add_column(time%t,'Time')
         call consfile%add_column(vol_total,'Volume')
         call consfile%add_column(mass,'Mass')
         call consfile%add_column(IE,'Internal energy')
         call consfile%add_column(KE,'Kinetic energy')
         call consfile%add_column(dkdt,'dkdt')
         call consfile%add_column(dedt,'dedt')
         call consfile%add_column(sdiss,'sdiss')
         call consfile%add_column(ddiss,'ddiss')
         call consfile%add_column(ddiss_sgsb,'ddiss_sgsb')
         call consfile%add_column(sdiss_sgss,'sdiss_sgss')
         call consfile%add_column(ddiss_sgss,'ddiss_sgss')
         call consfile%add_column(pdil,'pdil')
         call consfile%add_column(Etot,'Total energy')
         call consfile%add_column(dEtot,'dE/E0')
         call consfile%add_column(dmass,'dM/M0')
         call consfile%add_column(fs%rhoUint,'U Momentum')
         call consfile%add_column(fs%rhoVint,'V Momentum')
         call consfile%add_column(fs%rhoWint,'W Momentum')
         call consfile%write()
      end block create_monitors

   end subroutine simulation_init

   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time, ending exactly at Max time
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         call time%adjust_dt()
         time%dt=min(time%dt,time%tmax-time%t)
         call time%increment()

         ! Remember old state and energy integrals
         call fs%store_old()
         IEold=IE; KEold=KE

         ! Transport coefficients from the start-of-step state
         call get_viscosities()

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
         if (viz_evt%occurs()) call write_viz()

         ! Perform and output monitoring
         call fs%get_info()
         call analyze_conservation()
         call mfile%write()
         call consfile%write()
         call cflfile%write()

      end do

   end subroutine simulation_run

   !> Write a plotfile; the cell-centered velocity and the Mach number are only computed for it
   subroutine write_viz()
      implicit none
      call fs%get_cell_velocity()
      call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
      call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
      call viz%write(time=time%t)
   end subroutine write_viz

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Finalize time
      call time%finalize()
      ! Finalize solver and workspaces
      call fs%finalize()
      call dQdt%finalize(); call drhoUdt%finalize(); call drhoVdt%finalize(); call drhoWdt%finalize()
      call Umag%finalize(); call Mach%finalize(); call viscm%finalize()
      ! Finalize materials
      call liq%finalize()
      call gas%finalize()
      ! Finalize visualization
      call viz%finalize()
      call viz_evt%finalize()
      ! Finalize monitoring
      call mfile%finalize()
      call cflfile%finalize()
      call consfile%finalize()
      ! Finalize grid
      call amr%finalize()
   end subroutine simulation_final

end module simulation
