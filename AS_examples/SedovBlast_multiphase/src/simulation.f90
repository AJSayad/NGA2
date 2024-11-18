!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use mast_class,        only: mast
   use vfs_class,         only: vfs
   use matm_class,        only: matm,water !AS
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   use hypre_str_class,   only: hypre_str
   implicit none
   private

   !> Single two-phase flow solver, volume fraction solver, and material model set
   !> With corresponding time tracker
   type(mast),        public :: fs
   type(vfs),         public :: vf
   type(matm),        public :: matmod
   type(timetracker), public :: time
   type(hypre_str),   public :: ps
   type(hypre_str),   public :: vs

   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt

   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,cvgfile

   public :: simulation_init,simulation_run,simulation_final

   !> Problem definition
   real(WP) :: ddrop
   real(WP), dimension(3) :: dctr
   integer :: relax_model

contains

   !> Function that localizes the left (x-) of the domain
   function left_of_domain(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     isIn=.false.
     if (i.eq.pg%imin) isIn=.true.
   end function left_of_domain

   !> Function that localizes the right (x+) of the domain
   function right_of_domain(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     isIn=.false.
     if (i.eq.pg%imax+1) isIn=.true.
   end function right_of_domain

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         call param_read('Max steps',time%nmax)
         time%dt=time%dtmax
         time%itmax=2 
      end block initialize_timetracker


      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom, only: cube_refine_vol
         use vfs_class, only: r2p,lvira,elvira,VFhi,VFlo
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4

         real(WP) :: interface_loc,Lx,dx,gas_length !AS
         integer :: nx, spread, spread_ctr !AS

         call param_read('Gas Phase Spread',spread) !AS
         call param_read('Gas Phase Spread Center',spread_ctr)

         ! Create a VOF solver with lvira reconstruction
         ! vf=vfs(cfg=cfg,reconstruction_method=elvira,name='VOF')
         call vf%initialize(cfg=cfg,reconstruction_method=lvira,name='VOF')
         ! Initialize liquid at left
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  ! set  single phase cells
                  if ((i.le.spread_ctr + 0.5*spread).and.(i.gt.spread_ctr - 0.5*spread)) then
                     !gas
                     vf%VF(i,:,:)=0.0_WP
                  else
                     !liquid
                     vf%VF(i,:,:)=1.0_WP
                  end if
                  vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
               end do
            end do
         end do

         !AS set the phases
         call param_read('Lx',Lx)
         !print*, "Lx is read in: ",Lx
         call param_read('nx',nx)
         !print*, "nx is read in: ",nx
         dx = Lx/nx

         ! Boundary conditions on VF are built into the mast solver
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Set initial interface at the boundaries
         call vf%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call vf%polygonalize_interface()
         ! Calculate distance from polygons
         call vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call vf%subcell_vol()
         ! Calculate curvature
         call vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf%reset_volume_moments()
      end block create_and_initialize_vof

      ! Create a compressible two-phase flow solver
      create_and_initialize_flow_solver: block
         use mast_class, only: clipped_neumann,dirichlet,bc_scope,bcond,mech_egy_mech_hhz
         use hypre_str_class, only: pcg_pfmg
         use mathtools,  only: Pi
         use parallel,   only: amRoot
         integer :: i,j,k,n
         real(WP), dimension(3) :: xyz
         real(WP) :: gamm_l,Pref_l,gamm_g,visc_l,visc_g,Pref,cv_l0,cv_g0,kappa_l,kappa_g
         real(WP) :: Lrho0,LP0,Lu0,sigma
         real(WP) :: Grho0,GP0,Gu0 
         integer :: nx, spread, spread_ctr
         type(bcond), pointer :: mybc

         call param_read('Gas Phase Spread',spread) !AS
         print*, "Gas phase cell spread: ", spread
         call param_read('Gas Phase Spread Center',spread_ctr)
         print*, "Gas phase spread center (cell): ", spread_ctr

         
         ! Create material model class
         matmod=matm(cfg=cfg,name='Liquid-gas models')
         
         ! Get EOS parameters from input
         call param_read('Liquid Pref', Pref_l) 
         call param_read('Liquid gamma',gamm_l)
         call param_read('Gas gamma',gamm_g)
         call param_read('Liquid dynamic viscosity',visc_l)
         call param_read('Gas dynamic viscosity',visc_g)
         call param_read('Liquid thermal conductivity',kappa_l)
         call param_read('Gas thermal conductivity',kappa_g)
         call param_read('Liquid specific heat (constant vol)',cv_l0)
         call param_read('Gas specific heat (constant vol)',cv_g0)

         ! Register equations of state
         call matmod%register_stiffenedgas('liquid',gamm_l,Pref_l)
         call matmod%register_idealgas('gas',gamm_g)
         
         ! Create flow solver
         fs=mast(cfg=cfg,name='Two-phase All-Mach',vf=vf)

         ! Register flow solver variables with material models
         call matmod%register_thermoflow_variables('liquid',fs%Lrho,fs%Ui,fs%Vi,fs%Wi,fs%LrhoE,fs%LP)
         call matmod%register_thermoflow_variables('gas'   ,fs%Grho,fs%Ui,fs%Vi,fs%Wi,fs%GrhoE,fs%GP)
         call matmod%register_diffusion_thermo_models(viscconst_gas=visc_g, viscconst_liquid=visc_l,hdffconst_gas=kappa_g, hdffconst_liquid=kappa_l,sphtconst_gas=cv_g0,sphtconst_liquid=cv_l0)

         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs%sigma)
         
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg,nst=7)
         ps%maxlevel=10
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         
         ! Configure implicit velocity solver
         vs=hypre_str(cfg=cfg,name='Velocity',method=pcg_pfmg,nst=7)
         call param_read('Implicit iteration',vs%maxit)
         call param_read('Implicit tolerance',vs%rcvg)
         
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)

         !AS
         call param_read('Liquid density',Lrho0)
         call param_read('Liquid pressure',LP0)
         call param_read('Liquid initial velocity',Lu0)
         call param_read('Gas density',Grho0)
         call param_read('Gas pressure',GP0)
         call param_read('Gas initial velocity',Gu0)

         GP0 = GP0/(spread*cfg%dx(1))

         if (cfg%amRoot) then !AS
            print *, "===== Problem Setup Description ====="
            print *, "High Pressure Region Cell Spread: ",spread
            print *, "High Pressure Region Centered on Cell: ",spread_ctr
            print *, "Lrho0: ",Lrho0
            print *, "LP0: ",LP0
            print *, "Lu0: ",Lu0
            print *, "Grho0: ",Grho0
            print *, "GP0 (after division byspread*dx): ",GP0
            print *, "GP0: ",GP0
            print *, "Gu0: ",Gu0
         end if

         ! Liquid and gas density
         fs%Lrho = Lrho0; fs%Grho=Grho0         
         ! Initially 0 velocity
         fs%Ui = 0.0_WP; fs%Vi = 0.0_WP; fs%Wi = 0.0_WP
         ! Zero face velocities as well for the sake of dirichlet boundaries
         fs%U = 0.0_WP; fs%V = 0.0_WP; fs%W = 0.0_WP

         do i=vf%cfg%imino_,vf%cfg%imaxo_
            if ((i.le.spread_ctr + 0.5*spread).and.(i.gt.spread_ctr - 0.5*spread)) then
               !gas
               fs%GrhoE(i,:,:) = matmod%EOS_energy(GP0,Grho0,0.0_WP,0.0_WP,0.0_WP,'gas')
               fs%LrhoE(i,:,:) = matmod%EOS_energy(GP0,Lrho0,0.0_WP,0.0_WP,0.0_WP,'liquid')
            else
               !liquid  
               fs%GrhoE(i,:,:) = matmod%EOS_energy(LP0,Grho0,0.0_WP,0.0_WP,0.0_WP,'gas')
               fs%LrhoE(i,:,:) = matmod%EOS_energy(LP0,Lrho0,0.0_WP,0.0_WP,0.0_WP,'liquid')               
            end if
         end do

         fs%LP = LP0

         ! Initialize liquid energy, with surface tension
         fs%LrhoE = matmod%EOS_energy(LP0,Lrho0,0.0_WP,0.0_WP,0.0_WP,'liquid')

         ! Define boundary conditions - initialized values are intended dirichlet values too, for the cell centers
         call fs%add_bcond(name= 'inflow',type=dirichlet      ,locator=left_of_domain ,face='x',dir=-1)
         call fs%add_bcond(name='outflow',type=clipped_neumann,locator=right_of_domain,face='x',dir=+1)

         ! Calculate face velocities
         call fs%interp_vel_basic(vf,fs%Ui,fs%Vi,fs%Wi,fs%U,fs%V,fs%W)
         

         ! Apply face BC - inflow
         call fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%n_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%U(i,j,k)=0.0_WP!vshock
         end do
         ! Apply face BC - outflow
         bc_scope = 'velocity'
         call fs%apply_bcond(time%dt,bc_scope)

         ! Calculate mixture density and momenta
         fs%RHO   = (1.0_WP-vf%VF)*fs%Grho  + vf%VF*fs%Lrho
         fs%rhoUi = fs%RHO*fs%Ui; fs%rhoVi = fs%RHO*fs%Vi; fs%rhoWi = fs%RHO*fs%Wi

         ! Perform initial pressure relax
         relax_model = mech_egy_mech_hhz
         call fs%pressure_relax(vf,matmod,relax_model)
         ! Calculate initial phase and bulk moduli
         call fs%init_phase_bulkmod(vf,matmod)
         call fs%reinit_phase_pressure(vf,matmod)
         call fs%harmonize_advpressure_bulkmod(vf,matmod)

         !AS quick fix for pressure NaNs in IC
         do i=vf%cfg%imino_,vf%cfg%imaxo_
            if ((i.le.spread_ctr + 0.5*spread).and.(i.gt.spread_ctr - 0.5*spread)) then
               !gas phase
               fs%PA(i,:,:) = GP0
            else
               !liquid phase
               fs%PA(i,:,:) = LP0
            end if
         end do

         ! Set initial pressure to harmonized field based on internal energy
         fs%P = fs%PA

      end block create_and_initialize_flow_solver

      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='Sedov_Blastwave')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',fs%Ui,fs%Vi,fs%Wi)
         call ens_out%add_scalar('P',fs%P)
         call ens_out%add_scalar('PA',fs%PA)
         call ens_out%add_scalar('Grho',fs%Grho)
         call ens_out%add_scalar('Lrho',fs%Lrho)
         call ens_out%add_scalar('Density',fs%RHO)
         call ens_out%add_scalar('Bulkmod',fs%RHOSS2)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_scalar('Mach',fs%Mach)
         call ens_out%add_scalar('Vol_frac',cfg%VF)!AS
         call ens_out%add_scalar('T',fs%Tmptr) !AS
         call ens_out%add_scalar('SL_x',fs%sl_x) !AS
         call ens_out%add_scalar('SL_y',fs%sl_y) !AS
         call ens_out%add_scalar('SL_z',fs%sl_z) !AS
         call ens_out%add_scalar('LP',fs%LP) !AS
         call ens_out%add_scalar('GrhoE',fs%GrhoE) !AS
         call ens_out%add_scalar('LrhoE',fs%LrhoE) !AS
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight

      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call vf%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%RHOmin,'RHOmin')
         call mfile%add_column(fs%RHOmax,'RHOmax')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%Tmax,'Tmax')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLst,'STension CFL')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%add_column(fs%CFLa_x,'Acoustic xCFL')
         call cflfile%add_column(fs%CFLa_y,'Acoustic yCFL')
         call cflfile%add_column(fs%CFLa_z,'Acoustic zCFL')
         call cflfile%write()
         ! Create convergence monitor
         cvgfile=monitor(fs%cfg%amRoot,'cvg')
         call cvgfile%add_column(time%n,'Timestep number')
         call cvgfile%add_column(time%it,'Iteration')
         call cvgfile%add_column(time%t,'Time')
         call cvgfile%add_column(fs%impl_it_x,'Impl_x iteration')
         call cvgfile%add_column(fs%impl_rerr_x,'Impl_x error')
         call cvgfile%add_column(fs%impl_it_y,'Impl_y iteration')
         call cvgfile%add_column(fs%impl_rerr_y,'Impl_y error')
         call cvgfile%add_column(fs%implicit%it,'Impl_z iteration')
         call cvgfile%add_column(fs%implicit%rerr,'Impl_z error')
         call cvgfile%add_column(fs%psolv%it,'Pressure iteration')
         call cvgfile%add_column(fs%psolv%rerr,'Pressure error')
      end block create_monitor

   end subroutine simulation_init


   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      use messager, only: die
      implicit none

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Reinitialize phase pressure by syncing it with conserved phase energy
         call fs%reinit_phase_pressure(vf,matmod)
         fs%Uiold=fs%Ui; fs%Viold=fs%Vi; fs%Wiold=fs%Wi
         fs%RHOold = fs%RHO
         ! Remember old flow variables (phase)
         fs%Grhoold = fs%Grho; fs%Lrhoold = fs%Lrho
         fs%GrhoEold=fs%GrhoE; fs%LrhoEold=fs%LrhoE
         fs%GPold   =   fs%GP; fs%LPold   =   fs%LP

         ! Remember old interface, including VF and barycenters
         call vf%copy_interface_to_old()

         ! Create in-cell reconstruction
         call fs%flow_reconstruct(vf)

         ! Zero variables that will change during subiterations
         fs%P = 0.0_WP
         fs%Pjx = 0.0_WP; fs%Pjy = 0.0_WP; fs%Pjz = 0.0_WP
         fs%Hpjump = 0.0_WP

         ! Determine semi-Lagrangian advection flag
         call fs%flag_sl(time%dt,vf)

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

            ! Predictor step, involving advection and pressure terms
            call fs%advection_step(time%dt,vf,matmod)

            ! Viscous step
            call fs%diffusion_src_explicit_step(time%dt,vf,matmod)

            ! Prepare pressure projection
            call fs%pressureproj_prepare(time%dt,vf,matmod)
            ! Initialize and solve Helmholtz equation
            call fs%psolv%setup()
            fs%psolv%sol=fs%PA-fs%P
            call fs%psolv%solve()
            call fs%cfg%sync(fs%psolv%sol)
            ! Perform corrector step using solution
            fs%P=fs%P+fs%psolv%sol
            call fs%pressureproj_correct(time%dt,vf,fs%psolv%sol)

            ! Record convergence monitor
            call cvgfile%write()
            ! Increment sub-iteration counter
            time%it=time%it+1

         end do

         ! Pressure relaxation
         call fs%pressure_relax(vf,matmod,relax_model)

         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)

         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call fs%get_viz()
         call mfile%write()
         call cflfile%write()

      end do

   end subroutine simulation_run


   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none

   end subroutine simulation_final

end module simulation
