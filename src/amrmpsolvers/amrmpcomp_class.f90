!> AMR staggered compressible multiphase solver class
module amrmpcomp_class
   use iso_c_binding,     only: c_ptr,c_loc,c_f_pointer
   use precision,         only: WP
   use amrdata_class,     only: amrdata
   use amrmpflow_class,   only: amrmpflow
   use amrvof_class,      only: VFlo,VFhi,vol_eps,BC_LIQ,BC_GAS,BC_REFLECT,BC_USER
   use amrex_amr_module,  only: amrex_box,amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_multifab
   use material_class,    only: material
   use thermorelax_class, only: thermorelax,RELAX_OK
   implicit none
   private

   ! Expose type and constants
   public :: amrmpcomp,VFlo,VFhi,BC_LIQ,BC_GAS,BC_REFLECT,BC_USER

   ! Semi-Lagrangian face flux components: phasic Q (1:4) and CTU donor moments (5:16) sum(m*u_a(b)) and
   ! sum(m*u_a(b+e_a)) per phase and momentum axis
   integer, parameter :: nSLF=16

   ! Unit vectors
   integer, dimension(3,3), parameter :: eu=reshape([1,0,0,0,1,0,0,0,1],[3,3])

   !> Pointer to one FAB, for direction-indexed access
   type :: fab_ptr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p=>null()
   end type fab_ptr

   !> Diagnostic counts and bookings of the last get_dQdt, get_energy, and apply_relax calls
   type :: diag_counts
      real(WP) :: nvf_guard=0.0_WP,nsliver=0.0_WP,nrelax_fail=0.0_WP,ngate=0.0_WP
      real(WP) :: gate_mass=0.0_WP,gate_energy=0.0_WP,gate_rhoU=0.0_WP,gate_rhoV=0.0_WP,gate_rhoW=0.0_WP
   end type diag_counts

   !> AMR staggered compressible multiphase solver type
   type, extends(amrmpflow) :: amrmpcomp

      ! User-configurable callbacks
      procedure(stag_init_iface),    pointer, pass :: user_init   =>null()  !< User-defined initialization
      procedure(stag_tagging_iface), pointer, pass :: user_tagging=>null()  !< User-defined tagging
      procedure(stag_bc_iface),      pointer, pass :: user_bc     =>null()  !< User-defined boundary conditions
      procedure(stag_vofbc_iface),   pointer, pass :: user_vofbc  =>null()  !< User-defined VOF boundary conditions

      ! Liquid and gas materials
      class(material), pointer :: liq=>null()
      class(material), pointer :: gas=>null()

      ! Thermodynamic relaxation model
      class(thermorelax), pointer :: relax=>null()

      ! Surface tension coefficient
      real(WP) :: sigma=0.0_WP

      ! Staggered mixture momentum
      type(amrdata) :: rhoU,rhoUold
      type(amrdata) :: rhoV,rhoVold
      type(amrdata) :: rhoW,rhoWold

      ! Midpoint face velocities
      type(amrdata) :: Umid,Vmid,Wmid

      ! Cell-centered mixture velocity
      type(amrdata) :: UVW

      ! Mixture pressure
      type(amrdata) :: P

      ! Mixture speed of sound
      type(amrdata) :: C

      ! Current phasic primitive variables
      type(amrdata) :: RHOL,RHOLold,RHOG,RHOGold         !< Phasic densities
      type(amrdata) :: IL,IG                             !< Phasic internal energies
      type(amrdata) :: PL,PG,PLold,PGold                 !< Phasic pressures
      type(amrdata) :: TL,TG                             !< Phasic temperatures

      ! Physical properties
      type(amrdata) :: visc                              !< Dynamic viscosity
      type(amrdata) :: beta                              !< Bulk viscosity
      type(amrdata) :: diffL,diffG                       !< Phasic heat diffusivities

      ! Fluxes for get_dQdt
      type(amrex_multifab), dimension(:), allocatable :: Fx,Fy,Fz     !< Phasic mass and internal energy fluxes (faces)
      type(amrex_multifab), dimension(:), allocatable :: FUx,FVy,FWz  !< Diagonal momentum fluxes (cell centers)
      type(amrex_multifab), dimension(:), allocatable :: FUy,FVx      !< Momentum cross-fluxes (xy-edges)
      type(amrex_multifab), dimension(:), allocatable :: FUz,FWx      !< Momentum cross-fluxes (xz-edges)
      type(amrex_multifab), dimension(:), allocatable :: FVz,FWy      !< Momentum cross-fluxes (yz-edges)
      type(amrex_multifab), dimension(:), allocatable :: Ex,Ey,Ez     !< Mixture work fluxes of the energy pass (faces)
      type(amrex_multifab) :: band                                    !< Semi-Lagrangian band (finest level)
      type(amrex_multifab), dimension(3) :: SLF                       !< Band face fluxes of phasic Q and CTU donor moments (finest level)
      type(amrex_multifab), dimension(3) :: SLV                       !< Band face volume moment fluxes (finest level)
      type(amrex_multifab) :: curvold,SDold                           !< Start-of-step curvature and surface density (finest level)
      type(amrex_multifab), dimension(3) :: Pj                        !< Surface tension pressure jump to add to -grad(P) on faces (finest level)

      ! Band momentum flux and sliver parameters
      real(WP) :: alpha_ctu=0.5_WP       !< Staggered phase fraction below which CTU is used
      real(WP) :: sliver_frac=1.0e-2_WP  !< Transported gas fraction below which a mixed cell is equilibrated before the energy pass

      ! CFL numbers
      real(WP) :: CFLp=0.0_WP                                !< Pressure+convection
      real(WP) :: CFLa_x=0.0_WP,CFLa_y=0.0_WP,CFLa_z=0.0_WP  !< Acoustic
      real(WP) :: CFLv_x=0.0_WP,CFLv_y=0.0_WP,CFLv_z=0.0_WP  !< Viscous
      real(WP) :: CFLst=0.0_WP                               !< Capillary

      ! Monitoring quantities
      real(WP) :: RHOLmin=0.0_WP,RHOLmax=0.0_WP,RHOGmin=0.0_WP,RHOGmax=0.0_WP
      real(WP) :: ILmin=0.0_WP,ILmax=0.0_WP,IGmin=0.0_WP,IGmax=0.0_WP
      real(WP) :: PLmin=0.0_WP,PLmax=0.0_WP,PGmin=0.0_WP,PGmax=0.0_WP
      real(WP) :: TLmin=0.0_WP,TLmax=0.0_WP,TGmin=0.0_WP,TGmax=0.0_WP
      real(WP) :: Pmin=0.0_WP,Pmax=0.0_WP
      real(WP) :: Cmin=0.0_WP,Cmax=0.0_WP
      real(WP) :: dPmax=0.0_WP
      real(WP) :: rhoUint=0.0_WP,rhoVint=0.0_WP,rhoWint=0.0_WP
      real(WP) :: rhoKint=0.0_WP
      real(WP) :: nvf_guard=0.0_WP               !< Band cells whose transported phase volume sum vanished in the last get_dQdt call (given the majority phase)
      real(WP) :: nrelax_fail=0.0_WP             !< Mixed cells whose relaxation did not return RELAX_OK in the last call
      real(WP) :: nsliver=0.0_WP                 !< Mixed cells equilibrated as slivers in the last get_energy call
      real(WP) :: ngate=0.0_WP                   !< Mixed cells that lost a phase to the positivity gate in the last apply_relax call
      real(WP) :: gate_mass=0.0_WP,gate_energy=0.0_WP                 !< Mass and total energy the gate removed in the last apply_relax call
      real(WP) :: gate_rhoU=0.0_WP,gate_rhoV=0.0_WP,gate_rhoW=0.0_WP   !< Momentum change of the faces of gated cells (velocity kept) in that call
      type(diag_counts) :: cnt                   !< This rank's values of the counts above, which get_info sums across ranks

      ! Minimum density for stability
      real(WP) :: rho_floor=1.0e-10_WP

      ! Timing: per-rank accumulators (seconds), reduced to min/max across ranks and reset in get_info
      real(WP) :: wt_mid   =0.0_WP      !< Prepare_midpoint
      real(WP) :: wt_prim  =0.0_WP      !< Get_primitive cell loops (excluding face velocity)
      real(WP) :: wt_fvel  =0.0_WP      !< Get_face_velocity: momentum sync, face velocity, ghost fill, momentum BC
      real(WP) :: wt_dQdt  =0.0_WP      !< Full get_dQdt
      real(WP) :: wt_alloc =0.0_WP      !< Get_dQdt: band build (flux workspaces are built at regrid)
      real(WP) :: wt_slflux=0.0_WP      !< Get_dQdt: geometric band fluxes (flux polyhedra, tet cutting, sync)
      real(WP) :: wt_vmom  =0.0_WP      !< Get_dQdt: band VF/barycenter update
      real(WP) :: wt_fv    =0.0_WP      !< Get_dQdt: bulk and momentum flux loops
      real(WP) :: wt_div   =0.0_WP      !< Get_dQdt: C/F averaging, divergence, and sources
      real(WP) :: wt_relax =0.0_WP      !< Apply_relax
      real(WP) :: wt_visc  =0.0_WP      !< Viscosity models (artificial bulk viscosity and Vreman)
      real(WP) :: wtmax_mid   =0.0_WP, wtmin_mid   =0.0_WP
      real(WP) :: wtmax_prim  =0.0_WP, wtmin_prim  =0.0_WP
      real(WP) :: wtmax_fvel  =0.0_WP, wtmin_fvel  =0.0_WP
      real(WP) :: wtmax_dQdt  =0.0_WP, wtmin_dQdt  =0.0_WP
      real(WP) :: wtmax_alloc =0.0_WP, wtmin_alloc =0.0_WP
      real(WP) :: wtmax_slflux=0.0_WP, wtmin_slflux=0.0_WP
      real(WP) :: wtmax_vmom  =0.0_WP, wtmin_vmom  =0.0_WP
      real(WP) :: wtmax_fv    =0.0_WP, wtmin_fv    =0.0_WP
      real(WP) :: wtmax_div   =0.0_WP, wtmin_div   =0.0_WP
      real(WP) :: wtmax_relax =0.0_WP, wtmin_relax =0.0_WP
      real(WP) :: wtmax_visc  =0.0_WP, wtmin_visc  =0.0_WP

   contains
      ! Type-bound constructor/destructor
      procedure :: initialize
      procedure :: finalize
      ! Lifecycle callbacks
      procedure :: on_init
      procedure :: on_coarse
      procedure :: on_remake
      procedure :: on_clear
      procedure :: post_regrid
      ! BC overrides
      procedure :: apply_velbc=>stag_apply_velbc
      procedure :: apply_Qbc=>stag_apply_Qbc
      procedure :: apply_vofbc=>stag_apply_vofbc
      ! Utilities
      procedure :: store_old                 !< Store Q, momentum, velocity, VOF, and phasic primitives at the start of the step
      procedure :: prepare_midpoint          !< Set face velocities and phasic primitives to midpoint time
      procedure :: add_surface_tension_jump  !< Set the surface tension pressure jump on finest-level faces
      procedure :: get_face_velocity         !< Update face velocities from momentum and fill staggered ghosts
      procedure :: get_cell_velocity         !< Interpolate face velocities to cell centers in UVW
      procedure :: sync_momentum             !< Reconcile shared faces and restrict evolved momentum
      procedure :: apply_momentum_bc         !< Recover physical boundary momentum from filled Q and velocity
      procedure :: get_momentum              !< Set momentum from face velocity on initialization and regrid only
      procedure :: get_total_energy          !< Set the phasic total energies on initialization only
      ! Physics procedures
      procedure :: get_primitive             !< Get phasic and mixture primitive variables from conserved variables
      procedure :: get_dQdt                  !< Compute dQdt and drhoUVWdt over a step
      procedure :: get_energy                !< Set the phasic total energies of the new iterate from the transported internal energies, the phasic work, and the mixture work flux
      procedure :: apply_relax               !< Apply mechanical/thermal relaxation in mixed cells
      procedure :: add_viscartif             !< Add localized artificial bulk viscosity
      procedure :: add_vreman                !< Add Vreman SGS eddy viscosity
      procedure :: get_cfl                   !< Compute CFL numbers
      ! Print solver info
      procedure :: get_info
      procedure :: print=>amrmpcomp_print
      ! Checkpoint I/O
      procedure :: register_checkpoint
      procedure :: restore_checkpoint
   end type amrmpcomp

   !> Abstract interface for user-provided on_init callback
   abstract interface
      subroutine stag_init_iface(solver,lvl,time,ba,dm)
         import :: amrmpcomp,WP,amrex_boxarray,amrex_distromap
         class(amrmpcomp), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(amrex_boxarray), intent(in) :: ba
         type(amrex_distromap), intent(in) :: dm
      end subroutine stag_init_iface
   end interface

   !> Abstract interface for user-provided tagging callback
   abstract interface
      subroutine stag_tagging_iface(solver,lvl,time,tags)
         import :: amrmpcomp,c_ptr,WP
         class(amrmpcomp), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(c_ptr), intent(in) :: tags
      end subroutine stag_tagging_iface
   end interface

   !> Abstract interface for user-provided BC callback
   abstract interface
      subroutine stag_bc_iface(solver,lvl,time,face,bx,comp,p)
         import :: amrmpcomp,amrex_box,WP
         class(amrmpcomp), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         integer, intent(in) :: face                       !< 1=xlo,2=xhi,3=ylo,4=yhi,5=zlo,6=zhi
         type(amrex_box), intent(in) :: bx                 !< Boundary box to fill
         character(len=1), intent(in) :: comp              !< Can be 'U','V','W','Q'
         real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      end subroutine stag_bc_iface
   end interface

   !> Abstract interface for user-provided VOF BC callback
   abstract interface
      subroutine stag_vofbc_iface(solver,lvl,time,face,bx,pVF,pCL,pCG,pPLIC)
         import :: amrmpcomp,amrex_box,WP
         class(amrmpcomp), intent(inout) :: solver
         integer, intent(in) :: lvl,face
         real(WP), intent(in) :: time
         type(amrex_box), intent(in) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCL,pCG,pPLIC
      end subroutine stag_vofbc_iface
   end interface

contains

   !> Component of the CTU donor moment of phase (1=liquid, 2=gas), side (1: velocity at the lower
   !> a-face of the donor, 2: at its upper a-face), and momentum axis a in the SL face workspace
   pure integer function ictu(phase,side,a)
      implicit none
      integer, intent(in) :: phase,side,a
      ictu=4+6*(phase-1)+3*(side-1)+a
   end function ictu

   ! ============================================================================
   ! DISPATCHERS (module-level) - recover concrete amrmpcomp type
   ! ============================================================================

   !> Dispatch on_init: calls type-bound method then user callback
   subroutine amrmpcomp_on_init(ctx,lvl,time,ba,dm)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrmpcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_init(lvl,time,ba,dm)
      if (associated(this%user_init)) call this%user_init(lvl,time,ba,dm)
   end subroutine amrmpcomp_on_init

   !> Dispatch on_coarse: calls type-bound method
   subroutine amrmpcomp_on_coarse(ctx,lvl,time,ba,dm)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrmpcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_coarse(lvl,time,ba,dm)
   end subroutine amrmpcomp_on_coarse

   !> Dispatch on_remake: calls type-bound method
   subroutine amrmpcomp_on_remake(ctx,lvl,time,ba,dm)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrmpcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_remake(lvl,time,ba,dm)
   end subroutine amrmpcomp_on_remake

   !> Dispatch on_clear: calls type-bound method
   subroutine amrmpcomp_on_clear(ctx,lvl)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      type(amrmpcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%on_clear(lvl)
   end subroutine amrmpcomp_on_clear

   !> Dispatch post_regrid: calls type-bound method
   subroutine amrmpcomp_postregrid(ctx,lbase,time)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      type(amrmpcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%post_regrid(lbase,time)
   end subroutine amrmpcomp_postregrid

   !> Dispatch tagging: interface band tagging then user callback
   subroutine amrmpcomp_tagging(ctx,lvl,time,tags)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags
      type(amrmpcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%tagging(lvl,time,tags)
      if (associated(this%user_tagging)) call this%user_tagging(lvl,time,tags)
   end subroutine amrmpcomp_tagging

   !> Dispatch cost: calls type-bound method
   subroutine amrmpcomp_get_cost(ctx,lvl,nboxes,costs,ba)
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl,nboxes
      real(WP), intent(inout) :: costs(nboxes)
      type(amrex_boxarray), intent(in) :: ba
      type(amrmpcomp), pointer :: this
      call c_f_pointer(ctx,this)
      call this%get_cost(lvl,nboxes,costs,ba)
   end subroutine amrmpcomp_get_cost

   ! ============================================================================
   ! INITIALIZATION / FINALIZATION
   ! ============================================================================

   !> Initialize the staggered compressible multiphase solver
   subroutine initialize(this,amr,name)
      use amrgrid_class,    only: amrgrid
      use amrex_amr_module, only: amrex_bc_foextrap
      use amrdata_class,    only: interp_face_lin
      use messager,         only: die
      implicit none
      class(amrmpcomp), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), intent(in), optional :: name

      ! Materials must have been provided by the user before this call
      if (.not.associated(this%liq).or..not.associated(this%gas)) call die('[amrmpcomp initialize] liq and gas must be assigned before initialize')
      ! Assign number of arrays in Q, number of ghost cells
      this%nQ=4; this%nover=max(this%nover,3); this%nover_vel=max(this%nover_vel,2); this%nover_old=max(this%nover_old,2)
      ! Compressible face velocities are not divergence-free
      this%interp_vel=interp_face_lin
      call this%amrmpflow%initialize(amr,name); call this%set_parent()

      ! Initialize staggered momentum components
      call this%rhoU%initialize(amr,name='rhoU',ncomp=1,ng=0,nodal=[.true. ,.false.,.false.]); this%rhoU%parent=>this
      call this%rhoV%initialize(amr,name='rhoV',ncomp=1,ng=0,nodal=[.false.,.true. ,.false.]); this%rhoV%parent=>this
      call this%rhoW%initialize(amr,name='rhoW',ncomp=1,ng=0,nodal=[.false.,.false.,.true. ]); this%rhoW%parent=>this
      call this%rhoUold%initialize(amr,name='rhoUold',ncomp=1,ng=0,nodal=[.true. ,.false.,.false.]); this%rhoUold%parent=>this
      call this%rhoVold%initialize(amr,name='rhoVold',ncomp=1,ng=0,nodal=[.false.,.true. ,.false.]); this%rhoVold%parent=>this
      call this%rhoWold%initialize(amr,name='rhoWold',ncomp=1,ng=0,nodal=[.false.,.false.,.true. ]); this%rhoWold%parent=>this
      ! Initialize staggered midtime velocity components
      call this%Umid%initialize(amr,name='Umid',ncomp=1,ng=1,nodal=[.true. ,.false.,.false.]); this%Umid%parent=>this
      call this%Vmid%initialize(amr,name='Vmid',ncomp=1,ng=1,nodal=[.false.,.true. ,.false.]); this%Vmid%parent=>this
      call this%Wmid%initialize(amr,name='Wmid',ncomp=1,ng=1,nodal=[.false.,.false.,.true. ]); this%Wmid%parent=>this

      ! Initialize mixture variables
      call this%UVW%initialize(amr,name='UVW',ncomp=3,ng=0);                this%UVW%parent=>this
      call this%P%initialize   (amr,name='P'   ,ncomp=1,ng=this%nover_old); this%P%parent  =>this
      call this%C%initialize  (amr,name='C'  ,ncomp=1,ng=0);                this%C%parent  =>this

      ! Initialize phasic primitive variables
      call this%RHOL%initialize (amr,name='RHOL' ,ncomp=1,ng=this%nover); this%RHOL%parent =>this
      call this%RHOG%initialize (amr,name='RHOG' ,ncomp=1,ng=this%nover); this%RHOG%parent =>this
      call this%RHOLold%initialize(amr,name='RHOLold',ncomp=1,ng=this%nover); this%RHOLold%parent=>this
      call this%RHOGold%initialize(amr,name='RHOGold',ncomp=1,ng=this%nover); this%RHOGold%parent=>this
      call this%IL%initialize   (amr,name='IL'   ,ncomp=1,ng=0         ); this%IL%parent   =>this
      call this%IG%initialize   (amr,name='IG'   ,ncomp=1,ng=0         ); this%IG%parent   =>this
      call this%PL%initialize   (amr,name='PL'   ,ncomp=1,ng=this%nover_old); this%PL%parent=>this
      call this%PG%initialize   (amr,name='PG'   ,ncomp=1,ng=this%nover_old); this%PG%parent=>this
      call this%PLold%initialize  (amr,name='PLold'  ,ncomp=1,ng=this%nover_old); this%PLold%parent  =>this
      call this%PGold%initialize  (amr,name='PGold'  ,ncomp=1,ng=this%nover_old); this%PGold%parent  =>this
      call this%TL%initialize   (amr,name='TL'   ,ncomp=1,ng=1         ); this%TL%parent   =>this
      call this%TG%initialize   (amr,name='TG'   ,ncomp=1,ng=1         ); this%TG%parent   =>this

      ! Initialize physical properties (Neumann BCs on those)
      call this%visc%initialize (amr,name='visc' ,ncomp=1,ng=1); this%visc%parent =>this
      call this%beta%initialize (amr,name='beta' ,ncomp=1,ng=1); this%beta%parent =>this
      call this%diffL%initialize(amr,name='diffL',ncomp=1,ng=1); this%diffL%parent=>this
      call this%diffG%initialize(amr,name='diffG',ncomp=1,ng=1); this%diffG%parent=>this
      if (.not.amr%xper) then
         this%visc%lo_bc(1,1) =amrex_bc_foextrap; this%visc%hi_bc(1,1) =amrex_bc_foextrap
         this%beta%lo_bc(1,1) =amrex_bc_foextrap; this%beta%hi_bc(1,1) =amrex_bc_foextrap
         this%diffL%lo_bc(1,1)=amrex_bc_foextrap; this%diffL%hi_bc(1,1)=amrex_bc_foextrap
         this%diffG%lo_bc(1,1)=amrex_bc_foextrap; this%diffG%hi_bc(1,1)=amrex_bc_foextrap
      end if
      if (.not.amr%yper) then
         this%visc%lo_bc(2,1) =amrex_bc_foextrap; this%visc%hi_bc(2,1) =amrex_bc_foextrap
         this%beta%lo_bc(2,1) =amrex_bc_foextrap; this%beta%hi_bc(2,1) =amrex_bc_foextrap
         this%diffL%lo_bc(2,1)=amrex_bc_foextrap; this%diffL%hi_bc(2,1)=amrex_bc_foextrap
         this%diffG%lo_bc(2,1)=amrex_bc_foextrap; this%diffG%hi_bc(2,1)=amrex_bc_foextrap
      end if
      if (.not.amr%zper) then
         this%visc%lo_bc(3,1) =amrex_bc_foextrap; this%visc%hi_bc(3,1) =amrex_bc_foextrap
         this%beta%lo_bc(3,1) =amrex_bc_foextrap; this%beta%hi_bc(3,1) =amrex_bc_foextrap
         this%diffL%lo_bc(3,1)=amrex_bc_foextrap; this%diffL%hi_bc(3,1)=amrex_bc_foextrap
         this%diffG%lo_bc(3,1)=amrex_bc_foextrap; this%diffG%hi_bc(3,1)=amrex_bc_foextrap
      end if

      ! Allocate per-level flux workspaces (built in reset_fields when levels are created or remade)
      allocate(this%Fx (0:amr%maxlvl),this%Fy (0:amr%maxlvl),this%Fz (0:amr%maxlvl))
      allocate(this%FUx(0:amr%maxlvl),this%FVy(0:amr%maxlvl),this%FWz(0:amr%maxlvl))
      allocate(this%FUy(0:amr%maxlvl),this%FVx(0:amr%maxlvl))
      allocate(this%FUz(0:amr%maxlvl),this%FWx(0:amr%maxlvl))
      allocate(this%FVz(0:amr%maxlvl),this%FWy(0:amr%maxlvl))
      allocate(this%Ex (0:amr%maxlvl),this%Ey (0:amr%maxlvl),this%Ez (0:amr%maxlvl))

      ! Register callbacks with amrgrid using concrete dispatchers
      select type (this)
         type is (amrmpcomp)
            call this%amr%add_on_init   (amrmpcomp_on_init,   c_loc(this))
            call this%amr%add_on_coarse (amrmpcomp_on_coarse, c_loc(this))
            call this%amr%add_on_remake (amrmpcomp_on_remake, c_loc(this))
            call this%amr%add_on_clear  (amrmpcomp_on_clear,  c_loc(this))
            call this%amr%add_tagging   (amrmpcomp_tagging,   c_loc(this))
            call this%amr%add_postregrid(amrmpcomp_postregrid,c_loc(this))
            call this%amr%set_get_cost  (amrmpcomp_get_cost,  c_loc(this))
      end select

      ! Print solver info
      call this%print()

   end subroutine initialize

   !> Finalize the staggered compressible multiphase solver
   subroutine finalize(this)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      call this%rhoU%finalize(); call this%rhoUold%finalize()
      call this%rhoV%finalize(); call this%rhoVold%finalize()
      call this%rhoW%finalize(); call this%rhoWold%finalize()
      call this%Umid%finalize(); call this%Vmid%finalize(); call this%Wmid%finalize()
      call this%UVW%finalize()
      call this%P%finalize()
      call this%C%finalize()
      call this%RHOL%finalize(); call this%RHOG%finalize()
      call this%IL%finalize(); call this%IG%finalize()
      call this%PL%finalize(); call this%PG%finalize()
      call this%TL%finalize(); call this%TG%finalize()
      call this%RHOLold%finalize(); call this%RHOGold%finalize()
      call this%PLold%finalize(); call this%PGold%finalize()
      call this%visc%finalize(); call this%beta%finalize()
      call this%diffL%finalize(); call this%diffG%finalize()
      call destroy_fluxes(this)
      if (allocated(this%Fx)) deallocate(this%Fx,this%Fy,this%Fz,this%FUx,this%FVy,this%FWz,this%FUy,this%FVx,this%FUz,this%FWx,this%FVz,this%FWy,this%Ex,this%Ey,this%Ez)
      nullify(this%user_init); nullify(this%user_tagging); nullify(this%user_bc); nullify(this%user_vofbc)
      nullify(this%liq); nullify(this%gas); nullify(this%relax)
      call this%amrmpflow%finalize()
   end subroutine finalize

   ! ============================================================================
   ! INTERNAL CALLBACK OVERRIDES
   ! ============================================================================

   !> Override on_init: reset levels and set to zero
   subroutine on_init(this,lvl,time,ba,dm)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Parent handles face velocities, divergence, Q, and VOF data
      call this%amrmpflow%on_init(lvl,time,ba,dm)
      ! Reset level layouts
      call reset_fields(this,lvl,ba,dm)
      ! Set to zero
      call this%rhoU%setval(val=0.0_WP,lvl=lvl); call this%rhoUold%setval(val=0.0_WP,lvl=lvl)
      call this%rhoV%setval(val=0.0_WP,lvl=lvl); call this%rhoVold%setval(val=0.0_WP,lvl=lvl)
      call this%rhoW%setval(val=0.0_WP,lvl=lvl); call this%rhoWold%setval(val=0.0_WP,lvl=lvl)
      call this%Umid%setval(val=0.0_WP,lvl=lvl); call this%Vmid%setval(val=0.0_WP,lvl=lvl); call this%Wmid%setval(val=0.0_WP,lvl=lvl)
      call this%UVW%setval(val=0.0_WP,lvl=lvl)
      call this%P%setval(val=0.0_WP,lvl=lvl)
      call this%C%setval(val=0.0_WP,lvl=lvl)
      call this%RHOL%setval(val=0.0_WP,lvl=lvl); call this%RHOG%setval(val=0.0_WP,lvl=lvl)
      call this%IL%setval(val=0.0_WP,lvl=lvl); call this%IG%setval(val=0.0_WP,lvl=lvl)
      call this%PL%setval(val=0.0_WP,lvl=lvl); call this%PG%setval(val=0.0_WP,lvl=lvl)
      call this%TL%setval(val=0.0_WP,lvl=lvl); call this%TG%setval(val=0.0_WP,lvl=lvl)
      call this%RHOLold%setval(val=0.0_WP,lvl=lvl); call this%RHOGold%setval(val=0.0_WP,lvl=lvl)
      call this%PLold%setval(val=0.0_WP,lvl=lvl); call this%PGold%setval(val=0.0_WP,lvl=lvl)
      call this%visc%setval(val=0.0_WP,lvl=lvl); call this%beta%setval(val=0.0_WP,lvl=lvl)
      call this%diffL%setval(val=0.0_WP,lvl=lvl); call this%diffG%setval(val=0.0_WP,lvl=lvl)
   end subroutine on_init

   !> Override on_coarse: create new fine level from coarse
   !> Momentum is rebuilt from the interpolated face velocity in post_regrid
   subroutine on_coarse(this,lvl,time,ba,dm)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Parent handles face velocities, divergence, Q, and VOF data
      call this%amrmpflow%on_coarse(lvl,time,ba,dm)
      ! Momentum and derived variables are reset
      call reset_fields(this,lvl,ba,dm)
      ! No heat conduction on the new level until the user sets the phasic conductivities
      call this%diffL%setval(val=0.0_WP,lvl=lvl); call this%diffG%setval(val=0.0_WP,lvl=lvl)
   end subroutine on_coarse

   !> Override on_remake: migrate data on regrid
   subroutine on_remake(this,lvl,time,ba,dm)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      ! Parent handles face velocities, divergence, Q, and VOF data
      call this%amrmpflow%on_remake(lvl,time,ba,dm)
      ! Momentum and derived variables are reset
      call reset_fields(this,lvl,ba,dm)
      ! No heat conduction on the remade level until the user sets the phasic conductivities
      call this%diffL%setval(val=0.0_WP,lvl=lvl); call this%diffG%setval(val=0.0_WP,lvl=lvl)
   end subroutine on_remake

   !> Override on_clear: delete level
   subroutine on_clear(this,lvl)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lvl
      ! Parent handles face velocities, divergence, Q, and VOF data
      call this%amrmpflow%on_clear(lvl)
      ! Clear momentum and derived variables
      call this%rhoU%clear_level(lvl); call this%rhoUold%clear_level(lvl)
      call this%rhoV%clear_level(lvl); call this%rhoVold%clear_level(lvl)
      call this%rhoW%clear_level(lvl); call this%rhoWold%clear_level(lvl)
      call this%Umid%clear_level(lvl); call this%Vmid%clear_level(lvl); call this%Wmid%clear_level(lvl)
      call this%UVW%clear_level(lvl)
      call this%P%clear_level(lvl)
      call this%C%clear_level(lvl)
      call this%RHOL%clear_level(lvl); call this%RHOG%clear_level(lvl)
      call this%IL%clear_level(lvl); call this%IG%clear_level(lvl)
      call this%PL%clear_level(lvl); call this%PG%clear_level(lvl)
      call this%TL%clear_level(lvl); call this%TG%clear_level(lvl)
      call this%RHOLold%clear_level(lvl); call this%RHOGold%clear_level(lvl)
      call this%PLold%clear_level(lvl); call this%PGold%clear_level(lvl)
      call this%visc%clear_level(lvl); call this%beta%clear_level(lvl)
      call this%diffL%clear_level(lvl); call this%diffG%clear_level(lvl)
      call destroy_fluxes(this,lvl)
   end subroutine on_clear

   !> Reset the layouts of all solver-owned fields at a level
   subroutine reset_fields(this,lvl,ba,dm)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      call this%rhoU%reset_level(lvl,ba,dm); call this%rhoUold%reset_level(lvl,ba,dm)
      call this%rhoV%reset_level(lvl,ba,dm); call this%rhoVold%reset_level(lvl,ba,dm)
      call this%rhoW%reset_level(lvl,ba,dm); call this%rhoWold%reset_level(lvl,ba,dm)
      call this%Umid%reset_level(lvl,ba,dm); call this%Vmid%reset_level(lvl,ba,dm); call this%Wmid%reset_level(lvl,ba,dm)
      call this%UVW%reset_level(lvl,ba,dm)
      call this%P%reset_level(lvl,ba,dm)
      call this%C%reset_level(lvl,ba,dm)
      call this%RHOL%reset_level(lvl,ba,dm); call this%RHOG%reset_level(lvl,ba,dm)
      call this%IL%reset_level(lvl,ba,dm); call this%IG%reset_level(lvl,ba,dm)
      call this%PL%reset_level(lvl,ba,dm); call this%PG%reset_level(lvl,ba,dm)
      call this%TL%reset_level(lvl,ba,dm); call this%TG%reset_level(lvl,ba,dm)
      call this%RHOLold%reset_level(lvl,ba,dm); call this%RHOGold%reset_level(lvl,ba,dm)
      call this%PLold%reset_level(lvl,ba,dm); call this%PGold%reset_level(lvl,ba,dm)
      call this%visc%reset_level(lvl,ba,dm); call this%beta%reset_level(lvl,ba,dm)
      call this%diffL%reset_level(lvl,ba,dm); call this%diffG%reset_level(lvl,ba,dm)
      call build_fluxes(this,lvl,ba,dm)
   end subroutine reset_fields
   !> Build flux mfabs for get_dQdt and get_energy
   subroutine build_fluxes(this,lvl,ba,dm)
      use amrgrid_class, only: mfab_rebuild
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      integer :: a
      ! Fx, Fy, Fz: cell-face
      call mfab_rebuild(this%Fx(lvl),ba,dm,nc=this%nQ,ng=1,nodal=[.true. ,.false.,.false.])
      call mfab_rebuild(this%Fy(lvl),ba,dm,nc=this%nQ,ng=1,nodal=[.false.,.true. ,.false.])
      call mfab_rebuild(this%Fz(lvl),ba,dm,nc=this%nQ,ng=1,nodal=[.false.,.false.,.true. ])
      ! FUx, FVy, FWz: cell-centered (diagonal fluxes)
      call mfab_rebuild(this%FUx(lvl),ba,dm,nc=1,ng=1)
      call mfab_rebuild(this%FVy(lvl),ba,dm,nc=1,ng=1)
      call mfab_rebuild(this%FWz(lvl),ba,dm,nc=1,ng=1)
      ! FUy, FVx: xy-edge (cross-fluxes)
      call mfab_rebuild(this%FUy(lvl),ba,dm,nc=1,ng=0,nodal=[.true.,.true.,.false.])
      call mfab_rebuild(this%FVx(lvl),ba,dm,nc=1,ng=0,nodal=[.true.,.true.,.false.])
      ! FUz, FWx: xz-edge (cross-fluxes)
      call mfab_rebuild(this%FUz(lvl),ba,dm,nc=1,ng=0,nodal=[.true.,.false.,.true.])
      call mfab_rebuild(this%FWx(lvl),ba,dm,nc=1,ng=0,nodal=[.true.,.false.,.true.])
      ! FVz, FWy: yz-edge (cross-fluxes)
      call mfab_rebuild(this%FVz(lvl),ba,dm,nc=1,ng=0,nodal=[.false.,.true.,.true.])
      call mfab_rebuild(this%FWy(lvl),ba,dm,nc=1,ng=0,nodal=[.false.,.true.,.true.])
      ! Ex, Ey, Ez: cell-face (mixture work fluxes)
      call mfab_rebuild(this%Ex(lvl),ba,dm,nc=1,ng=0,nodal=[.true. ,.false.,.false.])
      call mfab_rebuild(this%Ey(lvl),ba,dm,nc=1,ng=0,nodal=[.false.,.true. ,.false.])
      call mfab_rebuild(this%Ez(lvl),ba,dm,nc=1,ng=0,nodal=[.false.,.false.,.true. ])
      ! The band, its face workspaces, and the surface tension data only exist on the finest level
      if (lvl.eq.this%amr%maxlvl) then
         call mfab_rebuild(this%band,ba,dm,nc=1,ng=2)
         call mfab_rebuild(this%curvold,ba,dm,nc=1,ng=1)
         call mfab_rebuild(this%SDold  ,ba,dm,nc=1,ng=1)
         do a=1,3
            call mfab_rebuild(this%SLF(a),ba,dm,nc=nSLF,ng=1,nodal=eu(:,a).eq.1)
            call mfab_rebuild(this%SLV(a),ba,dm,nc=8   ,ng=0,nodal=eu(:,a).eq.1)
            call mfab_rebuild(this%Pj(a) ,ba,dm,nc=1   ,ng=0,nodal=eu(:,a).eq.1)
         end do
      end if
   end subroutine build_fluxes

   !> Destroy the get_dQdt flux workspaces of a level, or of all levels if lvl is absent
   subroutine destroy_fluxes(this,lvl)
      use amrex_amr_module, only: amrex_multifab_destroy
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in), optional :: lvl
      if (.not.allocated(this%Fx)) return
      ! Per-level workspaces
      if (present(lvl)) then
         call amrex_multifab_destroy(this%Fx (lvl)); call amrex_multifab_destroy(this%Fy (lvl)); call amrex_multifab_destroy(this%Fz (lvl))
         call amrex_multifab_destroy(this%FUx(lvl)); call amrex_multifab_destroy(this%FVy(lvl)); call amrex_multifab_destroy(this%FWz(lvl))
         call amrex_multifab_destroy(this%FUy(lvl)); call amrex_multifab_destroy(this%FVx(lvl))
         call amrex_multifab_destroy(this%FUz(lvl)); call amrex_multifab_destroy(this%FWx(lvl))
         call amrex_multifab_destroy(this%FVz(lvl)); call amrex_multifab_destroy(this%FWy(lvl))
         call amrex_multifab_destroy(this%Ex (lvl)); call amrex_multifab_destroy(this%Ey (lvl)); call amrex_multifab_destroy(this%Ez (lvl))
         if (lvl.lt.this%amr%maxlvl) return
      else
         call amrex_multifab_destroy(this%Fx ); call amrex_multifab_destroy(this%Fy ); call amrex_multifab_destroy(this%Fz )
         call amrex_multifab_destroy(this%FUx); call amrex_multifab_destroy(this%FVy); call amrex_multifab_destroy(this%FWz)
         call amrex_multifab_destroy(this%FUy); call amrex_multifab_destroy(this%FVx)
         call amrex_multifab_destroy(this%FUz); call amrex_multifab_destroy(this%FWx)
         call amrex_multifab_destroy(this%FVz); call amrex_multifab_destroy(this%FWy)
         call amrex_multifab_destroy(this%Ex ); call amrex_multifab_destroy(this%Ey ); call amrex_multifab_destroy(this%Ez )
      end if
      ! Finest-level workspaces
      call amrex_multifab_destroy(this%band)
      call amrex_multifab_destroy(this%SLF)
      call amrex_multifab_destroy(this%SLV)
      call amrex_multifab_destroy(this%curvold); call amrex_multifab_destroy(this%SDold)
      call amrex_multifab_destroy(this%Pj)
   end subroutine destroy_fluxes

   !> Override post_regrid: average down and fill, rebuild momentum from the regridded face velocity, then rebuild primitive variables
   subroutine post_regrid(this,lbase,time)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      ! Parent handles face velocities, Q, and VOF data
      call this%amrmpflow%post_regrid(lbase,time)
      ! Momentum from face velocity and face density
      call this%get_momentum()
      ! Rebuild primitive variables
      call this%get_primitive(this%Q,this%rhoU,this%rhoV,this%rhoW,time)
   end subroutine post_regrid

   ! ============================================================================
   ! BOUNDARY CONDITIONS
   ! ============================================================================

   !> Velocity BC override: forward to user_bc with U/V/W component name
   subroutine stag_apply_velbc(this,lvl,time,face,bx,comp,p)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lvl,face
      real(WP), intent(in) :: time
      type(amrex_box), intent(in) :: bx
      character(len=1), intent(in) :: comp
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      if (associated(this%user_bc)) call this%user_bc(lvl=lvl,time=time,face=face,bx=bx,comp=comp,p=p)
   end subroutine stag_apply_velbc

   !> Q BC override: forward to user_bc with comp='Q'
   subroutine stag_apply_Qbc(this,lvl,time,face,bx,p)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lvl,face
      real(WP), intent(in) :: time
      type(amrex_box), intent(in) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      if (associated(this%user_bc)) call this%user_bc(lvl=lvl,time=time,face=face,bx=bx,comp='Q',p=p)
   end subroutine stag_apply_Qbc

   !> VOF BC override: forward to user_vofbc
   subroutine stag_apply_vofbc(this,lvl,time,face,bx,pVF,pCL,pCG,pPLIC)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer, intent(in) :: lvl,face
      real(WP), intent(in) :: time
      type(amrex_box), intent(in) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pCL,pCG,pPLIC
      if (associated(this%user_vofbc)) call this%user_vofbc(lvl=lvl,time=time,face=face,bx=bx,pVF=pVF,pCL=pCL,pCG=pCG,pPLIC=pPLIC)
   end subroutine stag_apply_vofbc

   ! ============================================================================
   ! UTILITIES
   ! ============================================================================

   !> Store the start-of-step state: Q, face velocity and VOF data (parent), face momenta,
   !> phasic primitives used by the semi-Lagrangian donors and the midpoint state, and the curvature and surface density
   subroutine store_old(this)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      call this%amrmpflow%store_old()
      call this%rhoUold%copy(src=this%rhoU)
      call this%rhoVold%copy(src=this%rhoV)
      call this%rhoWold%copy(src=this%rhoW)
      call this%RHOLold%copy(src=this%RHOL)
      call this%RHOGold%copy(src=this%RHOG)
      call this%PLold%copy(src=this%PL)
      call this%PGold%copy(src=this%PG)
      if (this%amr%clvl().lt.this%amr%maxlvl) return
      call this%curvold%copy(srcmf=this%curv,srccomp=1,dstcomp=1,nc=1,ng=1)
      call this%SDold%copy  (srcmf=this%SD  ,srccomp=1,dstcomp=1,nc=1,ng=1)
   end subroutine store_old

   !> Get midpoint for face velocities, phasic RHO, P, and T, and mixture pressure
   subroutine prepare_midpoint(this)
      use mpi_f08, only: MPI_Wtime
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP) :: tw0
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pRHOL,pPL,pTL,pRHOLold,pPLold
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pRHOG,pPG,pTG,pRHOGold,pPGold
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pP,pVFold
      real(WP) :: Pn,Iold,Inew
      logical :: has_old,has_new
      tw0=MPI_Wtime()
      ! Face velocities
      call this%U%lincomb(a=0.5_WP,src1=this%Uold,b=0.5_WP,src2=this%U)
      call this%V%lincomb(a=0.5_WP,src1=this%Vold,b=0.5_WP,src2=this%V)
      call this%W%lincomb(a=0.5_WP,src1=this%Wold,b=0.5_WP,src2=this%W)
      ! Phasic primitive variables
      do lvl=0,this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pRHOL=>this%RHOL%mf(lvl)%dataptr(mfi); pRHOLold=>this%RHOLold%mf(lvl)%dataptr(mfi)
            pPL  =>this%PL%mf(lvl)%dataptr(mfi);   pPLold  =>this%PLold%mf(lvl)%dataptr(mfi)
            pRHOG=>this%RHOG%mf(lvl)%dataptr(mfi); pRHOGold=>this%RHOGold%mf(lvl)%dataptr(mfi)
            pPG  =>this%PG%mf(lvl)%dataptr(mfi);   pPGold  =>this%PGold%mf(lvl)%dataptr(mfi)
            pTL  =>this%TL%mf(lvl)%dataptr(mfi)
            pTG  =>this%TG%mf(lvl)%dataptr(mfi)
            pP   =>this%P%mf(lvl)%dataptr(mfi);    pVFold  =>this%VFold%mf(lvl)%dataptr(mfi)
            ! Phasic temperature
            bx=mfi%growntilebox(min(this%TL%ng,this%TG%ng))
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Liquid
               has_old=pRHOLold(i,j,k,1).gt.0.0_WP; has_new=pRHOL(i,j,k,1).gt.0.0_WP
               Iold=0.0_WP; if (has_old) Iold=this%liq%get_e_from_p_rho(p=pPLold(i,j,k,1),rho=pRHOLold(i,j,k,1),y=[1.0_WP])
               Inew=0.0_WP; if (has_new) Inew=this%liq%get_e_from_p_rho(p=pPL   (i,j,k,1),rho=pRHOL   (i,j,k,1),y=[1.0_WP])
               pTL(i,j,k,1)=0.0_WP
               if (has_old.or.has_new) pTL(i,j,k,1)=this%liq%get_T_from_rho_e(rho=midval(pRHOLold(i,j,k,1),pRHOL(i,j,k,1),has_old,has_new),e=midval(Iold,Inew,has_old,has_new),y=[1.0_WP])
               ! Gas
               has_old=pRHOGold(i,j,k,1).gt.0.0_WP; has_new=pRHOG(i,j,k,1).gt.0.0_WP
               Iold=0.0_WP; if (has_old) Iold=this%gas%get_e_from_p_rho(p=pPGold(i,j,k,1),rho=pRHOGold(i,j,k,1),y=[1.0_WP])
               Inew=0.0_WP; if (has_new) Inew=this%gas%get_e_from_p_rho(p=pPG   (i,j,k,1),rho=pRHOG   (i,j,k,1),y=[1.0_WP])
               pTG(i,j,k,1)=0.0_WP
               if (has_old.or.has_new) pTG(i,j,k,1)=this%gas%get_T_from_rho_e(rho=midval(pRHOGold(i,j,k,1),pRHOG(i,j,k,1),has_old,has_new),e=midval(Iold,Inew,has_old,has_new),y=[1.0_WP])
            end do; end do; end do
            ! Phasic mixture pressure
            bx=mfi%growntilebox(min(this%PL%ng,this%PG%ng,this%P%ng,this%VFold%ng))
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               Pn=pVFold(i,j,k,1)*pPLold(i,j,k,1)+(1.0_WP-pVFold(i,j,k,1))*pPGold(i,j,k,1)
               ! Liquid
               has_old=pRHOLold(i,j,k,1).gt.0.0_WP; has_new=pRHOL(i,j,k,1).gt.0.0_WP
               pPL(i,j,k,1)=0.5_WP*(merge(pPLold(i,j,k,1),Pn,has_old)+merge(pPL(i,j,k,1),pP(i,j,k,1),has_new))
               ! Gas
               has_old=pRHOGold(i,j,k,1).gt.0.0_WP; has_new=pRHOG(i,j,k,1).gt.0.0_WP
               pPG(i,j,k,1)=0.5_WP*(merge(pPGold(i,j,k,1),Pn,has_old)+merge(pPG(i,j,k,1),pP(i,j,k,1),has_new))
               ! Mixture
               pP(i,j,k,1)=0.5_WP*(Pn+pP(i,j,k,1))
            end do; end do; end do
            ! Phasic densities
            bx=mfi%growntilebox(this%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pRHOL (i,j,k,1)=midval(pRHOLold (i,j,k,1),pRHOL (i,j,k,1),pRHOLold(i,j,k,1).gt.0.0_WP,pRHOL(i,j,k,1).gt.0.0_WP)
               pRHOG (i,j,k,1)=midval(pRHOGold (i,j,k,1),pRHOG (i,j,k,1),pRHOGold(i,j,k,1).gt.0.0_WP,pRHOG(i,j,k,1).gt.0.0_WP)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      this%wt_mid=this%wt_mid+(MPI_Wtime()-tw0)
   contains
      !> Average of the old and new values where available
      pure real(WP) function midval(vold,vnew,has_old,has_new)
         implicit none
         real(WP), intent(in) :: vold,vnew
         logical, intent(in) :: has_old,has_new
         if (has_old.and.has_new) then
            midval=0.5_WP*(vold+vnew)
         else if (has_old) then
            midval=vold
         else
            midval=vnew
         end if
      end function midval
   end subroutine prepare_midpoint

   !> Surface tension pressure jump on finest-level faces: mean of the start-of-step and current-iterate jumps
   subroutine add_surface_tension_jump(this)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer :: lvl,i,j,k,a
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pVFold,pCurv,pCurvold,pSD,pSDold,pPj
      real(WP), dimension(3) :: hinv
      if (this%amr%clvl().lt.this%amr%maxlvl) return
      lvl=this%amr%maxlvl
      hinv=1.0_WP/[this%amr%dx(lvl),this%amr%dy(lvl),this%amr%dz(lvl)]
      call this%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         pVF  =>this%VF%mf(lvl)%dataptr(mfi); pVFold  =>this%VFold%mf(lvl)%dataptr(mfi)
         pCurv=>this%curv%dataptr(mfi);       pCurvold=>this%curvold%dataptr(mfi)
         pSD  =>this%SD%dataptr(mfi);         pSDold  =>this%SDold%dataptr(mfi)
         do a=1,3
            pPj=>this%Pj(a)%dataptr(mfi)
            bx=mfi%nodaltilebox(a)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pPj(i,j,k,1)=0.5_WP*this%sigma*hinv(a)*(face_jump(a,i,j,k,pCurvold,pSDold,pVFold)+face_jump(a,i,j,k,pCurv,pSD,pVF))
            end do; end do; end do
         end do
      end do
      call this%amr%mfiter_destroy(mfi)
   contains
      !> Surface-density-weighted curvature times the VF jump across a-face (i,j,k), for one time level
      real(WP) function face_jump(a,i,j,k,kap,sd,vf)
         implicit none
         integer, intent(in) :: a,i,j,k
         real(WP), dimension(:,:,:,:), contiguous, pointer, intent(in) :: kap,sd,vf
         integer :: l,m,n
         l=i-eu(1,a); m=j-eu(2,a); n=k-eu(3,a)
         face_jump=0.0_WP
         if (sd(l,m,n,1)+sd(i,j,k,1).gt.0.0_WP) face_jump=(sd(l,m,n,1)*kap(l,m,n,1)+sd(i,j,k,1)*kap(i,j,k,1))/(sd(l,m,n,1)+sd(i,j,k,1))*(vf(i,j,k,1)-vf(l,m,n,1))
      end function face_jump
   end subroutine add_surface_tension_jump

   !> Synchronize/restrict momentum, recover valid face velocities, restrict them, and fill velocity BCs.
   subroutine get_face_velocity(this,Q,rhoU,rhoV,rhoW,time)
      use mpi_f08, only: MPI_Wtime
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP) :: tw0
      type(amrdata), intent(in) :: Q
      type(amrdata), intent(inout) :: rhoU,rhoV,rhoW
      real(WP), intent(in) :: time
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: fbx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pRHOU,pRHOV,pRHOW,pU,pV,pW
      tw0=MPI_Wtime()
      !  Sync momentum
      call this%sync_momentum(rhoU,rhoV,rhoW)
      ! Compute face velocities
      do lvl=0,this%amr%clvl()
         ! Loop over tiles
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pQ=>Q%mf(lvl)%dataptr(mfi)
            pRHOU=>rhoU%mf(lvl)%dataptr(mfi)
            pRHOV=>rhoV%mf(lvl)%dataptr(mfi)
            pRHOW=>rhoW%mf(lvl)%dataptr(mfi)
            pU=>this%U%mf(lvl)%dataptr(mfi)
            pV=>this%V%mf(lvl)%dataptr(mfi)
            pW=>this%W%mf(lvl)%dataptr(mfi)
            ! Get X-face velocity
            fbx=mfi%nodaltilebox(1)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               pU(i,j,k,1)=pRHOU(i,j,k,1)/max(0.5_WP*sum(pQ(i-1:i,j,k,1:2)),this%rho_floor)
            end do; end do; end do
            ! Get Y-face velocity
            fbx=mfi%nodaltilebox(2)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               pV(i,j,k,1)=pRHOV(i,j,k,1)/max(0.5_WP*sum(pQ(i,j-1:j,k,1:2)),this%rho_floor)
            end do; end do; end do
            ! Get Z-face velocity
            fbx=mfi%nodaltilebox(3)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               pW(i,j,k,1)=pRHOW(i,j,k,1)/max(0.5_WP*sum(pQ(i,j,k-1:k,1:2)),this%rho_floor)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      ! Replace covered coarse face velocities by the average of the fine velocities
      call this%average_down_velocity()
      ! Fill derived velocity ghosts
      call this%fill_velocity(time=time)
      call this%apply_momentum_bc(Q,rhoU,rhoV,rhoW)
      this%wt_fvel=this%wt_fvel+(MPI_Wtime()-tw0)
   end subroutine get_face_velocity

   !> Interpolate face velocities to cell centers in UVW, including at ghost cells
   subroutine get_cell_velocity(this)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW,pUVW
      do lvl=0,this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pU=>this%U%mf(lvl)%dataptr(mfi)
            pV=>this%V%mf(lvl)%dataptr(mfi)
            pW=>this%W%mf(lvl)%dataptr(mfi)
            pUVW=>this%UVW%mf(lvl)%dataptr(mfi)
            ! Loop over tiles grown by the UVW ghost cells
            bx=mfi%growntilebox(this%UVW%ng)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pUVW(i,j,k,1)=0.5_WP*sum(pU(i:i+1,j,k,1))
               pUVW(i,j,k,2)=0.5_WP*sum(pV(i,j:j+1,k,1))
               pUVW(i,j,k,3)=0.5_WP*sum(pW(i,j,k:k+1,1))
            end do; end do; end do
         end do
      call this%amr%mfiter_destroy(mfi)
      end do
end subroutine get_cell_velocity

   !> Reconcile duplicate valid faces, then restrict momentum from fine to coarse.
   subroutine sync_momentum(this,rhoU,rhoV,rhoW)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      type(amrdata), intent(inout) :: rhoU,rhoV,rhoW
      integer :: lvl
      do lvl=this%amr%clvl(),0,-1
         call rhoU%mf(lvl)%override_sync(this%amr%geom(lvl))
         call rhoV%mf(lvl)%override_sync(this%amr%geom(lvl))
         call rhoW%mf(lvl)%override_sync(this%amr%geom(lvl))
         if (lvl.eq.0) cycle
         call rhoU%average_downto(lvl-1)
         call rhoV%average_downto(lvl-1)
         call rhoW%average_downto(lvl-1)
      end do
   end subroutine sync_momentum

   !> Enforce physical momentum BCs after Q and face velocities have been filled.
   subroutine apply_momentum_bc(this,Q,rhoU,rhoV,rhoW)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      type(amrdata), intent(in) :: Q
      type(amrdata), intent(inout) :: rhoU,rhoV,rhoW
      integer :: lvl
      type(amrex_mfiter) :: mfi
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pRHOU,pRHOV,pRHOW,pU,pV,pW
      do lvl=0,this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi,tiling=.false.)
         do while (mfi%next())
            ! Get pointers to data
            pQ=>Q%mf(lvl)%dataptr(mfi)
            pRHOU=>rhoU%mf(lvl)%dataptr(mfi); pU=>this%U%mf(lvl)%dataptr(mfi)
            pRHOV=>rhoV%mf(lvl)%dataptr(mfi); pV=>this%V%mf(lvl)%dataptr(mfi)
            pRHOW=>rhoW%mf(lvl)%dataptr(mfi); pW=>this%W%mf(lvl)%dataptr(mfi)
            call set_momentum(pRHOU,pU,dir=1,fbx=mfi%nodaltilebox(1),lo_bc=this%U%lo_bc(1,1),hi_bc=this%U%hi_bc(1,1),per=this%amr%xper)
            call set_momentum(pRHOV,pV,dir=2,fbx=mfi%nodaltilebox(2),lo_bc=this%V%lo_bc(2,1),hi_bc=this%V%hi_bc(2,1),per=this%amr%yper)
            call set_momentum(pRHOW,pW,dir=3,fbx=mfi%nodaltilebox(3),lo_bc=this%W%lo_bc(3,1),hi_bc=this%W%hi_bc(3,1),per=this%amr%zper)
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      
   contains

      subroutine set_momentum(pM,pVel,dir,fbx,lo_bc,hi_bc,per)
         use amrex_amr_module, only: amrex_bc_ext_dir,amrex_bc_reflect_odd
         implicit none
         real(WP), dimension(:,:,:,:), contiguous, pointer, intent(in) :: pM,pVel
         integer, intent(in) :: dir,lo_bc,hi_bc
         type(amrex_box), intent(in) :: fbx
         logical, intent(in) :: per
         integer :: i,j,k
         integer, dimension(3) :: n,m,dlo,dhi,slo,shi
         logical :: lo_set,hi_set,physical_ghost
         logical, dimension(3) :: periodic
         dlo=this%amr%geom(lvl)%domain%lo
         dhi=this%amr%geom(lvl)%domain%hi; dhi(dir)=dhi(dir)+1
         periodic=[this%amr%xper,this%amr%yper,this%amr%zper]
         lo_set=.not.per.and.(lo_bc.eq.amrex_bc_ext_dir.or.lo_bc.eq.amrex_bc_reflect_odd)
         hi_set=.not.per.and.(hi_bc.eq.amrex_bc_ext_dir.or.hi_bc.eq.amrex_bc_reflect_odd)
         ! Intersect the allocated face boxes with faces having a complete Q stencil.
         slo=[lbound(pQ,1),lbound(pQ,2),lbound(pQ,3)]
         shi=[ubound(pQ,1),ubound(pQ,2),ubound(pQ,3)]
         slo(dir)=slo(dir)+1
         slo=max(slo,[lbound(pM,1),lbound(pM,2),lbound(pM,3)],[lbound(pVel,1),lbound(pVel,2),lbound(pVel,3)])
         shi=min(shi,[ubound(pM,1),ubound(pM,2),ubound(pM,3)],[ubound(pVel,1),ubound(pVel,2),ubound(pVel,3)])
         do k=slo(3),shi(3); do j=slo(2),shi(2); do i=slo(1),shi(1)
                  n=[i,j,k]
                  physical_ghost=any(((n.lt.dlo).or.(n.gt.dhi)).and..not.periodic)
                  if (.not.physical_ghost) then
                     ! Interior/shared/periodic faces are not physical BCs.
                     if (.not.(all(n.ge.fbx%lo).and.all(n.le.fbx%hi))) cycle
                     if (.not.((lo_set.and.n(dir).eq.dlo(dir)).or.(hi_set.and.n(dir).eq.dhi(dir)))) cycle
                  end if
                  m=n; m(dir)=n(dir)-1
                  pM(i,j,k,1)=max(0.5_WP*(sum(pQ(m(1),m(2),m(3),1:2))+sum(pQ(i,j,k,1:2))),this%rho_floor)*pVel(i,j,k,1)
               end do; end do; end do
      end subroutine set_momentum

   end subroutine apply_momentum_bc

   !> Set valid face momenta from face velocities and mixture face density for initialization and regrid only
   subroutine get_momentum(this)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: fbx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pRHOU,pRHOV,pRHOW,pU,pV,pW
      do lvl=0,this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pQ=>this%Q%mf(lvl)%dataptr(mfi)
            pRHOU=>this%rhoU%mf(lvl)%dataptr(mfi); pU=>this%U%mf(lvl)%dataptr(mfi)
            pRHOV=>this%rhoV%mf(lvl)%dataptr(mfi); pV=>this%V%mf(lvl)%dataptr(mfi)
            pRHOW=>this%rhoW%mf(lvl)%dataptr(mfi); pW=>this%W%mf(lvl)%dataptr(mfi)
            ! X-face momentum
            fbx=mfi%nodaltilebox(1)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               pRHOU(i,j,k,1)=max(0.5_WP*sum(pQ(i-1:i,j,k,1:2)),this%rho_floor)*pU(i,j,k,1)
            end do; end do; end do
            ! Y-face momentum
            fbx=mfi%nodaltilebox(2)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               pRHOV(i,j,k,1)=max(0.5_WP*sum(pQ(i,j-1:j,k,1:2)),this%rho_floor)*pV(i,j,k,1)
            end do; end do; end do
            ! Z-face momentum
            fbx=mfi%nodaltilebox(3)
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               pRHOW(i,j,k,1)=max(0.5_WP*sum(pQ(i,j,k-1:k,1:2)),this%rho_floor)*pW(i,j,k,1)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      ! Fill same-level/periodic ghosts
      call this%rhoU%sync(); call this%rhoV%sync(); call this%rhoW%sync()
   end subroutine get_momentum

   !> Set phasic total energies from internal+kinetic energy
   subroutine get_total_energy(this,time)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP), intent(in) :: time
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pU,pV,pW
      real(WP), dimension(2) :: mk
      call this%get_face_velocity(this%Q,this%rhoU,this%rhoV,this%rhoW,time)
      do lvl=0,this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pQ=>this%Q%mf(lvl)%dataptr(mfi)
            pU=>this%U%mf(lvl)%dataptr(mfi); pV=>this%V%mf(lvl)%dataptr(mfi); pW=>this%W%mf(lvl)%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               mk=max(pQ(i,j,k,1:2),0.0_WP)
               pQ(i,j,k,3:4)=pQ(i,j,k,3:4)+mk/max(sum(mk),this%rho_floor)*cell_ke(pQ,pU,pV,pW,i,j,k,this%rho_floor)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      call this%Q%average_down(); call this%Q%fill(time=time)
      call this%get_primitive(this%Q,this%rhoU,this%rhoV,this%rhoW,time)
   end subroutine get_total_energy

   !> Mixture kinetic energy of cell (i,j,k)
   pure real(WP) function cell_ke(pQ,pU,pV,pW,i,j,k,rho_floor) result(ke)
      implicit none
      real(WP), dimension(:,:,:,:), contiguous, pointer, intent(in) :: pQ,pU,pV,pW
      integer, intent(in) :: i,j,k
      real(WP), intent(in) :: rho_floor
      ke=0.25_WP*(max(0.5_WP*sum(pQ(i-1:i,j,k,1:2)),rho_floor)*pU(i,j,k,1)**2+max(0.5_WP*sum(pQ(i:i+1,j,k,1:2)),rho_floor)*pU(i+1,j,k,1)**2 &
      &          +max(0.5_WP*sum(pQ(i,j-1:j,k,1:2)),rho_floor)*pV(i,j,k,1)**2+max(0.5_WP*sum(pQ(i,j:j+1,k,1:2)),rho_floor)*pV(i,j+1,k,1)**2 &
      &          +max(0.5_WP*sum(pQ(i,j,k-1:k,1:2)),rho_floor)*pW(i,j,k,1)**2+max(0.5_WP*sum(pQ(i,j,k:k+1,1:2)),rho_floor)*pW(i,j,k+1,1)**2)
   end function cell_ke

   ! ============================================================================
   ! PHYSICS METHODS
   ! ============================================================================

   !> Calculate phasic and mixture primitive variables from conserved variables
   subroutine get_primitive(this,Q,rhoU,rhoV,rhoW,time)
      use messager, only: die
      use mpi_f08,  only: MPI_Wtime
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP) :: tw0
      type(amrdata), intent(in) :: Q
      type(amrdata), intent(inout) :: rhoU,rhoV,rhoW
      real(WP), intent(in) :: time
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF,pRHOL,pRHOG,pIL,pIG,pPL,pPG,pTL,pTG,pP,pC
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW
      real(WP), dimension(2) :: Y,mk
      real(WP) :: CL,CG
      ! Check passed Q is as expected
      if (Q%ncomp.ne.this%nQ) call die('[amrmpcomp get_primitive] Q has wrong number of components')
      if (Q%ng.lt.this%nover) call die('[amrmpcomp get_primitive] Q must have at least nover ghost cells')
      ! Face velocities from momentum
      call this%get_face_velocity(Q,rhoU,rhoV,rhoW,time)
      tw0=MPI_Wtime()
      ! Loop over levels
      do lvl=0,this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pQ  =>Q%mf(lvl)%dataptr(mfi)
            pU  =>this%U%mf(lvl)%dataptr(mfi); pV=>this%V%mf(lvl)%dataptr(mfi); pW=>this%W%mf(lvl)%dataptr(mfi)
            pVF =>this%VF%mf(lvl)%dataptr(mfi)
            pRHOL=>this%RHOL%mf(lvl)%dataptr(mfi); pRHOG=>this%RHOG%mf(lvl)%dataptr(mfi)
            pIL =>this%IL%mf(lvl)%dataptr(mfi);   pIG =>this%IG%mf(lvl)%dataptr(mfi)
            pPL =>this%PL%mf(lvl)%dataptr(mfi);   pPG =>this%PG%mf(lvl)%dataptr(mfi)
            pTL =>this%TL%mf(lvl)%dataptr(mfi);   pTG =>this%TG%mf(lvl)%dataptr(mfi)
            pP  =>this%P%mf(lvl)%dataptr(mfi)
            pC  =>this%C%mf(lvl)%dataptr(mfi)
            ! Phasic densities
            bx=mfi%growntilebox(this%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Liquid
               if (pVF(i,j,k,1).ge.VFlo.and.pQ(i,j,k,1).gt.0.0_WP) then
                  pRHOL(i,j,k,1)=pQ(i,j,k,1)/pVF(i,j,k,1)
               else
                  pRHOL(i,j,k,1)=0.0_WP
               end if
               ! Gas
               if (pVF(i,j,k,1).le.VFhi.and.pQ(i,j,k,2).gt.0.0_WP) then
                  pRHOG(i,j,k,1)=pQ(i,j,k,2)/(1.0_WP-pVF(i,j,k,1))
               else
                  pRHOG(i,j,k,1)=0.0_WP
               end if
            end do; end do; end do
            ! Phasic internal energies from total energies, then phasic and mixture pressures
            bx=mfi%growntilebox(min(this%PL%ng,this%PG%ng,this%P%ng))
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               mk=max(pQ(i,j,k,1:2),0.0_WP)
               Y=pQ(i,j,k,3:4)-mk/max(sum(mk),this%rho_floor)*cell_ke(pQ,pU,pV,pW,i,j,k,this%rho_floor)
               pPL(i,j,k,1)=0.0_WP
               if (pRHOL(i,j,k,1).gt.0.0_WP.and.Y(1).gt.0.0_WP) pPL(i,j,k,1)=this%liq%get_p_from_rho_e(rho=pRHOL(i,j,k,1),e=Y(1)/pQ(i,j,k,1),y=[1.0_WP])
               pPG(i,j,k,1)=0.0_WP
               if (pRHOG(i,j,k,1).gt.0.0_WP.and.Y(2).gt.0.0_WP) pPG(i,j,k,1)=this%gas%get_p_from_rho_e(rho=pRHOG(i,j,k,1),e=Y(2)/pQ(i,j,k,2),y=[1.0_WP])
               pP(i,j,k,1)=pVF(i,j,k,1)*pPL(i,j,k,1)+(1.0_WP-pVF(i,j,k,1))*pPG(i,j,k,1)
            end do; end do; end do
            ! Phasic temperatures
            bx=mfi%growntilebox(min(this%TL%ng,this%TG%ng))
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pTL(i,j,k,1)=0.0_WP
               if (pRHOL(i,j,k,1).gt.0.0_WP) pTL(i,j,k,1)=this%liq%get_T_from_p_rho(p=pPL(i,j,k,1),rho=pRHOL(i,j,k,1),y=[1.0_WP])
               pTG(i,j,k,1)=0.0_WP
               if (pRHOG(i,j,k,1).gt.0.0_WP) pTG(i,j,k,1)=this%gas%get_T_from_p_rho(p=pPG(i,j,k,1),rho=pRHOG(i,j,k,1),y=[1.0_WP])
            end do; end do; end do
            ! Phasic internal energies and mixture speed of sound 
            bx=mfi%growntilebox(min(this%IL%ng,this%IG%ng,this%C%ng))
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Liquid
               CL=0.0_WP; pIL(i,j,k,1)=0.0_WP
               if (pRHOL(i,j,k,1).gt.0.0_WP) then
                  pIL(i,j,k,1)=this%liq%get_e_from_p_rho(p=pPL(i,j,k,1),rho=pRHOL(i,j,k,1),y=[1.0_WP])
                  CL          =this%liq%get_c_from_p_rho(p=pPL(i,j,k,1),rho=pRHOL(i,j,k,1),y=[1.0_WP])
               end if
               ! Gas
               CG=0.0_WP; pIG(i,j,k,1)=0.0_WP
               if (pRHOG(i,j,k,1).gt.0.0_WP) then
                  pIG(i,j,k,1)=this%gas%get_e_from_p_rho(p=pPG(i,j,k,1),rho=pRHOG(i,j,k,1),y=[1.0_WP])
                  CG          =this%gas%get_c_from_p_rho(p=pPG(i,j,k,1),rho=pRHOG(i,j,k,1),y=[1.0_WP])
               end if
               ! Mixture speed of sound
               pC(i,j,k,1)=sqrt((pQ(i,j,k,1)*CL**2+pQ(i,j,k,2)*CG**2)/max(pQ(i,j,k,1)+pQ(i,j,k,2),this%rho_floor))
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      this%wt_prim=this%wt_prim+(MPI_Wtime()-tw0)
   end subroutine get_primitive

   !> Calculate dQdt(1:4) and drhoUVWdt. Geometric VOF/CTU is used in band region, dissipation-free scheme in bulk.
   subroutine get_dQdt(this,dQdt,drhoUdt,drhoVdt,drhoWdt,dt,time)
      use messager,         only: die
      use mpi_f08,          only: MPI_Wtime
      implicit none
      class(amrmpcomp), intent(inout) :: this
      type(amrdata), intent(inout) :: dQdt,drhoUdt,drhoVdt,drhoWdt
      real(WP) :: tw0,tw1
      real(WP), intent(in) :: dt
      real(WP), intent(in) :: time
      integer :: slvl,a
      real(WP) :: dx,dy,dz,dxi,dyi,dzi
      logical :: on_sl
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQold,pVFold,pPLICold,pPLold,pPGold,pVF,pBand,pRHOLold,pRHOGold
      type(fab_ptr), dimension(3) :: pVel,pVelold,pF,pM
      integer, dimension(3) :: halo_lo,halo_hi

      tw0=MPI_Wtime(); tw1=tw0
      ! Build the band on the finest level
      if (this%amr%clvl().lt.this%amr%maxlvl) call die('[amrmpcomp get_dQdt] the finest level must exist')
      slvl=this%amr%maxlvl
      band_cells: block
         integer :: i,j,k
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP), dimension(3,3,3) :: blk
         call this%band%setval(0.0_WP)
         call this%amr%mfiter_build(slvl,mfi)
         do while (mfi%next())
            pVFold=>this%VFold%mf(slvl)%dataptr(mfi)
            pBand =>this%band%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               blk=pVFold(i-1:i+1,j-1:j+1,k-1:k+1,1)
               if (any(blk.ge.VFlo.and.blk.le.VFhi).or.(any(blk.lt.VFlo).and.any(blk.gt.VFhi))) pBand(i,j,k,1)=1.0_WP
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         call this%band%fill_boundary(this%amr%geom(slvl))
      end block band_cells
      dx=this%amr%dx(slvl); dxi=1.0_WP/dx
      dy=this%amr%dy(slvl); dyi=1.0_WP/dy
      dz=this%amr%dz(slvl); dzi=1.0_WP/dz
      this%cnt%nvf_guard=0.0_WP
      this%wt_alloc=this%wt_alloc+(MPI_Wtime()-tw1)

      ! Geometric fluxes on band faces
      tw1=MPI_Wtime()
      geometric_fluxes: block
         use amrvof_geometry, only: tet_sign,tet_map,correct_flux_poly
         integer :: i,j,k,n,nn
         real(WP), dimension(3,9) :: face
         real(WP), dimension(3,4) :: tet
         integer , dimension(3,4) :: ijk
         integer , dimension(3,9) :: fijk
         real(WP), dimension(:,:,:,:), allocatable :: proj
         real(WP), dimension(8) :: Vflux,Vacc
         real(WP), dimension(nSLF) :: Qflux
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVx,pVy,pVz,pFx,pFy,pFz
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: fbx,nbx,vbx
         call this%amr%mfiter_build(slvl,mfi,tiling=.false.)
         do while (mfi%next())
            ! Get data pointers
            pQold   =>this%Qold%mf(slvl)%dataptr(mfi)
            pRHOLold=>this%RHOLold%mf(slvl)%dataptr(mfi); pRHOGold=>this%RHOGold%mf(slvl)%dataptr(mfi)
            pVFold  =>this%VFold%mf(slvl)%dataptr(mfi)
            pPLICold=>this%PLICold%dataptr(mfi)
            pPLold  =>this%PLold%mf(slvl)%dataptr(mfi)
            pPGold  =>this%PGold%mf(slvl)%dataptr(mfi)
            pBand   =>this%band%dataptr(mfi)
            pVel(1)%p=>this%U%mf(slvl)%dataptr(mfi); pVelold(1)%p=>this%Uold%mf(slvl)%dataptr(mfi)
            pVel(2)%p=>this%V%mf(slvl)%dataptr(mfi); pVelold(2)%p=>this%Vold%mf(slvl)%dataptr(mfi)
            pVel(3)%p=>this%W%mf(slvl)%dataptr(mfi); pVelold(3)%p=>this%Wold%mf(slvl)%dataptr(mfi)
            pVx=>this%SLV(1)%dataptr(mfi); pFx=>this%SLF(1)%dataptr(mfi)
            pVy=>this%SLV(2)%dataptr(mfi); pFy=>this%SLF(2)%dataptr(mfi)
            pVz=>this%SLV(3)%dataptr(mfi); pFz=>this%SLF(3)%dataptr(mfi)
            ! Donor range: cells covered by all start-of-step fields, including both time-n velocity faces of the donor
            do a=1,3
               halo_lo(a)=max(lbound(pQold,a),lbound(pVFold,a),lbound(pPLICold,a),lbound(pPLold,a),lbound(pPGold,a), &
               &              lbound(pVelold(1)%p,a),lbound(pVelold(2)%p,a),lbound(pVelold(3)%p,a))
               halo_hi(a)=min(ubound(pQold,a),ubound(pVFold,a),ubound(pPLICold,a),ubound(pPLold,a),ubound(pPGold,a), &
               &              ubound(pVelold(1)%p,a)-eu(a,1),ubound(pVelold(2)%p,a)-eu(a,2),ubound(pVelold(3)%p,a)-eu(a,3))
            end do
            ! Remap vertices of band faces back in time via RK2
            nbx=mfi%nodaltilebox(); nbx%lo=nbx%lo-1; nbx%hi=nbx%hi+1
            allocate(proj(3,nbx%lo(1):nbx%hi(1),nbx%lo(2):nbx%hi(2),nbx%lo(3):nbx%hi(3)))
            do k=nbx%lo(3),nbx%hi(3); do j=nbx%lo(2),nbx%hi(2); do i=nbx%lo(1),nbx%hi(1)
               if (maxval(pBand(i-1:i,j-1:j,k-1:k,1)).gt.0.0_WP) proj(:,i,j,k)=project([this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k,WP)*dz],-dt)
            end do; end do; end do
            ! X-fluxes
            vbx=mfi%nodaltilebox(1); fbx=vbx; fbx%lo=fbx%lo-1; fbx%hi=fbx%hi+1
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               ! Skip if outside band
               if (maxval(pBand(i-1:i,j,k,1)).eq.0.0_WP) cycle
               ! Build flux polyhedron
               face(:,1)=[this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]; face(:,5)=proj(:,i,j  ,k  )
               face(:,2)=[this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]; face(:,6)=proj(:,i,j  ,k+1)
               face(:,3)=[this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]; face(:,7)=proj(:,i,j+1,k+1)
               face(:,4)=[this%amr%xlo+real(i,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]; face(:,8)=proj(:,i,j+1,k  )
               face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
               call correct_flux_poly(poly=face,target_volume=dt*dy*dz*pVel(1)%p(i,j,k,1))
               ! Compute face indices
               do nn=1,9; fijk(:,nn)=floor([(face(1,nn)-this%amr%xlo)*dxi,(face(2,nn)-this%amr%ylo)*dyi,(face(3,nn)-this%amr%zlo)*dzi]); end do
               do nn=1,4; fijk(1,nn)=merge(i-1,i,pVel(1)%p(i,j,k,1).gt.0.0_WP); end do
               ! Decompose into tets, cut, and accumulate
               Vacc=0.0_WP
               pFx(i,j,k,:)=0.0_WP
               do n=1,8
                  do nn=1,4
                     tet(:,nn)=face(:,tet_map(nn,n))
                     ijk(:,nn)=fijk(:,tet_map(nn,n))
                  end do
                  call tet2flux(tet,ijk,Vflux,Qflux)
                  Vacc=Vacc+tet_sign(tet)*Vflux
                  pFx(i,j,k,:)=pFx(i,j,k,:)+tet_sign(tet)*Qflux
               end do
               ! Convert to flux rate
               pFx(i,j,k,:)=-pFx(i,j,k,:)/(dt*dy*dz)
               ! Keep the volume moment fluxes on valid faces
               if (all([i,j,k].ge.vbx%lo).and.all([i,j,k].le.vbx%hi)) pVx(i,j,k,:)=Vacc
            end do; end do; end do
            ! Y-fluxes
            vbx=mfi%nodaltilebox(2); fbx=vbx; fbx%lo=fbx%lo-1; fbx%hi=fbx%hi+1
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               ! Skip if outside band
               if (maxval(pBand(i,j-1:j,k,1)).eq.0.0_WP) cycle
               ! Build flux polyhedron
               face(:,1)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]; face(:,5)=proj(:,i+1,j,k+1)
               face(:,2)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]; face(:,6)=proj(:,i  ,j,k+1)
               face(:,3)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]; face(:,7)=proj(:,i  ,j,k  )
               face(:,4)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]; face(:,8)=proj(:,i+1,j,k  )
               face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
               call correct_flux_poly(poly=face,target_volume=dt*dz*dx*pVel(2)%p(i,j,k,1))
               ! Compute face indices
               do nn=1,9; fijk(:,nn)=floor([(face(1,nn)-this%amr%xlo)*dxi,(face(2,nn)-this%amr%ylo)*dyi,(face(3,nn)-this%amr%zlo)*dzi]); end do
               do nn=1,4; fijk(2,nn)=merge(j-1,j,pVel(2)%p(i,j,k,1).gt.0.0_WP); end do
               ! Decompose into tets, cut, and accumulate
               Vacc=0.0_WP
               pFy(i,j,k,:)=0.0_WP
               do n=1,8
                  do nn=1,4
                     tet(:,nn)=face(:,tet_map(nn,n))
                     ijk(:,nn)=fijk(:,tet_map(nn,n))
                  end do
                  call tet2flux(tet,ijk,Vflux,Qflux)
                  Vacc=Vacc+tet_sign(tet)*Vflux
                  pFy(i,j,k,:)=pFy(i,j,k,:)+tet_sign(tet)*Qflux
               end do
               ! Convert to flux rate
               pFy(i,j,k,:)=-pFy(i,j,k,:)/(dt*dz*dx)
               ! Keep the volume moment fluxes on valid faces
               if (all([i,j,k].ge.vbx%lo).and.all([i,j,k].le.vbx%hi)) pVy(i,j,k,:)=Vacc
            end do; end do; end do
            ! Z-fluxes
            vbx=mfi%nodaltilebox(3); fbx=vbx; fbx%lo=fbx%lo-1; fbx%hi=fbx%hi+1
            do k=fbx%lo(3),fbx%hi(3); do j=fbx%lo(2),fbx%hi(2); do i=fbx%lo(1),fbx%hi(1)
               ! Skip if outside band
               if (maxval(pBand(i,j,k-1:k,1)).eq.0.0_WP) cycle
               ! Build flux polyhedron
               face(:,1)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k,WP)*dz]; face(:,5)=proj(:,i+1,j  ,k)
               face(:,2)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k,WP)*dz]; face(:,6)=proj(:,i  ,j  ,k)
               face(:,3)=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k,WP)*dz]; face(:,7)=proj(:,i  ,j+1,k)
               face(:,4)=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k,WP)*dz]; face(:,8)=proj(:,i+1,j+1,k)
               face(:,9)=0.25_WP*(face(:,5)+face(:,6)+face(:,7)+face(:,8))
               call correct_flux_poly(poly=face,target_volume=dt*dx*dy*pVel(3)%p(i,j,k,1))
               ! Compute face indices
               do nn=1,9; fijk(:,nn)=floor([(face(1,nn)-this%amr%xlo)*dxi,(face(2,nn)-this%amr%ylo)*dyi,(face(3,nn)-this%amr%zlo)*dzi]); end do
               do nn=1,4; fijk(3,nn)=merge(k-1,k,pVel(3)%p(i,j,k,1).gt.0.0_WP); end do
               ! Decompose into tets, cut, and accumulate
               Vacc=0.0_WP
               pFz(i,j,k,:)=0.0_WP
               do n=1,8
                  do nn=1,4
                     tet(:,nn)=face(:,tet_map(nn,n))
                     ijk(:,nn)=fijk(:,tet_map(nn,n))
                  end do
                  call tet2flux(tet,ijk,Vflux,Qflux)
                  Vacc=Vacc+tet_sign(tet)*Vflux
                  pFz(i,j,k,:)=pFz(i,j,k,:)+tet_sign(tet)*Qflux
               end do
               ! Convert to flux rate
               pFz(i,j,k,:)=-pFz(i,j,k,:)/(dt*dx*dy)
               ! Keep the volume moment fluxes on valid faces
               if (all([i,j,k].ge.vbx%lo).and.all([i,j,k].le.vbx%hi)) pVz(i,j,k,:)=Vacc
            end do; end do; end do
            ! Deallocate proj for this tile
            deallocate(proj)
         end do
         call this%amr%mfiter_destroy(mfi)
         ! Handle periodicity
         if (this%amr%xper.or.this%amr%yper.or.this%amr%zper) then
            do a=1,3
               call this%SLF(a)%override_sync(this%amr%geom(slvl))
               call this%SLV(a)%override_sync(this%amr%geom(slvl))
            end do
         end if
      end block geometric_fluxes
      this%wt_slflux=this%wt_slflux+(MPI_Wtime()-tw1)

      ! Update VF and barycenters in the band from the volume moment fluxes
      tw1=MPI_Wtime()
      volume_moments: block
         integer :: i,j,k
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP) :: vol,VFo,Lvol_old,Gvol_old,Lvol_new,Gvol_new,Lvol_flux,Gvol_flux
         real(WP), dimension(3) :: Lbar_old,Gbar_old,Lbar_new,Gbar_new,Lbar_flux,Gbar_flux
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pCL,pCG,pCLold,pCGold,pVx,pVy,pVz
         vol=dx*dy*dz
         call this%amr%mfiter_build(slvl,mfi,tiling=.false.)
         do while (mfi%next())
            ! Get data pointers
            pBand =>this%band%dataptr(mfi)
            pVFold=>this%VFold%mf(slvl)%dataptr(mfi)
            pVF   =>this%VF%mf(slvl)%dataptr(mfi)
            pCL   =>this%CL%dataptr(mfi);    pCG   =>this%CG%dataptr(mfi)
            pCLold=>this%CLold%dataptr(mfi); pCGold=>this%CGold%dataptr(mfi)
            pVx=>this%SLV(1)%dataptr(mfi)
            pVy=>this%SLV(2)%dataptr(mfi)
            pVz=>this%SLV(3)%dataptr(mfi)
            pVel(1)%p=>this%U%mf(slvl)%dataptr(mfi)
            pVel(2)%p=>this%V%mf(slvl)%dataptr(mfi)
            pVel(3)%p=>this%W%mf(slvl)%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Only work in the band
               if (pBand(i,j,k,1).le.0.0_WP) cycle
               ! Old phasic moments
               VFo=pVFold(i,j,k,1)
               Lvol_old=(       VFo)*vol
               Gvol_old=(1.0_WP-VFo)*vol
               Lbar_old=pCLold(i,j,k,1:3)
               Gbar_old=pCGold(i,j,k,1:3)
               ! Net volume flux from SL volume moments
               Lvol_flux=pVx(i+1,j,k, 1 )-pVx(i,j,k, 1 )+pVy(i,j+1,k, 1 )-pVy(i,j,k, 1 )+pVz(i,j,k+1, 1 )-pVz(i,j,k, 1 )
               Gvol_flux=pVx(i+1,j,k, 2 )-pVx(i,j,k, 2 )+pVy(i,j+1,k, 2 )-pVy(i,j,k, 2 )+pVz(i,j,k+1, 2 )-pVz(i,j,k, 2 )
               Lbar_flux=pVx(i+1,j,k,3:5)-pVx(i,j,k,3:5)+pVy(i,j+1,k,3:5)-pVy(i,j,k,3:5)+pVz(i,j,k+1,3:5)-pVz(i,j,k,3:5)
               Gbar_flux=pVx(i+1,j,k,6:8)-pVx(i,j,k,6:8)+pVy(i,j+1,k,6:8)-pVy(i,j,k,6:8)+pVz(i,j,k+1,6:8)-pVz(i,j,k,6:8)
               ! New phasic volumes
               Lvol_new=Lvol_old-Lvol_flux
               Gvol_new=Gvol_old-Gvol_flux
               if (Lvol_new+Gvol_new.gt.vol_eps*vol) then
                  pVF(i,j,k,1)=Lvol_new/(Lvol_new+Gvol_new)
               else
                  pVF(i,j,k,1)=merge(1.0_WP,0.0_WP,Lvol_new.gt.Gvol_new)
                  this%cnt%nvf_guard=this%cnt%nvf_guard+1.0_WP
               end if
               pCL(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
               pCG(i,j,k,1:3)=pCL(i,j,k,1:3)
               ! Clip and update barycenters
               if (pVF(i,j,k,1).lt.VFlo) then
                  pVF(i,j,k,1)=0.0_WP
               else if (pVF(i,j,k,1).gt.VFhi) then
                  pVF(i,j,k,1)=1.0_WP
               else
                  ! Update barycenters from moment conservation and project forward with handling of degenerate cases. 
                  if (Lvol_new/(Lvol_new+Gvol_new).gt.vol_eps) then
                     Lbar_new=(Lbar_old*Lvol_old-Lbar_flux)/Lvol_new
                     if (all(abs(Lbar_new-pCL(i,j,k,1:3)).le.[dx,dy,dz])) pCL(i,j,k,1:3)=project(Lbar_new,dt)
                  end if
                  if (Gvol_new/(Lvol_new+Gvol_new).gt.vol_eps) then
                     Gbar_new=(Gbar_old*Gvol_old-Gbar_flux)/Gvol_new
                     if (all(abs(Gbar_new-pCG(i,j,k,1:3)).le.[dx,dy,dz])) pCG(i,j,k,1:3)=project(Gbar_new,dt)
                  end if
               end if
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         ! Synchronize: VF on all levels, barycenters on the finest
         call this%VF%average_down(); call this%fill(time=time)
      end block volume_moments
      this%wt_vmom=this%wt_vmom+(MPI_Wtime()-tw1)

      ! Compute fluxes for all levels
      tw1=MPI_Wtime()
      compute_fluxes: block
         integer :: lvl,i,j,k
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP) :: dxi,dyi,dzi,div,mu,lam,strain,stress,w,cent,conv,rhop,rhom
         real(WP), dimension(3) :: gradU,tii
         real(WP), dimension(-2: 0) :: wxp,wyp,wzp
         real(WP), dimension(-1:+1) :: wxm,wym,wzm
         real(WP), parameter :: eps=1.0e-15_WP
         logical :: bandx,bandy,bandz,efx,efy,efz
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW,pRHOL,pRHOG,pPL,pPG,pVisc,pBeta
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pTL,pTG,pDiffL,pDiffG
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pFx,pFy,pFz
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pFUx,pFUy,pFUz
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pFVx,pFVy,pFVz
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pFWx,pFWy,pFWz,pP
         ! Zero RHS: the flux divergences are added after C/F averaging
         call dQdt%setval(val=0.0_WP)
         ! Traverse levels
         do lvl=0,this%amr%clvl()
            ! Get mesh size
            dxi=1.0_WP/this%amr%dx(lvl)
            dyi=1.0_WP/this%amr%dy(lvl)
            dzi=1.0_WP/this%amr%dz(lvl)
            ! The band only exists on the finest level
            on_sl=(lvl.eq.slvl)
            ! Loop over all tiles
            call this%amr%mfiter_build(lvl=lvl,mfi=mfi)
            do while (mfi%next())
               ! Get data pointers
               pU=>this%U%mf(lvl)%dataptr(mfi)
               pV=>this%V%mf(lvl)%dataptr(mfi)
               pW=>this%W%mf(lvl)%dataptr(mfi)
               pVF=>this%VF%mf(lvl)%dataptr(mfi)
               pRHOL=>this%RHOL%mf(lvl)%dataptr(mfi); pRHOG=>this%RHOG%mf(lvl)%dataptr(mfi)
               pPL=>this%PL%mf(lvl)%dataptr(mfi); pPG=>this%PG%mf(lvl)%dataptr(mfi)
               pVisc=>this%visc%mf(lvl)%dataptr(mfi)
               pBeta=>this%beta%mf(lvl)%dataptr(mfi)
               pTL=>this%TL%mf(lvl)%dataptr(mfi); pTG=>this%TG%mf(lvl)%dataptr(mfi)
               pDiffL=>this%diffL%mf(lvl)%dataptr(mfi); pDiffG=>this%diffG%mf(lvl)%dataptr(mfi)
               pP=>this%P%mf(lvl)%dataptr(mfi)
               pFx=>this%Fx(lvl)%dataptr(mfi)
               pFy=>this%Fy(lvl)%dataptr(mfi)
               pFz=>this%Fz(lvl)%dataptr(mfi)
               pFUx=>this%FUx(lvl)%dataptr(mfi)
               pFUy=>this%FUy(lvl)%dataptr(mfi)
               pFUz=>this%FUz(lvl)%dataptr(mfi)
               pFVx=>this%FVx(lvl)%dataptr(mfi)
               pFVy=>this%FVy(lvl)%dataptr(mfi)
               pFVz=>this%FVz(lvl)%dataptr(mfi)
               pFWx=>this%FWx(lvl)%dataptr(mfi)
               pFWy=>this%FWy(lvl)%dataptr(mfi)
               pFWz=>this%FWz(lvl)%dataptr(mfi)
               ! Pointers shared with the band momentum flux
               pVel(1)%p=>pU;  pVel(2)%p=>pV;  pVel(3)%p=>pW
               pM(1)%p  =>pFx; pM(2)%p  =>pFy; pM(3)%p  =>pFz
               if (on_sl) then
                  pBand=>this%band%dataptr(mfi)
                  pVFold=>this%VFold%mf(lvl)%dataptr(mfi)
                  do a=1,3; pF(a)%p=>this%SLF(a)%dataptr(mfi); end do
                  pVelold(1)%p=>this%Uold%mf(lvl)%dataptr(mfi)
                  pVelold(2)%p=>this%Vold%mf(lvl)%dataptr(mfi)
                  pVelold(3)%p=>this%Wold%mf(lvl)%dataptr(mfi)
               end if
               bx=mfi%tilebox()
               ! Phasic mass and internal energy fluxes
               ! X-faces
               do k=bx%lo(3)-1,bx%hi(3)+1; do j=bx%lo(2)-1,bx%hi(2)+1; do i=bx%lo(1)-1,bx%hi(1)+2
                  bandx=.false.; if (on_sl) bandx=maxval(pBand(i-1:i,j,k,1)).gt.0.0_WP
                  efx=i.ge.bx%lo(1).and.i.le.bx%hi(1)+1.and.j.ge.bx%lo(2).and.j.le.bx%hi(2).and.k.ge.bx%lo(3).and.k.le.bx%hi(3)
                  if (bandx) then
                     pFx(i,j,k,1:4)=pF(1)%p(i,j,k,1:4)
                  else if (pVF(i,j,k,1).gt.0.5_WP) then
                     ! No gas fluxes
                     pFx(i,j,k,2)=0.0_WP; pFx(i,j,k,4)=0.0_WP
                     ! Liquid mass
                     w=weno_weight((abs(pRHOL(i-1,j,k,1)-pRHOL(i-2,j,k,1))+eps)/(abs(pRHOL(i,j,k,1)-pRHOL(i-1,j,k,1))+eps)); wxp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(pRHOL(i+1,j,k,1)-pRHOL(i  ,j,k,1))+eps)/(abs(pRHOL(i,j,k,1)-pRHOL(i-1,j,k,1))+eps)); wxm=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     rhop=sum(wxp*pRHOL(i-2:i  ,j,k,1)); rhom=sum(wxm*pRHOL(i-1:i+1,j,k,1))
                     pFx(i,j,k,1)=-0.5_WP*(pU(i,j,k,1)+abs(pU(i,j,k,1)))*rhop-0.5_WP*(pU(i,j,k,1)-abs(pU(i,j,k,1)))*rhom
                     ! Liquid internal energy from the reconstructed density (above) and pressure, through the equation of state
                     if (efx) then
                        w=weno_weight((abs(pPL(i-1,j,k,1)-pPL(i-2,j,k,1))+eps)/(abs(pPL(i,j,k,1)-pPL(i-1,j,k,1))+eps)); wxp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pPL(i+1,j,k,1)-pPL(i  ,j,k,1))+eps)/(abs(pPL(i,j,k,1)-pPL(i-1,j,k,1))+eps)); wxm=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                        pFx(i,j,k,3)=-0.5_WP*(pU(i,j,k,1)+abs(pU(i,j,k,1)))*this%liq%get_rhoe_from_p_rho(p=sum(wxp*pPL(i-2:i  ,j,k,1)),rho=rhop,y=[1.0_WP]) &
                        &            -0.5_WP*(pU(i,j,k,1)-abs(pU(i,j,k,1)))*this%liq%get_rhoe_from_p_rho(p=sum(wxm*pPL(i-1:i+1,j,k,1)),rho=rhom,y=[1.0_WP])
                     end if
                  else
                     ! No liquid fluxes
                     pFx(i,j,k,1)=0.0_WP; pFx(i,j,k,3)=0.0_WP
                     ! Gas mass
                     w=weno_weight((abs(pRHOG(i-1,j,k,1)-pRHOG(i-2,j,k,1))+eps)/(abs(pRHOG(i,j,k,1)-pRHOG(i-1,j,k,1))+eps)); wxp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(pRHOG(i+1,j,k,1)-pRHOG(i  ,j,k,1))+eps)/(abs(pRHOG(i,j,k,1)-pRHOG(i-1,j,k,1))+eps)); wxm=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     rhop=sum(wxp*pRHOG(i-2:i  ,j,k,1)); rhom=sum(wxm*pRHOG(i-1:i+1,j,k,1))
                     pFx(i,j,k,2)=-0.5_WP*(pU(i,j,k,1)+abs(pU(i,j,k,1)))*rhop-0.5_WP*(pU(i,j,k,1)-abs(pU(i,j,k,1)))*rhom
                     ! Gas internal energy from the reconstructed density (above) and pressure
                     if (efx) then
                        w=weno_weight((abs(pPG(i-1,j,k,1)-pPG(i-2,j,k,1))+eps)/(abs(pPG(i,j,k,1)-pPG(i-1,j,k,1))+eps)); wxp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pPG(i+1,j,k,1)-pPG(i  ,j,k,1))+eps)/(abs(pPG(i,j,k,1)-pPG(i-1,j,k,1))+eps)); wxm=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                        pFx(i,j,k,4)=-0.5_WP*(pU(i,j,k,1)+abs(pU(i,j,k,1)))*this%gas%get_rhoe_from_p_rho(p=sum(wxp*pPG(i-2:i  ,j,k,1)),rho=rhop,y=[1.0_WP]) &
                        &            -0.5_WP*(pU(i,j,k,1)-abs(pU(i,j,k,1)))*this%gas%get_rhoe_from_p_rho(p=sum(wxm*pPG(i-1:i+1,j,k,1)),rho=rhom,y=[1.0_WP])
                     end if
                  end if
                  if (.not.bandx.and.efx) then
                     if (all(pVF(i-1:i,j,k,1).gt.VFhi).and.all(pTL(i-1:i,j,k,1).gt.0.0_WP)) pFx(i,j,k,3)=pFx(i,j,k,3)+0.5_WP*sum(pDiffL(i-1:i,j,k,1))*dxi*(pTL(i,j,k,1)-pTL(i-1,j,k,1))
                     if (all(pVF(i-1:i,j,k,1).lt.VFlo).and.all(pTG(i-1:i,j,k,1).gt.0.0_WP)) pFx(i,j,k,4)=pFx(i,j,k,4)+0.5_WP*sum(pDiffG(i-1:i,j,k,1))*dxi*(pTG(i,j,k,1)-pTG(i-1,j,k,1))
                  end if
               end do; end do; end do
               ! Y-faces
               do k=bx%lo(3)-1,bx%hi(3)+1; do j=bx%lo(2)-1,bx%hi(2)+2; do i=bx%lo(1)-1,bx%hi(1)+1
                  bandy=.false.; if (on_sl) bandy=maxval(pBand(i,j-1:j,k,1)).gt.0.0_WP
                  efy=j.ge.bx%lo(2).and.j.le.bx%hi(2)+1.and.i.ge.bx%lo(1).and.i.le.bx%hi(1).and.k.ge.bx%lo(3).and.k.le.bx%hi(3)
                  if (bandy) then
                     pFy(i,j,k,1:4)=pF(2)%p(i,j,k,1:4)
                  else if (pVF(i,j,k,1).gt.0.5_WP) then
                     ! No gas fluxes
                     pFy(i,j,k,2)=0.0_WP; pFy(i,j,k,4)=0.0_WP
                     ! Liquid mass
                     w=weno_weight((abs(pRHOL(i,j-1,k,1)-pRHOL(i,j-2,k,1))+eps)/(abs(pRHOL(i,j,k,1)-pRHOL(i,j-1,k,1))+eps)); wyp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(pRHOL(i,j+1,k,1)-pRHOL(i,j  ,k,1))+eps)/(abs(pRHOL(i,j,k,1)-pRHOL(i,j-1,k,1))+eps)); wym=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     rhop=sum(wyp*pRHOL(i,j-2:j  ,k,1)); rhom=sum(wym*pRHOL(i,j-1:j+1,k,1))
                     pFy(i,j,k,1)=-0.5_WP*(pV(i,j,k,1)+abs(pV(i,j,k,1)))*rhop-0.5_WP*(pV(i,j,k,1)-abs(pV(i,j,k,1)))*rhom
                     ! Liquid internal energy from the reconstructed density (above) and pressure
                     if (efy) then
                        w=weno_weight((abs(pPL(i,j-1,k,1)-pPL(i,j-2,k,1))+eps)/(abs(pPL(i,j,k,1)-pPL(i,j-1,k,1))+eps)); wyp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pPL(i,j+1,k,1)-pPL(i,j  ,k,1))+eps)/(abs(pPL(i,j,k,1)-pPL(i,j-1,k,1))+eps)); wym=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                        pFy(i,j,k,3)=-0.5_WP*(pV(i,j,k,1)+abs(pV(i,j,k,1)))*this%liq%get_rhoe_from_p_rho(p=sum(wyp*pPL(i,j-2:j  ,k,1)),rho=rhop,y=[1.0_WP]) &
                        &            -0.5_WP*(pV(i,j,k,1)-abs(pV(i,j,k,1)))*this%liq%get_rhoe_from_p_rho(p=sum(wym*pPL(i,j-1:j+1,k,1)),rho=rhom,y=[1.0_WP])
                     end if
                  else
                     ! No liquid fluxes
                     pFy(i,j,k,1)=0.0_WP; pFy(i,j,k,3)=0.0_WP
                     ! Gas mass
                     w=weno_weight((abs(pRHOG(i,j-1,k,1)-pRHOG(i,j-2,k,1))+eps)/(abs(pRHOG(i,j,k,1)-pRHOG(i,j-1,k,1))+eps)); wyp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(pRHOG(i,j+1,k,1)-pRHOG(i,j  ,k,1))+eps)/(abs(pRHOG(i,j,k,1)-pRHOG(i,j-1,k,1))+eps)); wym=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     rhop=sum(wyp*pRHOG(i,j-2:j  ,k,1)); rhom=sum(wym*pRHOG(i,j-1:j+1,k,1))
                     pFy(i,j,k,2)=-0.5_WP*(pV(i,j,k,1)+abs(pV(i,j,k,1)))*rhop-0.5_WP*(pV(i,j,k,1)-abs(pV(i,j,k,1)))*rhom
                     ! Gas internal energy from the reconstructed density (above) and pressure
                     if (efy) then
                        w=weno_weight((abs(pPG(i,j-1,k,1)-pPG(i,j-2,k,1))+eps)/(abs(pPG(i,j,k,1)-pPG(i,j-1,k,1))+eps)); wyp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pPG(i,j+1,k,1)-pPG(i,j  ,k,1))+eps)/(abs(pPG(i,j,k,1)-pPG(i,j-1,k,1))+eps)); wym=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                        pFy(i,j,k,4)=-0.5_WP*(pV(i,j,k,1)+abs(pV(i,j,k,1)))*this%gas%get_rhoe_from_p_rho(p=sum(wyp*pPG(i,j-2:j  ,k,1)),rho=rhop,y=[1.0_WP]) &
                        &            -0.5_WP*(pV(i,j,k,1)-abs(pV(i,j,k,1)))*this%gas%get_rhoe_from_p_rho(p=sum(wym*pPG(i,j-1:j+1,k,1)),rho=rhom,y=[1.0_WP])
                     end if
                  end if
                  if (.not.bandy.and.efy) then
                     if (all(pVF(i,j-1:j,k,1).gt.VFhi).and.all(pTL(i,j-1:j,k,1).gt.0.0_WP)) pFy(i,j,k,3)=pFy(i,j,k,3)+0.5_WP*sum(pDiffL(i,j-1:j,k,1))*dyi*(pTL(i,j,k,1)-pTL(i,j-1,k,1))
                     if (all(pVF(i,j-1:j,k,1).lt.VFlo).and.all(pTG(i,j-1:j,k,1).gt.0.0_WP)) pFy(i,j,k,4)=pFy(i,j,k,4)+0.5_WP*sum(pDiffG(i,j-1:j,k,1))*dyi*(pTG(i,j,k,1)-pTG(i,j-1,k,1))
                  end if
               end do; end do; end do
               ! Z-faces
               do k=bx%lo(3)-1,bx%hi(3)+2; do j=bx%lo(2)-1,bx%hi(2)+1; do i=bx%lo(1)-1,bx%hi(1)+1
                  bandz=.false.; if (on_sl) bandz=maxval(pBand(i,j,k-1:k,1)).gt.0.0_WP
                  efz=k.ge.bx%lo(3).and.k.le.bx%hi(3)+1.and.i.ge.bx%lo(1).and.i.le.bx%hi(1).and.j.ge.bx%lo(2).and.j.le.bx%hi(2)
                  if (bandz) then
                     pFz(i,j,k,1:4)=pF(3)%p(i,j,k,1:4)
                  else if (pVF(i,j,k,1).gt.0.5_WP) then
                     ! No gas fluxes
                     pFz(i,j,k,2)=0.0_WP; pFz(i,j,k,4)=0.0_WP
                     ! Liquid mass
                     w=weno_weight((abs(pRHOL(i,j,k-1,1)-pRHOL(i,j,k-2,1))+eps)/(abs(pRHOL(i,j,k,1)-pRHOL(i,j,k-1,1))+eps)); wzp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(pRHOL(i,j,k+1,1)-pRHOL(i,j,k  ,1))+eps)/(abs(pRHOL(i,j,k,1)-pRHOL(i,j,k-1,1))+eps)); wzm=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     rhop=sum(wzp*pRHOL(i,j,k-2:k  ,1)); rhom=sum(wzm*pRHOL(i,j,k-1:k+1,1))
                     pFz(i,j,k,1)=-0.5_WP*(pW(i,j,k,1)+abs(pW(i,j,k,1)))*rhop-0.5_WP*(pW(i,j,k,1)-abs(pW(i,j,k,1)))*rhom
                     ! Liquid internal energy from the reconstructed density (above) and pressure
                     if (efz) then
                        w=weno_weight((abs(pPL(i,j,k-1,1)-pPL(i,j,k-2,1))+eps)/(abs(pPL(i,j,k,1)-pPL(i,j,k-1,1))+eps)); wzp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pPL(i,j,k+1,1)-pPL(i,j,k  ,1))+eps)/(abs(pPL(i,j,k,1)-pPL(i,j,k-1,1))+eps)); wzm=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                        pFz(i,j,k,3)=-0.5_WP*(pW(i,j,k,1)+abs(pW(i,j,k,1)))*this%liq%get_rhoe_from_p_rho(p=sum(wzp*pPL(i,j,k-2:k  ,1)),rho=rhop,y=[1.0_WP]) &
                        &            -0.5_WP*(pW(i,j,k,1)-abs(pW(i,j,k,1)))*this%liq%get_rhoe_from_p_rho(p=sum(wzm*pPL(i,j,k-1:k+1,1)),rho=rhom,y=[1.0_WP])
                     end if
                  else
                     ! No liquid fluxes
                     pFz(i,j,k,1)=0.0_WP; pFz(i,j,k,3)=0.0_WP
                     ! Gas mass
                     w=weno_weight((abs(pRHOG(i,j,k-1,1)-pRHOG(i,j,k-2,1))+eps)/(abs(pRHOG(i,j,k,1)-pRHOG(i,j,k-1,1))+eps)); wzp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(pRHOG(i,j,k+1,1)-pRHOG(i,j,k  ,1))+eps)/(abs(pRHOG(i,j,k,1)-pRHOG(i,j,k-1,1))+eps)); wzm=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     rhop=sum(wzp*pRHOG(i,j,k-2:k  ,1)); rhom=sum(wzm*pRHOG(i,j,k-1:k+1,1))
                     pFz(i,j,k,2)=-0.5_WP*(pW(i,j,k,1)+abs(pW(i,j,k,1)))*rhop-0.5_WP*(pW(i,j,k,1)-abs(pW(i,j,k,1)))*rhom
                     ! Gas internal energy from the reconstructed density (above) and pressure
                     if (efz) then
                        w=weno_weight((abs(pPG(i,j,k-1,1)-pPG(i,j,k-2,1))+eps)/(abs(pPG(i,j,k,1)-pPG(i,j,k-1,1))+eps)); wzp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                        w=weno_weight((abs(pPG(i,j,k+1,1)-pPG(i,j,k  ,1))+eps)/(abs(pPG(i,j,k,1)-pPG(i,j,k-1,1))+eps)); wzm=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                        pFz(i,j,k,4)=-0.5_WP*(pW(i,j,k,1)+abs(pW(i,j,k,1)))*this%gas%get_rhoe_from_p_rho(p=sum(wzp*pPG(i,j,k-2:k  ,1)),rho=rhop,y=[1.0_WP]) &
                        &            -0.5_WP*(pW(i,j,k,1)-abs(pW(i,j,k,1)))*this%gas%get_rhoe_from_p_rho(p=sum(wzm*pPG(i,j,k-1:k+1,1)),rho=rhom,y=[1.0_WP])
                     end if
                  end if
                  if (.not.bandz.and.efz) then
                     if (all(pVF(i,j,k-1:k,1).gt.VFhi).and.all(pTL(i,j,k-1:k,1).gt.0.0_WP)) pFz(i,j,k,3)=pFz(i,j,k,3)+0.5_WP*sum(pDiffL(i,j,k-1:k,1))*dzi*(pTL(i,j,k,1)-pTL(i,j,k-1,1))
                     if (all(pVF(i,j,k-1:k,1).lt.VFlo).and.all(pTG(i,j,k-1:k,1).gt.0.0_WP)) pFz(i,j,k,4)=pFz(i,j,k,4)+0.5_WP*sum(pDiffG(i,j,k-1:k,1))*dzi*(pTG(i,j,k,1)-pTG(i,j,k-1,1))
                  end if
               end do; end do; end do
               ! xy-edge (FUy, FVx): convection (centered, or band interpolant) and shear stress
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2)+1; do i=bx%lo(1),bx%hi(1)+1
                  strain=(pU(i,j,k,1)-pU(i,j-1,k,1))*dyi+(pV(i,j,k,1)-pV(i-1,j,k,1))*dxi
                  stress=0.25_WP*sum(pVisc(i-1:i,j-1:j,k,1))*strain
                  cent=0.25_WP*sum(pFy(i-1:i,j,k,1:2))*sum(pU(i,j-1:j,k,1)); conv=conv_flux(1,2,[i,j,k],cent)
                  pFUy(i,j,k,1)=conv+stress
                  cent=0.25_WP*sum(pFx(i,j-1:j,k,1:2))*sum(pV(i-1:i,j,k,1)); conv=conv_flux(2,1,[i,j,k],cent)
                  pFVx(i,j,k,1)=conv+stress
               end do; end do; end do
               ! yz-edge (FVz, FWy): the same
               do k=bx%lo(3),bx%hi(3)+1; do j=bx%lo(2),bx%hi(2)+1; do i=bx%lo(1),bx%hi(1)
                  strain=(pV(i,j,k,1)-pV(i,j,k-1,1))*dzi+(pW(i,j,k,1)-pW(i,j-1,k,1))*dyi
                  stress=0.25_WP*sum(pVisc(i,j-1:j,k-1:k,1))*strain
                  cent=0.25_WP*sum(pFz(i,j-1:j,k,1:2))*sum(pV(i,j,k-1:k,1)); conv=conv_flux(2,3,[i,j,k],cent)
                  pFVz(i,j,k,1)=conv+stress
                  cent=0.25_WP*sum(pFy(i,j,k-1:k,1:2))*sum(pW(i,j-1:j,k,1)); conv=conv_flux(3,2,[i,j,k],cent)
                  pFWy(i,j,k,1)=conv+stress
               end do; end do; end do
               ! zx-edge (FWx, FUz): the same
               do k=bx%lo(3),bx%hi(3)+1; do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)+1
                  strain=(pW(i,j,k,1)-pW(i-1,j,k,1))*dxi+(pU(i,j,k,1)-pU(i,j,k-1,1))*dzi
                  stress=0.25_WP*sum(pVisc(i-1:i,j,k-1:k,1))*strain
                  cent=0.25_WP*sum(pFx(i,j,k-1:k,1:2))*sum(pW(i-1:i,j,k,1)); conv=conv_flux(3,1,[i,j,k],cent)
                  pFWx(i,j,k,1)=conv+stress
                  cent=0.25_WP*sum(pFz(i-1:i,j,k,1:2))*sum(pU(i,j,k-1:k,1)); conv=conv_flux(1,3,[i,j,k],cent)
                  pFUz(i,j,k,1)=conv+stress
               end do; end do; end do
               ! Momentum fluxes from convection (centered, or band interpolant), normal viscous, stress, and pressure. 
               do k=bx%lo(3)-1,bx%hi(3)+1; do j=bx%lo(2)-1,bx%hi(2)+1; do i=bx%lo(1)-1,bx%hi(1)+1
                  gradU=[(pU(i+1,j,k,1)-pU(i,j,k,1))*dxi,(pV(i,j+1,k,1)-pV(i,j,k,1))*dyi,(pW(i,j,k+1,1)-pW(i,j,k,1))*dzi]
                  div=sum(gradU)
                  mu=pVisc(i,j,k,1); lam=pBeta(i,j,k,1)-2.0_WP*mu/3.0_WP
                  tii=2.0_WP*mu*gradU+lam*div
                  if (j.ge.bx%lo(2).and.j.le.bx%hi(2).and.k.ge.bx%lo(3).and.k.le.bx%hi(3)) then
                     cent=0.25_WP*sum(pFx(i:i+1,j,k,1:2))*sum(pU(i:i+1,j,k,1)); conv=conv_flux(1,1,[i,j,k],cent)
                     pFUx(i,j,k,1)=conv+tii(1)-pP(i,j,k,1)
                  end if
                  if (i.ge.bx%lo(1).and.i.le.bx%hi(1).and.k.ge.bx%lo(3).and.k.le.bx%hi(3)) then
                     cent=0.25_WP*sum(pFy(i,j:j+1,k,1:2))*sum(pV(i,j:j+1,k,1)); conv=conv_flux(2,2,[i,j,k],cent)
                     pFVy(i,j,k,1)=conv+tii(2)-pP(i,j,k,1)
                  end if
                  if (i.ge.bx%lo(1).and.i.le.bx%hi(1).and.j.ge.bx%lo(2).and.j.le.bx%hi(2)) then
                     cent=0.25_WP*sum(pFz(i,j,k:k+1,1:2))*sum(pW(i,j,k:k+1,1)); conv=conv_flux(3,3,[i,j,k],cent)
                     pFWz(i,j,k,1)=conv+tii(3)-pP(i,j,k,1)
                  end if
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
         end do
      end block compute_fluxes
      this%wt_fv=this%wt_fv+(MPI_Wtime()-tw1)

      ! Average down all fluxes for C/F consistency
      tw1=MPI_Wtime()
      c_f_consistency: block
         use amrex_interface, only: amrmfab_average_down_cell,amrmfab_average_down_edge,amrmfab_average_down_face
         integer :: lvl
         do lvl=this%amr%clvl(),1,-1
            ! Face-centered fluxes
            call amrmfab_average_down_face(fmf=this%Fx(lvl),cmf=this%Fx(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
            call amrmfab_average_down_face(fmf=this%Fy(lvl),cmf=this%Fy(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
            call amrmfab_average_down_face(fmf=this%Fz(lvl),cmf=this%Fz(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
            ! Cell-centered fluxes
            call amrmfab_average_down_cell(fmf=this%FUx(lvl),cmf=this%FUx(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1),ngcrse=0)
            call amrmfab_average_down_cell(fmf=this%FVy(lvl),cmf=this%FVy(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1),ngcrse=0)
            call amrmfab_average_down_cell(fmf=this%FWz(lvl),cmf=this%FWz(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1),ngcrse=0)
            ! Edge-centered fluxes
            call amrmfab_average_down_edge(fmf=this%FUy(lvl),cmf=this%FUy(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1),ngcrse=0)
            call amrmfab_average_down_edge(fmf=this%FVx(lvl),cmf=this%FVx(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1),ngcrse=0)
            call amrmfab_average_down_edge(fmf=this%FUz(lvl),cmf=this%FUz(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1),ngcrse=0)
            call amrmfab_average_down_edge(fmf=this%FWx(lvl),cmf=this%FWx(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1),ngcrse=0)
            call amrmfab_average_down_edge(fmf=this%FVz(lvl),cmf=this%FVz(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1),ngcrse=0)
            call amrmfab_average_down_edge(fmf=this%FWy(lvl),cmf=this%FWy(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1),ngcrse=0)
         end do
      end block c_f_consistency

      ! Compute divergence to get momentum RHS
      divergence_and_sources: block
         integer :: lvl,i,j,k
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP) :: dxi,dyi,dzi
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pFUx,pFUy,pFUz
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pFVx,pFVy,pFVz
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pFWx,pFWy,pFWz
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pdrhoUdt,pdrhoVdt,pdrhoWdt
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pFx,pFy,pFz,pdQdt
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pPjx,pPjy,pPjz
         logical :: finest
         ! Traverse levels
         do lvl=0,this%amr%clvl()
            ! Get mesh size
            dxi=1.0_WP/this%amr%dx(lvl)
            dyi=1.0_WP/this%amr%dy(lvl)
            dzi=1.0_WP/this%amr%dz(lvl)
            call this%amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               ! Get pointers to data
               pFUx=>this%FUx(lvl)%dataptr(mfi)
               pFUy=>this%FUy(lvl)%dataptr(mfi)
               pFUz=>this%FUz(lvl)%dataptr(mfi)
               pFVx=>this%FVx(lvl)%dataptr(mfi)
               pFVy=>this%FVy(lvl)%dataptr(mfi)
               pFVz=>this%FVz(lvl)%dataptr(mfi)
               pFWx=>this%FWx(lvl)%dataptr(mfi)
               pFWy=>this%FWy(lvl)%dataptr(mfi)
               pFWz=>this%FWz(lvl)%dataptr(mfi)
               pdrhoUdt=>drhoUdt%mf(lvl)%dataptr(mfi)
               pdrhoVdt=>drhoVdt%mf(lvl)%dataptr(mfi)
               pdrhoWdt=>drhoWdt%mf(lvl)%dataptr(mfi)
               pFx=>this%Fx(lvl)%dataptr(mfi)
               pFy=>this%Fy(lvl)%dataptr(mfi)
               pFz=>this%Fz(lvl)%dataptr(mfi)
               pdQdt=>dQdt%mf(lvl)%dataptr(mfi)
               ! Surface tension force, on the finest level only
               finest=lvl.eq.this%amr%maxlvl
               if (finest) then
                  pPjx=>this%Pj(1)%dataptr(mfi); pPjy=>this%Pj(2)%dataptr(mfi); pPjz=>this%Pj(3)%dataptr(mfi)
               end if
               ! U RHS at x-faces
               bx=mfi%nodaltilebox(1)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pdrhoUdt(i,j,k,1)=dxi*(pFUx(i,j,k,1)-pFUx(i-1,j,k,1))+dyi*(pFUy(i,j+1,k,1)-pFUy(i,j,k,1))+dzi*(pFUz(i,j,k+1,1)-pFUz(i,j,k,1))
                  if (finest) pdrhoUdt(i,j,k,1)=pdrhoUdt(i,j,k,1)+pPjx(i,j,k,1)
               end do; end do; end do
               ! V RHS at y-faces
               bx=mfi%nodaltilebox(2)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pdrhoVdt(i,j,k,1)=dxi*(pFVx(i+1,j,k,1)-pFVx(i,j,k,1))+dyi*(pFVy(i,j,k,1)-pFVy(i,j-1,k,1))+dzi*(pFVz(i,j,k+1,1)-pFVz(i,j,k,1))
                  if (finest) pdrhoVdt(i,j,k,1)=pdrhoVdt(i,j,k,1)+pPjy(i,j,k,1)
               end do; end do; end do
               ! W RHS at z-faces
               bx=mfi%nodaltilebox(3)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pdrhoWdt(i,j,k,1)=dxi*(pFWx(i+1,j,k,1)-pFWx(i,j,k,1))+dyi*(pFWy(i,j+1,k,1)-pFWy(i,j,k,1))+dzi*(pFWz(i,j,k,1)-pFWz(i,j,k-1,1))
                  if (finest) pdrhoWdt(i,j,k,1)=pdrhoWdt(i,j,k,1)+pPjz(i,j,k,1)
               end do; end do; end do
               ! Phasic masses and total energies at cell centers
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pdQdt(i,j,k,1:4)=pdQdt(i,j,k,1:4)+dxi*(pFx(i+1,j,k,1:4)-pFx(i,j,k,1:4))+dyi*(pFy(i,j+1,k,1:4)-pFy(i,j,k,1:4))+dzi*(pFz(i,j,k+1,1:4)-pFz(i,j,k,1:4))
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
         end do
      end block divergence_and_sources
      this%wt_div=this%wt_div+(MPI_Wtime()-tw1)

      this%wt_dQdt=this%wt_dQdt+(MPI_Wtime()-tw0)

   contains

      !> WENO switch function
      real(WP) function weno_weight(ratio)
         implicit none
         real(WP), intent(in) :: ratio
         real(WP), parameter :: lambda=0.13_WP
         real(WP), parameter :: delta=0.01_WP
         weno_weight=(1.0_WP-tanh((ratio-lambda)/delta))/3.0_WP+(1.0_WP-tanh((ratio-1.0_WP/lambda)/delta))/6.0_WP
      end function weno_weight

      ! Band momentum flux helpers
      ! a: momentum axis; 
      ! b: normal of the face crossed; 
      ! n: index triple. 
      ! A dual face is a face of a staggered momentum CV: the cell centre if a=b, an edge otherwise

      !> Checks whether b-face borders a band cell
      logical function band_face(b,n)
         implicit none
         integer, intent(in) :: b
         integer, dimension(3), intent(in) :: n
         integer, dimension(3) :: m
         m=n-eu(:,b)
         band_face=max(pBand(m(1),m(2),m(3),1),pBand(n(1),n(2),n(3),1)).gt.0.0_WP
      end function band_face

      !> Gets indices for the dual face normal to b
      subroutine dual_face(a,b,n,hlo,hhi,cv1,cv2)
         implicit none
         integer, intent(in) :: a,b
         integer, dimension(3), intent(in) :: n
         integer, dimension(3), intent(out) :: hlo,hhi,cv1,cv2
         if (a.eq.b) then
            hlo=n; hhi=n+eu(:,a); cv1=n; cv2=n+eu(:,a)
         else
            hlo=n-eu(:,a); hhi=n; cv1=n-eu(:,b); cv2=n
         end if
      end subroutine dual_face

      !> Convective a-momentum flux through the dual face normal to b at n
      !> Either uses centered flux or build up band flux
      real(WP) function conv_flux(a,b,n,centered) result(flux)
         implicit none
         integer, intent(in) :: a,b
         integer, dimension(3), intent(in) :: n
         real(WP), intent(in) :: centered
         integer, dimension(3) :: hlo,hhi,cv1,cv2
         logical :: sllo,slhi
         flux=centered
         if (.not.on_sl) return
         call dual_face(a,b,n,hlo,hhi,cv1,cv2)
         sllo=band_face(b,hlo); slhi=band_face(b,hhi)
         if (sllo.or.slhi) flux=dual_flux(a,b,n,sllo,slhi)
      end function conv_flux

      !> Band momentum flux through a dual face: CTU or MC per phase
      real(WP) function dual_flux(a,b,n,sllo,slhi) result(flux)
         implicit none
         integer, intent(in) :: a,b
         integer, dimension(3), intent(in) :: n
         logical, intent(in) :: sllo,slhi
         integer, dimension(3) :: hlo,hhi,cv1,cv2
         integer :: phase
         real(WP) :: phi,flo,fhi
         logical :: use_ctu
         ! Contributing b-faces and separated staggered CVs
         call dual_face(a,b,n,hlo,hhi,cv1,cv2)
         ! MC velocity, transported by the phasic mass fluxes where CTU is not used
         phi=mc_velocity(a,b,n,cv1,cv2)
         ! Per-phase MC or CTU flux
         flux=0.0_WP
         do phase=1,2
            use_ctu=min(alpha(phase,a,cv1),alpha(phase,a,cv2)).lt.this%alpha_ctu
            if (use_ctu) then
               ! Lower half-face: donor moment with the velocity of the donor's upper a-face, or, off the band, the whole
               ! flux from the upwind cell with the time-n velocity of its upper a-face
               if (sllo) then
                  flo=0.5_WP*flx(b,hlo,ictu(phase,2,a))
               else
                  flo=0.5_WP*mflx(b,hlo,phase)*velold(a,upwind(b,hlo)+eu(:,a))
               end if
               ! Upper half-face: the same with the lower a-face
               if (slhi) then
                  fhi=0.5_WP*flx(b,hhi,ictu(phase,1,a))
               else
                  fhi=0.5_WP*mflx(b,hhi,phase)*velold(a,upwind(b,hhi))
               end if
            else
               flo=0.5_WP*mflx(b,hlo,phase)*phi
               fhi=0.5_WP*mflx(b,hhi,phase)*phi
            end if
            flux=flux+flo+fhi
         end do
      end function dual_flux

      !> Construct MC flux interpolant velocity using one downwind and two upwind values along line traced back through face
      real(WP) function mc_velocity(a,b,n,cv1,cv2) result(phi)
         implicit none
         integer, intent(in) :: a,b
         integer, dimension(3), intent(in) :: n,cv1,cv2
         integer, dimension(3) :: fdown,fup,fupup
         real(WP), dimension(3) :: h,xg,ug,off
         real(WP) :: nu,uD,uU,uUU,r,psi
         ! Dual face centroid: cell center if a=b, otherwise nodal along a and b
         h=[dx,dy,dz]
         xg=[this%amr%xlo,this%amr%ylo,this%amr%zlo]+(real(n,WP)+0.5_WP)*h
         if (a.ne.b) then
            xg(a)=xg(a)-0.5_WP*h(a)
            xg(b)=xg(b)-0.5_WP*h(b)
         end if
         ! Side, CFL number, and line from the current velocity
         ug=interp_velocity(xg)
         nu=dt*abs(ug(b))/h(b)
         if (ug(b).ge.0.0_WP) then
            fdown=cv2; fup=cv1; fupup=cv1-eu(:,b)
         else
            fdown=cv1; fup=cv2; fupup=cv2+eu(:,b)
         end if
         ! Transverse offset of the line on the upwind plane, -dt*u_g/(1+nu), and twice that on the next plane
         off=-dt*ug/(1.0_WP+nu)
         ! Downwind and upwind values of the time-n velocity, combined with the MC limiter
         uD =velold(a,fdown)
         uU =plane_value(a,b,fup  ,       off)
         uUU=plane_value(a,b,fupup,2.0_WP*off)
         phi=uU
         if (uD.eq.uU) return
         r=(uU-uUU)/(uD-uU)
         psi=max(0.0_WP,min(2.0_WP*r,0.5_WP*(1.0_WP+r),2.0_WP))
         phi=uU+0.5_WP*(1.0_WP-nu)*psi*(uD-uU)
      end function mc_velocity
   
      !>Time-n a-velocity on the plane normal to b through a-face f, bilinearly interpolated at the transverse offset off
      real(WP) function plane_value(a,b,f,off) result(val)
         implicit none
         integer, intent(in) :: a,b
         integer, dimension(3), intent(in) :: f
         real(WP), dimension(3), intent(in) :: off
         integer :: t1,t2,i1,i2
         integer, dimension(3) :: m
         real(WP) :: q1,q2,w1,w2
         real(WP), dimension(3) :: h
         h=[dx,dy,dz]
         t1=1+mod(b,3); t2=1+mod(b+1,3)
         q1=off(t1)/h(t1); i1=floor(q1); w1=q1-real(i1,WP)
         q2=off(t2)/h(t2); i2=floor(q2); w2=q2-real(i2,WP)
         m=f+i1*eu(:,t1)+i2*eu(:,t2)
         val=(1.0_WP-w1)*(1.0_WP-w2)*velold(a,m)+w1*(1.0_WP-w2)*velold(a,m+eu(:,t1)) &
         &  +(1.0_WP-w1)*       w2 *velold(a,m+eu(:,t2))+w1*w2*velold(a,m+eu(:,t1)+eu(:,t2))
      end function plane_value

      !> Semi-Lagrangian band flux component comp for CTU donor moments
      real(WP) function flx(b,n,comp)
         implicit none
         integer, intent(in) :: b,comp
         integer, dimension(3), intent(in) :: n
         flx=pF(b)%p(n(1),n(2),n(3),comp)
      end function flx

      !> Mass flux of phase at b-face n (semi-Lagrangian in the band, WENO elsewhere)
      real(WP) function mflx(b,n,phase)
         implicit none
         integer, intent(in) :: b,phase
         integer, dimension(3), intent(in) :: n
         mflx=pM(b)%p(n(1),n(2),n(3),phase)
      end function mflx

      !> Time-n face velocity of axis a at a-face
      real(WP) function velold(a,n)
         implicit none
         integer, intent(in) :: a
         integer, dimension(3), intent(in) :: n
         velold=pVelold(a)%p(n(1),n(2),n(3),1)
      end function velold

      !> Upwind cell of b-face n, from the sign of the current velocity
      function upwind(b,n) result(c)
         implicit none
         integer, intent(in) :: b
         integer, dimension(3), intent(in) :: n
         integer, dimension(3) :: c
         c=n
         if (pVel(b)%p(n(1),n(2),n(3),1).gt.0.0_WP) c=n-eu(:,b)
      end function upwind

      !> Time-n phase fraction of the staggered CV of a-face n, compared with alpha_ctu by the CTU switch
      real(WP) function alpha(phase,a,n)
         implicit none
         integer, intent(in) :: phase,a
         integer, dimension(3), intent(in) :: n
         integer, dimension(3) :: m
         m=n-eu(:,a)
         alpha=0.5_WP*(pVFold(m(1),m(2),m(3),1)+pVFold(n(1),n(2),n(3),1))
         if (phase.eq.2) alpha=1.0_WP-alpha
      end function alpha

      !> Recursive subroutine that cuts a tet by grid planes to compute volume and phasic fluxes
      recursive subroutine tet2flux(mytet,myind,myVflux,myQflux)
         use amrvof_geometry, only: cut_side,cut_v1,cut_v2,cut_vtet,cut_ntets,cut_nvert
         implicit none
         real(WP), dimension(3,4), intent(in) :: mytet
         integer,  dimension(3,4), intent(in) :: myind
         real(WP), dimension(8),    intent(out) :: myVflux
         real(WP), dimension(nSLF), intent(out) :: myQflux
         integer :: dir,cut_ind,icase,n1,n2,v1,v2
         real(WP), dimension(4) :: dd
         real(WP), dimension(3,8) :: vert
         integer,  dimension(3,8,2) :: vert_ind
         real(WP) :: mu,my_vol
         real(WP), dimension(3,4) :: newtet
         integer,  dimension(3,4) :: newind
         real(WP), dimension(3) :: a1,b1,c1
         real(WP), dimension(8) :: subVflux
         real(WP), dimension(nSLF) :: subQflux
         real(WP) :: xcut,ycut,zcut

         myVflux=0.0_WP
         myQflux=0.0_WP

         ! Determine if tet spans multiple cells and needs cutting
         if (maxval(myind(1,:))-minval(myind(1,:)).gt.0) then
            dir=1; cut_ind=maxval(myind(1,:))
            xcut=this%amr%xlo+real(cut_ind,WP)*dx
            dd(:)=mytet(1,:)-xcut
         else if (maxval(myind(2,:))-minval(myind(2,:)).gt.0) then
            dir=2; cut_ind=maxval(myind(2,:))
            ycut=this%amr%ylo+real(cut_ind,WP)*dy
            dd(:)=mytet(2,:)-ycut
         else if (maxval(myind(3,:))-minval(myind(3,:)).gt.0) then
            dir=3; cut_ind=maxval(myind(3,:))
            zcut=this%amr%zlo+real(cut_ind,WP)*dz
            dd(:)=mytet(3,:)-zcut
         else
            ! All vertices in same cell - cut by PLIC and return
            call tet2flux_plic(mytet,myind(1,1),myind(2,1),myind(3,1),myVflux,myQflux)
            return
         end if

         ! Find cut case (1-indexed: 1-16)
         icase=1+int(0.5_WP+sign(0.5_WP,dd(1))) &
         &    +2*int(0.5_WP+sign(0.5_WP,dd(2))) &
         &    +4*int(0.5_WP+sign(0.5_WP,dd(3))) &
         &    +8*int(0.5_WP+sign(0.5_WP,dd(4)))

         ! Copy vertices and indices
         do n1=1,4
            vert(:,n1)=mytet(:,n1)
            vert_ind(:,n1,1)=myind(:,n1)
            vert_ind(:,n1,2)=myind(:,n1)
            vert_ind(dir,n1,1)=min(vert_ind(dir,n1,1),cut_ind-1)
            vert_ind(dir,n1,2)=max(vert_ind(dir,n1,1),cut_ind)
         end do

         ! Create interpolated vertices on cut plane
         do n1=1,cut_nvert(icase)
            v1=cut_v1(n1,icase); v2=cut_v2(n1,icase)
            mu=min(1.0_WP,max(0.0_WP,-dd(v1)/(sign(abs(dd(v2)-dd(v1))+epsilon(1.0_WP),dd(v2)-dd(v1)))))
            vert(:,4+n1)=(1.0_WP-mu)*vert(:,v1)+mu*vert(:,v2)
            vert_ind(1,4+n1,1)=floor((vert(1,4+n1)-this%amr%xlo)*dxi)
            vert_ind(2,4+n1,1)=floor((vert(2,4+n1)-this%amr%ylo)*dyi)
            vert_ind(3,4+n1,1)=floor((vert(3,4+n1)-this%amr%zlo)*dzi)
            vert_ind(:,4+n1,1)=max(vert_ind(:,4+n1,1),min(vert_ind(:,v1,1),vert_ind(:,v2,1)))
            vert_ind(:,4+n1,1)=min(vert_ind(:,4+n1,1),max(vert_ind(:,v1,1),vert_ind(:,v2,1)))
            vert_ind(:,4+n1,2)=vert_ind(:,4+n1,1)
            vert_ind(dir,4+n1,1)=cut_ind-1
            vert_ind(dir,4+n1,2)=cut_ind
         end do

         ! Create and process sub-tets
         do n1=1,cut_ntets(icase)
            do n2=1,4
               newtet(:,n2)=vert(:,cut_vtet(n2,n1,icase))
               newind(:,n2)=vert_ind(:,cut_vtet(n2,n1,icase),cut_side(n1,icase))
            end do
            a1=newtet(:,1)-newtet(:,4)
            b1=newtet(:,2)-newtet(:,4)
            c1=newtet(:,3)-newtet(:,4)
            my_vol=abs(a1(1)*(b1(2)*c1(3)-c1(2)*b1(3))-a1(2)*(b1(1)*c1(3)-c1(1)*b1(3))+a1(3)*(b1(1)*c1(2)-c1(1)*b1(2)))/6.0_WP
            if (my_vol.lt.VFlo*dx*dy*dz) cycle
            call tet2flux(newtet,newind,subVflux,subQflux)
            myVflux=myVflux+subVflux
            myQflux=myQflux+subQflux
         end do

      end subroutine tet2flux

      !> Cut a tet inside donor cell (i0,j0,k0) by its old PLIC into phasic fluxes and CTU donor moments
      subroutine tet2flux_plic(mytet,i0,j0,k0,myVflux,myQflux)
         use amrvof_geometry, only: cut_v1,cut_v2,cut_vtet,cut_ntets,cut_nvert,cut_nntet,tet_vol
         implicit none
         real(WP), dimension(3,4), intent(in) :: mytet
         integer,  intent(in) :: i0,j0,k0
         real(WP), dimension(8),    intent(out) :: myVflux
         real(WP), dimension(nSLF), intent(out) :: myQflux
         integer :: icase,n1,v1,v2
         real(WP), dimension(4) :: dd
         real(WP), dimension(3,8) :: vert
         real(WP), dimension(3) :: a1,b1,c1,bary,normal,bary_tot
         real(WP) :: mu,my_vol,dist,VF0,vol_tot

         ! Zero out flux arrays
         myVflux=0.0_WP
         myQflux=0.0_WP

         ! Donor and its time-n velocity faces must be within the old-state halo
         if (any([i0,j0,k0].lt.halo_lo).or.any([i0,j0,k0].gt.halo_hi)) call die('[amrmpcomp get_dQdt] donor outside halo: reduce dt')

         ! Get old VF for this cell
         VF0=pVFold(i0,j0,k0,1)

         ! Tet volume and barycenter
         vol_tot=abs(tet_vol(mytet))
         bary_tot=0.25_WP*(mytet(:,1)+mytet(:,2)+mytet(:,3)+mytet(:,4))

         ! Pure cell shortcut
         if (pPLICold(i0,j0,k0,4).gt.+1.0e9_WP) then
            ! Pure liquid
            myVflux( 1 )=vol_tot
            myVflux(3:5)=vol_tot*bary_tot
            myQflux(1:2)=vol_tot*pQold(i0,j0,k0,1:2)
            myQflux(3:4)=vol_tot*[rhoe_old(1,i0,j0,k0),rhoe_old(2,i0,j0,k0)]
            call donor_moments(i0,j0,k0,myQflux)
            return
         else if (pPLICold(i0,j0,k0,4).lt.-1.0e9_WP) then
            ! Pure gas
            myVflux( 2 )=vol_tot
            myVflux(6:8)=vol_tot*bary_tot
            myQflux(1:2)=vol_tot*pQold(i0,j0,k0,1:2)
            myQflux(3:4)=vol_tot*[rhoe_old(1,i0,j0,k0),rhoe_old(2,i0,j0,k0)]
            call donor_moments(i0,j0,k0,myQflux)
            return
         end if

         ! Get PLIC from this cell
         normal=pPLICold(i0,j0,k0,1:3)
         dist=pPLICold(i0,j0,k0,4)

         ! Compute signed distance to plane for each vertex
         dd(1)=normal(1)*mytet(1,1)+normal(2)*mytet(2,1)+normal(3)*mytet(3,1)-dist
         dd(2)=normal(1)*mytet(1,2)+normal(2)*mytet(2,2)+normal(3)*mytet(3,2)-dist
         dd(3)=normal(1)*mytet(1,3)+normal(2)*mytet(2,3)+normal(3)*mytet(3,3)-dist
         dd(4)=normal(1)*mytet(1,4)+normal(2)*mytet(2,4)+normal(3)*mytet(3,4)-dist

         ! Find cut case
         icase=1+int(0.5_WP+sign(0.5_WP,dd(1))) &
         &    +2*int(0.5_WP+sign(0.5_WP,dd(2))) &
         &    +4*int(0.5_WP+sign(0.5_WP,dd(3))) &
         &    +8*int(0.5_WP+sign(0.5_WP,dd(4)))

         ! Copy vertices
         vert(:,1:4)=mytet(:,1:4)

         ! Create interpolated vertices on cut plane
         do n1=1,cut_nvert(icase)
            v1=cut_v1(n1,icase); v2=cut_v2(n1,icase)
            mu=min(1.0_WP,max(0.0_WP,-dd(v1)/(sign(abs(dd(v2)-dd(v1))+epsilon(1.0_WP),dd(v2)-dd(v1)))))
            vert(:,4+n1)=(1.0_WP-mu)*vert(:,v1)+mu*vert(:,v2)
         end do

         ! Cut the minority phase (safer as we subtract small from large)
         if (VF0.gt.0.5_WP) then
            ! Liquid is dominant: compute gas directly
            do n1=1,cut_nntet(icase)-1
               a1=vert(:,cut_vtet(1,n1,icase))-vert(:,cut_vtet(4,n1,icase))
               b1=vert(:,cut_vtet(2,n1,icase))-vert(:,cut_vtet(4,n1,icase))
               c1=vert(:,cut_vtet(3,n1,icase))-vert(:,cut_vtet(4,n1,icase))
               my_vol=abs(a1(1)*(b1(2)*c1(3)-c1(2)*b1(3))-a1(2)*(b1(1)*c1(3)-c1(1)*b1(3))+a1(3)*(b1(1)*c1(2)-c1(1)*b1(2)))/6.0_WP
               bary=0.25_WP*(vert(:,cut_vtet(1,n1,icase))+vert(:,cut_vtet(2,n1,icase)) &
               &            +vert(:,cut_vtet(3,n1,icase))+vert(:,cut_vtet(4,n1,icase)))
               myVflux( 2 )=myVflux( 2 )+my_vol
               myVflux(6:8)=myVflux(6:8)+my_vol*bary
            end do
            ! Liquid = total - gas
            myVflux( 1 )=vol_tot-myVflux( 2 )
            myVflux(3:5)=vol_tot*bary_tot-myVflux(6:8)
         else
            ! Gas is dominant: compute liquid directly
            do n1=cut_ntets(icase),cut_nntet(icase),-1
               a1=vert(:,cut_vtet(1,n1,icase))-vert(:,cut_vtet(4,n1,icase))
               b1=vert(:,cut_vtet(2,n1,icase))-vert(:,cut_vtet(4,n1,icase))
               c1=vert(:,cut_vtet(3,n1,icase))-vert(:,cut_vtet(4,n1,icase))
               my_vol=abs(a1(1)*(b1(2)*c1(3)-c1(2)*b1(3))-a1(2)*(b1(1)*c1(3)-c1(1)*b1(3))+a1(3)*(b1(1)*c1(2)-c1(1)*b1(2)))/6.0_WP
               bary=0.25_WP*(vert(:,cut_vtet(1,n1,icase))+vert(:,cut_vtet(2,n1,icase)) &
               &            +vert(:,cut_vtet(3,n1,icase))+vert(:,cut_vtet(4,n1,icase)))
               myVflux( 1 )=myVflux( 1 )+my_vol
               myVflux(3:5)=myVflux(3:5)+my_vol*bary
            end do
            ! Gas = total - liquid
            myVflux( 2 )=vol_tot-myVflux( 1 )
            myVflux(6:8)=vol_tot*bary_tot-myVflux(3:5)
         end if

         ! Phasic masses and internal energies from the donor's old phasic state
         if (VF0.ge.VFlo) then
            myQflux(1)=myVflux(1)*pQold(i0,j0,k0,1)/VF0
            myQflux(3)=myVflux(1)*rhoe_old(1,i0,j0,k0)
         end if
         if (VF0.le.VFhi) then
            myQflux(2)=myVflux(2)*pQold(i0,j0,k0,2)/(1.0_WP-VF0)
            myQflux(4)=myVflux(2)*rhoe_old(2,i0,j0,k0)
         end if
         call donor_moments(i0,j0,k0,myQflux)

      end subroutine tet2flux_plic

      !> Start-of-step internal energy per unit volume of phase n (1 liquid, 2 gas) at cell (i,j,k), from its density and
      !> pressure; zero where the phase is absent
      real(WP) function rhoe_old(n,i,j,k)
         implicit none
         integer, intent(in) :: n,i,j,k
         rhoe_old=0.0_WP
         if (n.eq.1) then
            if (pRHOLold(i,j,k,1).gt.0.0_WP) rhoe_old=this%liq%get_rhoe_from_p_rho(p=pPLold(i,j,k,1),rho=pRHOLold(i,j,k,1),y=[1.0_WP])
         else
            if (pRHOGold(i,j,k,1).gt.0.0_WP) rhoe_old=this%gas%get_rhoe_from_p_rho(p=pPGold(i,j,k,1),rho=pRHOGold(i,j,k,1),y=[1.0_WP])
         end if
      end function rhoe_old

      !> CTU donor moments of donor (i0,j0,k0): per phase and axis a, the phasic mass it gives times the time-n velocity of its lower and upper a-faces
      subroutine donor_moments(i0,j0,k0,myQflux)
         implicit none
         integer, intent(in) :: i0,j0,k0
         real(WP), dimension(nSLF), intent(inout) :: myQflux
         integer :: phase,ax
         integer, dimension(3) :: up
         do phase=1,2
            do ax=1,3
               up=[i0,j0,k0]+eu(:,ax)
               myQflux(ictu(phase,1,ax))=myQflux(phase)*pVelold(ax)%p(i0,j0,k0,1)
               myQflux(ictu(phase,2,ax))=myQflux(phase)*pVelold(ax)%p(up(1),up(2),up(3),1)
            end do
         end do
      end subroutine donor_moments

      !> RK2 vertex projection with the current face velocity 
      function project(p1,mydt) result(p2)
         implicit none
         real(WP), dimension(3), intent(in) :: p1
         real(WP), dimension(3)             :: p2
         real(WP),               intent(in) :: mydt
         p2=p1+mydt*interp_velocity(        p1    )
         p2=p1+mydt*interp_velocity(0.5_WP*(p1+p2))
      end function project

      !> Trilinear interpolation of the current staggered velocity at pos
      function interp_velocity(pos) result(vel3)
         implicit none
         real(WP), dimension(3), intent(in) :: pos
         real(WP), dimension(3) :: vel3
         integer  :: ipc,jpc,kpc   ! Cell-centered indices
         integer  :: ipu,jpv,kpw   ! Face-centered indices
         real(WP) :: wxc1,wyc1,wzc1,wxc2,wyc2,wzc2  ! Cell-centered weights
         real(WP) :: wxu1,wyv1,wzw1,wxu2,wyv2,wzw2  ! Face-centered weights
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW
         pU=>pVel(1)%p; pV=>pVel(2)%p; pW=>pVel(3)%p
         ! Compute raw indices
         ipc=floor((pos(1)-this%amr%xlo)*dxi-0.5_WP)
         jpc=floor((pos(2)-this%amr%ylo)*dyi-0.5_WP)
         kpc=floor((pos(3)-this%amr%zlo)*dzi-0.5_WP)
         ipu=floor((pos(1)-this%amr%xlo)*dxi)
         jpv=floor((pos(2)-this%amr%ylo)*dyi)
         kpw=floor((pos(3)-this%amr%zlo)*dzi)
         ! Cell-centered weights
         wxc1=(pos(1)-(this%amr%xlo+(real(ipc,WP)+0.5_WP)*dx))*dxi
         wyc1=(pos(2)-(this%amr%ylo+(real(jpc,WP)+0.5_WP)*dy))*dyi
         wzc1=(pos(3)-(this%amr%zlo+(real(kpc,WP)+0.5_WP)*dz))*dzi
         wxc1=max(0.0_WP,min(1.0_WP,wxc1)); wxc2=1.0_WP-wxc1
         wyc1=max(0.0_WP,min(1.0_WP,wyc1)); wyc2=1.0_WP-wyc1
         wzc1=max(0.0_WP,min(1.0_WP,wzc1)); wzc2=1.0_WP-wzc1
         ! Face-centered weights
         wxu1=(pos(1)-(this%amr%xlo+real(ipu,WP)*dx))*dxi
         wyv1=(pos(2)-(this%amr%ylo+real(jpv,WP)*dy))*dyi
         wzw1=(pos(3)-(this%amr%zlo+real(kpw,WP)*dz))*dzi
         wxu1=max(0.0_WP,min(1.0_WP,wxu1)); wxu2=1.0_WP-wxu1
         wyv1=max(0.0_WP,min(1.0_WP,wyv1)); wyv2=1.0_WP-wyv1
         wzw1=max(0.0_WP,min(1.0_WP,wzw1)); wzw2=1.0_WP-wzw1
         ! U at x-faces: face-centered in x, cell-centered in y,z
         vel3(1)=wzc1*(wyc1*(wxu1*pU(ipu+1,jpc+1,kpc+1,1)+wxu2*pU(ipu,jpc+1,kpc+1,1)) +&
         &             wyc2*(wxu1*pU(ipu+1,jpc  ,kpc+1,1)+wxu2*pU(ipu,jpc  ,kpc+1,1)))+&
         &       wzc2*(wyc1*(wxu1*pU(ipu+1,jpc+1,kpc  ,1)+wxu2*pU(ipu,jpc+1,kpc  ,1)) +&
         &             wyc2*(wxu1*pU(ipu+1,jpc  ,kpc  ,1)+wxu2*pU(ipu,jpc  ,kpc  ,1)))
         ! V at y-faces: cell-centered in x, face-centered in y, cell-centered in z
         vel3(2)=wzc1*(wyv1*(wxc1*pV(ipc+1,jpv+1,kpc+1,1)+wxc2*pV(ipc,jpv+1,kpc+1,1)) +&
         &             wyv2*(wxc1*pV(ipc+1,jpv  ,kpc+1,1)+wxc2*pV(ipc,jpv  ,kpc+1,1)))+&
         &       wzc2*(wyv1*(wxc1*pV(ipc+1,jpv+1,kpc  ,1)+wxc2*pV(ipc,jpv+1,kpc  ,1)) +&
         &             wyv2*(wxc1*pV(ipc+1,jpv  ,kpc  ,1)+wxc2*pV(ipc,jpv  ,kpc  ,1)))
         ! W at z-faces: cell-centered in x,y, face-centered in z
         vel3(3)=wzw1*(wyc1*(wxc1*pW(ipc+1,jpc+1,kpw+1,1)+wxc2*pW(ipc,jpc+1,kpw+1,1)) +&
         &             wyc2*(wxc1*pW(ipc+1,jpc  ,kpw+1,1)+wxc2*pW(ipc,jpc  ,kpw+1,1)))+&
         &       wzw2*(wyc1*(wxc1*pW(ipc+1,jpc+1,kpw  ,1)+wxc2*pW(ipc,jpc+1,kpw  ,1)) +&
         &             wyc2*(wxc1*pW(ipc+1,jpc  ,kpw  ,1)+wxc2*pW(ipc,jpc  ,kpw  ,1)))
      end function interp_velocity

   end subroutine get_dQdt

   !> Build up the rest of phasic total energy using the updated face velocity, mass and momentum fluxes, and pressure and viscous work
   !> Bulk phase simplifies to adding centered KE contribution alongside the pressure and viscous work directly
   !> Mixture region requires redistributing the KE contribution by mass while viscous and pressure work is weighted by VF
   subroutine get_energy(this,dt,time)
      use mpi_f08,         only: MPI_Wtime
      use amrex_interface, only: amrmfab_average_down_face
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(in) :: time
      real(WP) :: tw0
      integer :: lvl,i,j,k,n1,n2,rerr
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP) :: dxi,dyi,dzi,div,dss,dW,Kn,Cap
      real(WP), dimension(3) :: gradU
      real(WP), dimension(2) :: S,mk,mkn,Yb
      real(WP), dimension(4) :: q
      logical :: inband
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pQold,pVF,pPL,pPG,pP,pU,pV,pW,pUold,pVold,pWold,pBand
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pFx,pFy,pFz,pEx,pEy,pEz,pUm,pVm,pWm,pVisc,pBeta
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pFUx,pFUy,pFUz,pFVx,pFVy,pFVz,pFWx,pFWy,pFWz
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pPjx,pPjy,pPjz
      tw0=MPI_Wtime()
      ! Midpoint face velocities of the iterate that get_dQdt transported with
      call this%Umid%copy(src=this%U); call this%Vmid%copy(src=this%V); call this%Wmid%copy(src=this%W)
      ! Face velocities of the new momenta, then the midpoint velocities of the energy pass
      call this%get_face_velocity(this%Q,this%rhoU,this%rhoV,this%rhoW,time)
      call this%U%lincomb(a=0.5_WP,src1=this%Uold,b=0.5_WP,src2=this%U)
      call this%V%lincomb(a=0.5_WP,src1=this%Vold,b=0.5_WP,src2=this%V)
      call this%W%lincomb(a=0.5_WP,src1=this%Wold,b=0.5_WP,src2=this%W)
      ! Sliver equilibration (torrens): a mixed cell of the finest level whose transported gas fraction is below sliver_frac,
      ! with both phases present, is relaxed by the relaxation model (with relax_isentropic, the gas keeps its entropy). The
      ! model acts on the transported internal contents Y^tr_k, and their change is added to Q(3:4). This
      ! changes VF and the transported internal contents, and so the phasic work below
      this%cnt%nsliver=0.0_WP
      if (associated(this%relax)) then
         lvl=this%amr%maxlvl
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pQ=>this%Q%mf(lvl)%dataptr(mfi); pVF=>this%VF%mf(lvl)%dataptr(mfi); pQold=>this%Qold%mf(lvl)%dataptr(mfi)
            pUold=>this%Uold%mf(lvl)%dataptr(mfi); pVold=>this%Vold%mf(lvl)%dataptr(mfi); pWold=>this%Wold%mf(lvl)%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               if (pVF(i,j,k,1).lt.VFlo.or.pVF(i,j,k,1).gt.VFhi.or.1.0_WP-pVF(i,j,k,1).ge.this%sliver_frac) cycle
               if (any(pQ(i,j,k,1:2).le.0.0_WP)) cycle
               mkn=max(pQold(i,j,k,1:2),0.0_WP)
               q(1:2)=pQ(i,j,k,1:2)
               q(3:4)=pQ(i,j,k,3:4)-mkn/max(sum(mkn),this%rho_floor)*cell_ke(pQold,pUold,pVold,pWold,i,j,k,this%rho_floor)
               Yb=q(3:4)
               call this%relax%apply(dt=dt,VF=pVF(i,j,k,1),Q=q,Pjump=0.0_WP,ierr=rerr)
               pQ(i,j,k,3:4)=pQ(i,j,k,3:4)+q(3:4)-Yb
               if (rerr.eq.RELAX_OK) this%cnt%nsliver=this%cnt%nsliver+1.0_WP
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         call this%VF%average_down(); call this%fill(lvl=this%amr%maxlvl,time=time)
      end if
      ! Mixture work flux on the valid faces: per direction, the momentum through the face at the two cells (its own axis)
      ! and at the two edges of each other axis, then the pressure work and the normal viscous work
      do lvl=0,this%amr%clvl()
         dxi=1.0_WP/this%amr%dx(lvl)
         dyi=1.0_WP/this%amr%dy(lvl)
         dzi=1.0_WP/this%amr%dz(lvl)
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pP=>this%P%mf(lvl)%dataptr(mfi)
            pU=>this%U%mf(lvl)%dataptr(mfi); pV=>this%V%mf(lvl)%dataptr(mfi); pW=>this%W%mf(lvl)%dataptr(mfi)
            pUm=>this%Umid%mf(lvl)%dataptr(mfi); pVm=>this%Vmid%mf(lvl)%dataptr(mfi); pWm=>this%Wmid%mf(lvl)%dataptr(mfi)
            pVisc=>this%visc%mf(lvl)%dataptr(mfi); pBeta=>this%beta%mf(lvl)%dataptr(mfi)
            pFx=>this%Fx(lvl)%dataptr(mfi); pFy=>this%Fy(lvl)%dataptr(mfi); pFz=>this%Fz(lvl)%dataptr(mfi)
            pFUx=>this%FUx(lvl)%dataptr(mfi); pFUy=>this%FUy(lvl)%dataptr(mfi); pFUz=>this%FUz(lvl)%dataptr(mfi)
            pFVx=>this%FVx(lvl)%dataptr(mfi); pFVy=>this%FVy(lvl)%dataptr(mfi); pFVz=>this%FVz(lvl)%dataptr(mfi)
            pFWx=>this%FWx(lvl)%dataptr(mfi); pFWy=>this%FWy(lvl)%dataptr(mfi); pFWz=>this%FWz(lvl)%dataptr(mfi)
            pEx=>this%Ex(lvl)%dataptr(mfi); pEy=>this%Ey(lvl)%dataptr(mfi); pEz=>this%Ez(lvl)%dataptr(mfi)
            if (lvl.eq.this%amr%maxlvl) pBand=>this%band%dataptr(mfi)
            ! X-faces: U at cells i-1 and i, V at xy-edges j and j+1, W at zx-edges k and k+1
            bx=mfi%nodaltilebox(1)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pEx(i,j,k,1)=0.5_WP*(site(0.5_WP*sum(pFx(i-1:i,j,k,1:2)),pU(i-1,j,k,1),pU(i  ,j,k,1),dev(1,1,i-1,j,k)) &
               &                   +site(0.5_WP*sum(pFx(i:i+1,j,k,1:2)),pU(i  ,j,k,1),pU(i+1,j,k,1),dev(1,1,i  ,j,k))) &
               &           +0.5_WP*(site(0.5_WP*sum(pFx(i,j-1:j,k,1:2)),pV(i-1,j  ,k,1),pV(i,j  ,k,1),dev(2,1,i,j  ,k)) &
               &                   +site(0.5_WP*sum(pFx(i,j:j+1,k,1:2)),pV(i-1,j+1,k,1),pV(i,j+1,k,1),dev(2,1,i,j+1,k))) &
               &           +0.5_WP*(site(0.5_WP*sum(pFx(i,j,k-1:k,1:2)),pW(i-1,j,k  ,1),pW(i,j,k  ,1),dev(3,1,i,j,k  )) &
               &                   +site(0.5_WP*sum(pFx(i,j,k:k+1,1:2)),pW(i-1,j,k+1,1),pW(i,j,k+1,1),dev(3,1,i,j,k+1))) &
               &           +0.5_WP*(stress(1,1,i-1,j,k)+stress(1,1,i,j,k)-pP(i-1,j,k,1)-pP(i,j,k,1))*pU(i,j,k,1)
            end do; end do; end do
            ! Y-faces: V at cells j-1 and j, U at xy-edges i and i+1, W at yz-edges k and k+1
            bx=mfi%nodaltilebox(2)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pEy(i,j,k,1)=0.5_WP*(site(0.5_WP*sum(pFy(i,j-1:j,k,1:2)),pV(i,j-1,k,1),pV(i,j  ,k,1),dev(2,2,i,j-1,k)) &
               &                   +site(0.5_WP*sum(pFy(i,j:j+1,k,1:2)),pV(i,j  ,k,1),pV(i,j+1,k,1),dev(2,2,i,j  ,k))) &
               &           +0.5_WP*(site(0.5_WP*sum(pFy(i-1:i,j,k,1:2)),pU(i  ,j-1,k,1),pU(i  ,j,k,1),dev(1,2,i  ,j,k)) &
               &                   +site(0.5_WP*sum(pFy(i:i+1,j,k,1:2)),pU(i+1,j-1,k,1),pU(i+1,j,k,1),dev(1,2,i+1,j,k))) &
               &           +0.5_WP*(site(0.5_WP*sum(pFy(i,j,k-1:k,1:2)),pW(i,j-1,k  ,1),pW(i,j,k  ,1),dev(3,2,i,j,k  )) &
               &                   +site(0.5_WP*sum(pFy(i,j,k:k+1,1:2)),pW(i,j-1,k+1,1),pW(i,j,k+1,1),dev(3,2,i,j,k+1))) &
               &           +0.5_WP*(stress(2,2,i,j-1,k)+stress(2,2,i,j,k)-pP(i,j-1,k,1)-pP(i,j,k,1))*pV(i,j,k,1)
            end do; end do; end do
            ! Z-faces: W at cells k-1 and k, U at zx-edges i and i+1, V at yz-edges j and j+1
            bx=mfi%nodaltilebox(3)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pEz(i,j,k,1)=0.5_WP*(site(0.5_WP*sum(pFz(i,j,k-1:k,1:2)),pW(i,j,k-1,1),pW(i,j,k  ,1),dev(3,3,i,j,k-1)) &
               &                   +site(0.5_WP*sum(pFz(i,j,k:k+1,1:2)),pW(i,j,k  ,1),pW(i,j,k+1,1),dev(3,3,i,j,k  ))) &
               &           +0.5_WP*(site(0.5_WP*sum(pFz(i-1:i,j,k,1:2)),pU(i  ,j,k-1,1),pU(i  ,j,k,1),dev(1,3,i  ,j,k)) &
               &                   +site(0.5_WP*sum(pFz(i:i+1,j,k,1:2)),pU(i+1,j,k-1,1),pU(i+1,j,k,1),dev(1,3,i+1,j,k))) &
               &           +0.5_WP*(site(0.5_WP*sum(pFz(i,j-1:j,k,1:2)),pV(i,j  ,k-1,1),pV(i,j  ,k,1),dev(2,3,i,j  ,k)) &
               &                   +site(0.5_WP*sum(pFz(i,j:j+1,k,1:2)),pV(i,j+1,k-1,1),pV(i,j+1,k,1),dev(2,3,i,j+1,k))) &
               &           +0.5_WP*(stress(3,3,i,j,k-1)+stress(3,3,i,j,k)-pP(i,j,k-1,1)-pP(i,j,k,1))*pW(i,j,k,1)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      ! Average down the work fluxes for C/F consistency
      do lvl=this%amr%clvl(),1,-1
         call amrmfab_average_down_face(fmf=this%Ex(lvl),cmf=this%Ex(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
         call amrmfab_average_down_face(fmf=this%Ey(lvl),cmf=this%Ey(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
         call amrmfab_average_down_face(fmf=this%Ez(lvl),cmf=this%Ez(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
      end do
      ! Phasic total energies on the valid cells: flux form off the band, direct treatment in band cells
      do lvl=0,this%amr%clvl()
         dxi=1.0_WP/this%amr%dx(lvl)
         dyi=1.0_WP/this%amr%dy(lvl)
         dzi=1.0_WP/this%amr%dz(lvl)
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pQ=>this%Q%mf(lvl)%dataptr(mfi); pQold=>this%Qold%mf(lvl)%dataptr(mfi)
            pVF=>this%VF%mf(lvl)%dataptr(mfi)
            pPL=>this%PL%mf(lvl)%dataptr(mfi); pPG=>this%PG%mf(lvl)%dataptr(mfi)
            pU=>this%U%mf(lvl)%dataptr(mfi); pUold=>this%Uold%mf(lvl)%dataptr(mfi)
            pV=>this%V%mf(lvl)%dataptr(mfi); pVold=>this%Vold%mf(lvl)%dataptr(mfi)
            pW=>this%W%mf(lvl)%dataptr(mfi); pWold=>this%Wold%mf(lvl)%dataptr(mfi)
            pUm=>this%Umid%mf(lvl)%dataptr(mfi); pVm=>this%Vmid%mf(lvl)%dataptr(mfi); pWm=>this%Wmid%mf(lvl)%dataptr(mfi)
            pVisc=>this%visc%mf(lvl)%dataptr(mfi); pBeta=>this%beta%mf(lvl)%dataptr(mfi)
            pEx=>this%Ex(lvl)%dataptr(mfi); pEy=>this%Ey(lvl)%dataptr(mfi); pEz=>this%Ez(lvl)%dataptr(mfi)
            if (lvl.eq.this%amr%maxlvl) then
               pBand=>this%band%dataptr(mfi)
               pPjx=>this%Pj(1)%dataptr(mfi); pPjy=>this%Pj(2)%dataptr(mfi); pPjz=>this%Pj(3)%dataptr(mfi)
            end if
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Rate of the mixture work flux (negative convention) and the new mass shares
               dW=dxi*(pEx(i+1,j,k,1)-pEx(i,j,k,1))+dyi*(pEy(i,j+1,k,1)-pEy(i,j,k,1))+dzi*(pEz(i,j,k+1,1)-pEz(i,j,k,1))
               mk=max(pQ(i,j,k,1:2),0.0_WP); mk=mk/max(sum(mk),this%rho_floor)
               ! Off the band: flux form
               inband=.false.; if (lvl.eq.this%amr%maxlvl) inband=pBand(i,j,k,1).gt.0.0_WP
               if (.not.inband) then
                  pQ(i,j,k,3:4)=pQ(i,j,k,3:4)+dt*mk*dW
                  cycle
               end if
               ! Viscous dissipation: normal stresses, and the four surrounding edges of each orientation
               gradU=[(pU(i+1,j,k,1)-pU(i,j,k,1))*dxi,(pV(i,j+1,k,1)-pV(i,j,k,1))*dyi,(pW(i,j,k+1,1)-pW(i,j,k,1))*dzi]
               div=sum(gradU)
               dss=stress(1,1,i,j,k)*gradU(1)+stress(2,2,i,j,k)*gradU(2)+stress(3,3,i,j,k)*gradU(3)
               do n2=0,1; do n1=0,1
                  dss=dss+0.25_WP*(stress(1,2,i+n1,j+n2,k)*((pU(i+n1,j+n2,k,1)-pU(i+n1,j+n2-1,k,1))*dyi+(pV(i+n1,j+n2,k,1)-pV(i+n1-1,j+n2,k,1))*dxi) &
                  &               +stress(2,3,i,j+n1,k+n2)*((pV(i,j+n1,k+n2,1)-pV(i,j+n1,k+n2-1,1))*dzi+(pW(i,j+n1,k+n2,1)-pW(i,j+n1-1,k+n2,1))*dyi) &
                  &               +stress(3,1,i+n1,j,k+n2)*((pW(i+n1,j,k+n2,1)-pW(i+n1-1,j,k+n2,1))*dxi+(pU(i+n1,j,k+n2,1)-pU(i+n1,j,k+n2-1,1))*dzi))
               end do; end do
               ! Phasic work against the phasic midpoint pressures, with the transported VF
               S=[pVF(i,j,k,1)*(dss-pPL(i,j,k,1)*div),(1.0_WP-pVF(i,j,k,1))*(dss-pPG(i,j,k,1)*div)]
               ! Capillary work: the surface tension force times the midpoint velocity, averaged over the cell's faces
               Cap=0.5_WP*(pU(i,j,k,1)*pPjx(i,j,k,1)+pU(i+1,j,k,1)*pPjx(i+1,j,k,1) &
               &          +pV(i,j,k,1)*pPjy(i,j,k,1)+pV(i,j+1,k,1)*pPjy(i,j+1,k,1) &
               &          +pW(i,j,k,1)*pPjz(i,j,k,1)+pW(i,j,k+1,1)*pPjz(i,j,k+1,1))
               ! Transported internal content, phasic work, and the new mass share of the cell residual
               Kn=cell_ke(pQold,pUold,pVold,pWold,i,j,k,this%rho_floor)
               mkn=max(pQold(i,j,k,1:2),0.0_WP)
               pQ(i,j,k,3:4)=pQ(i,j,k,3:4)-mkn/max(sum(mkn),this%rho_floor)*Kn+dt*S+mk*(Kn+dt*(dW+Cap-sum(S)))
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      ! Restrict and fill
      call this%Q%average_down(); call this%Q%fill(time=time)
      this%wt_div=this%wt_div+(MPI_Wtime()-tw0); this%wt_dQdt=this%wt_dQdt+(MPI_Wtime()-tw0)
   contains
      !> Kinetic energy flux through a dual face: the centered mass flux g times half the product of
      !> the velocities ua and ub of the two staggered CVs, plus the deviation of the convective flux from the centered
      !> one, and the shear stress at edges
      pure real(WP) function site(g,ua,ub,dev)
         implicit none
         real(WP), intent(in) :: g,ua,ub,dev
         site=0.5_WP*(g*ua*ub+dev*(ua+ub))
      end function site
      !> Viscous stress tau_ab of Umid, as get_dQdt formed it
      real(WP) function stress(a,b,i,j,k)
         implicit none
         integer, intent(in) :: a,b,i,j,k
         real(WP), dimension(3) :: g,t
         real(WP) :: mu,lam
         select case (10*a+b)
         case (11,22,33)
            g=[(pUm(i+1,j,k,1)-pUm(i,j,k,1))*dxi,(pVm(i,j+1,k,1)-pVm(i,j,k,1))*dyi,(pWm(i,j,k+1,1)-pWm(i,j,k,1))*dzi]
            mu=pVisc(i,j,k,1); lam=pBeta(i,j,k,1)-2.0_WP*mu/3.0_WP
            t=2.0_WP*mu*g+lam*sum(g)
            stress=t(a)
         case (12,21); stress=0.25_WP*sum(pVisc(i-1:i,j-1:j,k,1))*((pUm(i,j,k,1)-pUm(i,j-1,k,1))*dyi+(pVm(i,j,k,1)-pVm(i-1,j,k,1))*dxi)
         case (23,32); stress=0.25_WP*sum(pVisc(i,j-1:j,k-1:k,1))*((pVm(i,j,k,1)-pVm(i,j,k-1,1))*dzi+(pWm(i,j,k,1)-pWm(i,j-1,k,1))*dyi)
         case default; stress=0.25_WP*sum(pVisc(i-1:i,j,k-1:k,1))*((pWm(i,j,k,1)-pWm(i-1,j,k,1))*dxi+(pUm(i,j,k,1)-pUm(i,j,k-1,1))*dzi)
         end select
      end function stress
      !> Deviation of the stored a-momentum flux through the b-faces from the centered flux of Umid, at cell (i,j,k) 
      real(WP) function dev(a,b,i,j,k)
         implicit none
         integer, intent(in) :: a,b,i,j,k
         integer, dimension(3) :: lo,hi
         logical :: onband
         onband=.false.
         if (lvl.eq.this%amr%maxlvl) then
            if (a.eq.b) then
               lo=[i,j,k]-eu(:,a); hi=[i,j,k]+eu(:,a)
            else
               lo=[i,j,k]-eu(:,a)-eu(:,b); hi=[i,j,k]
            end if
            onband=maxval(pBand(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1)).gt.0.0_WP
         end if
         if (.not.onband) then
            dev=0.0_WP; if (a.ne.b) dev=stress(a,b,i,j,k)
            return
         end if
         select case (10*a+b)
         case (11); dev=pFUx(i,j,k,1)-stress(1,1,i,j,k)+pP(i,j,k,1)-0.25_WP*sum(pFx(i:i+1,j,k,1:2))*sum(pUm(i:i+1,j,k,1))
         case (22); dev=pFVy(i,j,k,1)-stress(2,2,i,j,k)+pP(i,j,k,1)-0.25_WP*sum(pFy(i,j:j+1,k,1:2))*sum(pVm(i,j:j+1,k,1))
         case (33); dev=pFWz(i,j,k,1)-stress(3,3,i,j,k)+pP(i,j,k,1)-0.25_WP*sum(pFz(i,j,k:k+1,1:2))*sum(pWm(i,j,k:k+1,1))
         case (12); dev=pFUy(i,j,k,1)-0.25_WP*sum(pFy(i-1:i,j,k,1:2))*sum(pUm(i,j-1:j,k,1))
         case (21); dev=pFVx(i,j,k,1)-0.25_WP*sum(pFx(i,j-1:j,k,1:2))*sum(pVm(i-1:i,j,k,1))
         case (23); dev=pFVz(i,j,k,1)-0.25_WP*sum(pFz(i,j-1:j,k,1:2))*sum(pVm(i,j,k-1:k,1))
         case (32); dev=pFWy(i,j,k,1)-0.25_WP*sum(pFy(i,j,k-1:k,1:2))*sum(pWm(i,j-1:j,k,1))
         case (31); dev=pFWx(i,j,k,1)-0.25_WP*sum(pFx(i,j,k-1:k,1:2))*sum(pWm(i-1:i,j,k,1))
         case default; dev=pFUz(i,j,k,1)-0.25_WP*sum(pFz(i-1:i,j,k,1:2))*sum(pUm(i,j,k-1:k,1))
         end select
      end function dev
   end subroutine get_energy

   !> Apply the relaxation model in mixed cells of the finest level, then re-match the PLIC
   !> plane to the relaxed VF. Q(1:2) and Q(3)+Q(4) are conserved 
   subroutine apply_relax(this,dt,time)
      use amrvof_geometry, only: get_plane_dist,cut_hex_vol
      use mpi_f08,         only: MPI_Wtime
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(in) :: time
      real(WP) :: tw0
      integer :: lvl,i,j,k,rerr
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pQ,pCL,pCG,pPLIC,pU,pV,pW,pRHOU,pRHOV,pRHOW,pCurv
      logical :: oldmix,newmix,gated
      real(WP) :: dx,dy,dz,cell_vol,vol_liq,vol_gas
      real(WP), dimension(3) :: lo,hi,bary_liq,bary_gas
      real(WP), dimension(3,8) :: hex
      real(WP), dimension(4) :: plane,q
      real(WP), dimension(2) :: mk,Yb
      ! Reset failure and gate counters
      this%cnt%nrelax_fail=0.0_WP
      this%cnt%ngate=0.0_WP; this%cnt%gate_mass=0.0_WP; this%cnt%gate_energy=0.0_WP
      this%cnt%gate_rhoU=0.0_WP; this%cnt%gate_rhoV=0.0_WP; this%cnt%gate_rhoW=0.0_WP
      ! If no relaxation model was provided, return
      if (.not.associated(this%relax)) return
      ! Return if clvl<maxlvl
      if (this%amr%clvl().lt.this%amr%maxlvl) return
      ! Face velocities of the new iterate, for its kinetic energy
      call this%get_face_velocity(this%Q,this%rhoU,this%rhoV,this%rhoW,time)
      tw0=MPI_Wtime()
      ! Apply relaxation on finest level only (mixture cells are always at finest)
      lvl=this%amr%maxlvl
      dx=this%amr%dx(lvl); dy=this%amr%dy(lvl); dz=this%amr%dz(lvl)
      cell_vol=this%amr%cell_vol(lvl)
      call this%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         ! Get pointers to data
         pVF  =>this%VF%mf(lvl)%dataptr(mfi)
         pQ   =>this%Q%mf(lvl)%dataptr(mfi)
         pCL  =>this%CL%dataptr(mfi)
         pCG  =>this%CG%dataptr(mfi)
         pPLIC=>this%PLIC%dataptr(mfi)
         pCurv=>this%curv%dataptr(mfi)
         pU=>this%U%mf(lvl)%dataptr(mfi); pRHOU=>this%rhoU%mf(lvl)%dataptr(mfi)
         pV=>this%V%mf(lvl)%dataptr(mfi); pRHOV=>this%rhoV%mf(lvl)%dataptr(mfi)
         pW=>this%W%mf(lvl)%dataptr(mfi); pRHOW=>this%rhoW%mf(lvl)%dataptr(mfi)
         ! Loop over valid cells
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Check if mixture cell prior to relaxation
            oldmix=(pVF(i,j,k,1).ge.VFlo.and.pVF(i,j,k,1).le.VFhi)
            gated=.false.; rerr=RELAX_OK
            if (oldmix) then
               ! Internal contents of the mixed cell, Y_k=E_k-m_k/rho*K
               mk=max(pQ(i,j,k,1:2),0.0_WP)
               q(1:2)=pQ(i,j,k,1:2)
               q(3:4)=pQ(i,j,k,3:4)-mk/max(sum(mk),this%rho_floor)*cell_ke(pQ,pU,pV,pW,i,j,k,this%rho_floor)
               ! Positivity gate before the relaxation
               call gate_cell(i,j,k,q,gated)
               ! User-provided relaxation model to the capillary pressure jump (modifies VF and q), and the internal energy
               ! it exchanges added to E_k
               if (.not.gated) then
                  Yb=q(3:4)
                  call this%relax%apply(dt=dt,VF=pVF(i,j,k,1),Q=q,Pjump=this%sigma*pCurv(i,j,k,1),ierr=rerr)
                  pQ(i,j,k,3:4)=pQ(i,j,k,3:4)+q(3:4)-Yb
               end if
               if (rerr.ne.RELAX_OK) this%cnt%nrelax_fail=this%cnt%nrelax_fail+1.0_WP
            end if
            ! Check if mixture cell after relaxation
            newmix=(pVF(i,j,k,1).ge.VFlo.and.pVF(i,j,k,1).le.VFhi)
            ! If not mixture cells, clean up and cycle
            if (.not.newmix) then
               if (pVF(i,j,k,1).lt.VFlo) then
                  ! Pure gas
                  pVF(i,j,k,1)=0.0_WP
                  pCL(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
                  pCG(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
                  pQ(i,j,k,1)=0.0_WP
                  pQ(i,j,k,4)=pQ(i,j,k,4)+pQ(i,j,k,3); pQ(i,j,k,3)=0.0_WP
                  pPLIC(i,j,k,:)=[0.0_WP,0.0_WP,0.0_WP,-1.0e10_WP]
               else if (pVF(i,j,k,1).gt.VFhi) then
                  ! Pure liquid
                  pVF(i,j,k,1)=1.0_WP
                  pCL(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
                  pCG(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
                  pQ(i,j,k,2)=0.0_WP
                  pQ(i,j,k,3)=pQ(i,j,k,3)+pQ(i,j,k,4); pQ(i,j,k,4)=0.0_WP
                  pPLIC(i,j,k,:)=[0.0_WP,0.0_WP,0.0_WP,+1.0e10_WP]
               end if
               cycle
            end if
            ! If mixture cell, post-process PLIC and barycenters
            if (.not.oldmix) pPLIC(i,j,k,1:3)=[1.0_WP,0.0_WP,0.0_WP]
            ! Reposition plane: keep normal, adjust distance for new VF
            lo=[this%amr%xlo+real(i  ,WP)*dx,this%amr%ylo+real(j  ,WP)*dy,this%amr%zlo+real(k  ,WP)*dz]
            hi=[this%amr%xlo+real(i+1,WP)*dx,this%amr%ylo+real(j+1,WP)*dy,this%amr%zlo+real(k+1,WP)*dz]
            pPLIC(i,j,k,4)=get_plane_dist(pPLIC(i,j,k,1:3),lo,hi,pVF(i,j,k,1))
            ! Recompute barycenters from adjusted PLIC
            hex(:,1)=[lo(1),lo(2),lo(3)]
            hex(:,2)=[hi(1),lo(2),lo(3)]
            hex(:,3)=[hi(1),hi(2),lo(3)]
            hex(:,4)=[lo(1),hi(2),lo(3)]
            hex(:,5)=[lo(1),lo(2),hi(3)]
            hex(:,6)=[hi(1),lo(2),hi(3)]
            hex(:,7)=[hi(1),hi(2),hi(3)]
            hex(:,8)=[lo(1),hi(2),hi(3)]
            plane=pPLIC(i,j,k,:)
            call cut_hex_vol(hex,plane,vol_liq,vol_gas,bary_liq,bary_gas)
            if (this%reset_vf) pVF(i,j,k,1)=vol_liq/cell_vol
            pCL(i,j,k,1:3)=bary_liq
            pCG(i,j,k,1:3)=bary_gas
         end do; end do; end do
      end do
      call this%amr%mfiter_destroy(mfi)
      ! Sync and apply BC
      call this%VF%average_down(); call this%fill(lvl=this%amr%maxlvl,time=time)
      call this%Q%average_down(); call this%Q%fill(time=time)
      this%wt_relax=this%wt_relax+(MPI_Wtime()-tw0)
   contains
      !> Positivity gate of mixed cell (i,j,k) (torrens), with q=[m_1,m_2,Y_1,Y_2] its masses and internal contents: a phase
      !> without positive mass, internal content, or sound speed (p>-pinf) is removed, momentum and energy is rebuilt
      subroutine gate_cell(i,j,k,q,gated)
         implicit none
         integer, intent(in) :: i,j,k
         real(WP), dimension(4), intent(in) :: q
         logical, intent(out) :: gated
         logical :: badL,badG
         integer :: n,o,s
         real(WP) :: E0,old
         badL=q(1).le.0.0_WP.or.q(3).le.0.0_WP
         if (.not.badL) badL=this%liq%get_c_from_rho_e(rho=q(1)/pVF(i,j,k,1),e=q(3)/q(1),y=[1.0_WP]).le.0.0_WP
         badG=q(2).le.0.0_WP.or.q(4).le.0.0_WP
         if (.not.badG) badG=this%gas%get_c_from_rho_e(rho=q(2)/(1.0_WP-pVF(i,j,k,1)),e=q(4)/q(2),y=[1.0_WP]).le.0.0_WP
         gated=badL.or.badG
         if (.not.gated) return
         ! Remove phase n
         n=2; if (badL.and.(.not.badG.or.q(1).le.q(2))) n=1
         o=3-n
         E0=pQ(i,j,k,3)+pQ(i,j,k,4)
         this%cnt%ngate=this%cnt%ngate+1.0_WP
         this%cnt%gate_mass=this%cnt%gate_mass+pQ(i,j,k,n)*cell_vol
         pQ(i,j,k,n)=0.0_WP; pQ(i,j,k,2+n)=0.0_WP
         pVF(i,j,k,1)=merge(0.0_WP,1.0_WP,n.eq.1)
         ! Face momenta from the new face masses at the current face velocities
         do s=0,1
            old=pRHOU(i+s,j,k,1); pRHOU(i+s,j,k,1)=max(0.5_WP*sum(pQ(i+s-1:i+s,j,k,1:2)),this%rho_floor)*pU(i+s,j,k,1)
            this%cnt%gate_rhoU=this%cnt%gate_rhoU+(pRHOU(i+s,j,k,1)-old)*cell_vol
            old=pRHOV(i,j+s,k,1); pRHOV(i,j+s,k,1)=max(0.5_WP*sum(pQ(i,j+s-1:j+s,k,1:2)),this%rho_floor)*pV(i,j+s,k,1)
            this%cnt%gate_rhoV=this%cnt%gate_rhoV+(pRHOV(i,j+s,k,1)-old)*cell_vol
            old=pRHOW(i,j,k+s,1); pRHOW(i,j,k+s,1)=max(0.5_WP*sum(pQ(i,j,k+s-1:k+s,1:2)),this%rho_floor)*pW(i,j,k+s,1)
            this%cnt%gate_rhoW=this%cnt%gate_rhoW+(pRHOW(i,j,k+s,1)-old)*cell_vol
         end do
         ! The remaining phase owns the kinetic energy of the cell
         pQ(i,j,k,2+o)=q(2+o)+cell_ke(pQ,pU,pV,pW,i,j,k,this%rho_floor)
         this%cnt%gate_energy=this%cnt%gate_energy+(E0-pQ(i,j,k,2+o))*cell_vol
      end subroutine gate_cell
   end subroutine apply_relax

   !> Add artificial bulk viscosity to this%beta (and optionally this%visc). The kinematic
   !> viscosity is converted with the harmonic phasic density
   subroutine add_viscartif(this,dt,Cartif,Cvisc,max_cfl)
      use amrsgs,  only: get_viscartif
      use mpi_f08, only: MPI_Wtime
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP) :: tw0
      real(WP), intent(in) :: dt
      real(WP), intent(in), optional :: Cartif
      real(WP), intent(in), optional :: Cvisc
      real(WP), intent(in), optional :: max_cfl
      real(WP) :: myCvisc,my_beta
      type(amrdata) :: beta_t
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pBeta_t,pBeta,pVisc,pVF,pRHOL,pRHOG
      tw0=MPI_Wtime()
      ! Set shear viscosity constant
      if (present(Cvisc)) then; myCvisc=Cvisc; else; myCvisc=0.0_WP; end if
      ! Create temp amrdata
      call beta_t%initialize(this%amr,name='beta_t',ncomp=1,ng=this%beta%ng); call beta_t%reset()
      ! Compute kinematic artificial bulk viscosity into temp from staggered velocity
      call get_viscartif(dt=dt,visc=beta_t,U=this%U,V=this%V,W=this%W,C=this%C,Cartif=Cartif,max_cfl=max_cfl)
      ! Convert kinematic to dynamic via harmonic phasic density, add to beta and visc
      do lvl=0,this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            pBeta_t=>beta_t%mf(lvl)%dataptr(mfi)
            pBeta=>this%beta%mf(lvl)%dataptr(mfi)
            pVisc=>this%visc%mf(lvl)%dataptr(mfi)
            pRHOL=>this%RHOL%mf(lvl)%dataptr(mfi)
            pRHOG=>this%RHOG%mf(lvl)%dataptr(mfi)
            pVF  =>this%VF%mf(lvl)%dataptr(mfi)
            bx=mfi%growntilebox(min(this%beta%ng,this%visc%ng))
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               my_beta=pBeta_t(i,j,k,1)/(pVF(i,j,k,1)/max(pRHOL(i,j,k,1),this%rho_floor)+(1.0_WP-pVF(i,j,k,1))/max(pRHOG(i,j,k,1),this%rho_floor))
               pBeta(i,j,k,1)=pBeta(i,j,k,1)+my_beta
               pVisc(i,j,k,1)=pVisc(i,j,k,1)+my_beta*myCvisc
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      ! Destroy temp amrdata
      call beta_t%finalize()
      this%wt_visc=this%wt_visc+(MPI_Wtime()-tw0)
   end subroutine add_viscartif

   !> Add Vreman SGS eddy viscosity to this%visc from the staggered velocity
   subroutine add_vreman(this,dt,Cs)
      use amrsgs,  only: get_vreman
      use mpi_f08, only: MPI_Wtime
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP) :: tw0
      real(WP), intent(in) :: dt
      real(WP), intent(in), optional :: Cs
      type(amrdata) :: visc_t
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVisc_t,pVisc,pVF,pRHOL,pRHOG
      tw0=MPI_Wtime()
      ! Create temp amrdata
      call visc_t%initialize(amr=this%amr,name='visc_t',ncomp=1,ng=this%visc%ng); call visc_t%reset()
      ! Compute kinematic eddy viscosity into temp from staggered velocity
      call get_vreman(dt=dt,visc=visc_t,U=this%U,V=this%V,W=this%W,Cs=Cs)
      ! Convert kinematic to dynamic via harmonic phasic density, add to visc
      do lvl=0,this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            pVisc_t=>visc_t%mf(lvl)%dataptr(mfi)
            pVisc=>this%visc%mf(lvl)%dataptr(mfi)
            pRHOL=>this%RHOL%mf(lvl)%dataptr(mfi)
            pRHOG=>this%RHOG%mf(lvl)%dataptr(mfi)
            pVF  =>this%VF%mf(lvl)%dataptr(mfi)
            bx=mfi%growntilebox(this%visc%ng)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                     pVisc(i,j,k,1)=pVisc(i,j,k,1)+pVisc_t(i,j,k,1)/(pVF(i,j,k,1)/max(pRHOL(i,j,k,1),this%rho_floor)+(1.0_WP-pVF(i,j,k,1))/max(pRHOG(i,j,k,1),this%rho_floor))
                  end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      ! Destroy temp amrdata
      call visc_t%finalize()
      this%wt_visc=this%wt_visc+(MPI_Wtime()-tw0)
   end subroutine add_vreman

   !> Compute CFL numbers (convective, acoustic, viscous/heat diffusion, and capillary)
   subroutine get_cfl(this,dt,cfl)
      use mathtools, only: Pi
      use parallel, only: MPI_REAL_WP
      use mpi_f08, only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_MAX
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(out) :: cfl
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      integer :: lvl,i,j,k,ierr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pP,pU,pV,pW,pVisc,pBeta,pC
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pRHOL,pRHOG,pIL,pIG,pDiffL,pDiffG
      real(WP) :: dxi,dyi,dzi,rho,conv,pgrad,viscmax,alpha_heat
      ! Get convective CFL from parent
      call this%amrmpflow%get_cflc(dt=dt)
      ! Reset CFLs
      this%CFLp=0.0_WP
      this%CFLa_x=0.0_WP; this%CFLa_y=0.0_WP; this%CFLa_z=0.0_WP
      this%CFLv_x=0.0_WP; this%CFLv_y=0.0_WP; this%CFLv_z=0.0_WP
      this%CFLst=0.0_WP
      ! Compute CFL at each level (finest level determines dt)
      do lvl=0,this%amr%clvl()
         ! Get mesh spacing
         dxi=1.0_WP/this%amr%dx(lvl)
         dyi=1.0_WP/this%amr%dy(lvl)
         dzi=1.0_WP/this%amr%dz(lvl)
         ! Loop over tiles
         call this%amr%mfiter_build(lvl,mfi)
         do while(mfi%next())
            ! Get data pointers
            pQ=>this%Q%mf(lvl)%dataptr(mfi)
            pVisc=>this%visc%mf(lvl)%dataptr(mfi)
            pBeta=>this%beta%mf(lvl)%dataptr(mfi)
            pP=>this%P%mf(lvl)%dataptr(mfi)
            pU=>this%U%mf(lvl)%dataptr(mfi)
            pV=>this%V%mf(lvl)%dataptr(mfi)
            pW=>this%W%mf(lvl)%dataptr(mfi)
            pC=>this%C%mf(lvl)%dataptr(mfi)
            pVF=>this%VF%mf(lvl)%dataptr(mfi)
            pRHOL=>this%RHOL%mf(lvl)%dataptr(mfi); pRHOG=>this%RHOG%mf(lvl)%dataptr(mfi)
            pIL=>this%IL%mf(lvl)%dataptr(mfi); pIG=>this%IG%mf(lvl)%dataptr(mfi)
            pDiffL=>this%diffL%mf(lvl)%dataptr(mfi); pDiffG=>this%diffG%mf(lvl)%dataptr(mfi)
            ! Loop over interior tiles
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               rho=max(pQ(i,j,k,1)+pQ(i,j,k,2),this%rho_floor)
               ! Heat-diffusion CFL
               alpha_heat=0.0_WP
               if (pVF(i,j,k,1).gt.VFhi.and.pRHOL(i,j,k,1).gt.0.0_WP) then
                  alpha_heat=pDiffL(i,j,k,1)/max(pRHOL(i,j,k,1)*this%liq%get_cv_from_rho_e(rho=pRHOL(i,j,k,1),e=pIL(i,j,k,1),y=[1.0_WP]),tiny(1.0_WP))
               else if (pVF(i,j,k,1).lt.VFlo.and.pRHOG(i,j,k,1).gt.0.0_WP) then
                  alpha_heat=pDiffG(i,j,k,1)/max(pRHOG(i,j,k,1)*this%gas%get_cv_from_rho_e(rho=pRHOG(i,j,k,1),e=pIG(i,j,k,1),y=[1.0_WP]),tiny(1.0_WP))
               end if
               ! Viscous CFL with the longitudinal coefficient (beta+4*visc/3)/rho, or the thermal diffusivity if larger
               viscmax=max((pBeta(i,j,k,1)+4.0_WP*pVisc(i,j,k,1)/3.0_WP)/rho,alpha_heat)
               if (this%amr%nx.gt.1) this%CFLv_x=max(this%CFLv_x,4.0_WP*viscmax*dt*dxi**2)
               if (this%amr%ny.gt.1) this%CFLv_y=max(this%CFLv_y,4.0_WP*viscmax*dt*dyi**2)
               if (this%amr%nz.gt.1) this%CFLv_z=max(this%CFLv_z,4.0_WP*viscmax*dt*dzi**2)
               ! Convective+pressure CFL
               conv=abs(pU(i,j,k,1))*dxi+abs(pV(i,j,k,1))*dyi+abs(pW(i,j,k,1))*dzi
               pgrad=(abs(pP(i,j,k,1)-pP(i-1,j,k,1))*dxi**2+abs(pP(i,j,k,1)-pP(i,j-1,k,1))*dyi**2+abs(pP(i,j,k,1)-pP(i,j,k-1,1))*dzi**2)/rho
               this%CFLp=max(this%CFLp,0.5_WP*dt*(conv+sqrt(conv**2+4.0_WP*pgrad)))
               ! Acoustic CFL
               if (this%amr%nx.gt.1) this%CFLa_x=max(this%CFLa_x,(abs(0.5_WP*(pU(i,j,k,1)+pU(i+1,j,k,1)))+pC(i,j,k,1))*dt*dxi)
               if (this%amr%ny.gt.1) this%CFLa_y=max(this%CFLa_y,(abs(0.5_WP*(pV(i,j,k,1)+pV(i,j+1,k,1)))+pC(i,j,k,1))*dt*dyi)
               if (this%amr%nz.gt.1) this%CFLa_z=max(this%CFLa_z,(abs(0.5_WP*(pW(i,j,k,1)+pW(i,j,k+1,1)))+pC(i,j,k,1))*dt*dzi)
               ! Capillary CFL in finest-level mixed cells
               if (lvl.eq.this%amr%maxlvl.and.pVF(i,j,k,1).ge.VFlo.and.pVF(i,j,k,1).le.VFhi) this%CFLst=max(this%CFLst,dt*sqrt(4.0_WP*Pi*this%sigma/(rho*this%amr%min_meshsize(lvl)**3)))
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      ! Reduce across ranks
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLp  ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLa_x,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLa_y,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLa_z,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLv_x,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLv_y,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLv_z,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLst ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      ! Return max CFL
      cfl=max(this%CFLc_x,this%CFLc_y,this%CFLc_z,this%CFLv_x,this%CFLv_y,this%CFLv_z)
      cfl=max(cfl,this%CFLa_x,this%CFLa_y,this%CFLa_z,this%CFLst)
   end subroutine get_cfl

   ! ============================================================================
   ! SOLVER INFO
   ! ============================================================================

   !> Get solver information: phasic extrema, mixture pressure and sound speed, momentum, kinetic energy, and the
   !> diagnostic counts of the last calls summed across ranks
   subroutine get_info(this)
      use amrex_amr_module, only: amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
      use amrex_interface,  only: amrmask_make_fine
      use parallel,         only: MPI_REAL_WP
      use mpi_f08,          only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM,MPI_MAX,MPI_MIN
      implicit none
      class(amrmpcomp), intent(inout) :: this
      integer :: lvl,i,j,k,ierr
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_imultifab) :: mask
      integer, dimension(:,:,:,:), contiguous, pointer :: pMask
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pRHOL,pRHOG,pIL,pIG,pPL,pPG,pTL,pTG
      real(WP), dimension(:,:,:,:), contiguous, pointer :: prhoU,prhoV,prhoW,pU,pV,pW
      real(WP), dimension(9) :: buf

      ! Use parent's method first (velocity, Q, and VF statistics)
      call this%amrmpflow%get_info()

      ! Mixture pressure and speed of sound extrema
      this%Pmin=+huge(1.0_WP); this%Pmax=-huge(1.0_WP)
      this%Cmin=+huge(1.0_WP); this%Cmax=-huge(1.0_WP)
      do lvl=0,this%amr%clvl()
         this%Pmin=min(this%Pmin,this%P%get_min(lvl=lvl)); this%Pmax=max(this%Pmax,this%P%get_max(lvl=lvl))
         this%Cmin=min(this%Cmin,this%C%get_min(lvl=lvl)); this%Cmax=max(this%Cmax,this%C%get_max(lvl=lvl))
      end do

      ! Momentum integrals (rho * U * dV, summed over faces at level 0)
      this%rhoUint=this%rhoU%get_sum(lvl=0)*this%amr%cell_vol(0)
      this%rhoVint=this%rhoV%get_sum(lvl=0)*this%amr%cell_vol(0)
      this%rhoWint=this%rhoW%get_sum(lvl=0)*this%amr%cell_vol(0)

      ! Phasic extrema where the phase exists, and kinetic energy, over uncovered cells
      this%RHOLmin=huge(1.0_WP); this%RHOLmax=-huge(1.0_WP); this%RHOGmin=huge(1.0_WP); this%RHOGmax=-huge(1.0_WP)
      this%ILmin=huge(1.0_WP); this%ILmax=-huge(1.0_WP); this%IGmin=huge(1.0_WP); this%IGmax=-huge(1.0_WP)
      this%PLmin=huge(1.0_WP); this%PLmax=-huge(1.0_WP); this%PGmin=huge(1.0_WP); this%PGmax=-huge(1.0_WP)
      this%TLmin=huge(1.0_WP); this%TLmax=-huge(1.0_WP); this%TGmin=huge(1.0_WP); this%TGmax=-huge(1.0_WP)
      this%dPmax=0.0_WP
      this%rhoKint=0.0_WP
      do lvl=0,this%amr%clvl()
         ! Build fine mask for this level (if not finest)
         if (lvl.lt.this%amr%clvl()) then
            call amrex_imultifab_build(mask,this%amr%ba(lvl),this%amr%dm(lvl),1,0)
            call amrmask_make_fine(mask,this%amr%ba(lvl+1),[this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],0,1)
         end if
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get data pointers
            pVF=>this%VF%mf(lvl)%dataptr(mfi)
            pRHOL=>this%RHOL%mf(lvl)%dataptr(mfi); pRHOG=>this%RHOG%mf(lvl)%dataptr(mfi)
            pIL=>this%IL%mf(lvl)%dataptr(mfi); pIG=>this%IG%mf(lvl)%dataptr(mfi)
            pPL=>this%PL%mf(lvl)%dataptr(mfi); pPG=>this%PG%mf(lvl)%dataptr(mfi)
            pTL=>this%TL%mf(lvl)%dataptr(mfi); pTG=>this%TG%mf(lvl)%dataptr(mfi)
            prhoU=>this%rhoU%mf(lvl)%dataptr(mfi); pU=>this%U%mf(lvl)%dataptr(mfi)
            prhoV=>this%rhoV%mf(lvl)%dataptr(mfi); pV=>this%V%mf(lvl)%dataptr(mfi)
            prhoW=>this%rhoW%mf(lvl)%dataptr(mfi); pW=>this%W%mf(lvl)%dataptr(mfi)
            if (lvl.lt.this%amr%clvl()) pMask=>mask%dataptr(mfi)
            ! Loop over interior tiles
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Skip cells covered by finer level
               if (lvl.lt.this%amr%clvl()) then; if (pMask(i,j,k,1).eq.0) cycle; end if
               ! Liquid variables
               if (pVF(i,j,k,1).ge.VFlo) then
                  this%RHOLmin=min(this%RHOLmin,pRHOL(i,j,k,1)); this%RHOLmax=max(this%RHOLmax,pRHOL(i,j,k,1))
                  this%ILmin  =min(this%ILmin  ,pIL  (i,j,k,1)); this%ILmax  =max(this%ILmax  ,pIL  (i,j,k,1))
                  this%PLmin  =min(this%PLmin  ,pPL  (i,j,k,1)); this%PLmax  =max(this%PLmax  ,pPL  (i,j,k,1))
                  this%TLmin  =min(this%TLmin  ,pTL  (i,j,k,1)); this%TLmax  =max(this%TLmax  ,pTL  (i,j,k,1))
               end if
               ! Gas variables
               if (pVF(i,j,k,1).le.VFhi) then
                  this%RHOGmin=min(this%RHOGmin,pRHOG(i,j,k,1)); this%RHOGmax=max(this%RHOGmax,pRHOG(i,j,k,1))
                  this%IGmin  =min(this%IGmin  ,pIG  (i,j,k,1)); this%IGmax  =max(this%IGmax  ,pIG  (i,j,k,1))
                  this%PGmin  =min(this%PGmin  ,pPG  (i,j,k,1)); this%PGmax  =max(this%PGmax  ,pPG  (i,j,k,1))
                  this%TGmin  =min(this%TGmin  ,pTG  (i,j,k,1)); this%TGmax  =max(this%TGmax  ,pTG  (i,j,k,1))
               end if
               ! Pressure gap in mixed cells
               if (pVF(i,j,k,1).ge.VFlo.and.pVF(i,j,k,1).le.VFhi) this%dPmax=max(this%dPmax,abs(pPL(i,j,k,1)-pPG(i,j,k,1)))
               ! Kinetic energy from face momenta and velocities
               this%rhoKint=this%rhoKint+0.25_WP*(sum(prhoU(i:i+1,j,k,1)*pU(i:i+1,j,k,1))+sum(prhoV(i,j:j+1,k,1)*pV(i,j:j+1,k,1))+sum(prhoW(i,j,k:k+1,1)*pW(i,j,k:k+1,1)))*this%amr%cell_vol(lvl)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         if (lvl.lt.this%amr%clvl()) call amrex_imultifab_destroy(mask)
      end do
      ! Reduce across MPI ranks
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%RHOLmin,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr); call MPI_ALLREDUCE(MPI_IN_PLACE,this%RHOLmax,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%ILmin  ,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr); call MPI_ALLREDUCE(MPI_IN_PLACE,this%ILmax  ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%PLmin  ,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr); call MPI_ALLREDUCE(MPI_IN_PLACE,this%PLmax  ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%TLmin  ,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr); call MPI_ALLREDUCE(MPI_IN_PLACE,this%TLmax  ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%RHOGmin,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr); call MPI_ALLREDUCE(MPI_IN_PLACE,this%RHOGmax,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%IGmin  ,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr); call MPI_ALLREDUCE(MPI_IN_PLACE,this%IGmax  ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%PGmin  ,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr); call MPI_ALLREDUCE(MPI_IN_PLACE,this%PGmax  ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%TGmin  ,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr); call MPI_ALLREDUCE(MPI_IN_PLACE,this%TGmax  ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%dPmax  ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%rhoKint,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)

      ! Diagnostic counts of the last calls, summed across ranks from this rank's values
      buf=[this%cnt%nvf_guard,this%cnt%nsliver,this%cnt%nrelax_fail,this%cnt%ngate,this%cnt%gate_mass,this%cnt%gate_energy, &
      &    this%cnt%gate_rhoU,this%cnt%gate_rhoV,this%cnt%gate_rhoW]
      call MPI_ALLREDUCE(MPI_IN_PLACE,buf,9,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
      this%nvf_guard=buf(1); this%nsliver=buf(2); this%nrelax_fail=buf(3); this%ngate=buf(4); this%gate_mass=buf(5)
      this%gate_energy=buf(6); this%gate_rhoU=buf(7); this%gate_rhoV=buf(8); this%gate_rhoW=buf(9)

      ! Reduce per-rank timing to min/max across ranks
      call MPI_ALLREDUCE(this%wt_mid   ,this%wtmax_mid   ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr); call MPI_ALLREDUCE(this%wt_mid   ,this%wtmin_mid   ,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_prim  ,this%wtmax_prim  ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr); call MPI_ALLREDUCE(this%wt_prim  ,this%wtmin_prim  ,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_fvel  ,this%wtmax_fvel  ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr); call MPI_ALLREDUCE(this%wt_fvel  ,this%wtmin_fvel  ,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_dQdt  ,this%wtmax_dQdt  ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr); call MPI_ALLREDUCE(this%wt_dQdt  ,this%wtmin_dQdt  ,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_alloc ,this%wtmax_alloc ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr); call MPI_ALLREDUCE(this%wt_alloc ,this%wtmin_alloc ,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_slflux,this%wtmax_slflux,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr); call MPI_ALLREDUCE(this%wt_slflux,this%wtmin_slflux,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_vmom  ,this%wtmax_vmom  ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr); call MPI_ALLREDUCE(this%wt_vmom  ,this%wtmin_vmom  ,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_fv    ,this%wtmax_fv    ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr); call MPI_ALLREDUCE(this%wt_fv    ,this%wtmin_fv    ,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_div   ,this%wtmax_div   ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr); call MPI_ALLREDUCE(this%wt_div   ,this%wtmin_div   ,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_relax ,this%wtmax_relax ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr); call MPI_ALLREDUCE(this%wt_relax ,this%wtmin_relax ,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(this%wt_visc  ,this%wtmax_visc  ,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr); call MPI_ALLREDUCE(this%wt_visc  ,this%wtmin_visc  ,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
      ! Reset per-rank timing accumulators for the next interval
      this%wt_mid=0.0_WP; this%wt_prim=0.0_WP; this%wt_fvel=0.0_WP; this%wt_dQdt=0.0_WP; this%wt_alloc=0.0_WP; this%wt_slflux=0.0_WP
      this%wt_vmom=0.0_WP; this%wt_fv=0.0_WP; this%wt_div=0.0_WP; this%wt_relax=0.0_WP; this%wt_visc=0.0_WP

   end subroutine get_info

   !> Print solver info to screen
   subroutine amrmpcomp_print(this)
      use messager, only: log
      implicit none
      class(amrmpcomp), intent(in) :: this
      call log("Staggered compressible multiphase solver (phasic total energies): "//trim(this%name))
      call log("  Grid: "//trim(this%amr%name))
   end subroutine amrmpcomp_print

   ! ============================================================================
   ! CHECKPOINT IO
   ! ============================================================================

   !> Register solver data for checkpoint: VOF data, face velocities, Q, and face momenta
   subroutine register_checkpoint(this,io)
      use amrio_class, only: amrio
      implicit none
      class(amrmpcomp), intent(inout) :: this
      class(amrio), intent(inout) :: io
      call this%amrmpflow%register_checkpoint(io)
      call io%add_data(this%rhoU,'rhoU')
      call io%add_data(this%rhoV,'rhoV')
      call io%add_data(this%rhoW,'rhoW')
   end subroutine register_checkpoint

   !> Restore solver data from checkpoint, then rebuild primitive variables from the restored momentum
   subroutine restore_checkpoint(this,io,dirname,time)
      use amrio_class, only: amrio
      implicit none
      class(amrmpcomp), intent(inout) :: this
      class(amrio), intent(inout) :: io
      character(len=*), intent(in) :: dirname
      real(WP), intent(in) :: time
      ! VOF data, face velocities, and Q 
      call this%amrmpflow%restore_checkpoint(io,dirname,time)
      ! Face momenta
      call io%read_data(dirname,this%rhoU,'rhoU')
      call io%read_data(dirname,this%rhoV,'rhoV')
      call io%read_data(dirname,this%rhoW,'rhoW')
      ! Rebuild primitive variables
      call this%get_primitive(this%Q,this%rhoU,this%rhoV,this%rhoW,time)
   end subroutine restore_checkpoint

end module amrmpcomp_class
