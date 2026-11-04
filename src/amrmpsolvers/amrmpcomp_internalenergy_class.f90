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

   ! Momentum flux interpolant on band dual faces that do not use CTU (modified CTU, paper eq. 55)
   integer, parameter, public :: flux_lw      =1  !< Lax-Wendroff, eq. (51)
   integer, parameter, public :: flux_centered=2  !< Centered current velocity

   ! Semi-Lagrangian face flux components: phasic Q (1:4), phasic pressure-volume (5:6), and CTU donor
   ! moments (7:18) sum(m*u_a(b)) and sum(m*u_a(b+e_a)) per phase and momentum axis, see ictu
   integer, parameter :: nSLF=18

   ! Unit vectors
   integer, dimension(3,3), parameter :: eu=reshape([1,0,0,0,1,0,0,0,1],[3,3])

   !> Pointer to one FAB, for direction-indexed access
   type :: fab_ptr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p=>null()
   end type fab_ptr

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

      ! Staggered mixture momentum
      type(amrdata) :: rhoU,rhoUold
      type(amrdata) :: rhoV,rhoVold
      type(amrdata) :: rhoW,rhoWold

      ! Cell-centered mixture velocity
      type(amrdata) :: UVW

      ! Mixture pressure
      type(amrdata) :: P

      ! Mixture speed of sound
      type(amrdata) :: C

      ! Current phasic primitive variables
      type(amrdata) :: RHOL,RHOLold,RHOG,RHOGold         !< Phasic densities
      type(amrdata) :: IL,IG                             !< Phasic internal energies
      type(amrdata) :: RHOIL,RHOIG,RHOILold,RHOIGold     !< Phasic internal energies per unit volume (rho*I)
      type(amrdata) :: PL,PG,PLold,PGold                 !< Phasic pressures
      type(amrdata) :: TL,TG,TLold,TGold                 !< Phasic temperatures

      ! Physical properties
      type(amrdata) :: visc              !< Dynamic viscosity
      type(amrdata) :: beta              !< Bulk viscosity
      type(amrdata) :: diffL,diffG       !< Phasic heat diffusivities (thermal conductivities, liquid and gas), zero after regrid

      ! Flux workspaces of get_dQdt: built with the level layouts in reset_fields (at regrid), not at every call
      type(amrex_multifab), dimension(:), allocatable :: Fx,Fy,Fz     !< Phasic mass and internal energy fluxes (faces)
      type(amrex_multifab), dimension(:), allocatable :: FUx,FVy,FWz  !< Diagonal momentum fluxes (cell centers)
      type(amrex_multifab), dimension(:), allocatable :: FUy,FVx      !< Momentum cross-fluxes (xy-edges)
      type(amrex_multifab), dimension(:), allocatable :: FUz,FWx      !< Momentum cross-fluxes (xz-edges)
      type(amrex_multifab), dimension(:), allocatable :: FVz,FWy      !< Momentum cross-fluxes (yz-edges)
      type(amrex_multifab) :: band                                    !< Semi-Lagrangian band (finest level)
      type(amrex_multifab), dimension(3) :: SLF                       !< Band face fluxes of phasic Q, pressure-volume, and CTU donor moments (finest level)
      type(amrex_multifab), dimension(3) :: SLV                       !< Band face volume moment fluxes (finest level)

      ! Momentum flux interpolant in the band: per phase, CTU where a neighboring staggered phase fraction<alpha_ctu,
      ! flux_interp otherwise (modified CTU, eq. 55)
      real(WP) :: alpha_ctu=0.5_WP       !< Staggered phase fraction below which CTU is used (must be <=1/2)
      integer  :: flux_interp=flux_lw    !< flux_lw or flux_centered

      ! CFL numbers
      real(WP) :: CFLp=0.0_WP                                !< Pressure+convection
      real(WP) :: CFLa_x=0.0_WP,CFLa_y=0.0_WP,CFLa_z=0.0_WP  !< Acoustic
      real(WP) :: CFLv_x=0.0_WP,CFLv_y=0.0_WP,CFLv_z=0.0_WP  !< Viscous

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
      real(WP) :: nctuL=0.0_WP,nctuG=0.0_WP      !< Band dual faces that used CTU in the last get_dQdt call (liquid, gas)
      real(WP) :: nrelax_fail=0.0_WP             !< Mixed cells whose relaxation did not return RELAX_OK in the last call

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
      procedure :: prepare_midpoint          !< Set face velocities and phasic primitives (incl. temperatures) to the midpoint of the old state and current iterate
      procedure :: get_face_velocity         !< Update face velocities from momentum and fill staggered ghosts
      procedure :: get_cell_velocity         !< Interpolate face velocities to cell centers in UVW
      procedure :: sync_momentum             !< Reconcile shared faces and restrict evolved momentum
      procedure :: apply_momentum_bc         !< Recover physical boundary momentum from filled Q and velocity
      procedure :: get_momentum              !< Set momentum from face velocity (initialization and regrid only)
      ! Physics procedures
      procedure :: get_primitive             !< Get phasic and mixture primitive variables from conserved variables
      procedure :: get_dQdt                  !< Compute time derivative of Q and face momenta over a step: SL band transport and bulk fluxes, with pressure and heat conduction
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
      ictu=6+6*(phase-1)+3*(side-1)+a
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

      ! Q layout: 4 phasic components, momentum is staggered. Ghost cell counts follow the 3D stencil reach:
      ! - nover=3 (Q, VF/PLIC, phasic densities): the momentum flux at the cell beyond the last valid face needs the WENO3
      !   mass flux one face further out, whose upwind stencil reaches 3 cells
      ! - nover_vel=2 (face velocities): Lax-Wendroff interpolation, transverse offsets around the outer dual faces, and
      !   RK2 vertex projection of the ghost band faces (CFL<1)
      ! - nover_old=2 (start-of-step state read by the SL donors): the band is also swept on one ghost face layer, whose
      !   donors lie within 2 cells at CFL<1; a donor beyond the halo stops the run. For headroom at CFL near 1, set both
      !   nover_old and nover_vel to 3 before initialize (the halo also covers the time-n velocity faces of the donors)
      this%nQ=4; this%nover=max(this%nover,3); this%nover_vel=max(this%nover_vel,2); this%nover_old=max(this%nover_old,2)
      ! Compressible face velocities are not divergence-free
      this%interp_vel=interp_face_lin
      call this%amrmpflow%initialize(amr,name); call this%set_parent()

      ! Initialize staggered momentum components: valid faces only, as the RHS reads Q and velocity ghosts, not momentum ghosts
      call this%rhoU%initialize(amr,name='rhoU',ncomp=1,ng=0,nodal=[.true. ,.false.,.false.]); this%rhoU%parent=>this
      call this%rhoV%initialize(amr,name='rhoV',ncomp=1,ng=0,nodal=[.false.,.true. ,.false.]); this%rhoV%parent=>this
      call this%rhoW%initialize(amr,name='rhoW',ncomp=1,ng=0,nodal=[.false.,.false.,.true. ]); this%rhoW%parent=>this
      call this%rhoUold%initialize(amr,name='rhoUold',ncomp=1,ng=0,nodal=[.true. ,.false.,.false.]); this%rhoUold%parent=>this
      call this%rhoVold%initialize(amr,name='rhoVold',ncomp=1,ng=0,nodal=[.false.,.true. ,.false.]); this%rhoVold%parent=>this
      call this%rhoWold%initialize(amr,name='rhoWold',ncomp=1,ng=0,nodal=[.false.,.false.,.true. ]); this%rhoWold%parent=>this

      ! Initialize mixture variables: UVW and C are only read in valid cells, P also in one ghost layer (get_cfl)
      call this%UVW%initialize(amr,name='UVW',ncomp=3,ng=0); this%UVW%parent=>this
      call this%P%initialize  (amr,name='P'  ,ncomp=1,ng=1); this%P%parent  =>this
      call this%C%initialize  (amr,name='C'  ,ncomp=1,ng=0); this%C%parent  =>this

      ! Initialize phasic primitive variables: densities over the WENO3 reach, pressures over the old state reach (they are
      ! stored into PLold/PGold for the SL donors, and the momentum flux needs one layer), temperatures over one layer (heat
      ! conduction), diagnostics in valid cells only
      call this%RHOL%initialize (amr,name='RHOL' ,ncomp=1,ng=this%nover); this%RHOL%parent =>this
      call this%RHOG%initialize (amr,name='RHOG' ,ncomp=1,ng=this%nover); this%RHOG%parent =>this
      call this%IL%initialize   (amr,name='IL'   ,ncomp=1,ng=0         ); this%IL%parent   =>this
      call this%IG%initialize   (amr,name='IG'   ,ncomp=1,ng=0         ); this%IG%parent   =>this
      call this%RHOIL%initialize(amr,name='RHOIL',ncomp=1,ng=this%nover); this%RHOIL%parent=>this
      call this%RHOIG%initialize(amr,name='RHOIG',ncomp=1,ng=this%nover); this%RHOIG%parent=>this
      call this%PL%initialize   (amr,name='PL'   ,ncomp=1,ng=this%nover_old); this%PL%parent=>this
      call this%PG%initialize   (amr,name='PG'   ,ncomp=1,ng=this%nover_old); this%PG%parent=>this
      call this%TL%initialize   (amr,name='TL'   ,ncomp=1,ng=1         ); this%TL%parent   =>this
      call this%TG%initialize   (amr,name='TG'   ,ncomp=1,ng=1         ); this%TG%parent   =>this

      ! Initialize start-of-step phasic primitive variables: densities feed the midpoint densities, pressures the SL donors,
      ! temperatures the midpoint temperatures
      call this%RHOLold%initialize(amr,name='RHOLold',ncomp=1,ng=this%nover); this%RHOLold%parent=>this
      call this%RHOGold%initialize(amr,name='RHOGold',ncomp=1,ng=this%nover); this%RHOGold%parent=>this
      call this%RHOILold%initialize(amr,name='RHOILold',ncomp=1,ng=this%nover); this%RHOILold%parent=>this
      call this%RHOIGold%initialize(amr,name='RHOIGold',ncomp=1,ng=this%nover); this%RHOIGold%parent=>this
      call this%PLold%initialize  (amr,name='PLold'  ,ncomp=1,ng=this%nover_old); this%PLold%parent  =>this
      call this%PGold%initialize  (amr,name='PGold'  ,ncomp=1,ng=this%nover_old); this%PGold%parent  =>this
      call this%TLold%initialize  (amr,name='TLold'  ,ncomp=1,ng=this%TL%ng    ); this%TLold%parent  =>this
      call this%TGold%initialize  (amr,name='TGold'  ,ncomp=1,ng=this%TG%ng    ); this%TGold%parent  =>this

      ! Initialize physical properties (Neumann BCs on those): one ghost layer for the edge and cell momentum fluxes, and
      ! for the conduction face fluxes
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
      call this%UVW%finalize()
      call this%P%finalize()
      call this%C%finalize()
      call this%RHOL%finalize(); call this%RHOG%finalize()
      call this%IL%finalize(); call this%IG%finalize()
      call this%RHOIL%finalize(); call this%RHOIG%finalize()
      call this%PL%finalize(); call this%PG%finalize()
      call this%TL%finalize(); call this%TG%finalize()
      call this%RHOLold%finalize(); call this%RHOGold%finalize()
      call this%RHOILold%finalize(); call this%RHOIGold%finalize()
      call this%PLold%finalize(); call this%PGold%finalize()
      call this%TLold%finalize(); call this%TGold%finalize()
      call this%visc%finalize(); call this%beta%finalize()
      call this%diffL%finalize(); call this%diffG%finalize()
      call destroy_fluxes(this)
      if (allocated(this%Fx)) deallocate(this%Fx,this%Fy,this%Fz,this%FUx,this%FVy,this%FWz,this%FUy,this%FVx,this%FUz,this%FWx,this%FVz,this%FWy)
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
      call this%UVW%setval(val=0.0_WP,lvl=lvl)
      call this%P%setval(val=0.0_WP,lvl=lvl)
      call this%C%setval(val=0.0_WP,lvl=lvl)
      call this%RHOL%setval(val=0.0_WP,lvl=lvl); call this%RHOG%setval(val=0.0_WP,lvl=lvl)
      call this%IL%setval(val=0.0_WP,lvl=lvl); call this%IG%setval(val=0.0_WP,lvl=lvl)
      call this%RHOIL%setval(val=0.0_WP,lvl=lvl); call this%RHOIG%setval(val=0.0_WP,lvl=lvl)
      call this%PL%setval(val=0.0_WP,lvl=lvl); call this%PG%setval(val=0.0_WP,lvl=lvl)
      call this%TL%setval(val=0.0_WP,lvl=lvl); call this%TG%setval(val=0.0_WP,lvl=lvl)
      call this%RHOLold%setval(val=0.0_WP,lvl=lvl); call this%RHOGold%setval(val=0.0_WP,lvl=lvl)
      call this%RHOILold%setval(val=0.0_WP,lvl=lvl); call this%RHOIGold%setval(val=0.0_WP,lvl=lvl)
      call this%PLold%setval(val=0.0_WP,lvl=lvl); call this%PGold%setval(val=0.0_WP,lvl=lvl)
      call this%TLold%setval(val=0.0_WP,lvl=lvl); call this%TGold%setval(val=0.0_WP,lvl=lvl)
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
   !> Momentum is rebuilt from the remade face velocity in post_regrid
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
      call this%UVW%clear_level(lvl)
      call this%P%clear_level(lvl)
      call this%C%clear_level(lvl)
      call this%RHOL%clear_level(lvl); call this%RHOG%clear_level(lvl)
      call this%IL%clear_level(lvl); call this%IG%clear_level(lvl)
      call this%RHOIL%clear_level(lvl); call this%RHOIG%clear_level(lvl)
      call this%PL%clear_level(lvl); call this%PG%clear_level(lvl)
      call this%TL%clear_level(lvl); call this%TG%clear_level(lvl)
      call this%RHOLold%clear_level(lvl); call this%RHOGold%clear_level(lvl)
      call this%RHOILold%clear_level(lvl); call this%RHOIGold%clear_level(lvl)
      call this%PLold%clear_level(lvl); call this%PGold%clear_level(lvl)
      call this%TLold%clear_level(lvl); call this%TGold%clear_level(lvl)
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
      call this%UVW%reset_level(lvl,ba,dm)
      call this%P%reset_level(lvl,ba,dm)
      call this%C%reset_level(lvl,ba,dm)
      call this%RHOL%reset_level(lvl,ba,dm); call this%RHOG%reset_level(lvl,ba,dm)
      call this%IL%reset_level(lvl,ba,dm); call this%IG%reset_level(lvl,ba,dm)
      call this%RHOIL%reset_level(lvl,ba,dm); call this%RHOIG%reset_level(lvl,ba,dm)
      call this%PL%reset_level(lvl,ba,dm); call this%PG%reset_level(lvl,ba,dm)
      call this%TL%reset_level(lvl,ba,dm); call this%TG%reset_level(lvl,ba,dm)
      call this%RHOLold%reset_level(lvl,ba,dm); call this%RHOGold%reset_level(lvl,ba,dm)
      call this%RHOILold%reset_level(lvl,ba,dm); call this%RHOIGold%reset_level(lvl,ba,dm)
      call this%PLold%reset_level(lvl,ba,dm); call this%PGold%reset_level(lvl,ba,dm)
      call this%TLold%reset_level(lvl,ba,dm); call this%TGold%reset_level(lvl,ba,dm)
      call this%visc%reset_level(lvl,ba,dm); call this%beta%reset_level(lvl,ba,dm)
      call this%diffL%reset_level(lvl,ba,dm); call this%diffG%reset_level(lvl,ba,dm)
      call build_fluxes(this,lvl,ba,dm)
   end subroutine reset_fields

   !> Build the get_dQdt flux workspaces of a level with its new layout: phasic fluxes and SLF with one ghost face layer
   !> (the faces [lo-1,hi+2] along the normal and [lo-1,hi+1] across that the momentum fluxes read), diagonal momentum
   !> fluxes with one ghost cell, edge fluxes and SLV on valid faces only
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
      ! The band (2 ghost layers for the flux stencils) and its face workspaces only exist on the finest level
      if (lvl.eq.this%amr%maxlvl) then
         call mfab_rebuild(this%band,ba,dm,nc=1,ng=2)
         do a=1,3
            call mfab_rebuild(this%SLF(a),ba,dm,nc=nSLF,ng=1,nodal=eu(:,a).eq.1)
            call mfab_rebuild(this%SLV(a),ba,dm,nc=8   ,ng=0,nodal=eu(:,a).eq.1)
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
         if (lvl.lt.this%amr%maxlvl) return
      else
         call amrex_multifab_destroy(this%Fx ); call amrex_multifab_destroy(this%Fy ); call amrex_multifab_destroy(this%Fz )
         call amrex_multifab_destroy(this%FUx); call amrex_multifab_destroy(this%FVy); call amrex_multifab_destroy(this%FWz)
         call amrex_multifab_destroy(this%FUy); call amrex_multifab_destroy(this%FVx)
         call amrex_multifab_destroy(this%FUz); call amrex_multifab_destroy(this%FWx)
         call amrex_multifab_destroy(this%FVz); call amrex_multifab_destroy(this%FWy)
      end if
      ! Finest-level workspaces
      call amrex_multifab_destroy(this%band)
      call amrex_multifab_destroy(this%SLF)
      call amrex_multifab_destroy(this%SLV)
   end subroutine destroy_fluxes

   !> Override post_regrid: average down and fill, rebuild momentum from the regridded face
   !> velocity (not conservative across refinement changes), then rebuild primitive variables
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
   !> and phasic primitives used by the semi-Lagrangian donors and the midpoint state
   subroutine store_old(this)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      call this%amrmpflow%store_old()
      call this%rhoUold%copy(src=this%rhoU)
      call this%rhoVold%copy(src=this%rhoV)
      call this%rhoWold%copy(src=this%rhoW)
      call this%RHOLold%copy(src=this%RHOL)
      call this%RHOGold%copy(src=this%RHOG)
      call this%RHOILold%copy(src=this%RHOIL)
      call this%RHOIGold%copy(src=this%RHOIG)
      call this%PLold%copy(src=this%PL)
      call this%PGold%copy(src=this%PG)
      call this%TLold%copy(src=this%TL)
      call this%TGold%copy(src=this%TG)
   end subroutine store_old

   !> Get midpoint for primitive quantities: face velocities, and the phasic RHO, RHOI, P, and T
   subroutine prepare_midpoint(this)
      use mpi_f08, only: MPI_Wtime
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP) :: tw0
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pRHOL,pRHOIL,pPL,pTL,pRHOLold,pRHOILold,pPLold,pTLold
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pRHOG,pRHOIG,pPG,pTG,pRHOGold,pRHOIGold,pPGold,pTGold
      logical :: Iold,Inew
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
            pRHOIL=>this%RHOIL%mf(lvl)%dataptr(mfi); pRHOILold=>this%RHOILold%mf(lvl)%dataptr(mfi)
            pPL  =>this%PL%mf(lvl)%dataptr(mfi);   pPLold  =>this%PLold%mf(lvl)%dataptr(mfi)
            pRHOG=>this%RHOG%mf(lvl)%dataptr(mfi); pRHOGold=>this%RHOGold%mf(lvl)%dataptr(mfi)
            pRHOIG=>this%RHOIG%mf(lvl)%dataptr(mfi); pRHOIGold=>this%RHOIGold%mf(lvl)%dataptr(mfi)
            pPG  =>this%PG%mf(lvl)%dataptr(mfi);   pPGold  =>this%PGold%mf(lvl)%dataptr(mfi)
            pTL  =>this%TL%mf(lvl)%dataptr(mfi);   pTLold  =>this%TLold%mf(lvl)%dataptr(mfi)
            pTG  =>this%TG%mf(lvl)%dataptr(mfi);   pTGold  =>this%TGold%mf(lvl)%dataptr(mfi)
            ! Phasic pressures over their ghost cells, first since phase presence is read from the unmodified RHOIL/RHOIG
            bx=mfi%growntilebox(min(this%PL%ng,this%PG%ng))
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Liquid
               Iold=pRHOILold(i,j,k,1).gt.0.0_WP; Inew=pRHOIL(i,j,k,1).gt.0.0_WP
               pPL   (i,j,k,1)=midval(pPLold   (i,j,k,1),pPL   (i,j,k,1),Iold,Inew)
               ! Gas
               Iold=pRHOIGold(i,j,k,1).gt.0.0_WP; Inew=pRHOIG(i,j,k,1).gt.0.0_WP
               pPG   (i,j,k,1)=midval(pPGold   (i,j,k,1),pPG   (i,j,k,1),Iold,Inew)
            end do; end do; end do
            ! Phasic temperatures over their ghost cells (heat conduction), also before RHOIL/RHOIG are modified
            bx=mfi%growntilebox(min(this%TL%ng,this%TG%ng))
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Liquid
               Iold=pRHOILold(i,j,k,1).gt.0.0_WP; Inew=pRHOIL(i,j,k,1).gt.0.0_WP
               pTL   (i,j,k,1)=midval(pTLold   (i,j,k,1),pTL   (i,j,k,1),Iold,Inew)
               ! Gas
               Iold=pRHOIGold(i,j,k,1).gt.0.0_WP; Inew=pRHOIG(i,j,k,1).gt.0.0_WP
               pTG   (i,j,k,1)=midval(pTGold   (i,j,k,1),pTG   (i,j,k,1),Iold,Inew)
            end do; end do; end do
            ! Phasic densities and internal energies per unit volume over the WENO3 reach
            bx=mfi%growntilebox(this%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Liquid
               Iold=pRHOILold(i,j,k,1).gt.0.0_WP; Inew=pRHOIL(i,j,k,1).gt.0.0_WP
               pRHOL (i,j,k,1)=midval(pRHOLold (i,j,k,1),pRHOL (i,j,k,1),pRHOLold(i,j,k,1).gt.0.0_WP,pRHOL(i,j,k,1).gt.0.0_WP)
               pRHOIL(i,j,k,1)=midval(pRHOILold(i,j,k,1),pRHOIL(i,j,k,1),Iold,Inew)
               ! Gas
               Iold=pRHOIGold(i,j,k,1).gt.0.0_WP; Inew=pRHOIG(i,j,k,1).gt.0.0_WP
               pRHOG (i,j,k,1)=midval(pRHOGold (i,j,k,1),pRHOG (i,j,k,1),pRHOGold(i,j,k,1).gt.0.0_WP,pRHOG(i,j,k,1).gt.0.0_WP)
               pRHOIG(i,j,k,1)=midval(pRHOIGold(i,j,k,1),pRHOIG(i,j,k,1),Iold,Inew)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      this%wt_mid=this%wt_mid+(MPI_Wtime()-tw0)
   contains
      !> Average of the old and new values where both exist, otherwise the one that exists (zero if neither)
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
      ! Momentum, rather than derived velocity, is the field to synchronize/restrict.
      call this%sync_momentum(rhoU,rhoV,rhoW)
      ! Compute face velocities on interior faces from mixture face density
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
      ! Replace covered coarse face velocities by the average of the fine ones: the coarse face density averages a wider
      ! stencil than the fine face densities, so restricted momentum over it overshoots across steep density gradients.
      ! Covered momentum is left as restricted by sync_momentum.
      call this%average_down_velocity()
      ! Fill derived velocity ghosts
      call this%fill_velocity(time=time)
      call this%apply_momentum_bc(Q,rhoU,rhoV,rhoW)
      this%wt_fvel=this%wt_fvel+(MPI_Wtime()-tw0)
   end subroutine get_face_velocity

   !> Interpolate face velocities to cell centers in UVW, over its ghost cells if any (face velocity ghosts must be filled)
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
   !> Q must be restricted separately before recovering velocity. This does not
   !> implement coarse/fine momentum conservation.
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
   !> Valid faces are overwritten only for prescribed normal velocity (ext_dir,
   !> reflect_odd). Extrapolation outlets remain evolved unknowns. Physical ghosts
   !> use rho_face*velocity only where both adjacent density cells exist; the RHS
   !> reads Q/velocity ghosts, not the unused outermost momentum ghosts.
   !> Same-level/periodic momentum ghosts are exchanged, not reconstructed.
   subroutine apply_momentum_bc(this,Q,rhoU,rhoV,rhoW)
      implicit none
      class(amrmpcomp), intent(inout) :: this
      type(amrdata), intent(in) :: Q
      type(amrdata), intent(inout) :: rhoU,rhoV,rhoW
      integer :: lvl
      type(amrex_mfiter) :: mfi
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pRHOU,pRHOV,pRHOW,pU,pV,pW
      do lvl=0,this%amr%clvl()
         ! Whole FABs including ghosts are processed, so disable tiling
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
      ! Propagate prescribed boundary values and fill same-level/periodic ghosts.
      call rhoU%sync(); call rhoV%sync(); call rhoW%sync()
      
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
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF,pRHOL,pRHOG,pIL,pIG,pRHOIL,pRHOIG,pPL,pPG,pTL,pTG,pP,pC
      real(WP) :: CL,CG
      ! Check passed Q is as expected
      if (Q%ncomp.ne.this%nQ) call die('[amrmpcomp get_primitive] Q has wrong number of components')
      if (Q%ng.lt.this%nover) call die('[amrmpcomp get_primitive] Q must have at least nover ghost cells')
      tw0=MPI_Wtime()
      ! Loop over levels
      do lvl=0,this%amr%clvl()
         call this%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pQ  =>Q%mf(lvl)%dataptr(mfi)
            pVF =>this%VF%mf(lvl)%dataptr(mfi)
            pRHOL=>this%RHOL%mf(lvl)%dataptr(mfi); pRHOG=>this%RHOG%mf(lvl)%dataptr(mfi)
            pIL =>this%IL%mf(lvl)%dataptr(mfi);   pIG =>this%IG%mf(lvl)%dataptr(mfi)
            pRHOIL=>this%RHOIL%mf(lvl)%dataptr(mfi); pRHOIG=>this%RHOIG%mf(lvl)%dataptr(mfi)
            pPL =>this%PL%mf(lvl)%dataptr(mfi);   pPG =>this%PG%mf(lvl)%dataptr(mfi)
            pTL =>this%TL%mf(lvl)%dataptr(mfi);   pTG =>this%TG%mf(lvl)%dataptr(mfi)
            pP  =>this%P%mf(lvl)%dataptr(mfi)
            pC  =>this%C%mf(lvl)%dataptr(mfi)
            ! Each quantity is computed over its own ghost cells, from the widest (densities) to the narrowest (diagnostics);
            ! a phasic internal energy exists (RHOI>0) where the phase is present with positive mass and energy
            ! Phasic densities and internal energies per unit volume over the WENO3 reach
            bx=mfi%growntilebox(this%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Liquid
               if (pVF(i,j,k,1).ge.VFlo.and.pQ(i,j,k,1).gt.0.0_WP) then
                  pRHOL(i,j,k,1)=pQ(i,j,k,1)/pVF(i,j,k,1)
                  if (pQ(i,j,k,3).gt.0.0_WP) then
                     pRHOIL(i,j,k,1)=pQ(i,j,k,3)/pVF(i,j,k,1)
                  else
                     pRHOIL(i,j,k,1)=0.0_WP
                  end if
               else
                  pRHOL(i,j,k,1)=0.0_WP; pRHOIL(i,j,k,1)=0.0_WP
               end if
               ! Gas
               if (pVF(i,j,k,1).le.VFhi.and.pQ(i,j,k,2).gt.0.0_WP) then
                  pRHOG(i,j,k,1)=pQ(i,j,k,2)/(1.0_WP-pVF(i,j,k,1))
                  if (pQ(i,j,k,4).gt.0.0_WP) then
                     pRHOIG(i,j,k,1)=pQ(i,j,k,4)/(1.0_WP-pVF(i,j,k,1))
                  else
                     pRHOIG(i,j,k,1)=0.0_WP
                  end if
               else
                  pRHOG(i,j,k,1)=0.0_WP; pRHOIG(i,j,k,1)=0.0_WP
               end if
            end do; end do; end do
            ! Phasic pressures over their ghost cells (SL donor reach)
            bx=mfi%growntilebox(min(this%PL%ng,this%PG%ng))
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pPL(i,j,k,1)=0.0_WP
               if (pRHOIL(i,j,k,1).gt.0.0_WP) pPL(i,j,k,1)=this%liq%get_p_from_rho_e(rho=pRHOL(i,j,k,1),e=pQ(i,j,k,3)/pQ(i,j,k,1),y=[1.0_WP])
               pPG(i,j,k,1)=0.0_WP
               if (pRHOIG(i,j,k,1).gt.0.0_WP) pPG(i,j,k,1)=this%gas%get_p_from_rho_e(rho=pRHOG(i,j,k,1),e=pQ(i,j,k,4)/pQ(i,j,k,2),y=[1.0_WP])
            end do; end do; end do
            ! Mixture pressure over its ghost cells
            bx=mfi%growntilebox(min(this%P%ng,this%PL%ng,this%PG%ng))
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pP(i,j,k,1)=pVF(i,j,k,1)*pPL(i,j,k,1)+(1.0_WP-pVF(i,j,k,1))*pPG(i,j,k,1)
            end do; end do; end do
            ! Phasic temperatures over their ghost cells (heat conduction), zero where the phase does not exist
            bx=mfi%growntilebox(min(this%TL%ng,this%TG%ng))
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pTL(i,j,k,1)=0.0_WP
               if (pRHOIL(i,j,k,1).gt.0.0_WP) pTL(i,j,k,1)=this%liq%get_T_from_rho_e(rho=pRHOL(i,j,k,1),e=pQ(i,j,k,3)/pQ(i,j,k,1),y=[1.0_WP])
               pTG(i,j,k,1)=0.0_WP
               if (pRHOIG(i,j,k,1).gt.0.0_WP) pTG(i,j,k,1)=this%gas%get_T_from_rho_e(rho=pRHOG(i,j,k,1),e=pQ(i,j,k,4)/pQ(i,j,k,2),y=[1.0_WP])
            end do; end do; end do
            ! Phasic internal energies and mixture speed of sound over their ghost cells (valid cells by default)
            bx=mfi%growntilebox(min(this%IL%ng,this%IG%ng,this%C%ng))
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Liquid
               CL=0.0_WP; pIL(i,j,k,1)=0.0_WP
               if (pRHOIL(i,j,k,1).gt.0.0_WP) then
                  pIL(i,j,k,1)=pQ(i,j,k,3)/pQ(i,j,k,1)
                  CL          =this%liq%get_c_from_rho_e(rho=pRHOL(i,j,k,1),e=pIL(i,j,k,1),y=[1.0_WP])
               end if
               ! Gas
               CG=0.0_WP; pIG(i,j,k,1)=0.0_WP
               if (pRHOIG(i,j,k,1).gt.0.0_WP) then
                  pIG(i,j,k,1)=pQ(i,j,k,4)/pQ(i,j,k,2)
                  CG          =this%gas%get_c_from_rho_e(rho=pRHOG(i,j,k,1),e=pIG(i,j,k,1),y=[1.0_WP])
               end if
               ! Mixture speed of sound
               pC(i,j,k,1)=sqrt((pQ(i,j,k,1)*CL**2+pQ(i,j,k,2)*CG**2)/max(pQ(i,j,k,1)+pQ(i,j,k,2),this%rho_floor))
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      this%wt_prim=this%wt_prim+(MPI_Wtime()-tw0)
      ! Face velocities from momentum, with staggered ghosts
      call this%get_face_velocity(Q,rhoU,rhoV,rhoW,time)
   end subroutine get_primitive

   !> Calculate the time derivative of Q and face momenta over a step of size dt from the current (midpoint)
   !> state, including pressure, as rates: Q=Qold+dt*dQdt. On the finest level, faces touching the
   !> semi-Lagrangian band carry geometric fluxes: the start-of-step state (Qold, VFold, CLold, CGold,
   !> PLICold, PLold, PGold) is swept over dt by flux polyhedra built from the current face velocity, whose
   !> ghosts must be filled, and VF and barycenters are updated in the band (as mpcomp). Other faces carry
   !> bulk fluxes of the phase present. Band dual faces use, per phase, the CTU momentum flux interpolant or
   !> flux_interp (Lax-Wendroff or centered) with the modified switch (Remmerswaal & Veldman, arXiv:2209.14934,
   !> eqs. (46), (51), (53)-(55), and (B.9)); CTU and Lax-Wendroff transport the time-n velocity (Uold/Vold/Wold).
   !> The mixture pressure, weighted by the geometrically transported VF of this sweep (no midpoint VF), acts as a
   !> normal stress, and the phasic pressure dilatation is added to the phasic energies with the same VF (as mpcomp).
   !> Phasic heat conduction from the midpoint phasic temperatures and diffL/diffG is added to the phasic energy
   !> fluxes on all faces, band included (as amrmpcomp_col); interphase heat exchange is left to the relaxation model.
   subroutine get_dQdt(this,dQdt,drhoUdt,drhoVdt,drhoWdt,dt,time)
      use messager,         only: die
      use parallel,         only: MPI_REAL_WP
      use mpi_f08,          only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM,MPI_Wtime
      implicit none
      class(amrmpcomp), intent(inout) :: this
      type(amrdata), intent(inout) :: dQdt,drhoUdt,drhoVdt,drhoWdt
      real(WP) :: tw0,tw1
      real(WP), intent(in) :: dt
      real(WP), intent(in) :: time
      ! Flux workspaces are solver components (built at regrid): this%Fx..this%FWy per level, and on the finest level
      ! this%band, this%SLF (band face fluxes of phasic Q, phasic pressure-volume, and CTU donor moments: rates, negative
      ! convention) and this%SLV (band face volume moment fluxes Lvol, Gvol, Lvol*Lbar, Gvol*Gbar: integrated over dt,
      ! positive along the axis). SLF and SLV are not zeroed: they are only read on band faces, which are recomputed here
      ! Finest level and its mesh info, used by the semi-Lagrangian procedures
      integer :: slvl,a,ierr
      real(WP) :: dx,dy,dz,dxi,dyi,dzi
      ! Whether the level being processed is the finest level
      logical :: on_sl
      ! Data shared with internal procedures
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQold,pVFold,pPLICold,pPLold,pPGold,pVF,pBand
      type(fab_ptr), dimension(3) :: pVel,pVelold,pF,pM
      type(amrex_box) :: tbx
      ! Cell range of the FAB where a semi-Lagrangian donor has all its start-of-step data
      integer, dimension(3) :: halo_lo,halo_hi

      tw0=MPI_Wtime(); tw1=tw0
      ! Build the band on the finest level from the start-of-step interface (2 ghost layers for the flux stencils)
      if (this%amr%clvl().lt.this%amr%maxlvl) call die('[amrmpcomp get_dQdt] the finest level must exist')
      slvl=this%amr%maxlvl
      call this%build_band(lvl=slvl,VF=this%VFold%mf(slvl),band=this%band,nband=2)
      dx=this%amr%dx(slvl); dxi=1.0_WP/dx
      dy=this%amr%dy(slvl); dyi=1.0_WP/dy
      dz=this%amr%dz(slvl); dzi=1.0_WP/dz
      this%nctuL=0.0_WP; this%nctuG=0.0_WP
      this%wt_alloc=this%wt_alloc+(MPI_Wtime()-tw1)

      ! Geometric fluxes on band faces from flux polyhedra cut by the old PLIC (as amrmpcomp_col)
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
         ! Band faces are swept on the valid faces and on the ghost face layer read by the flux loops ([lo-1,hi+2] along
         ! the normal, [lo-1,hi+1] across), so SLF needs no ghost exchange: a face shared by two boxes is computed from the
         ! same data by both. SLV (volume moments) is only used on valid faces and has no ghosts
         call this%amr%mfiter_build(slvl,mfi,tiling=.false.)
         do while (mfi%next())
            ! Get data pointers
            pQold   =>this%Qold%mf(slvl)%dataptr(mfi)
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
            ! Remap vertices of band faces back in time via RK2, including the vertices of the ghost face layer. The loop is
            ! untiled (tile box = valid box), so boxes are grown explicitly by one: AMReX's Fortran grownnodaltilebox passes
            ! dir to C++ without converting it to 0-based (see examples/amrvof_tester), so it is not used here
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
         ! Ghost faces were computed locally, so no ghost exchange. Shared valid faces are computed identically by both
         ! boxes, except periodic images, which use coordinates shifted by the domain length: the owner override is only
         ! needed with periodicity (it keeps mass, energy, and volume conservation exact there)
         if (this%amr%xper.or.this%amr%yper.or.this%amr%zper) then
            do a=1,3
               call this%SLF(a)%override_sync(this%amr%geom(slvl))
               call this%SLV(a)%override_sync(this%amr%geom(slvl))
            end do
         end if
      end block geometric_fluxes
      this%wt_slflux=this%wt_slflux+(MPI_Wtime()-tw1)

      ! Update VF and barycenters in the band from the volume moment fluxes, and set the phasic pressure of a phase
      ! appearing in a cell without a midpoint pressure from its transported pressure-volume
      tw1=MPI_Wtime()
      volume_moments: block
         integer :: i,j,k
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP) :: vol,VFo,Lvol_old,Gvol_old,Lvol_new,Gvol_new,Lvol_flux,Gvol_flux
         real(WP), dimension(3) :: Lbar_old,Gbar_old,Lbar_new,Gbar_new,Lbar_flux,Gbar_flux
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pCL,pCG,pCLold,pCGold,pRHOIL,pRHOIG,pPL,pPG,pVx,pVy,pVz,pFx,pFy,pFz
         vol=dx*dy*dz
         call this%amr%mfiter_build(slvl,mfi,tiling=.false.)
         do while (mfi%next())
            ! Get data pointers
            pBand =>this%band%dataptr(mfi)
            pVFold=>this%VFold%mf(slvl)%dataptr(mfi)
            pVF   =>this%VF%mf(slvl)%dataptr(mfi)
            pCL   =>this%CL%dataptr(mfi);    pCG   =>this%CG%dataptr(mfi)
            pCLold=>this%CLold%dataptr(mfi); pCGold=>this%CGold%dataptr(mfi)
            pRHOIL=>this%RHOIL%mf(slvl)%dataptr(mfi); pRHOIG=>this%RHOIG%mf(slvl)%dataptr(mfi)
            pPL   =>this%PL%mf(slvl)%dataptr(mfi);    pPG   =>this%PG%mf(slvl)%dataptr(mfi)
            pVx=>this%SLV(1)%dataptr(mfi); pFx=>this%SLF(1)%dataptr(mfi)
            pVy=>this%SLV(2)%dataptr(mfi); pFy=>this%SLF(2)%dataptr(mfi)
            pVz=>this%SLV(3)%dataptr(mfi); pFz=>this%SLF(3)%dataptr(mfi)
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
               ! Net volume flux (outflow positive) from SL volume moments
               Lvol_flux=pVx(i+1,j,k, 1 )-pVx(i,j,k, 1 )+pVy(i,j+1,k, 1 )-pVy(i,j,k, 1 )+pVz(i,j,k+1, 1 )-pVz(i,j,k, 1 )
               Gvol_flux=pVx(i+1,j,k, 2 )-pVx(i,j,k, 2 )+pVy(i,j+1,k, 2 )-pVy(i,j,k, 2 )+pVz(i,j,k+1, 2 )-pVz(i,j,k, 2 )
               Lbar_flux=pVx(i+1,j,k,3:5)-pVx(i,j,k,3:5)+pVy(i,j+1,k,3:5)-pVy(i,j,k,3:5)+pVz(i,j,k+1,3:5)-pVz(i,j,k,3:5)
               Gbar_flux=pVx(i+1,j,k,6:8)-pVx(i,j,k,6:8)+pVy(i,j+1,k,6:8)-pVy(i,j,k,6:8)+pVz(i,j,k+1,6:8)-pVz(i,j,k,6:8)
               ! New phasic volumes, normalized by their sum (compressible VF update)
               Lvol_new=Lvol_old-Lvol_flux
               Gvol_new=Gvol_old-Gvol_flux
               pVF(i,j,k,1)=Lvol_new/(Lvol_new+Gvol_new)
               pCL(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
               pCG(i,j,k,1:3)=pCL(i,j,k,1:3)
               ! Clip and update barycenters
               if (pVF(i,j,k,1).lt.VFlo) then
                  pVF(i,j,k,1)=0.0_WP
               else if (pVF(i,j,k,1).gt.VFhi) then
                  pVF(i,j,k,1)=1.0_WP
               else
                  ! Update barycenters from moment conservation and project forward
                  if (Lvol_new/(Lvol_new+Gvol_new).gt.vol_eps) then; Lbar_new=(Lbar_old*Lvol_old-Lbar_flux)/Lvol_new; pCL(i,j,k,1:3)=project(Lbar_new,dt); end if
                  if (Gvol_new/(Lvol_new+Gvol_new).gt.vol_eps) then; Gbar_new=(Gbar_old*Gvol_old-Gbar_flux)/Gvol_new; pCG(i,j,k,1:3)=project(Gbar_new,dt); end if
               end if
               ! Phasic pressure in cells where a phase appears and has no midpoint pressure (absent at both levels, which
               ! happens in the first sub-iteration): transported pressure-volume over new phase volume (mpcomp)
               if (VFo.lt.VFlo.and.pVF(i,j,k,1).ge.VFlo.and.pRHOIL(i,j,k,1).le.0.0_WP) pPL(i,j,k,1)=dt*((pFx(i+1,j,k,5)-pFx(i,j,k,5))*dy*dz+(pFy(i,j+1,k,5)-pFy(i,j,k,5))*dz*dx+(pFz(i,j,k+1,5)-pFz(i,j,k,5))*dx*dy)/Lvol_new
               if (VFo.gt.VFhi.and.pVF(i,j,k,1).le.VFhi.and.pRHOIG(i,j,k,1).le.0.0_WP) pPG(i,j,k,1)=dt*((pFx(i+1,j,k,6)-pFx(i,j,k,6))*dy*dz+(pFy(i,j+1,k,6)-pFy(i,j,k,6))*dz*dx+(pFz(i,j,k+1,6)-pFz(i,j,k,6))*dx*dy)/Gvol_new
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         ! Synchronize: VF on all levels, barycenters on the finest, and phasic pressures
         call this%VF%average_down(); call this%fill(time=time)
         call this%PL%sync(); call this%PG%sync()
      end block volume_moments
      this%wt_vmom=this%wt_vmom+(MPI_Wtime()-tw1)

      ! Compute fluxes for all levels
      tw1=MPI_Wtime()
      compute_fluxes: block
         integer :: lvl,i,j,k
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP) :: dxi,dyi,dzi,heat,pres,div,mu,lam,strain,stress,w,condL,condG
         real(WP), dimension(3) :: gradU,tii
         real(WP), dimension(-2: 0) :: wxp,wyp,wzp
         real(WP), dimension(-1:+1) :: wxm,wym,wzm
         real(WP), parameter :: eps=1.0e-15_WP
         logical :: bandx,bandy,bandz
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW,pRHOL,pRHOG,pRHOIL,pRHOIG,pPL,pPG,pVisc,pBeta
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pTL,pTG,pDiffL,pDiffG
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pFx,pFy,pFz,pdQdt
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pFUx,pFUy,pFUz
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pFVx,pFVy,pFVz
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pFWx,pFWy,pFWz
         ! Viscous shear dissipation stress*strain at the xy-, yz-, and zx-edges of one tile
         real(WP), dimension(:,:,:), allocatable :: dxy,dyz,dzx
         ! Zero RHS: sources are accumulated here, flux divergences are added after C/F averaging
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
               pRHOIL=>this%RHOIL%mf(lvl)%dataptr(mfi); pRHOIG=>this%RHOIG%mf(lvl)%dataptr(mfi)
               pPL=>this%PL%mf(lvl)%dataptr(mfi); pPG=>this%PG%mf(lvl)%dataptr(mfi)
               pVisc=>this%visc%mf(lvl)%dataptr(mfi)
               pBeta=>this%beta%mf(lvl)%dataptr(mfi)
               pTL=>this%TL%mf(lvl)%dataptr(mfi); pTG=>this%TG%mf(lvl)%dataptr(mfi)
               pDiffL=>this%diffL%mf(lvl)%dataptr(mfi); pDiffG=>this%diffG%mf(lvl)%dataptr(mfi)
               pQold=>this%Qold%mf(lvl)%dataptr(mfi)
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
               pdQdt=>dQdt%mf(lvl)%dataptr(mfi)
               ! Pointers shared with the band momentum flux
               pVel(1)%p=>pU;  pVel(2)%p=>pV;  pVel(3)%p=>pW
               pM(1)%p  =>pFx; pM(2)%p  =>pFy; pM(3)%p  =>pFz
               if (on_sl) then
                  pBand=>this%band%dataptr(mfi)
                  do a=1,3; pF(a)%p=>this%SLF(a)%dataptr(mfi); end do
                  pVelold(1)%p=>this%Uold%mf(lvl)%dataptr(mfi)
                  pVelold(2)%p=>this%Vold%mf(lvl)%dataptr(mfi)
                  pVelold(3)%p=>this%Wold%mf(lvl)%dataptr(mfi)
               end if
               bx=mfi%tilebox(); tbx=bx
               ! Phasic mass and internal energy fluxes: geometric on band faces, and upwind-biased WENO3 of the phasic
               ! density and internal energy per unit volume of the phase present elsewhere (mpcomp); all four components
               ! are written (the absent phase gets zero) since the persistent flux workspaces are not zeroed per call.
               ! Each direction covers the faces read by this tile: [lo-1,hi+2] along the normal (momentum fluxes at the
               ! cells lo-1..hi+1) and [lo-1,hi+1] across (edge and cell momentum fluxes), i.e., one ghost face layer.
               ! Band faces touch a band cell (the band only exists on the finest level)
               ! Phasic heat conduction (as amrmpcomp_col) is added to the energy fluxes in the same pass, on the faces whose
               ! two cells have a temperature (all but the outermost ghost faces of the FAB), so that a face rewritten by a
               ! neighboring tile keeps it: own-phase conductance through the smaller phasic aperture of the two cells, zero
               ! where the phase is absent on either side, and content-limited (limit_cond)
               ! X-faces
               do k=bx%lo(3)-1,bx%hi(3)+1; do j=bx%lo(2)-1,bx%hi(2)+1; do i=bx%lo(1)-1,bx%hi(1)+2
                  bandx=.false.; if (on_sl) bandx=maxval(pBand(i-1:i,j,k,1)).gt.0.0_WP
                  if (bandx) then
                     pFx(i,j,k,1:4)=pF(1)%p(i,j,k,1:4)
                  else if (pVF(i,j,k,1).gt.0.5_WP) then
                     ! No gas fluxes
                     pFx(i,j,k,2)=0.0_WP; pFx(i,j,k,4)=0.0_WP
                     ! Liquid mass
                     w=weno_weight((abs(pRHOL(i-1,j,k,1)-pRHOL(i-2,j,k,1))+eps)/(abs(pRHOL(i,j,k,1)-pRHOL(i-1,j,k,1))+eps)); wxp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(pRHOL(i+1,j,k,1)-pRHOL(i  ,j,k,1))+eps)/(abs(pRHOL(i,j,k,1)-pRHOL(i-1,j,k,1))+eps)); wxm=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     pFx(i,j,k,1)=-0.5_WP*(pU(i,j,k,1)+abs(pU(i,j,k,1)))*sum(wxp*pRHOL(i-2:i  ,j,k,1)) &
                     &            -0.5_WP*(pU(i,j,k,1)-abs(pU(i,j,k,1)))*sum(wxm*pRHOL(i-1:i+1,j,k,1))
                     ! Liquid internal energy
                     w=weno_weight((abs(pRHOIL(i-1,j,k,1)-pRHOIL(i-2,j,k,1))+eps)/(abs(pRHOIL(i,j,k,1)-pRHOIL(i-1,j,k,1))+eps)); wxp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(pRHOIL(i+1,j,k,1)-pRHOIL(i  ,j,k,1))+eps)/(abs(pRHOIL(i,j,k,1)-pRHOIL(i-1,j,k,1))+eps)); wxm=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     pFx(i,j,k,3)=-0.5_WP*(pU(i,j,k,1)+abs(pU(i,j,k,1)))*sum(wxp*pRHOIL(i-2:i  ,j,k,1)) &
                     &            -0.5_WP*(pU(i,j,k,1)-abs(pU(i,j,k,1)))*sum(wxm*pRHOIL(i-1:i+1,j,k,1))
                  else
                     ! No liquid fluxes
                     pFx(i,j,k,1)=0.0_WP; pFx(i,j,k,3)=0.0_WP
                     ! Gas mass
                     w=weno_weight((abs(pRHOG(i-1,j,k,1)-pRHOG(i-2,j,k,1))+eps)/(abs(pRHOG(i,j,k,1)-pRHOG(i-1,j,k,1))+eps)); wxp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(pRHOG(i+1,j,k,1)-pRHOG(i  ,j,k,1))+eps)/(abs(pRHOG(i,j,k,1)-pRHOG(i-1,j,k,1))+eps)); wxm=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     pFx(i,j,k,2)=-0.5_WP*(pU(i,j,k,1)+abs(pU(i,j,k,1)))*sum(wxp*pRHOG(i-2:i  ,j,k,1)) &
                     &            -0.5_WP*(pU(i,j,k,1)-abs(pU(i,j,k,1)))*sum(wxm*pRHOG(i-1:i+1,j,k,1))
                     ! Gas internal energy
                     w=weno_weight((abs(pRHOIG(i-1,j,k,1)-pRHOIG(i-2,j,k,1))+eps)/(abs(pRHOIG(i,j,k,1)-pRHOIG(i-1,j,k,1))+eps)); wxp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(pRHOIG(i+1,j,k,1)-pRHOIG(i  ,j,k,1))+eps)/(abs(pRHOIG(i,j,k,1)-pRHOIG(i-1,j,k,1))+eps)); wxm=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     pFx(i,j,k,4)=-0.5_WP*(pU(i,j,k,1)+abs(pU(i,j,k,1)))*sum(wxp*pRHOIG(i-2:i  ,j,k,1)) &
                     &            -0.5_WP*(pU(i,j,k,1)-abs(pU(i,j,k,1)))*sum(wxm*pRHOIG(i-1:i+1,j,k,1))
                  end if
                  if (i.gt.lbound(pTL,1).and.i.le.ubound(pTL,1)) then
                     condL=0.0_WP; if (all(pVF(i-1:i,j,k,1).ge.VFlo).and.all(pTL(i-1:i,j,k,1).gt.0.0_WP)) condL=(       minval(pVF(i-1:i,j,k,1)))*0.5_WP*sum(pDiffL(i-1:i,j,k,1))*dxi
                     condG=0.0_WP; if (all(pVF(i-1:i,j,k,1).le.VFhi).and.all(pTG(i-1:i,j,k,1).gt.0.0_WP)) condG=(1.0_WP-maxval(pVF(i-1:i,j,k,1)))*0.5_WP*sum(pDiffG(i-1:i,j,k,1))*dxi
                     pFx(i,j,k,3)=pFx(i,j,k,3)+limit_cond(condL*(pTL(i,j,k,1)-pTL(i-1,j,k,1)),pQold(i,j,k,3),pQold(i-1,j,k,3),dxi)
                     pFx(i,j,k,4)=pFx(i,j,k,4)+limit_cond(condG*(pTG(i,j,k,1)-pTG(i-1,j,k,1)),pQold(i,j,k,4),pQold(i-1,j,k,4),dxi)
                  end if
               end do; end do; end do
               ! Y-faces
               do k=bx%lo(3)-1,bx%hi(3)+1; do j=bx%lo(2)-1,bx%hi(2)+2; do i=bx%lo(1)-1,bx%hi(1)+1
                  bandy=.false.; if (on_sl) bandy=maxval(pBand(i,j-1:j,k,1)).gt.0.0_WP
                  if (bandy) then
                     pFy(i,j,k,1:4)=pF(2)%p(i,j,k,1:4)
                  else if (pVF(i,j,k,1).gt.0.5_WP) then
                     ! No gas fluxes
                     pFy(i,j,k,2)=0.0_WP; pFy(i,j,k,4)=0.0_WP
                     ! Liquid mass
                     w=weno_weight((abs(pRHOL(i,j-1,k,1)-pRHOL(i,j-2,k,1))+eps)/(abs(pRHOL(i,j,k,1)-pRHOL(i,j-1,k,1))+eps)); wyp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(pRHOL(i,j+1,k,1)-pRHOL(i,j  ,k,1))+eps)/(abs(pRHOL(i,j,k,1)-pRHOL(i,j-1,k,1))+eps)); wym=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     pFy(i,j,k,1)=-0.5_WP*(pV(i,j,k,1)+abs(pV(i,j,k,1)))*sum(wyp*pRHOL(i,j-2:j  ,k,1)) &
                     &            -0.5_WP*(pV(i,j,k,1)-abs(pV(i,j,k,1)))*sum(wym*pRHOL(i,j-1:j+1,k,1))
                     ! Liquid internal energy
                     w=weno_weight((abs(pRHOIL(i,j-1,k,1)-pRHOIL(i,j-2,k,1))+eps)/(abs(pRHOIL(i,j,k,1)-pRHOIL(i,j-1,k,1))+eps)); wyp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(pRHOIL(i,j+1,k,1)-pRHOIL(i,j  ,k,1))+eps)/(abs(pRHOIL(i,j,k,1)-pRHOIL(i,j-1,k,1))+eps)); wym=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     pFy(i,j,k,3)=-0.5_WP*(pV(i,j,k,1)+abs(pV(i,j,k,1)))*sum(wyp*pRHOIL(i,j-2:j  ,k,1)) &
                     &            -0.5_WP*(pV(i,j,k,1)-abs(pV(i,j,k,1)))*sum(wym*pRHOIL(i,j-1:j+1,k,1))
                  else
                     ! No liquid fluxes
                     pFy(i,j,k,1)=0.0_WP; pFy(i,j,k,3)=0.0_WP
                     ! Gas mass
                     w=weno_weight((abs(pRHOG(i,j-1,k,1)-pRHOG(i,j-2,k,1))+eps)/(abs(pRHOG(i,j,k,1)-pRHOG(i,j-1,k,1))+eps)); wyp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(pRHOG(i,j+1,k,1)-pRHOG(i,j  ,k,1))+eps)/(abs(pRHOG(i,j,k,1)-pRHOG(i,j-1,k,1))+eps)); wym=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     pFy(i,j,k,2)=-0.5_WP*(pV(i,j,k,1)+abs(pV(i,j,k,1)))*sum(wyp*pRHOG(i,j-2:j  ,k,1)) &
                     &            -0.5_WP*(pV(i,j,k,1)-abs(pV(i,j,k,1)))*sum(wym*pRHOG(i,j-1:j+1,k,1))
                     ! Gas internal energy
                     w=weno_weight((abs(pRHOIG(i,j-1,k,1)-pRHOIG(i,j-2,k,1))+eps)/(abs(pRHOIG(i,j,k,1)-pRHOIG(i,j-1,k,1))+eps)); wyp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(pRHOIG(i,j+1,k,1)-pRHOIG(i,j  ,k,1))+eps)/(abs(pRHOIG(i,j,k,1)-pRHOIG(i,j-1,k,1))+eps)); wym=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     pFy(i,j,k,4)=-0.5_WP*(pV(i,j,k,1)+abs(pV(i,j,k,1)))*sum(wyp*pRHOIG(i,j-2:j  ,k,1)) &
                     &            -0.5_WP*(pV(i,j,k,1)-abs(pV(i,j,k,1)))*sum(wym*pRHOIG(i,j-1:j+1,k,1))
                  end if
                  if (j.gt.lbound(pTL,2).and.j.le.ubound(pTL,2)) then
                     condL=0.0_WP; if (all(pVF(i,j-1:j,k,1).ge.VFlo).and.all(pTL(i,j-1:j,k,1).gt.0.0_WP)) condL=(       minval(pVF(i,j-1:j,k,1)))*0.5_WP*sum(pDiffL(i,j-1:j,k,1))*dyi
                     condG=0.0_WP; if (all(pVF(i,j-1:j,k,1).le.VFhi).and.all(pTG(i,j-1:j,k,1).gt.0.0_WP)) condG=(1.0_WP-maxval(pVF(i,j-1:j,k,1)))*0.5_WP*sum(pDiffG(i,j-1:j,k,1))*dyi
                     pFy(i,j,k,3)=pFy(i,j,k,3)+limit_cond(condL*(pTL(i,j,k,1)-pTL(i,j-1,k,1)),pQold(i,j,k,3),pQold(i,j-1,k,3),dyi)
                     pFy(i,j,k,4)=pFy(i,j,k,4)+limit_cond(condG*(pTG(i,j,k,1)-pTG(i,j-1,k,1)),pQold(i,j,k,4),pQold(i,j-1,k,4),dyi)
                  end if
               end do; end do; end do
               ! Z-faces
               do k=bx%lo(3)-1,bx%hi(3)+2; do j=bx%lo(2)-1,bx%hi(2)+1; do i=bx%lo(1)-1,bx%hi(1)+1
                  bandz=.false.; if (on_sl) bandz=maxval(pBand(i,j,k-1:k,1)).gt.0.0_WP
                  if (bandz) then
                     pFz(i,j,k,1:4)=pF(3)%p(i,j,k,1:4)
                  else if (pVF(i,j,k,1).gt.0.5_WP) then
                     ! No gas fluxes
                     pFz(i,j,k,2)=0.0_WP; pFz(i,j,k,4)=0.0_WP
                     ! Liquid mass
                     w=weno_weight((abs(pRHOL(i,j,k-1,1)-pRHOL(i,j,k-2,1))+eps)/(abs(pRHOL(i,j,k,1)-pRHOL(i,j,k-1,1))+eps)); wzp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(pRHOL(i,j,k+1,1)-pRHOL(i,j,k  ,1))+eps)/(abs(pRHOL(i,j,k,1)-pRHOL(i,j,k-1,1))+eps)); wzm=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     pFz(i,j,k,1)=-0.5_WP*(pW(i,j,k,1)+abs(pW(i,j,k,1)))*sum(wzp*pRHOL(i,j,k-2:k  ,1)) &
                     &            -0.5_WP*(pW(i,j,k,1)-abs(pW(i,j,k,1)))*sum(wzm*pRHOL(i,j,k-1:k+1,1))
                     ! Liquid internal energy
                     w=weno_weight((abs(pRHOIL(i,j,k-1,1)-pRHOIL(i,j,k-2,1))+eps)/(abs(pRHOIL(i,j,k,1)-pRHOIL(i,j,k-1,1))+eps)); wzp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(pRHOIL(i,j,k+1,1)-pRHOIL(i,j,k  ,1))+eps)/(abs(pRHOIL(i,j,k,1)-pRHOIL(i,j,k-1,1))+eps)); wzm=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     pFz(i,j,k,3)=-0.5_WP*(pW(i,j,k,1)+abs(pW(i,j,k,1)))*sum(wzp*pRHOIL(i,j,k-2:k  ,1)) &
                     &            -0.5_WP*(pW(i,j,k,1)-abs(pW(i,j,k,1)))*sum(wzm*pRHOIL(i,j,k-1:k+1,1))
                  else
                     ! No liquid fluxes
                     pFz(i,j,k,1)=0.0_WP; pFz(i,j,k,3)=0.0_WP
                     ! Gas mass
                     w=weno_weight((abs(pRHOG(i,j,k-1,1)-pRHOG(i,j,k-2,1))+eps)/(abs(pRHOG(i,j,k,1)-pRHOG(i,j,k-1,1))+eps)); wzp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(pRHOG(i,j,k+1,1)-pRHOG(i,j,k  ,1))+eps)/(abs(pRHOG(i,j,k,1)-pRHOG(i,j,k-1,1))+eps)); wzm=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     pFz(i,j,k,2)=-0.5_WP*(pW(i,j,k,1)+abs(pW(i,j,k,1)))*sum(wzp*pRHOG(i,j,k-2:k  ,1)) &
                     &            -0.5_WP*(pW(i,j,k,1)-abs(pW(i,j,k,1)))*sum(wzm*pRHOG(i,j,k-1:k+1,1))
                     ! Gas internal energy
                     w=weno_weight((abs(pRHOIG(i,j,k-1,1)-pRHOIG(i,j,k-2,1))+eps)/(abs(pRHOIG(i,j,k,1)-pRHOIG(i,j,k-1,1))+eps)); wzp=0.5_WP*[      -w,1.0_WP+2.0_WP*w,1.0_WP-w]
                     w=weno_weight((abs(pRHOIG(i,j,k+1,1)-pRHOIG(i,j,k  ,1))+eps)/(abs(pRHOIG(i,j,k,1)-pRHOIG(i,j,k-1,1))+eps)); wzm=0.5_WP*[1.0_WP-w,1.0_WP+2.0_WP*w,      -w]
                     pFz(i,j,k,4)=-0.5_WP*(pW(i,j,k,1)+abs(pW(i,j,k,1)))*sum(wzp*pRHOIG(i,j,k-2:k  ,1)) &
                     &            -0.5_WP*(pW(i,j,k,1)-abs(pW(i,j,k,1)))*sum(wzm*pRHOIG(i,j,k-1:k+1,1))
                  end if
                  if (k.gt.lbound(pTL,3).and.k.le.ubound(pTL,3)) then
                     condL=0.0_WP; if (all(pVF(i,j,k-1:k,1).ge.VFlo).and.all(pTL(i,j,k-1:k,1).gt.0.0_WP)) condL=(       minval(pVF(i,j,k-1:k,1)))*0.5_WP*sum(pDiffL(i,j,k-1:k,1))*dzi
                     condG=0.0_WP; if (all(pVF(i,j,k-1:k,1).le.VFhi).and.all(pTG(i,j,k-1:k,1).gt.0.0_WP)) condG=(1.0_WP-maxval(pVF(i,j,k-1:k,1)))*0.5_WP*sum(pDiffG(i,j,k-1:k,1))*dzi
                     pFz(i,j,k,3)=pFz(i,j,k,3)+limit_cond(condL*(pTL(i,j,k,1)-pTL(i,j,k-1,1)),pQold(i,j,k,3),pQold(i,j,k-1,3),dzi)
                     pFz(i,j,k,4)=pFz(i,j,k,4)+limit_cond(condG*(pTG(i,j,k,1)-pTG(i,j,k-1,1)),pQold(i,j,k,4),pQold(i,j,k-1,4),dzi)
                  end if
               end do; end do; end do
               ! Edge shear dissipation storage: [lo,hi+1] along the two nodal directions, [lo,hi] along the third
               allocate(dxy(bx%lo(1):bx%hi(1)+1,bx%lo(2):bx%hi(2)+1,bx%lo(3):bx%hi(3)  ))
               allocate(dyz(bx%lo(1):bx%hi(1)  ,bx%lo(2):bx%hi(2)+1,bx%lo(3):bx%hi(3)+1))
               allocate(dzx(bx%lo(1):bx%hi(1)+1,bx%lo(2):bx%hi(2)  ,bx%lo(3):bx%hi(3)+1))
               ! xy-edge (FUy, FVx): convection (centered, or band interpolant) and shear stress
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2)+1; do i=bx%lo(1),bx%hi(1)+1
                  strain=(pU(i,j,k,1)-pU(i,j-1,k,1))*dyi+(pV(i,j,k,1)-pV(i-1,j,k,1))*dxi
                  stress=0.25_WP*sum(pVisc(i-1:i,j-1:j,k,1))*strain
                  dxy(i,j,k)=stress*strain
                  pFUy(i,j,k,1)=conv_flux(1,2,[i,j,k],0.25_WP*sum(pFy(i-1:i,j,k,1:2))*sum(pU(i,j-1:j,k,1)))+stress
                  pFVx(i,j,k,1)=conv_flux(2,1,[i,j,k],0.25_WP*sum(pFx(i,j-1:j,k,1:2))*sum(pV(i-1:i,j,k,1)))+stress
               end do; end do; end do
               ! yz-edge (FVz, FWy): convection (centered, or band interpolant) and shear stress
               do k=bx%lo(3),bx%hi(3)+1; do j=bx%lo(2),bx%hi(2)+1; do i=bx%lo(1),bx%hi(1)
                  strain=(pV(i,j,k,1)-pV(i,j,k-1,1))*dzi+(pW(i,j,k,1)-pW(i,j-1,k,1))*dyi
                  stress=0.25_WP*sum(pVisc(i,j-1:j,k-1:k,1))*strain
                  dyz(i,j,k)=stress*strain
                  pFVz(i,j,k,1)=conv_flux(2,3,[i,j,k],0.25_WP*sum(pFz(i,j-1:j,k,1:2))*sum(pV(i,j,k-1:k,1)))+stress
                  pFWy(i,j,k,1)=conv_flux(3,2,[i,j,k],0.25_WP*sum(pFy(i,j,k-1:k,1:2))*sum(pW(i,j-1:j,k,1)))+stress
               end do; end do; end do
               ! zx-edge (FWx, FUz): convection (centered, or band interpolant) and shear stress
               do k=bx%lo(3),bx%hi(3)+1; do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)+1
                  strain=(pW(i,j,k,1)-pW(i-1,j,k,1))*dxi+(pU(i,j,k,1)-pU(i,j,k-1,1))*dzi
                  stress=0.25_WP*sum(pVisc(i-1:i,j,k-1:k,1))*strain
                  dzx(i,j,k)=stress*strain
                  pFWx(i,j,k,1)=conv_flux(3,1,[i,j,k],0.25_WP*sum(pFx(i,j,k-1:k,1:2))*sum(pW(i-1:i,j,k,1)))+stress
                  pFUz(i,j,k,1)=conv_flux(1,3,[i,j,k],0.25_WP*sum(pFz(i-1:i,j,k,1:2))*sum(pU(i,j,k-1:k,1)))+stress
               end do; end do; end do
               ! Cell centers: momentum fluxes from convection (centered, or band interpolant), normal viscous stress, and
               ! pressure over [lo-1,hi+1]; in the valid cells [lo,hi], viscous dissipation (normal, and the four
               ! surrounding edges of each orientation) split by VF between phasic energies, and phasic pressure dilatation
               do k=bx%lo(3)-1,bx%hi(3)+1; do j=bx%lo(2)-1,bx%hi(2)+1; do i=bx%lo(1)-1,bx%hi(1)+1
                  gradU=[(pU(i+1,j,k,1)-pU(i,j,k,1))*dxi,(pV(i,j+1,k,1)-pV(i,j,k,1))*dyi,(pW(i,j,k+1,1)-pW(i,j,k,1))*dzi]
                  div=sum(gradU)
                  mu=pVisc(i,j,k,1); lam=pBeta(i,j,k,1)-2.0_WP*mu/3.0_WP
                  tii=2.0_WP*mu*gradU+lam*div
                  pres=pVF(i,j,k,1)*pPL(i,j,k,1)+(1.0_WP-pVF(i,j,k,1))*pPG(i,j,k,1)
                  pFUx(i,j,k,1)=conv_flux(1,1,[i,j,k],0.25_WP*sum(pFx(i:i+1,j,k,1:2))*sum(pU(i:i+1,j,k,1)))+tii(1)-pres
                  pFVy(i,j,k,1)=conv_flux(2,2,[i,j,k],0.25_WP*sum(pFy(i,j:j+1,k,1:2))*sum(pV(i,j:j+1,k,1)))+tii(2)-pres
                  pFWz(i,j,k,1)=conv_flux(3,3,[i,j,k],0.25_WP*sum(pFz(i,j,k:k+1,1:2))*sum(pW(i,j,k:k+1,1)))+tii(3)-pres
                  if (any([i,j,k].lt.bx%lo).or.any([i,j,k].gt.bx%hi)) cycle
                  heat=sum(tii*gradU)+0.25_WP*(sum(dxy(i:i+1,j:j+1,k))+sum(dyz(i,j:j+1,k:k+1))+sum(dzx(i:i+1,j,k:k+1)))
                  pdQdt(i,j,k,3)=pdQdt(i,j,k,3)+(       pVF(i,j,k,1))*(heat-pPL(i,j,k,1)*div)
                  pdQdt(i,j,k,4)=pdQdt(i,j,k,4)+(1.0_WP-pVF(i,j,k,1))*(heat-pPG(i,j,k,1)*div)
               end do; end do; end do
               ! Free edge dissipation storage
               deallocate(dxy,dyz,dzx)
            end do
            call this%amr%mfiter_destroy(mfi)
         end do
         ! Reduce CTU face counts
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%nctuL,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%nctuG,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
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
               ! U RHS at x-faces
               bx=mfi%nodaltilebox(1)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pdrhoUdt(i,j,k,1)=dxi*(pFUx(i,j,k,1)-pFUx(i-1,j,k,1))+dyi*(pFUy(i,j+1,k,1)-pFUy(i,j,k,1))+dzi*(pFUz(i,j,k+1,1)-pFUz(i,j,k,1))
               end do; end do; end do
               ! V RHS at y-faces
               bx=mfi%nodaltilebox(2)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pdrhoVdt(i,j,k,1)=dxi*(pFVx(i+1,j,k,1)-pFVx(i,j,k,1))+dyi*(pFVy(i,j,k,1)-pFVy(i,j-1,k,1))+dzi*(pFVz(i,j,k+1,1)-pFVz(i,j,k,1))
               end do; end do; end do
               ! W RHS at z-faces
               bx=mfi%nodaltilebox(3)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pdrhoWdt(i,j,k,1)=dxi*(pFWx(i+1,j,k,1)-pFWx(i,j,k,1))+dyi*(pFWy(i,j+1,k,1)-pFWy(i,j,k,1))+dzi*(pFWz(i,j,k,1)-pFWz(i,j,k-1,1))
               end do; end do; end do
               ! Phasic mass and energy RHS at cell centers: add flux divergence to sources from compute_fluxes
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

      !> Content-limited phasic conduction flux through a face of inverse spacing hi (as amrmpcomp_col), from the unlimited
      !> flux (positive when the upper cell is hotter) and the start-of-step phasic energies Eup and Elo of the upper and
      !> lower cells: one face may drain at most 10% of the hotter cell's phasic energy over the step (Q=Qold+dt*dQdt),
      !> so 2*ndim faces never exceed ~60%. Extraction otherwise scales with dT*aperture, not content, and a thin-mass hot
      !> phase can lose more than it holds.
      real(WP) function limit_cond(flux,Eup,Elo,hi) result(lflux)
         implicit none
         real(WP), intent(in) :: flux,Eup,Elo,hi
         if (flux.gt.0.0_WP) then
            lflux=min(flux,+0.1_WP*max(Eup,0.0_WP)/(hi*dt))
         else
            lflux=max(flux,-0.1_WP*max(Elo,0.0_WP)/(hi*dt))
         end if
      end function limit_cond

      !> Whether b-face n of the finest level touches the semi-Lagrangian band
      logical function band_face(b,n)
         implicit none
         integer, intent(in) :: b
         integer, dimension(3), intent(in) :: n
         integer, dimension(3) :: m
         m=n-eu(:,b)
         band_face=max(pBand(m(1),m(2),m(3),1),pBand(n(1),n(2),n(3),1)).gt.0.0_WP
      end function band_face

      !> Dual face normal to b at index n for momentum axis a (the cell center if a=b, the edge otherwise):
      !> it is made of halves of the b-faces hlo and hhi, and separates the staggered CVs of the a-faces cv1 and cv2
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

      !> Convective momentum flux of axis a through the dual face normal to b at index n: the centered flux
      !> passed in, unless a half-face of the dual face is a band face (dual_flux)
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
         if (sllo.or.slhi) flux=dual_flux(a,b,n,sllo,slhi,owned_dual(n))
      end function conv_flux

      !> Band momentum flux of axis a through the dual face normal to b at index n, in the negative flux
      !> convention; sllo and slhi tell whether its half-faces hlo and hhi are band faces. Per phase and
      !> half-face h, the flux is the phasic mass flux of h (eq. 46) times the flux_interp velocity: Lax-Wendroff
      !> (eq. 51) or the centered current velocity. Where h is a band face and the modified switch selects CTU
      !> (eq. 55), it is the donor moment of h instead (eqs. 53-54), whose donors are mapped to the staggered
      !> donors with (B.9): the upper a-face of the donor for hlo, the lower a-face for hhi. CTU uses are counted
      !> when cnt is true.
      real(WP) function dual_flux(a,b,n,sllo,slhi,cnt) result(flux)
         implicit none
         integer, intent(in) :: a,b
         integer, dimension(3), intent(in) :: n
         logical, intent(in) :: sllo,slhi,cnt
         integer, dimension(3) :: hlo,hhi,cv1,cv2
         integer :: phase
         real(WP) :: phi
         logical :: use_ctu
         ! Contributing b-faces and separated staggered CVs
         call dual_face(a,b,n,hlo,hhi,cv1,cv2)
         ! Velocity transported by the phasic mass fluxes where CTU is not used
         if (this%flux_interp.eq.flux_lw) then
            phi=lw_velocity(a,b,n,cv1,cv2)
         else
            phi=0.5_WP*(vel(a,cv1)+vel(a,cv2))
         end if
         ! Per-phase flux_interp or CTU flux
         flux=0.0_WP
         do phase=1,2
            use_ctu=min(alpha(phase,a,cv1),alpha(phase,a,cv2)).lt.this%alpha_ctu
            ! Lower half-face: donor moment with the velocity of the donor's upper a-face
            if (use_ctu.and.sllo) then
               flux=flux+0.5_WP*flx(b,hlo,ictu(phase,2,a))
            else
               flux=flux+0.5_WP*mflx(b,hlo,phase)*phi
            end if
            ! Upper half-face: donor moment with the velocity of the donor's lower a-face
            if (use_ctu.and.slhi) then
               flux=flux+0.5_WP*flx(b,hhi,ictu(phase,1,a))
            else
               flux=flux+0.5_WP*mflx(b,hhi,phase)*phi
            end if
            ! Count CTU uses on band half-faces carrying this phase
            if (cnt.and.use_ctu) then
               if ((sllo.and.mflx(b,hlo,phase).ne.0.0_WP).or.(slhi.and.mflx(b,hhi,phase).ne.0.0_WP)) then
                  if (phase.eq.1) this%nctuL=this%nctuL+1.0_WP
                  if (phase.eq.2) this%nctuG=this%nctuG+1.0_WP
               end if
            end if
         end do
      end function dual_flux

      !> Lax-Wendroff value (eq. 51) of the time-n velocity component a at the dual face normal to b at index n,
      !> which separates the staggered CVs of a-faces cv1 (lower along b) and cv2 (upper), paper section 4.2.1.
      !> With the interpolated time-n velocity u_g at the dual face centroid x_g, the LW reconstruction centroid is
      !> x_g-dt/2*u_g and the face CFL number is nu=dt*|u_g,b|/dx_b. U0 is the time-n velocity of the downwind CV,
      !> and U1 the bilinear interpolation of the time-n velocity on the plane of the upwind CV, where the line
      !> from the downwind node through the donor-region centroid crosses it (fig. 11).
      real(WP) function lw_velocity(a,b,n,cv1,cv2) result(phi)
         implicit none
         integer, intent(in) :: a,b
         integer, dimension(3), intent(in) :: n,cv1,cv2
         integer, dimension(3) :: fdown,fup
         real(WP), dimension(3) :: h,xg,ug,off
         real(WP) :: nu,U0,U1
         ! Dual face centroid: cell center if a=b, otherwise nodal along a and b
         h=[dx,dy,dz]
         xg=[this%amr%xlo,this%amr%ylo,this%amr%zlo]+(real(n,WP)+0.5_WP)*h
         if (a.ne.b) then
            xg(a)=xg(a)-0.5_WP*h(a)
            xg(b)=xg(b)-0.5_WP*h(b)
         end if
         ! LW uses time-n velocity for its direction, CFL, and reconstruction offset. The geometric sweep
         ! and VF/barycenter projection continue to use the current midpoint velocity.
         ug=interp_velocity(xg,old_velocity=.true.)
         nu=dt*abs(ug(b))/h(b)
         if (ug(b).ge.0.0_WP) then
            fdown=cv2; fup=cv1
         else
            fdown=cv1; fup=cv2
         end if
         ! Transverse offset of the interpolation line on the upwind plane: -dt*u_g/(1+nu)
         off=-dt*ug/(1.0_WP+nu)
         ! Downwind and upwind values of the time-n velocity, combined with eq. (51)
         U0=velold(a,fdown)
         U1=plane_value(a,b,fup,off)
         phi=U1+0.5_WP*(1.0_WP-nu)*(U0-U1)
      end function lw_velocity
   
      !> Bilinear interpolation of the time-n velocity component a on the plane normal to b through the a-face f,
      !> at the transverse offset off from f (the b component of off is ignored)
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

      !> Whether dual face index n is counted by the current tile: cells and edges indexed within the
      !> tile box, so that shared edges are counted once by the tile owning their lower index
      logical function owned_dual(n)
         implicit none
         integer, dimension(3), intent(in) :: n
         owned_dual=all(n.ge.tbx%lo).and.all(n.le.tbx%hi)
      end function owned_dual

      !> Band face flux component comp at b-face n
      real(WP) function flx(b,n,comp)
         implicit none
         integer, intent(in) :: b,comp
         integer, dimension(3), intent(in) :: n
         flx=pF(b)%p(n(1),n(2),n(3),comp)
      end function flx

      !> Phasic mass or internal energy flux component comp at b-face n (band or bulk)
      real(WP) function mflx(b,n,comp)
         implicit none
         integer, intent(in) :: b,comp
         integer, dimension(3), intent(in) :: n
         mflx=pM(b)%p(n(1),n(2),n(3),comp)
      end function mflx

      !> Current face velocity of axis a at a-face n
      real(WP) function vel(a,n)
         implicit none
         integer, intent(in) :: a
         integer, dimension(3), intent(in) :: n
         vel=pVel(a)%p(n(1),n(2),n(3),1)
      end function vel

      !> Time-n face velocity of axis a at a-face n
      real(WP) function velold(a,n)
         implicit none
         integer, intent(in) :: a
         integer, dimension(3), intent(in) :: n
         velold=pVelold(a)%p(n(1),n(2),n(3),1)
      end function velold

      !> Post-transport phase fraction of the staggered CV of a-face n: arithmetic mean of its two cells, eq. (37)
      real(WP) function alpha(phase,a,n)
         implicit none
         integer, intent(in) :: phase,a
         integer, dimension(3), intent(in) :: n
         integer, dimension(3) :: m
         m=n-eu(:,a)
         alpha=0.5_WP*(pVF(m(1),m(2),m(3),1)+pVF(n(1),n(2),n(3),1))
         if (phase.eq.2) alpha=1.0_WP-alpha
      end function alpha

      !> Recursive subroutine that cuts a tet by grid planes to compute volume and Q fluxes (as amrmpcomp_col)
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

      !> Cut a tet inside donor cell (i0,j0,k0) by its old PLIC: phase volumes and first moments, phasic Q and
      !> pressure-volume from the donor's old state (as amrmpcomp_col and mpcomp), and CTU donor moments
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
            myQflux(1:4)=vol_tot*pQold(i0,j0,k0,1:4)
            myQflux( 5 )=vol_tot*pPLold(i0,j0,k0,1)
            call donor_moments(i0,j0,k0,myQflux)
            return
         else if (pPLICold(i0,j0,k0,4).lt.-1.0e9_WP) then
            ! Pure gas
            myVflux( 2 )=vol_tot
            myVflux(6:8)=vol_tot*bary_tot
            myQflux(1:4)=vol_tot*pQold(i0,j0,k0,1:4)
            myQflux( 6 )=vol_tot*pPGold(i0,j0,k0,1)
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

         ! Phasic Q and pressure-volume from the donor's old phasic state
         if (VF0.ge.VFlo) then
            myQflux(1)=myVflux(1)*pQold(i0,j0,k0,1)/VF0
            myQflux(3)=myVflux(1)*pQold(i0,j0,k0,3)/VF0
            myQflux(5)=myVflux(1)*pPLold(i0,j0,k0,1)
         end if
         if (VF0.le.VFhi) then
            myQflux(2)=myVflux(2)*pQold(i0,j0,k0,2)/(1.0_WP-VF0)
            myQflux(4)=myVflux(2)*pQold(i0,j0,k0,4)/(1.0_WP-VF0)
            myQflux(6)=myVflux(2)*pPGold(i0,j0,k0,1)
         end if
         call donor_moments(i0,j0,k0,myQflux)

      end subroutine tet2flux_plic

      !> CTU donor moments of the phasic masses taken from donor b=(i0,j0,k0), eqs. (53) and (B.9): for
      !> each axis a, the phasic mass times the time-n velocity of the lower and upper a-faces of b
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

      !> RK2 vertex projection with the current face velocity (as amrmpcomp_col)
      function project(p1,mydt) result(p2)
         implicit none
         real(WP), dimension(3), intent(in) :: p1
         real(WP), dimension(3)             :: p2
         real(WP),               intent(in) :: mydt
         p2=p1+mydt*interp_velocity(        p1    )
         p2=p1+mydt*interp_velocity(0.5_WP*(p1+p2))
      end function project

      !> Trilinear interpolation of current staggered velocity, or time-n velocity for LW reconstruction
      function interp_velocity(pos,old_velocity) result(vel3)
         implicit none
         real(WP), dimension(3), intent(in) :: pos
         logical, intent(in), optional :: old_velocity
         real(WP), dimension(3) :: vel3
         integer  :: ipc,jpc,kpc   ! Cell-centered indices
         integer  :: ipu,jpv,kpw   ! Face-centered indices
         real(WP) :: wxc1,wyc1,wzc1,wxc2,wyc2,wzc2  ! Cell-centered weights
         real(WP) :: wxu1,wyv1,wzw1,wxu2,wyv2,wzw2  ! Face-centered weights
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW
         pU=>pVel(1)%p; pV=>pVel(2)%p; pW=>pVel(3)%p
         if (present(old_velocity)) then
            if (old_velocity) then
               pU=>pVelold(1)%p; pV=>pVelold(2)%p; pW=>pVelold(3)%p
            end if
         end if
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

   !> Apply the relaxation model in mixed cells of the finest level, then re-match the PLIC
   !> plane to the relaxed VF and recompute barycenters (as amrmpcomp_col). Q(1:2) and Q(3)+Q(4)
   !> are conserved by the model; face momenta are untouched. No surface tension (Pjump=0).
   subroutine apply_relax(this,dt,time)
      use amrvof_geometry, only: get_plane_dist,cut_hex_vol
      use parallel,        only: MPI_REAL_WP
      use mpi_f08,         only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM,MPI_Wtime
      implicit none
      class(amrmpcomp), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), intent(in) :: time
      real(WP) :: tw0
      integer :: lvl,i,j,k,ierr,rerr
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pQ,pCL,pCG,pPLIC
      logical :: oldmix,newmix
      real(WP) :: dx,dy,dz,cell_vol,vol_liq,vol_gas
      real(WP), dimension(3) :: lo,hi,bary_liq,bary_gas
      real(WP), dimension(3,8) :: hex
      real(WP), dimension(4) :: plane
      ! Reset failure counter
      this%nrelax_fail=0.0_WP
      ! If no relaxation model was provided, return
      if (.not.associated(this%relax)) return
      ! Return if clvl<maxlvl
      if (this%amr%clvl().lt.this%amr%maxlvl) return
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
         ! Loop over valid cells
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Check if mixture cell prior to relaxation
            oldmix=(pVF(i,j,k,1).ge.VFlo.and.pVF(i,j,k,1).le.VFhi)
            ! Apply user-provided relaxation model (modifies VF and Q)
            call this%relax%apply(dt=dt,VF=pVF(i,j,k,1),Q=pQ(i,j,k,:),Pjump=0.0_WP,ierr=rerr)
            if (oldmix.and.rerr.ne.RELAX_OK) this%nrelax_fail=this%nrelax_fail+1.0_WP
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
                  pQ(i,j,k,3)=0.0_WP
                  pPLIC(i,j,k,:)=[0.0_WP,0.0_WP,0.0_WP,-1.0e10_WP]
               else if (pVF(i,j,k,1).gt.VFhi) then
                  ! Pure liquid
                  pVF(i,j,k,1)=1.0_WP
                  pCL(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
                  pCG(i,j,k,1:3)=[this%amr%xlo+(real(i,WP)+0.5_WP)*dx,this%amr%ylo+(real(j,WP)+0.5_WP)*dy,this%amr%zlo+(real(k,WP)+0.5_WP)*dz]
                  pQ(i,j,k,2)=0.0_WP
                  pQ(i,j,k,4)=0.0_WP
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
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%nrelax_fail,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
      ! Sync and apply BC
      call this%VF%average_down(); call this%fill(lvl=this%amr%maxlvl,time=time)
      call this%Q%average_down(); call this%Q%fill(time=time)
      this%wt_relax=this%wt_relax+(MPI_Wtime()-tw0)
   end subroutine apply_relax

   !> Add artificial bulk viscosity to this%beta (and optionally this%visc). The kinematic
   !> viscosity is converted with the harmonic phasic density, as in amrmpcomp_col.
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

   !> Add Vreman SGS eddy viscosity to this%visc from the staggered velocity: assumes face velocity ghosts are filled.
   !> User must reset visc to molecular value before calling this routine.
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

   !> Compute CFL numbers (convective, acoustic, and viscous/heat diffusion)
   subroutine get_cfl(this,dt,cfl)
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
               ! Heat-diffusion CFL: thermal diffusivity alpha=lambda/(rho*cv), phasic in pure cells (as amrmpcomp_col)
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
      ! Return max CFL
      cfl=max(this%CFLc_x,this%CFLc_y,this%CFLc_z,this%CFLv_x,this%CFLv_y,this%CFLv_z)
      cfl=max(cfl,this%CFLa_x,this%CFLa_y,this%CFLa_z)
   end subroutine get_cfl

   ! ============================================================================
   ! SOLVER INFO
   ! ============================================================================

   !> Get solver information: phasic extrema, mixture pressure and sound speed, momentum, kinetic energy
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
      call log("Staggered compressible multiphase solver: "//trim(this%name))
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
      ! VOF data, face velocities, and Q via parent
      call this%amrmpflow%restore_checkpoint(io,dirname,time)
      ! Face momenta
      call io%read_data(dirname,this%rhoU,'rhoU')
      call io%read_data(dirname,this%rhoV,'rhoV')
      call io%read_data(dirname,this%rhoW,'rhoW')
      ! Rebuild primitive variables
      call this%get_primitive(this%Q,this%rhoU,this%rhoV,this%rhoW,time)
   end subroutine restore_checkpoint

end module amrmpcomp_class
