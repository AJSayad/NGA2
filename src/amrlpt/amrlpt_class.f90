!> AMR-aware Lagrangian particle tracking solver
!> Mirrors lpt_class capabilities on an AMReX AMR hierarchy.
module amrlpt_class
   use material_class, only: material
   use precision,     only: WP,I8
   use string,        only: str_medium
   use amrgrid_class, only: amrgrid
   use amrdata_class, only: amrdata
   use timer_class,   only: timer
   use iso_c_binding
   implicit none
   private

   ! Public exports
   public :: amrlpt,part,particle_finite

   ! Boundary condition flags for collision detection
   integer, parameter, public :: AMRLPT_OPEN=0  !< No wall collision (open or periodic boundary)
   integer, parameter, public :: AMRLPT_WALL=1  !< Hard wall collision at domain boundary

   ! Particle struct layout constants (must match #define in amrlpt_wrapper.cpp)
   integer, parameter, public :: AMRLPT_NREAL=16  !< extra reals per particle
   integer, parameter, public :: AMRLPT_NINT =1   !< extra ints  per particle

   ! Particle flags
   integer(c_int), parameter, public :: PART_MOVES    =1
   integer(c_int), parameter, public :: PART_COLLIDES =2
   integer(c_int), parameter, public :: PART_EXCHANGES=4
   integer(c_int), parameter, public :: PART_IS_DEAD  =0

   ! Particle struct -- must match C++ Particle<16,1> memory layout exactly:
   ! - pos[3]    (pos, managed by AMReX)
   ! - rdata[16] (d, vel[3], angVel[3], Acol[3], Tcol[3], dt, mass, eint)
   ! - idcpu     (packed id+cpu, private)
   ! - idata[1]  (flag)
   type, bind(C), public :: part
      !> AMReX position (physical coordinates)
      real(c_double) :: pos(3)
      !> Extra reals (NStructReal=16)
      real(c_double) :: d               !< Particle diameter
      real(c_double) :: vel(3)          !< Particle velocity
      real(c_double) :: angVel(3)       !< Angular velocity
      real(c_double) :: Acol(3)         !< Collision acceleration
      real(c_double) :: Tcol(3)         !< Collision torque
      real(c_double) :: dt              !< Particle sub-timestep
      real(c_double) :: mass            !< Constant particle mass
      real(c_double) :: eint            !< Specific internal energy; zero for nonthermal callers
      !> Packed id (39 bits, sign=valid) + cpu (24 bits)  [AMReX internal]
      integer(c_int64_t), private :: idcpu
      !> Extra int (NStructInt=1)
      integer(c_int) :: flag            !< 0=active, 1=mark for removal
   end type part

   !> C interface: bind(C) declarations to amrlpt_wrapper.cpp
   interface
      function amrlpt_particle_size() result(nbytes) bind(c)
         import :: c_size_t
         integer(c_size_t) :: nbytes
      end function
      subroutine amrlpt_new_pc(pc, amrcore) bind(c)
         import :: c_ptr
         type(c_ptr) :: pc
         type(c_ptr), value :: amrcore
      end subroutine

      subroutine amrlpt_delete_pc(pc) bind(c)
         import :: c_ptr
         type(c_ptr), value :: pc
      end subroutine

      subroutine amrlpt_redistribute(pc,lev_min,lev_max,ng) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: lev_min,lev_max,ng
      end subroutine

      subroutine amrlpt_enforce_max_level(pc,lev_max,number_demoted) bind(c)
         import :: c_ptr,c_int,c_int64_t
         type(c_ptr), value :: pc
         integer(c_int), value :: lev_max
         integer(c_int64_t) :: number_demoted
      end subroutine

      subroutine amrlpt_level_counts(pc,counts,nlevels) bind(c)
         import :: c_ptr,c_int,c_int64_t
         type(c_ptr), value :: pc
         integer(c_int64_t), dimension(*) :: counts
         integer(c_int), value :: nlevels
      end subroutine

      subroutine amrlpt_fill_neighbors(pc,ngrow) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: ngrow
      end subroutine

      subroutine amrlpt_fill_neighbors_radius(pc,search_radius) bind(c)
         import :: c_ptr,c_double
         type(c_ptr), value :: pc
         real(c_double), value :: search_radius
      end subroutine

      subroutine amrlpt_clear_neighbors(pc) bind(c)
         import :: c_ptr
         type(c_ptr), value :: pc
      end subroutine

      subroutine amrlpt_get_neighbor_particles_mfi(pc,lev,mfi,dp,np) bind(c)
         import :: c_ptr,c_int,c_int64_t
         type(c_ptr), value :: pc,mfi
         integer(c_int), value :: lev
         type(c_ptr) :: dp
         integer(c_int64_t) :: np
      end subroutine

      subroutine amrlpt_build_neighbor_list(pc,rcrit) bind(c)
         import :: c_ptr,c_double
         type(c_ptr), value :: pc
         real(c_double), value :: rcrit
      end subroutine

      subroutine amrlpt_get_neighbor_list_mfi(pc,lev,mfi,offsets,list,np,ntot) bind(c)
         import :: c_ptr,c_int,c_int64_t
         type(c_ptr), value :: pc,mfi
         integer(c_int), value :: lev
         type(c_ptr) :: offsets
         type(c_ptr) :: list
         integer(c_int64_t) :: np
         integer(c_int64_t) :: ntot
      end subroutine

      subroutine amrlpt_get_particles_mfi(pc,lev,mfi,dp,np) bind(c)
         import :: c_ptr,c_int,c_int64_t
         type(c_ptr), value :: pc,mfi
         integer(c_int), value :: lev
         type(c_ptr) :: dp
         integer(c_int64_t) :: np
      end subroutine

      subroutine amrlpt_get_all_particles_mfi(pc,lev,mfi,dp,np_total,np_valid) bind(c)
         import :: c_ptr,c_int,c_int64_t
         type(c_ptr), value :: pc,mfi
         integer(c_int), value :: lev
         type(c_ptr) :: dp
         integer(c_int64_t) :: np_total
         integer(c_int64_t) :: np_valid
      end subroutine

      subroutine amrlpt_add_particle_i(pc,lev,grid,tile,p) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc,p
         integer(c_int), value :: lev,grid,tile
      end subroutine

      subroutine amrlpt_get_next_id(id) bind(c)
         import :: c_int64_t
         integer(c_int64_t) :: id
      end subroutine

      subroutine amrlpt_set_next_id(id) bind(c)
         import :: c_int64_t
         integer(c_int64_t), value :: id
      end subroutine

      subroutine amrlpt_get_cpu(cpu) bind(c)
         import :: c_int
         integer(c_int) :: cpu
      end subroutine

      subroutine amrlpt_set_particle_id(id,p) bind(c)
         import :: c_int64_t,c_ptr
         integer(c_int64_t), value :: id
         type(c_ptr), value :: p
      end subroutine

      subroutine amrlpt_set_particle_cpu(cpu,p) bind(c)
         import :: c_int,c_ptr
         integer(c_int), value :: cpu
         type(c_ptr), value :: p
      end subroutine

      subroutine amrlpt_particle_is_valid(valid,p) bind(c)
         import :: c_int,c_ptr
         integer(c_int) :: valid
         type(c_ptr), value :: p
      end subroutine

      subroutine amrlpt_total_np(pc,np) bind(c)
         import :: c_ptr,c_int64_t
         type(c_ptr), value :: pc
         integer(c_int64_t) :: np
      end subroutine

      subroutine amrlpt_write(pc,path,is_chk) bind(c)
         import :: c_ptr,c_char,c_int
         type(c_ptr), value :: pc
         character(kind=c_char), dimension(*) :: path
         integer(c_int), value :: is_chk
      end subroutine

      subroutine amrlpt_read(pc,path) bind(c)
         import :: c_ptr,c_char
         type(c_ptr), value :: pc
         character(kind=c_char), dimension(*) :: path
      end subroutine

      subroutine amrlpt_append_particles(pc,raw,n) bind(c)
         import :: c_ptr,c_int64_t
         type(c_ptr), value :: pc
         type(c_ptr), value :: raw
         integer(c_int64_t), value :: n
      end subroutine

      subroutine amrlpt_get_particle_boxarray(pc,lev,ba) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: lev
         type(c_ptr) :: ba
      end subroutine

      subroutine amrlpt_get_particle_distromap(pc,lev,dm) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: lev
         type(c_ptr) :: dm
      end subroutine

      subroutine amrlpt_set_particle_boxarray(pc,lev,ba) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: lev
         type(c_ptr), value :: ba
      end subroutine

      subroutine amrlpt_set_particle_distromap(pc,lev,dm) bind(c)
         import :: c_ptr,c_int
         type(c_ptr), value :: pc
         integer(c_int), value :: lev
         type(c_ptr), value :: dm
      end subroutine

   end interface

   !> AMR LPT solver type
   type :: amrlpt

      !> Associated AMR grid
      class(amrgrid), pointer :: amr=>null()

      !> AMReX NeighborParticleContainer<16,1> opaque handle
      type(c_ptr) :: pc=c_null_ptr

      !> Solver name
      character(len=str_medium) :: name='UNNAMED_AMRLPT'

      !> User-provided injection callback
      procedure(lpt_inject_iface ), pointer, pass :: inject=>null()
      procedure(lpt_tagging_iface), pointer, pass :: user_lpt_tagging=>null()

      !> Global particle count
      integer(I8) :: np=0
      integer(I8) :: np_demoted=0 !< Particles moved below maxlvl by the last explicit enforcement
      integer :: occupied_minlvl=-1,occupied_maxlvl=-1
      integer(I8) :: np_above_maxlvl=0

      !> Physics parameters
      real(WP) :: rho                                   !< Reference density for injection and immersed carrier fallback
      real(WP) :: gas_gamma=1.4_WP                     !< Gas specific-heat ratio used by compressible drag laws
      real(WP), dimension(3) :: gravity=0.0_WP          !< Acceleration of gravity
      integer :: nstep=1                                !< Substeps per timestep
      character(len=str_medium) :: drag_model='Tenneti' !< Drag model
      !> Optional prescribed solids fraction used only by the drag closure. A negative value keeps
      !> the normal locally deposited particle VF. This is useful for isolated-particle/dilute-limit
      !> verification without changing VF deposition or two-way source transfer.
      real(WP) :: drag_vf_override=-1.0_WP
      !> Select total-energy reaction; default retains legacy drag-dissipation source semantics.
      logical :: total_energy_coupling=.false.
      !> Associated EOS enables heat exchange. Diameter/mass are fixed, so thermal capacity is m*cv.
      class(material), pointer :: heat_material=>null()
      real(WP) :: heat_prandtl=0.71_WP
      !> Cumulative heat returned to fluid (negative for particle heating).
      real(WP) :: heat_exchange=0.0_WP,heat_exchange_loc=0.0_WP
      !> Cumulative fluid impulse(3), mechanical energy, and drag-only dissipation.
      real(WP) :: exchange(5)=0.0_WP,exchange_loc(5)=0.0_WP

      !> Collision parameters
      real(WP) :: tau_col=1.0e-3_WP             !< Collision time scale
      real(WP) :: e_n=1.0_WP                    !< Normal restitution
      real(WP) :: e_w=1.0_WP                    !< Wall restitution
      real(WP) :: mu_f=0.0_WP                   !< Friction coefficient
      real(WP) :: clip_col=0.2_WP               !< Max overlap fraction

      !> Drag-regime diagnostics: max particle Reynolds, Mach and Knudsen numbers seen during the
      !> last advance. Worth watching - every drag model here has a limited range of validity.
      real(WP) :: Remax=0.0_WP,Mamax=0.0_WP,Knmax=0.0_WP
      !> Optional co-located records for the global Osnes Re/Ma/Kn extrema, disabled by default.
      logical :: drag_diagnostics=.false.
      logical :: drag_diagnostic_pending=.false.
      real(WP), dimension(30,3) :: drag_diagnostic_record=0.0_WP

      !> CFL numbers
      real(WP) :: CFLp_x=0.0_WP,CFLp_y=0.0_WP,CFLp_z=0.0_WP
      real(WP) :: CFL_col=0.0_WP

      !> Monitoring data
      real(WP) :: dmin,dmax,dmean,dvar
      real(WP) :: Umin,Umax,Umean,Uvar
      real(WP) :: Vmin,Vmax,Vmean,Vvar
      real(WP) :: Wmin,Wmax,Wmean,Wvar
      real(WP) :: VFmin,VFmax,VFmean,VFvar
      integer  :: np_new=0,np_out=0
      real(WP) :: Vp_new=0.0_WP,Vp_out=0.0_WP,Vp_tot=0.0_WP
      integer  :: np_new_loc=0,np_out_loc=0
      real(WP) :: Vp_new_loc=0.0_WP,Vp_out_loc=0.0_WP
      !> Physical particle inventory, summed over active particles by get_info.
      real(WP) :: mass_tot=0.0_WP,kinetic_tot=0.0_WP,internal_tot=0.0_WP
      real(WP), dimension(3) :: momentum_tot=0.0_WP
      !> Cumulative open-boundary losses: mass, momentum(3), total energy.
      real(WP), dimension(5) :: outflow=0.0_WP,outflow_loc=0.0_WP
      !> Case-selected gas-free removal; legacy callers retain immersed-liquid fallback.
      logical :: cull_gas_free=.false.
      integer :: np_culled=0,np_culled_loc=0
      real(WP) :: np_culled_tot=0.0_WP ! real storage for live scalar checkpoint registration
      real(WP) :: culled(5)=0.0_WP,culled_loc(5)=0.0_WP ! mass, momentum(3), total energy
      integer  :: ncol=0
      
      !> Engulfed-particle fallback: a particle whose whole sampling stencil is pure liquid is
      !> treated as a neutrally-buoyant tracer in the liquid, using liquid properties in a
      !> Schiller-Naumann law. Set visc_immersed to the liquid dynamic viscosity (required).
      real(WP) :: visc_immersed=-1.0_WP
      integer  :: nsub_max=1000              !< Cap on per-particle sub-steps within one step
      !> Counters. np_immersed is PARTICLES treated in the most recent advance; np_immersed_tot
      !> accumulates particle-advances over the run. E_immersed is the drag dissipation that was
      !> NOT deposited, because the containing cell holds no gas to receive it.
      integer  :: np_immersed=0,       np_immersed_loc=0
      integer  :: np_immersed_tot=0
      integer  :: np_substep_capped=0, np_substep_capped_loc=0
      real(WP) :: E_immersed=0.0_WP,   E_immersed_loc=0.0_WP
      !> First engulfed event of the current advance on this rank, for the log:
      !> [lvl, pos(3), d, |vrel|, gas fraction min, gas fraction max]
      real(WP), dimension(8) :: immersed_sample=0.0_WP

      !> Particle distribution across ranks
      integer(I8) :: np_loc=0           !< This rank's particle count
      integer(I8) :: np_min=0           !< Min particle count across ranks
      integer(I8) :: np_max=0           !< Max particle count across ranks
      real(WP)    :: np_eff=0.0_WP      !< Load efficiency (mean/max)

      !> Overlap size
      integer :: nover=2

      !> Maximum AMR level particles are allowed to occupy. Particles use the
      !> fluid AMR hierarchy but are placed on the finest covering level in
      !> [0,maxlvl]. Fluid sampling, VF deposition and two-way source deposition
      !> consequently all occur on this capped particle level. Defaults to the
      !> fluid maximum; applications may lower it to satisfy point-particle scale
      !> separation. Explicit enforcement recovers particles placed above the cap
      !> by AMReX's unrestricted checkpoint Restart().
      integer :: maxlvl=0

      !> Filter width
      real(WP) :: filter_width=0.0_WP
      !> Cells across the filter width on the level the deposit is actually diffused on. The cost
      !> of the explicit filter is ~0.54*(filter_width/dx)^2 sub-steps, so resolving the Gaussian
      !> on the finest level is needlessly expensive; 4 cells is plenty for a smoothing kernel.
      integer  :: filter_ncell=4

      !> Particle volume fraction
      type(amrdata) :: VF

      !> VF threshold for tagging; disabled if <=0
      real(WP) :: VF_tag=-1.0_WP

      !> Two-way coupling source terms, 4 components. SIGN CONVENTION: src is what the FLUID
      !> GAINS - 1:3 react to drag + optional resolved stress
      type(amrdata) :: src

      !> Boundary condition flags for collision detection
      integer, dimension(3) :: lo_bc=AMRLPT_OPEN
      integer, dimension(3) :: hi_bc=AMRLPT_OPEN

      !> Particle sub-stepping state
      real(WP) :: t=0.0_WP            !< Current particle time
      real(WP) :: dt=huge(1.0_WP)     !< Persistent CFL-limited sub-step size
      real(WP) :: dtmax=huge(1.0_WP)  !< Maximum allowed particle sub-step
      real(WP) :: cflmax=huge(1.0_WP) !< CFL limit for particle sub-stepping

      !> Knapsack load balancing
      logical :: rebalance=.false.    !< Enable knapsack DM for particles

      !> Per-phase wall-clock timers
      type(timer) :: tmr_coll         !< Collision phase total
      type(timer) :: tmr_coll_        !< Collision phase w/o MPI
      type(timer) :: tmr_fill         !< Ghost fill time
      type(timer) :: tmr_nbl          !< Neighbor list build time
      type(timer) :: tmr_step         !< Particle stepping phase
      type(timer) :: tmr_step_        !< Particle stepping phase w/o MPI
      type(timer) :: tmr_vf           !< Volume fraction computation
      type(timer) :: tmr_src          !< Source post-processing

   contains
      ! Lifecycle callbacks
      procedure :: post_regrid            !< Post-regrid callback
      procedure :: tagging                !< Tag cells for refinement based on VF
      procedure :: get_cost               !< Compute per-box particle cost for load balancing
      ! Type-bound constructor/destructor
      procedure :: initialize
      procedure :: finalize
      ! Physics procedures
      procedure :: collide                !< Soft-sphere collision model
      procedure :: advance                !< Advance particle ODEs one timestep
      procedure :: advance_to             !< CFL-substepped advance to a target time
      procedure, private :: step          !< Single sub-step: MFIter loop + redistribute
      procedure :: update_VF              !< Compute particle volume fraction field
      procedure :: get_cfl                !< Compute particle CFL numbers
      ! Utilities
      procedure :: redistribute           !< Call AMReX redistribute
      procedure :: enforce_max_level      !< Recover particles already stored above maxlvl
      procedure :: get_level_counts       !< Global real-particle count on every AMR level
      procedure :: fill_ghosts            !< Fill ghost particle buffer
      procedure :: fill_ghosts_radius     !< Fill ghosts with position-based filtering
      procedure :: clear_ghosts           !< Release ghost particle buffer
      procedure :: build_neighbor_list    !< Build explicit pair list within rcrit
      procedure :: get_np                 !< Update global particle count
      procedure, private :: get_particles     !< Get valid particle array for MFIter tile
      procedure, private :: get_ghosts        !< Get ghost particle array for MFIter tile
      procedure, private :: get_all_particles !< Get combined valid+ghost array for MFIter tile
      procedure, private :: get_neighbor_list !< Get CSR neighbor list for MFIter tile
      procedure :: interp                 !< Trilinear cell-centered interpolation
      procedure :: interp_face_velocities !< Trilinear face-centered interpolation
      procedure, private :: process_deposit   !< Post-process deposit: intensive conversion + F↔C transfers
      procedure, private :: filter        !< Explicit diffusion filter
      procedure :: gather_region          !< Allgather particles within a bounding box
      procedure :: append                 !< Append a Fortran part array
      procedure :: mfiter_build           !< Build MFIter for particle grid
      procedure :: mfiter_destroy         !< Destroy MFIter
      procedure :: set_particle_ba        !< Set particle BoxArray for a level
      procedure :: set_particle_dm        !< Set particle DistributionMapping for a level
      procedure :: get_particle_dm        !< Get particle DistributionMapping for a level
      ! Print solver info
      procedure :: get_info
      procedure :: print
      ! Checkpoint I/O
      procedure :: read
      procedure :: write
   end type amrlpt

   !> Abstract interface for user-provided particle injection callback
   abstract interface
      subroutine lpt_inject_iface(this,dt)
         import :: amrlpt,WP
         class(amrlpt), intent(inout) :: this
         real(WP), intent(in) :: dt
      end subroutine lpt_inject_iface
   end interface

   !> Abstract interface for user-overridable tagging callback
   abstract interface
      subroutine lpt_tagging_iface(solver,lvl,time,tags)
         import :: amrlpt,c_ptr,WP
         class(amrlpt), intent(inout) :: solver
         integer, intent(in) :: lvl
         real(WP), intent(in) :: time
         type(c_ptr), intent(in) :: tags
      end subroutine lpt_tagging_iface
   end interface

contains

   ! IEEE binary64 exponent check remains valid under the production -ffast-math flags.
   elemental logical function particle_finite(x) result(finite)
      real(WP), intent(in) :: x
      finite=iand(transfer(x,0_I8),int(z'7ff0000000000000',I8)).ne.int(z'7ff0000000000000',I8)
   end function particle_finite

   ! ============================================================================
   ! DISPATCHERS (module-level) - recover concrete amrlpt type
   ! ============================================================================

   !> Dispatch post_regrid: calls type-bound method
   subroutine amrlpt_postregrid(ctx,lbase,time)
      use iso_c_binding, only: c_ptr,c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      type(amrlpt), pointer :: this
      call c_f_pointer(ctx,this)
      call this%post_regrid(lbase,time)
   end subroutine amrlpt_postregrid

   !> Dispatch tagging: calls type-bound method then user callback
   subroutine amrlpt_tagging(ctx,lvl,time,tags)
      use iso_c_binding, only: c_ptr,c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags
      type(amrlpt), pointer :: this
      call c_f_pointer(ctx,this)
      call this%tagging(lvl,time,tags)
      if (associated(this%user_lpt_tagging)) call this%user_lpt_tagging(lvl,time,tags)
   end subroutine amrlpt_tagging

   !> Dispatch get_cost: calls type-bound method
   subroutine amrlpt_get_cost(ctx,lvl,nboxes,costs,ba)
      use amrex_amr_module, only: amrex_boxarray
      use iso_c_binding, only: c_ptr,c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lvl,nboxes
      real(WP), intent(inout) :: costs(nboxes)
      type(amrex_boxarray), intent(in) :: ba
      type(amrlpt), pointer :: this
      call c_f_pointer(ctx,this)
      call this%get_cost(lvl,nboxes,costs,ba)
   end subroutine amrlpt_get_cost

   ! ============================================================================
   ! LIFECYCLE CALLBACKS
   ! ============================================================================

   !> Post-regrid callback: redistribute particles
   subroutine post_regrid(this,lbase,time)
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time

      ! Force resync particle container with new ba/dm
      resync: block
         integer :: lvl
         ! NOTE: this must cover EVERY fluid level, not just [0,maxlvl]. Capping the BA/DM sync
         ! leaves the container's level structure inconsistent with AmrCore and the next regrid can
         ! fail in MPI. enforce_max_level() below imposes the particle ceiling after this full sync.
         do lvl=0,this%amr%clvl()
            call this%set_particle_ba(lvl,this%amr%get_boxarray(lvl))
            call this%set_particle_dm(lvl,this%amr%get_distromap(lvl))
         end do
      end block resync

      ! Redistribute onto the new grid while also recovering any particle that
      ! an earlier unrestricted operation may have placed above maxlvl.
      call this%enforce_max_level()

      ! Optionally rebalance particles with knapsack
      if (this%rebalance) then
         rebalance_particles: block
            use amrex_amr_module, only: amrex_distromap,amrex_distromap_destroy,amrex_mfiter,amrex_boxarray
            use amrex_interface,  only: amrdm_make_knapsack
            use parallel, only: MPI_REAL_WP
            use mpi_f08, only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
            integer :: lvl,nboxes,ierr
            real(WP), dimension(:), allocatable :: costs
            type(amrex_boxarray) :: ba
            type(amrex_distromap) :: new_dm
            type(amrex_mfiter) :: mfi
            type(part), dimension(:), pointer :: p
            integer(I8) :: np_
            ! Loop over levels
            do lvl=0,this%amr%clvl()
               ! Count number of boxes at that level and prepare costs array
               ba=this%amr%get_boxarray(lvl); nboxes=int(ba%nboxes())
               allocate(costs(nboxes));costs=0.0_WP
               ! Count particles per box via MFIter
               call this%mfiter_build(lvl,mfi)
               do while (mfi%next())
                  call this%get_particles(lvl,mfi,p,np_)
                  costs(mfi%grid_index()+1)=real(np_,WP)
               end do
               call this%mfiter_destroy(mfi)
               ! Global sum
               call MPI_ALLREDUCE(MPI_IN_PLACE,costs,nboxes,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
               ! Build knapsack DM from costs
               call amrdm_make_knapsack(new_dm,costs,nboxes)
               ! Apply to particle container
               call this%set_particle_dm(lvl,new_dm)
               ! Clean up
               call amrex_distromap_destroy(new_dm)
               deallocate(costs)
            end do
            ! Redistribute particles onto new DM
            call this%redistribute()
         end block rebalance_particles
      end if

      ! Recompute VF from particle positions
      call this%update_VF()

   end subroutine post_regrid

   !> Tag cells for refinement where particle VF exceeds VF_tag
   subroutine tagging(this,lvl,time,tags)
      use amrex_amr_module, only: amrex_tagboxarray,amrex_box,amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy
      use amrgrid_class, only: SETtag
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags
      type(amrex_tagboxarray) :: tba
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF
      integer :: i,j,k
      ! Skip if VF tagging is disabled
      if (this%VF_tag.le.0.0_WP) return
      ! Resolve tagboxarray pointer
      tba=tags
      ! Loop over tiles and tag
      call amrex_mfiter_build(mfi,this%VF%mf(lvl))
      do while (mfi%next())
         tagarr=>tba%dataPtr(mfi)
         pVF=>this%VF%mf(lvl)%dataptr(mfi)
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            if (pVF(i,j,k,1).gt.this%VF_tag) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine tagging

   !> Estimate per-box cost for load balancing based on particle count
   !> This assumes that particles don't change level much...
   subroutine get_cost(this,lvl,nboxes,costs,ba)
      use amrex_amr_module, only: amrex_boxarray,amrex_box,amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy
      use parallel, only: MPI_REAL_WP
      use mpi_f08, only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in) :: lvl,nboxes
      real(WP), intent(inout) :: costs(nboxes)
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: old_bx,new_bx
      type(part), dimension(:), pointer :: p
      integer(I8) :: np_
      integer :: n,m,ierr
      integer, dimension(3) :: ijk
      real(WP) :: dxi,dyi,dzi
      ! Inverse cell size at this level
      dxi=1.0_WP/this%amr%dx(lvl); dyi=1.0_WP/this%amr%dy(lvl); dzi=1.0_WP/this%amr%dz(lvl)
      ! Zero costs, then accumulate particle counts
      costs=0.0_WP
      call this%mfiter_build(lvl,mfi)
      do while (mfi%next())
         old_bx=mfi%tilebox()
         call this%get_particles(lvl,mfi,p,np_)
         do m=1,np_
            if (p(m)%flag.eq.PART_IS_DEAD) cycle
            ! Compute cell index
            ijk(1)=floor((p(m)%pos(1)-this%amr%xlo)*dxi)
            ijk(2)=floor((p(m)%pos(2)-this%amr%ylo)*dyi)
            ijk(3)=floor((p(m)%pos(3)-this%amr%zlo)*dzi)
            ! Find which proposed new box owns this particle
            do n=1,nboxes
               new_bx=ba%get_box(n-1)
               if (.not.old_bx%intersects(new_bx)) cycle
               if (new_bx%contains(ijk)) then
                  costs(n)=costs(n)+1.0_WP; exit
               end if
            end do
         end do
      end do
      call this%mfiter_destroy(mfi)
      ! Global sum
      call MPI_ALLREDUCE(MPI_IN_PLACE,costs,nboxes,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
   end subroutine get_cost

   ! ============================================================================
   ! INITIALIZATION / FINALIZATION
   ! ============================================================================

   !> Initialize amrlpt solver
   subroutine initialize(this,amr,name,gas_gamma)
      use amrex_amr_module, only: amrex_bc_foextrap
      use messager, only: die
      implicit none
      class(amrlpt), intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), optional :: name
      real(WP), intent(in), optional :: gas_gamma
      type(part) :: layout_probe
      if (c_sizeof(layout_probe).ne.amrlpt_particle_size()) call die('[amrlpt initialize] Fortran/C++ particle layout size mismatch')
      ! Assign solver name
      if (present(name)) this%name=trim(adjustl(name))
      ! Store the thermodynamic property used by compressible drag correlations. Keep the
      ! historical air value when callers do not provide it, but reject values that would make
      ! the Osnes Mach/Knudsen corrections undefined.
      if (present(gas_gamma)) this%gas_gamma=gas_gamma
      if (this%gas_gamma.le.0.0_WP) call die('[amrlpt initialize] gas_gamma must be positive')
      ! Point to amr
      this%amr=>amr
      ! Default level cap: allow particles up to the AMR grid's max refinement
      this%maxlvl=amr%maxlvl
      ! Create particle container
      call amrlpt_new_pc(this%pc,this%amr%amrcore)
      ! Initialize VF field
      call this%VF%initialize(amr=amr,name='VF',ncomp=1,ng=this%nover); call this%VF%register()
      if (.not.this%amr%xper) then; this%VF%lo_bc(1,1)=amrex_bc_foextrap; this%VF%hi_bc(1,1)=amrex_bc_foextrap; end if
      if (.not.this%amr%yper) then; this%VF%lo_bc(2,1)=amrex_bc_foextrap; this%VF%hi_bc(2,1)=amrex_bc_foextrap; end if
      if (.not.this%amr%zper) then; this%VF%lo_bc(3,1)=amrex_bc_foextrap; this%VF%hi_bc(3,1)=amrex_bc_foextrap; end if
      ! Initialize source terms
      call this%src%initialize(amr=amr,name='src',ncomp=4,ng=this%nover); call this%src%register()
      if (.not.this%amr%xper) then; this%src%lo_bc(1,:)=amrex_bc_foextrap; this%src%hi_bc(1,:)=amrex_bc_foextrap; end if
      if (.not.this%amr%yper) then; this%src%lo_bc(2,:)=amrex_bc_foextrap; this%src%hi_bc(2,:)=amrex_bc_foextrap; end if
      if (.not.this%amr%zper) then; this%src%lo_bc(3,:)=amrex_bc_foextrap; this%src%hi_bc(3,:)=amrex_bc_foextrap; end if
      ! Register post-regrid callback
      select type (this)
       type is (amrlpt)
         call this%amr%add_postregrid(amrlpt_postregrid,c_loc(this))
         call this%amr%add_tagging   (amrlpt_tagging,   c_loc(this))
         call this%amr%set_get_cost  (amrlpt_get_cost,  c_loc(this))
      end select
      ! Initialize per-phase timers (unsynchronized for imbalance detection)
      this%tmr_coll =timer(comm=this%amr%comm,name='coll', sync=.false.)
      this%tmr_coll_=timer(comm=this%amr%comm,name='coll_',sync=.false.)
      this%tmr_fill =timer(comm=this%amr%comm,name='fill', sync=.false.)
      this%tmr_nbl  =timer(comm=this%amr%comm,name='nbl',  sync=.false.)
      this%tmr_step =timer(comm=this%amr%comm,name='step', sync=.false.)
      this%tmr_step_=timer(comm=this%amr%comm,name='step_',sync=.false.)
      this%tmr_vf   =timer(comm=this%amr%comm,name='VF',   sync=.false.)
      this%tmr_src  =timer(comm=this%amr%comm,name='src',  sync=.false.)
      ! Print out info
      call this%print()
   end subroutine initialize

   !> Finalize: destroy particle container and release grid pointer
   subroutine finalize(this)
      implicit none
      class(amrlpt), intent(inout) :: this
      call this%VF%finalize()
      call this%src%finalize()
      nullify(this%inject)
      nullify(this%user_lpt_tagging)
      call amrlpt_delete_pc(this%pc)
      this%pc=c_null_ptr
      nullify(this%amr)
      call this%tmr_coll%finalize()
      call this%tmr_coll_%finalize()
      call this%tmr_fill%finalize()
      call this%tmr_nbl%finalize()
      call this%tmr_step%finalize()
      call this%tmr_step_%finalize()
      call this%tmr_vf%finalize()
      call this%tmr_src%finalize()
   end subroutine finalize

   ! ============================================================================
   ! PHYSICS METHODS
   ! ============================================================================

   !> Soft-sphere collision model: computes Acol, Tcol on every valid particle with flag PART_COLLIDES set
   subroutine collide(this,dt,Gib,Gibcomp)
      use amrex_amr_module, only: amrex_mfiter,amrex_multifab,amrex_multifab_build,amrex_multifab_destroy,amrex_distromap
      use amrex_distromap_module, only: operator(.eq.)
      use amrdata_class, only: amrdata
      use mathtools, only: Pi,cross_product
      implicit none
      class(amrlpt), intent(inout) :: this
      real(WP), intent(in) :: dt
      type(amrdata), intent(in), optional :: Gib
      integer, intent(in), optional :: Gibcomp
      type(amrex_mfiter) :: mfi
      type(part), dimension(:), pointer :: p
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pG
      integer(c_int32_t), dimension(:), pointer :: nbr_off,nbr_lst
      integer(I8) :: np_v,np_all
      integer :: lvl,no,i1,k
      integer(c_int32_t) :: j_c
      real(WP) :: d1,m1,d2,m2
      real(WP), dimension(3) :: r1,v1,w1,r2,v2,w2
      real(WP) :: k_coeff,eta_coeff,k_coeff_w,eta_coeff_w
      real(WP) :: dx,dy,dz
      logical :: hit,dual_dm
      integer :: gc
      type(amrex_multifab) :: tmpGib

      ! Start collision timer
      call this%tmr_coll%start()

      ! Precompute spring/damping coefficients
      k_coeff=(Pi**2+log(this%e_n)**2)/this%tau_col**2
      eta_coeff=-2.0_WP*log(this%e_n)/this%tau_col
      k_coeff_w=(Pi**2+log(this%e_w)**2)/this%tau_col**2
      eta_coeff_w=-2.0_WP*log(this%e_w)/this%tau_col

      ! Resolve Gib component
      gc=1; if (present(Gibcomp)) gc=Gibcomp

      ! Reset collision counter
      this%ncol=0

      ! Ghost fill and neighbor list radius: dmax+r_influ_max=1.2*dmax
      call this%tmr_fill%start()
      call this%fill_ghosts_radius(1.2_WP*this%dmax)
      call this%tmr_fill%stop()
      call this%tmr_nbl%start()
      call this%build_neighbor_list(rcrit=1.2_WP*this%dmax)
      call this%tmr_nbl%stop()

      ! Loop over levels
      call this%tmr_coll_%start()
      do lvl=0,this%amr%clvl()

         ! Get mesh size
         dx=this%amr%dx(lvl); dy=this%amr%dy(lvl); dz=this%amr%dz(lvl)

         ! Check if we need dual-grid bridge for Gib
         if (present(Gib)) then
            dual_dm=(.not.(Gib%mf(lvl)%dm.eq.this%get_particle_dm(lvl)))
            if (dual_dm) then
               call amrex_multifab_build(mf=tmpGib,ba=this%amr%ba(lvl),dm=this%get_particle_dm(lvl),nc=Gib%mf(lvl)%ncomp(),ng=Gib%mf(lvl)%nghost(),nodal=Gib%nodal)
               call tmpGib%parallel_copy(Gib%mf(lvl),1,1,Gib%mf(lvl)%ncomp(),Gib%mf(lvl)%nghost(),Gib%mf(lvl)%nghost(),this%amr%geom(lvl))
            end if
         end if

         ! Loop over tiles
         call this%mfiter_build(lvl,mfi)
         do while (mfi%next())

            ! Get combined valid+ghost particle array and neighbor list
            call this%get_all_particles(lvl,mfi,p,np_v,np_all)
            call this%get_neighbor_list(lvl,mfi,nbr_off,nbr_lst)

            ! Gib array pointer for IB collision
            if (present(Gib)) then
               if (dual_dm) then
                  pG=>tmpGib%dataptr(mfi)
               else
                  pG=>Gib%mf(lvl)%dataptr(mfi)
               end if
            end if

            ! Loop over valid particles
            do i1=1,np_v

               ! Only particles with PART_COLLIDES
               if (IAND(p(i1)%flag,PART_COLLIDES).eq.0) cycle

               ! Zero collision acceleration and torque
               p(i1)%Acol=0.0_WP; p(i1)%Tcol=0.0_WP

               ! Set first particle properties
               r1=p(i1)%pos; v1=p(i1)%vel; w1=p(i1)%angVel
               d1=p(i1)%d;   m1=p(i1)%mass

               ! Wall and IB collisions: r2=contact point, d2=0, v2=w2=0, and m2=huge so that m_eff=m1
               v2=0.0_WP; w2=0.0_WP; d2=0.0_WP; m2=huge(m1)
               if (this%lo_bc(1).eq.AMRLPT_WALL) then; r2=[this%amr%xlo,r1(2),r1(3)]; call col_force(k_coeff_w,eta_coeff_w); end if
               if (this%hi_bc(1).eq.AMRLPT_WALL) then; r2=[this%amr%xhi,r1(2),r1(3)]; call col_force(k_coeff_w,eta_coeff_w); end if
               if (this%lo_bc(2).eq.AMRLPT_WALL) then; r2=[r1(1),this%amr%ylo,r1(3)]; call col_force(k_coeff_w,eta_coeff_w); end if
               if (this%hi_bc(2).eq.AMRLPT_WALL) then; r2=[r1(1),this%amr%yhi,r1(3)]; call col_force(k_coeff_w,eta_coeff_w); end if
               if (this%lo_bc(3).eq.AMRLPT_WALL) then; r2=[r1(1),r1(2),this%amr%zlo]; call col_force(k_coeff_w,eta_coeff_w); end if
               if (this%hi_bc(3).eq.AMRLPT_WALL) then; r2=[r1(1),r1(2),this%amr%zhi]; call col_force(k_coeff_w,eta_coeff_w); end if
               if (present(Gib)) then
                  ib_col: block
                     real(WP) :: d_ib,buf
                     real(WP), dimension(3) :: pos_p,pos_m,n12_out
                     ! Signed distance and outward IB normal nabla G/|nabla G| (from IB into fluid)
                     d_ib=this%interp(lvl,r1,pG,gc)
                     pos_p=[r1(1)+0.5_WP*dx,r1(2),r1(3)]; pos_m=[r1(1)-0.5_WP*dx,r1(2),r1(3)]; n12_out(1)=(this%interp(lvl,pos_p,pG,gc)-this%interp(lvl,pos_m,pG,gc))/dx
                     pos_p=[r1(1),r1(2)+0.5_WP*dy,r1(3)]; pos_m=[r1(1),r1(2)-0.5_WP*dy,r1(3)]; n12_out(2)=(this%interp(lvl,pos_p,pG,gc)-this%interp(lvl,pos_m,pG,gc))/dy
                     pos_p=[r1(1),r1(2),r1(3)+0.5_WP*dz]; pos_m=[r1(1),r1(2),r1(3)-0.5_WP*dz]; n12_out(3)=(this%interp(lvl,pos_p,pG,gc)-this%interp(lvl,pos_m,pG,gc))/dz
                     buf=norm2(n12_out)+epsilon(1.0_WP); n12_out=n12_out/buf; r2=r1-d_ib*n12_out
                     call col_force(k_coeff_w,eta_coeff_w)
                  end block ib_col
               end if

               ! Particle-particle collisions
               if (associated(nbr_off)) then
                  do k=nbr_off(i1)+1,nbr_off(i1+1)
                     ! Get neighbor index (note the 1-shift)
                     j_c=nbr_lst(k)+1
                     ! Skip self-collision
                     if (j_c.eq.i1) cycle
                     ! Skip if neighbor doesn't collide
                     if (IAND(p(j_c)%flag,PART_COLLIDES).eq.0) cycle
                     ! Set neighbor properties
                     r2=p(j_c)%pos; v2=p(j_c)%vel; w2=p(j_c)%angVel
                     d2=p(j_c)%d;   m2=p(j_c)%mass
                     ! Compute collision force
                     call col_force(k_coeff,eta_coeff,hit)
                     if (hit) this%ncol=this%ncol+1
                  end do
               end if

               ! Suppress forces in collapsed directions
               if (this%amr%nx.eq.1) then; p(i1)%Acol(1)=0.0_WP; p(i1)%Tcol(2)=0.0_WP; p(i1)%Tcol(3)=0.0_WP; end if
               if (this%amr%ny.eq.1) then; p(i1)%Acol(2)=0.0_WP; p(i1)%Tcol(1)=0.0_WP; p(i1)%Tcol(3)=0.0_WP; end if
               if (this%amr%nz.eq.1) then; p(i1)%Acol(3)=0.0_WP; p(i1)%Tcol(1)=0.0_WP; p(i1)%Tcol(2)=0.0_WP; end if

            end do

         end do
         call this%mfiter_destroy(mfi)

         ! Clean up temporary Gib multifab if used
         if (dual_dm.and.present(Gib)) call amrex_multifab_destroy(tmpGib)

      end do
      call this%tmr_coll_%stop()

      ! Clear ghosts
      call this%clear_ghosts()

      ! Sum collision count (each pair counted once per partner => divide by 2)
      reduce_collision_count: block
         use mpi_f08, only: MPI_IN_PLACE,MPI_INTEGER,MPI_SUM
         integer :: ierr
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%ncol,1,MPI_INTEGER,MPI_SUM,this%amr%comm,ierr)
         this%ncol=this%ncol/2
      end block reduce_collision_count

      ! Stop collision timer
      call this%tmr_coll%stop()

   contains

      !> Soft-sphere collision force between (r1,v1,w1,d1,m1) and virtual partner (r2,v2,w2,d2,m2).
      !> n12=(r2-r1)/|r2-r1| computed internally => points from particle toward partner.
      !> Wall/IB: set r2=contact point, d2=0, v2=w2=0, m2=huge => m_eff->m1, d_eff=0.5*d1.
      subroutine col_force(kk,ee,collided)
         real(WP), intent(in)  :: kk,ee
         logical,  intent(out), optional :: collided
         real(WP) :: d12,d_eff,rnv,r_influ,delta_n,rtv,m_eff
         real(WP), dimension(3) :: n12,v12,t12,f_n,f_t
         real(WP), parameter :: aclipnorm=1.0e-6_WP,acliptan=1.0e-9_WP,rcliptan=0.05_WP
         ! No collision yet
         if (present(collided)) collided=.false.
         ! Get distance
         d12=norm2(r2-r1)
         ! Skip if particles are too close - likely self-collision
         if (d12.lt.10.0_WP*epsilon(d12)) return
         ! Get normal
         n12=(r2-r1)/d12
         ! Get effective diameter, relative velocity, and relative normal velocity
         d_eff=0.5_WP*(d1+d2); v12=v1-v2; rnv=dot_product(v12,n12)
         ! Get influence radius
         r_influ=min(abs(rnv)*dt,0.2_WP*d_eff)
         ! Get overlap
         delta_n=min(d_eff+r_influ-d12,this%clip_col*d_eff)
         ! Done if no overlap
         if (delta_n.le.0.0_WP) return
         ! Collision detected
         if (present(collided)) collided=.true.
         ! Get effective mass
         m_eff=m1*m2/(m1+m2)
         ! Get tangential velocity
         t12=v12-rnv*n12+cross_product(0.5_WP*(d1*w1+d2*w2),n12)
         ! Get normal force
         f_n=(-m_eff*kk*delta_n-m_eff*ee*rnv)*n12
         ! Get tangential force
         f_t=0.0_WP
         if (this%mu_f.gt.0.0_WP) then
            rtv=sqrt(sum(t12*t12))
            if (rnv*dt/d_eff.gt.aclipnorm) then; if (   rtv/rnv  .lt.rcliptan) rtv=0.0_WP
            else;                                if (dt*rtv/d_eff.lt.acliptan) rtv=0.0_WP; end if
            if (rtv.gt.0.0_WP) f_t=-this%mu_f*norm2(f_n)*t12/rtv
         end if
         ! Increment accelerations on p(i1)
         p(i1)%Acol=p(i1)%Acol+(f_n+f_t)/m1
         p(i1)%Tcol=p(i1)%Tcol+cross_product(0.5_WP*d1*n12,f_t/m1)
      end subroutine col_force

   end subroutine collide

   !> Advance all particles by dt using Euler-midpoint substepping
   !> U,V,W: velocity amrdata (staggered MAC or collocated, with ghosts filled)
   !> rho,visc: cell-centered amrdata (with ghosts filled)
   !> Ucomp/Vcomp/Wcomp: component to use from each velocity field (optional, default 1)
   !> rhocomp/visccomp: component to use from rho/visc fields (optional, default 1)
   !> cst_rho/cst_visc: constant scalar alternatives to rho/visc amrdata
   subroutine advance(this,dt,U,Ucomp,V,Vcomp,W,Wcomp,rho,rhocomp,visc,visccomp,snd,sndcomp,cst_rho,cst_visc,cst_snd,stress_div,thermal)
      use amrex_amr_module, only: amrex_distromap
      use amrex_distromap_module, only: operator(.eq.)
      use amrex_interface, only: amrmfab_parallel_add
      use messager, only: die
      implicit none
      class(amrlpt), intent(inout) :: this
      real(WP), intent(in) :: dt
      type(amrdata), intent(in) :: U,V,W
      type(amrdata), intent(in), optional :: rho,visc,snd
      ! Resolved stress divergence (force/volume), cell-centered components x,y,z.
      type(amrdata), intent(in), optional :: stress_div,thermal
      real(WP),      intent(in), optional :: cst_rho,cst_visc,cst_snd
      integer,       intent(in), optional :: Ucomp,Vcomp,Wcomp,rhocomp,visccomp,sndcomp
      type(amrex_distromap), dimension(:), allocatable :: pdms
      logical :: any_dual_dm
      type(amrdata), allocatable :: stress_pdm,thermal_pdm
      type(amrdata) :: U_pdm,V_pdm,W_pdm,rho_pdm,visc_pdm,VF_pdm,src_pdm
      integer :: lvl

      ! Build per-level particle DM array and check for dual-grid
      allocate(pdms(0:this%amr%clvl()))
      any_dual_dm=.false.
      do lvl=0,this%amr%clvl()
         pdms(lvl)=this%get_particle_dm(lvl)
         if (.not.(this%amr%dm(lvl).eq.pdms(lvl))) any_dual_dm=.true.
      end do

      ! Clone fields onto particle DMs if needed
      if (any_dual_dm) then
         if (present(thermal)) then
            allocate(thermal_pdm)
            call thermal%clone(thermal_pdm,dm=pdms)
         end if
         if (present(stress_div)) then
            allocate(stress_pdm)
            call stress_div%clone(stress_pdm,dm=pdms)
         end if
         call U%clone(U_pdm,dm=pdms)
         call V%clone(V_pdm,dm=pdms)
         call W%clone(W_pdm,dm=pdms)
         if (present(rho))  call rho%clone(rho_pdm,dm=pdms)
         if (present(visc)) call visc%clone(visc_pdm,dm=pdms)
         call this%VF%clone(VF_pdm,dm=pdms)
         call this%src%clone(src_pdm,dm=pdms)
      end if

      ! Reset sources and drag diagnostics
      call this%src%setval(0.0_WP)
      if (any_dual_dm) call src_pdm%setval(0.0_WP)
      this%Remax=0.0_WP; this%Mamax=0.0_WP; this%Knmax=0.0_WP
      this%drag_diagnostic_record=0.0_WP
      this%drag_diagnostic_pending=this%drag_diagnostics
      if (any_dual_dm.and.present(snd)) call die('[amrlpt advance] a gridded sound speed is not supported alongside particle rebalancing (dual grid) - use cst_snd or set rebalance=.false.')

      ! Perform one step
      if (any_dual_dm) then
         if (present(rho).and.present(visc)) then
            call this%step(dt=dt,U=U_pdm,V=V_pdm,W=W_pdm,VF=VF_pdm,src=src_pdm,rho=rho_pdm,visc=visc_pdm,cst_rho=cst_rho,cst_visc=cst_visc,Ucomp=Ucomp,Vcomp=Vcomp,Wcomp=Wcomp,rhocomp=rhocomp,visccomp=visccomp,stress_div=stress_pdm,thermal=thermal_pdm)
         else if (present(rho)) then
            call this%step(dt=dt,U=U_pdm,V=V_pdm,W=W_pdm,VF=VF_pdm,src=src_pdm,rho=rho_pdm,cst_visc=cst_visc,Ucomp=Ucomp,Vcomp=Vcomp,Wcomp=Wcomp,rhocomp=rhocomp,stress_div=stress_pdm,thermal=thermal_pdm)
         else if (present(visc)) then
            call this%step(dt=dt,U=U_pdm,V=V_pdm,W=W_pdm,VF=VF_pdm,src=src_pdm,visc=visc_pdm,cst_rho=cst_rho,Ucomp=Ucomp,Vcomp=Vcomp,Wcomp=Wcomp,visccomp=visccomp,stress_div=stress_pdm,thermal=thermal_pdm)
         else
            call this%step(dt=dt,U=U_pdm,V=V_pdm,W=W_pdm,VF=VF_pdm,src=src_pdm,cst_rho=cst_rho,cst_visc=cst_visc,Ucomp=Ucomp,Vcomp=Vcomp,Wcomp=Wcomp,stress_div=stress_pdm,thermal=thermal_pdm)
         end if
      else
         call this%step(dt=dt,U=U,V=V,W=W,VF=this%VF,src=this%src,rho=rho,visc=visc,snd=snd,cst_rho=cst_rho,cst_visc=cst_visc,cst_snd=cst_snd,Ucomp=Ucomp,Vcomp=Vcomp,Wcomp=Wcomp,rhocomp=rhocomp,visccomp=visccomp,sndcomp=sndcomp,stress_div=stress_div,thermal=thermal)
      end if

      ! Writeback src_pdm to src and clean up clones
      if (any_dual_dm) then
         do lvl=0,this%amr%clvl()
            call amrmfab_parallel_add(this%src%mf(lvl),src_pdm%mf(lvl),1,1,this%src%mf(lvl)%ncomp(),0,0,this%amr%geom(lvl))
         end do
         if (allocated(thermal_pdm)) then
            call thermal_pdm%finalize()
            deallocate(thermal_pdm)
         end if
         if (allocated(stress_pdm)) then
            call stress_pdm%finalize()
            deallocate(stress_pdm)
         end if
         call U_pdm%finalize()
         call V_pdm%finalize()
         call W_pdm%finalize()
         if (present(rho))  call rho_pdm%finalize()
         if (present(visc)) call visc_pdm%finalize()
         call VF_pdm%finalize()
         call src_pdm%finalize()
      end if
      deallocate(pdms)

      ! Recompute particle volume fraction
      call this%update_VF()

      ! Process accumulated sources
      call this%tmr_src%start()
      process_sources: block
         ! Process deposit
         call this%process_deposit(this%src)
         ! Fill ghost cells
         call this%src%fill(time=0.0_WP)
         ! Filter
         call this%filter(this%src)
      end block process_sources
      call this%tmr_src%stop()

      ! Log particle advance
      log_particle_advance: block
         use iso_fortran_env, only: output_unit
         use string,          only: str_long
         use param,           only: verbose
         use messager,        only: log
         character(len=str_long) :: message
         if (verbose.gt.0) then
            write(message,'(" [",a,"] LPT advance | Np=",I0," | dt=",ES9.3)') trim(this%name),this%np,dt
            if (verbose.gt.1.and.this%amr%amRoot) write(output_unit,'(a)') trim(message)
            call log(message)
         end if
      end block log_particle_advance
      
   end subroutine advance

   !> Advance particles to a target time using CFL-limited sub-steps:
   !> Sub-step size is adapted from this%dt (persistent) using this%cflmax
   !> On exit, this%t=t_target and this%dt holds the last adapted sub-step
   subroutine advance_to(this,time,do_collide,U,V,W,rho,visc,snd,cst_rho,cst_visc,cst_snd,Ucomp,Vcomp,Wcomp,rhocomp,visccomp,sndcomp,Gib,Gibcomp,stress_div,gas_fraction,thermal)
      use amrex_amr_module, only: amrex_distromap
      use amrex_distromap_module, only: operator(.eq.)
      use amrex_interface, only: amrmfab_parallel_add
      use messager, only: die
      implicit none
      class(amrlpt), intent(inout) :: this
      real(WP), intent(in) :: time
      logical,  intent(in) :: do_collide
      type(amrdata), intent(in) :: U,V,W
      type(amrdata), intent(in), optional :: rho,visc,snd
      ! Resolved stress divergence (force/volume), cell-centered components x,y,z.
      type(amrdata), intent(in), optional :: stress_div,thermal
      real(WP),      intent(in), optional :: cst_rho,cst_visc,cst_snd
      integer,       intent(in), optional :: Ucomp,Vcomp,Wcomp,rhocomp,visccomp,sndcomp
      type(amrdata), intent(in), optional :: Gib
      ! Optional Eulerian gas fraction for normalized gas-property interpolation.
      type(amrdata), intent(in), optional :: gas_fraction
      integer,       intent(in), optional :: Gibcomp
      real(WP) :: dt_done,mydt,cfl
      integer  :: n_sub,lvl
      type(amrex_distromap), dimension(:), allocatable :: pdms
      logical :: any_dual_dm
      type(amrdata), allocatable :: stress_pdm,thermal_pdm
      type(amrdata) :: U_pdm,V_pdm,W_pdm,rho_pdm,visc_pdm,Gib_pdm,VF_pdm,src_pdm

      ! Validate time target
      if (time.le.this%t) return

      ! CFL-adapt the persistent sub-step size
      call this%get_cfl(this%dt,cflc=cfl,cfl=cfl)
      if (cfl.gt.0.0_WP) this%dt=min(this%dt*this%cflmax/cfl,this%dtmax)

      ! Build per-level particle DM array and check for dual-grid
      allocate(pdms(0:this%amr%clvl()))
      any_dual_dm=.false.
      do lvl=0,this%amr%clvl()
         pdms(lvl)=this%get_particle_dm(lvl)
         if (.not.(this%amr%dm(lvl).eq.pdms(lvl))) any_dual_dm=.true.
      end do

      ! Clone fields onto particle DMs if needed (ONE TIME before sub-step loop)
      if (any_dual_dm.and.present(gas_fraction)) &
      & call die('[amrlpt advance_to] gas-fraction sampling requires the fluid particle DistributionMapping')
      if (any_dual_dm) then
         if (present(thermal)) then
            allocate(thermal_pdm)
            call thermal%clone(thermal_pdm,dm=pdms)
         end if
         if (present(stress_div)) then
            allocate(stress_pdm)
            call stress_div%clone(stress_pdm,dm=pdms)
         end if
         call U%clone(U_pdm,dm=pdms)
         call V%clone(V_pdm,dm=pdms)
         call W%clone(W_pdm,dm=pdms)
         if (present(rho))  call rho%clone(rho_pdm,dm=pdms)
         if (present(visc)) call visc%clone(visc_pdm,dm=pdms)
         if (present(Gib))  call Gib%clone(Gib_pdm,dm=pdms)
         call this%VF%clone(VF_pdm,dm=pdms)
         call this%src%clone(src_pdm,dm=pdms)
      end if

      ! Reset sources and drag diagnostics
      call this%src%setval(0.0_WP)
      if (any_dual_dm) call src_pdm%setval(0.0_WP)
      this%Remax=0.0_WP; this%Mamax=0.0_WP; this%Knmax=0.0_WP
      this%np_immersed_loc=0; this%np_substep_capped_loc=0
      this%E_immersed_loc=0.0_WP; this%immersed_sample=0.0_WP
      this%drag_diagnostic_record=0.0_WP
      this%drag_diagnostic_pending=this%drag_diagnostics
      if (any_dual_dm.and.present(snd)) call die('[amrlpt advance_to] a gridded sound speed is not supported alongside particle rebalancing (dual grid) - use cst_snd or set rebalance=.false.')

      ! Sub-step loop: inject, collide, step, accumulate src
      dt_done=0.0_WP; n_sub=0
      do while (dt_done.lt.(time-this%t)-epsilon(1.0_WP))
         ! Select mydt
         mydt=min(this%dt,time-this%t-dt_done)
         ! Inject if needed
         if (associated(this%inject)) call this%inject(dt=mydt)
         ! Collide if needed
         if (do_collide) then
            if (any_dual_dm.and.present(Gib)) then
               call this%collide(dt=mydt,Gib=Gib_pdm,Gibcomp=Gibcomp)
            else
               call this%collide(dt=mydt,Gib=Gib,Gibcomp=Gibcomp)
            end if
         end if
         ! Step particles
         if (any_dual_dm) then
            if (present(rho).and.present(visc)) then
               call this%step(dt=mydt,U=U_pdm,V=V_pdm,W=W_pdm,VF=VF_pdm,src=src_pdm,rho=rho_pdm,visc=visc_pdm,cst_rho=cst_rho,cst_visc=cst_visc,Ucomp=Ucomp,Vcomp=Vcomp,Wcomp=Wcomp,rhocomp=rhocomp,visccomp=visccomp,stress_div=stress_pdm,thermal=thermal_pdm)
            else if (present(rho)) then
               call this%step(dt=mydt,U=U_pdm,V=V_pdm,W=W_pdm,VF=VF_pdm,src=src_pdm,rho=rho_pdm,cst_visc=cst_visc,Ucomp=Ucomp,Vcomp=Vcomp,Wcomp=Wcomp,rhocomp=rhocomp,stress_div=stress_pdm,thermal=thermal_pdm)
            else if (present(visc)) then
               call this%step(dt=mydt,U=U_pdm,V=V_pdm,W=W_pdm,VF=VF_pdm,src=src_pdm,visc=visc_pdm,cst_rho=cst_rho,Ucomp=Ucomp,Vcomp=Vcomp,Wcomp=Wcomp,visccomp=visccomp,stress_div=stress_pdm,thermal=thermal_pdm)
            else
               call this%step(dt=mydt,U=U_pdm,V=V_pdm,W=W_pdm,VF=VF_pdm,src=src_pdm,cst_rho=cst_rho,cst_visc=cst_visc,Ucomp=Ucomp,Vcomp=Vcomp,Wcomp=Wcomp,stress_div=stress_pdm,thermal=thermal_pdm)
            end if
         else
            call this%step(dt=mydt,U=U,V=V,W=W,VF=this%VF,src=this%src,rho=rho,visc=visc,snd=snd,cst_rho=cst_rho,cst_visc=cst_visc,cst_snd=cst_snd,Ucomp=Ucomp,Vcomp=Vcomp,Wcomp=Wcomp,rhocomp=rhocomp,visccomp=visccomp,sndcomp=sndcomp,stress_div=stress_div,gas_fraction=gas_fraction,thermal=thermal)
         end if
         ! Increment time and sub-step counter
         dt_done=dt_done+mydt; n_sub=n_sub+1
      end do

      ! Writeback src_pdm to src and clean up clones
      if (any_dual_dm) then
         do lvl=0,this%amr%clvl()
            call amrmfab_parallel_add(this%src%mf(lvl),src_pdm%mf(lvl),1,1,this%src%mf(lvl)%ncomp(),0,0,this%amr%geom(lvl))
         end do
         if (allocated(thermal_pdm)) then
            call thermal_pdm%finalize()
            deallocate(thermal_pdm)
         end if
         if (allocated(stress_pdm)) then
            call stress_pdm%finalize()
            deallocate(stress_pdm)
         end if
         call U_pdm%finalize()
         call V_pdm%finalize()
         call W_pdm%finalize()
         if (present(rho))  call rho_pdm%finalize()
         if (present(visc)) call visc_pdm%finalize()
         if (present(Gib))  call Gib_pdm%finalize()
         call VF_pdm%finalize()
         call src_pdm%finalize()
      end if
      deallocate(pdms)

      ! Recompute particle volume fraction
      call this%update_VF()

      ! Process accumulated sources
      call this%tmr_src%start()
      process_sources: block
         ! Process deposit
         call this%process_deposit(this%src)
         ! Fill ghost cells
         call this%src%fill(time=0.0_WP)
         ! Filter
         call this%filter(this%src)
      end block process_sources
      call this%tmr_src%stop()

      ! Advance particle time
      this%t=time

      ! Reduce the engulfed-particle counters once per advance. Doing it here rather than in
      ! get_info keeps the numbers meaningful when get_info is called more than once per step.
      reduce_immersed: block
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM
         use parallel, only: MPI_REAL_WP
         real(WP), dimension(3) :: ibuf
         integer :: ierr
         ibuf=[real(this%np_immersed_loc,WP),real(this%np_substep_capped_loc,WP),this%E_immersed_loc]
         call MPI_ALLREDUCE(MPI_IN_PLACE,ibuf,3,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
         this%np_immersed      =nint(ibuf(1))
         this%np_immersed_tot  =this%np_immersed_tot+nint(ibuf(1))
         this%np_substep_capped=this%np_substep_capped+nint(ibuf(2))
         this%E_immersed       =this%E_immersed+ibuf(3)
         this%np_immersed_loc=0; this%np_substep_capped_loc=0; this%E_immersed_loc=0.0_WP
      end block reduce_immersed

      ! Log particle advance
      log_particle_advance: block
         use iso_fortran_env, only: output_unit
         use string,          only: str_long
         use param,           only: verbose
         use messager,        only: log
         character(len=str_long) :: message
         if (verbose.gt.0) then
            write(message,'(" [",a,"] ",i0," LPT sub-steps | Np=",I0," | dt=",ES9.3," | CFL=",F6.4)') trim(this%name),n_sub,this%np,this%dt,cfl
            if (verbose.gt.1.and.this%amr%amRoot) write(output_unit,'(a)') trim(message)
            call log(message)
         end if
      end block log_particle_advance

   end subroutine advance_to

   !> Step all particles on all AMR levels by dt, deposit momentum source
   !> All amrdata arguments must be on the particle DM (caller ensures via clone)
   subroutine step(this,dt,U,V,W,VF,src,rho,visc,snd,cst_rho,cst_visc,cst_snd,Ucomp,Vcomp,Wcomp,rhocomp,visccomp,sndcomp,stress_div,gas_fraction,thermal)
      use amrex_amr_module, only: amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy
      use mathtools, only: Pi
      use messager,  only: die
      implicit none
      class(amrlpt), intent(inout) :: this
      real(WP), intent(in) :: dt
      type(amrdata), intent(in) :: U,V,W,VF
      type(amrdata), intent(in), optional :: gas_fraction,thermal
      type(amrdata), intent(inout) :: src
      type(amrdata), intent(in), optional :: rho,visc,snd
      ! Resolved stress divergence (force/volume), cell-centered components x,y,z.
      type(amrdata), intent(in), optional :: stress_div
      real(WP),      intent(in), optional :: cst_rho,cst_visc,cst_snd
      integer,       intent(in), optional :: Ucomp,Vcomp,Wcomp,rhocomp,visccomp,sndcomp

      real(WP), dimension(:,:,:,:), contiguous, pointer :: pU,pV,pW
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pRho,pVisc,pSnd,pVolFrac,pSrc,pStress,pGasFraction,pThermal
      type(amrex_mfiter) :: mfi
      type(part), dimension(:), pointer :: p
      integer(I8) :: np_
      integer :: lvl,i,uc,vc,wc,rhoc,viscc,sndc
      real(WP) :: dx,dy,dz,dxi,dyi,dzi,mydt,dt_done,Ip
      real(WP), dimension(3) :: acc,dmom,acc_drag,dmom_drag
      !> Slip velocity and acc_drag are retained from the midpoint RHS evaluation so
      !> heating uses drag alone even when acc also includes resolved stress.
      real(WP), dimension(3) :: vrel
      real(WP), dimension(4) :: dsrc
      type(part), target :: myp
      type(part) :: pold
      logical :: is_stag
      logical :: immersed_this_particle,immersed_now,no_gas
      integer :: nsub_p,nretry
      real(WP) :: dedt,thermal_dt,dheat,Tcheck
      logical :: heat_enabled

      if (this%cull_gas_free.and..not.present(gas_fraction)) call die('[amrlpt step] gas-free culling requires gas_fraction')
      heat_enabled=associated(this%heat_material)
      if (heat_enabled.neqv.present(thermal)) call die('[amrlpt heat] supply both heat_material and thermal field')
      if (heat_enabled) then
         if (.not.this%total_energy_coupling) call die('[amrlpt heat] requires total-energy coupling')
         if (.not.particle_finite(this%heat_prandtl).or.this%heat_prandtl.le.0.0_WP) call die('[amrlpt heat] invalid Prandtl number')
         if (thermal%ncomp.ne.2.or.any(thermal%nodal).or.thermal%ng.lt.this%nover) &
            call die('[amrlpt heat] thermal field requires cell-centered temperature/conductivity and particle ghosts')
      end if

      ! Start step timer
      call this%tmr_step%start()

      ! Resolve optional component indices
      uc=1; if (present(Ucomp)) uc=Ucomp
      vc=1; if (present(Vcomp)) vc=Vcomp
      wc=1; if (present(Wcomp)) wc=Wcomp
      rhoc =1; if (present(rhocomp))  rhoc =rhocomp
      viscc=1; if (present(visccomp)) viscc=visccomp
      sndc =1; if (present(sndcomp))  sndc =sndcomp
      ! Osnes needs a particle Mach number, hence a speed of sound
      if (trim(this%drag_model).eq.'Osnes'.and..not.present(snd).and..not.present(cst_snd)) &
      &   call die("[amrlpt step] drag_model='Osnes' requires snd or cst_snd")

      ! Validate: each of rho and visc must be provided in exactly one form
      if (.not.present(rho) .and..not.present(cst_rho))  call die('[amrlpt step] rho or cst_rho required')
      if (.not.present(visc).and..not.present(cst_visc)) call die('[amrlpt step] visc or cst_visc required')

      if (present(stress_div)) then
         if (any(stress_div%nodal).or.stress_div%ncomp.ne.3.or.stress_div%ng.lt.1) &
         & call die('[amrlpt step] stress_div requires 3 cell-centered components and filled ghosts')
      end if

      ! Check velocity nodal locations
      if (present(gas_fraction)) then
         if (any(gas_fraction%nodal).or.gas_fraction%ncomp.ne.1.or.gas_fraction%ng.lt.this%nover) &
         & call die('[amrlpt step] gas_fraction requires one cell-centered component and particle ghost width')
      end if
      check_velocity: block
         logical, dimension(3) :: nU,nV,nW
         nU=U%nodal; nV=V%nodal; nW=W%nodal
         if (all(nU.eqv.[.true.,.false.,.false.]).and. &
         &   all(nV.eqv.[.false.,.true.,.false.]).and. &
         &   all(nW.eqv.[.false.,.false.,.true.])) then
            is_stag=.true.
         else if (.not.any(nU).and..not.any(nV).and..not.any(nW)) then
            is_stag=.false.
         else
            call die('[amrlpt step] U/V/W must be staggered (face-centered) or collocated (cell-centered)')
         end if
      end block check_velocity

      ! Loop over all AMR levels
      call this%tmr_step_%start()
      do lvl=0,this%amr%clvl()

         ! Get mesh size
         dx=this%amr%dx(lvl); dxi=1.0_WP/dx
         dy=this%amr%dy(lvl); dyi=1.0_WP/dy
         dz=this%amr%dz(lvl); dzi=1.0_WP/dz

         ! MFIter over the level
         call this%mfiter_build(lvl,mfi)
         do while (mfi%next())

            ! Get pointers to data
            pU      =>U%mf(lvl)%dataptr(mfi)
            pV      =>V%mf(lvl)%dataptr(mfi)
            pW      =>W%mf(lvl)%dataptr(mfi)
            if (present(rho))  pRho =>rho%mf(lvl)%dataptr(mfi)
            if (present(visc)) pVisc=>visc%mf(lvl)%dataptr(mfi)
            if (present(snd))  pSnd =>snd%mf(lvl)%dataptr(mfi)
            if (present(stress_div)) pStress=>stress_div%mf(lvl)%dataptr(mfi)
            if (present(gas_fraction)) pGasFraction=>gas_fraction%mf(lvl)%dataptr(mfi)
            if (heat_enabled) pThermal=>thermal%mf(lvl)%dataptr(mfi)
            pVolFrac=>VF%mf(lvl)%dataptr(mfi)
            pSrc    =>src%mf(lvl)%dataptr(mfi)

            ! Get particles on this tile
            call this%get_particles(lvl=lvl,mfi=mfi,p=p,np=np_)

            ! Loop over local particles
            do i=1,np_
               ! Skip particles that are not moving or exchanging
               if (IAND(p(i)%flag,PART_MOVES+PART_EXCHANGES).eq.0) cycle
               ! Create copy of particle
               myp=p(i)
               immersed_this_particle=.false.; no_gas=.false.
               ! Time-integrate until dt_done=dt
               dt_done=0.0_WP; nsub_p=0
               do while (dt_done.lt.dt)
                  mydt=min(myp%dt,dt-dt_done)
                  if (mydt.le.0.0_WP) mydt=dt-dt_done
                  nsub_p=nsub_p+1
                  if (heat_enabled.and.nsub_p.ge.this%nsub_max) &
                     call die('[amrlpt heat] particle substep limit reached; increase nsub_max or reduce flow timestep')
                  if (nsub_p.ge.this%nsub_max) then                ! <-- add: backstop
                     mydt=dt-dt_done
                     this%np_substep_capped_loc=this%np_substep_capped_loc+1
                  end if
                  immersed_now=.false.
                  ! Remember the particle
                  pold=myp
                  ! Precompute moment of inertia for a sphere
                  Ip=0.1_WP*myp%d**2
                  ! Existing explicit midpoint stages, extended to internal energy. Retry only
                  ! an unaccepted trial if the midpoint thermal time becomes more restrictive.
                  nretry=0
                  do
                     nretry=nretry+1
                     if (nretry.gt.50) call die('[amrlpt heat] failed to find a valid thermal substep')
                     myp=pold
                     acc=get_rhs()
                     if (no_gas) exit
                     if (heat_enabled) mydt=min(mydt,0.5_WP*thermal_dt)
                     if (dt_done+mydt.le.dt_done) call die('[amrlpt step] substep underflow')
                     if (heat_enabled) myp%eint=pold%eint+0.5_WP*mydt*dedt
                     if (IAND(myp%flag,PART_MOVES).ne.0) then
                        myp%pos=pold%pos+0.5_WP*mydt*myp%vel
                        myp%vel=pold%vel+0.5_WP*mydt*(acc+this%gravity+myp%Acol)
                        myp%angVel=pold%angVel+0.5_WP*mydt*myp%Tcol/Ip
                     end if
                     acc=get_rhs()
                     if (no_gas) exit
                     if (heat_enabled) then
                        if (mydt.gt.0.5_WP*thermal_dt*(1.0_WP+1.0e-12_WP)) then
                           mydt=0.5_WP*thermal_dt
                           cycle
                        end if
                        myp%eint=pold%eint+mydt*dedt
                        Tcheck=this%heat_material%get_T_from_rho_e(rho=myp%mass/(Pi/6.0_WP*myp%d**3),e=myp%eint,y=[1.0_WP])
                        if (.not.particle_finite(Tcheck).or.Tcheck.le.0.0_WP.or..not.particle_finite(myp%eint)) then
                           mydt=0.5_WP*mydt
                           cycle
                        end if
                     end if
                     if (IAND(myp%flag,PART_MOVES).ne.0) then
                        myp%pos=pold%pos+mydt*myp%vel
                        myp%vel=pold%vel+mydt*(acc+this%gravity+myp%Acol)
                        myp%angVel=pold%angVel+mydt*myp%Tcol/Ip
                     end if
                     exit
                  end do

                  if (no_gas) then
                     ! Discard the trial at either stage; retain all previously accepted deposits.
                     myp=pold
                     this%np_culled_loc=this%np_culled_loc+1
                     this%culled_loc=this%culled_loc+[myp%mass,myp%mass*myp%vel, &
                        myp%mass*(myp%eint+0.5_WP*dot_product(myp%vel,myp%vel))]
                     call amrlpt_set_particle_id(-1_c_int64_t,c_loc(myp))
                     exit
                  end if

                  ! Catch a particle that has left its deposit stencil (or gone NaN) within one
                  ! sub-step - almost always a cflmax/dtmax that is too loose for the local flow
                  if (floor((myp%pos(1)-this%amr%xlo)*dxi-0.5_WP).lt.lbound(pSrc,1) .or. &
                  &   floor((myp%pos(1)-this%amr%xlo)*dxi-0.5_WP)+1.gt.ubound(pSrc,1) .or. &
                  &   myp%pos(1).ne.myp%pos(1)) then
                     write(*,'(A,I0,A,I0)') 'lvl=',lvl,' i=',i
                     write(*,'(A,3ES18.10)') 'pold%pos=',pold%pos
                     write(*,'(A,3ES18.10)') ' myp%pos=',myp%pos
                     write(*,'(A,3ES18.10)') ' myp%vel=',myp%vel
                     write(*,'(A,3ES18.10)') '     acc=',acc
                     write(*,'(A,3ES18.10)') '    Acol=',myp%Acol
                     write(*,'(A,ES18.10)')  '      dt=',mydt
                     write(*,'(A,I0,A,I0)')  'pSrc x bounds: ',lbound(pSrc,1),' to ',ubound(pSrc,1)
                     call die('[amrlpt step] particle left its deposit stencil in one sub-step - tighten cflmax/dtmax')
                  end if

                  ! Hydrodynamic reaction only: gravity and contact forces are external to the fluid.
                  if (IAND(myp%flag,PART_EXCHANGES).ne.0) then
                     dmom=mydt*acc*myp%mass
                     dsrc(1:3)=-dmom
                     dmom_drag=mydt*acc_drag*myp%mass
                     dsrc(4)=dot_product(dmom_drag,vrel)
                     this%exchange_loc(1:3)=this%exchange_loc(1:3)+dsrc(1:3)
                     this%exchange_loc(5)=this%exchange_loc(5)+dsrc(4)
                     if (this%total_energy_coupling) then
                        ! Endpoint-average velocity gives exact hydrodynamic work in the discrete
                        ! kinetic-energy balance, including when gravity/contact also change velocity.
                        ! A prescribed nonmoving particle keeps its velocity: this is support work,
                        ! not a change of particle KE (zero work for a stationary supported particle).
                        dsrc(4)=-dot_product(0.5_WP*(pold%vel+myp%vel),dmom)
                        this%exchange_loc(4)=this%exchange_loc(4)+dsrc(4)
                     else if (immersed_now) then
                        ! Legacy gas-internal-energy callers cannot heat an absent gas phase.
                        this%E_immersed_loc=this%E_immersed_loc+dsrc(4)
                        dsrc(4)=0.0_WP
                     end if
                     ! Book the actual accepted internal-energy increment, including roundoff.
                     dheat=myp%mass*(myp%eint-pold%eint)
                     dsrc(4)=dsrc(4)-dheat
                     this%heat_exchange_loc=this%heat_exchange_loc-dheat
                     call deposit(val=dsrc)
                  end if
                  ! Increment
                  dt_done=dt_done+mydt
               end do
               ! Count the particle once if it was engulfed at any point this step   <-- add
               if (immersed_this_particle) this%np_immersed_loc=this%np_immersed_loc+1
               ! Track escape from non-periodic boundaries
               if (.not.no_gas) then
                  if ((.not.this%amr%xper.and.(myp%pos(1).lt.this%amr%xlo.or.myp%pos(1).gt.this%amr%xhi)).or. &
                  &   (.not.this%amr%yper.and.(myp%pos(2).lt.this%amr%ylo.or.myp%pos(2).gt.this%amr%yhi)).or. &
                  &   (.not.this%amr%zper.and.(myp%pos(3).lt.this%amr%zlo.or.myp%pos(3).gt.this%amr%zhi))) then
                     this%np_out_loc=this%np_out_loc+1
                     this%Vp_out_loc=this%Vp_out_loc+Pi/6.0_WP*myp%d**3
                     this%outflow_loc=this%outflow_loc+[myp%mass,myp%mass*myp%vel, &
                        myp%mass*(myp%eint+0.5_WP*dot_product(myp%vel,myp%vel))]
                  end if
               end if
               ! Write back
               p(i)=myp
            end do
         end do
         call this%mfiter_destroy(mfi)

      end do
      call this%tmr_step_%stop()

      ! Redistribute particles
      call this%redistribute()

      ! Stop step timer
      call this%tmr_step%stop()

   contains

      !> Calculate rhs of particle equations of motion
      function get_rhs() result(acc)
         use mathtools, only: Pi
         real(WP), dimension(3) :: acc,fvel
         integer :: c
         real(WP), dimension(30) :: diag_record
         real(WP) :: pVF,fVF,frho,fvisc,fC,Re,tau,corr,b1,b2,diag_Re,diag_Ma,diag_Kn,rho_p
         real(WP) :: Tg,kg,Tp,cvp,Nu,H,Reheat
         real(WP) :: wx,wy,wz,fx,fy,fz
         real(WP) :: gas_weights(0:1,0:1,0:1),gas_sum,gas_scale
         integer :: a,b,cg
         integer  :: ic,jc,kc,ix,iy,iz
         dedt=0.0_WP; thermal_dt=huge(1.0_WP); no_gas=.false.
         rho_p=myp%mass/(Pi/6.0_WP*myp%d**3)
         ! Cell-centered stencil indices and weights
         ic=floor((myp%pos(1)-this%amr%xlo)*dxi-0.5_WP); wx=(myp%pos(1)-this%amr%xlo)*dxi-0.5_WP-real(ic,WP)
         jc=floor((myp%pos(2)-this%amr%ylo)*dyi-0.5_WP); wy=(myp%pos(2)-this%amr%ylo)*dyi-0.5_WP-real(jc,WP)
         kc=floor((myp%pos(3)-this%amr%zlo)*dzi-0.5_WP); wz=(myp%pos(3)-this%amr%zlo)*dzi-0.5_WP-real(kc,WP)
         ! Velocity interpolation
         if (is_stag) then
            ! Face-centered indices and weights
            ix=floor((myp%pos(1)-this%amr%xlo)*dxi); fx=(myp%pos(1)-this%amr%xlo)*dxi-real(ix,WP)
            iy=floor((myp%pos(2)-this%amr%ylo)*dyi); fy=(myp%pos(2)-this%amr%ylo)*dyi-real(iy,WP)
            iz=floor((myp%pos(3)-this%amr%zlo)*dzi); fz=(myp%pos(3)-this%amr%zlo)*dzi-real(iz,WP)
            fvel(1)=(1.0_WP-fx)*(1.0_WP-wy)*(1.0_WP-wz)*pU(ix,jc,kc,uc)+fx*(1.0_WP-wy)*(1.0_WP-wz)*pU(ix+1,jc,kc,uc)+(1.0_WP-fx)*wy*(1.0_WP-wz)*pU(ix,jc+1,kc,uc)+fx*wy*(1.0_WP-wz)*pU(ix+1,jc+1,kc,uc)+(1.0_WP-fx)*(1.0_WP-wy)*wz*pU(ix,jc,kc+1,uc)+fx*(1.0_WP-wy)*wz*pU(ix+1,jc,kc+1,uc)+(1.0_WP-fx)*wy*wz*pU(ix,jc+1,kc+1,uc)+fx*wy*wz*pU(ix+1,jc+1,kc+1,uc)
            fvel(2)=(1.0_WP-wx)*(1.0_WP-fy)*(1.0_WP-wz)*pV(ic,iy,kc,vc)+wx*(1.0_WP-fy)*(1.0_WP-wz)*pV(ic+1,iy,kc,vc)+(1.0_WP-wx)*fy*(1.0_WP-wz)*pV(ic,iy+1,kc,vc)+wx*fy*(1.0_WP-wz)*pV(ic+1,iy+1,kc,vc)+(1.0_WP-wx)*(1.0_WP-fy)*wz*pV(ic,iy,kc+1,vc)+wx*(1.0_WP-fy)*wz*pV(ic+1,iy,kc+1,vc)+(1.0_WP-wx)*fy*wz*pV(ic,iy+1,kc+1,vc)+wx*fy*wz*pV(ic+1,iy+1,kc+1,vc)
            fvel(3)=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-fz)*pW(ic,jc,iz,wc)+wx*(1.0_WP-wy)*(1.0_WP-fz)*pW(ic+1,jc,iz,wc)+(1.0_WP-wx)*wy*(1.0_WP-fz)*pW(ic,jc+1,iz,wc)+wx*wy*(1.0_WP-fz)*pW(ic+1,jc+1,iz,wc)+(1.0_WP-wx)*(1.0_WP-wy)*fz*pW(ic,jc,iz+1,wc)+wx*(1.0_WP-wy)*fz*pW(ic+1,jc,iz+1,wc)+(1.0_WP-wx)*wy*fz*pW(ic,jc+1,iz+1,wc)+wx*wy*fz*pW(ic+1,jc+1,iz+1,wc)
         else
            fvel(1)=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz)*pU(ic,jc,kc,uc)+wx*(1.0_WP-wy)*(1.0_WP-wz)*pU(ic+1,jc,kc,uc)+(1.0_WP-wx)*wy*(1.0_WP-wz)*pU(ic,jc+1,kc,uc)+wx*wy*(1.0_WP-wz)*pU(ic+1,jc+1,kc,uc)+(1.0_WP-wx)*(1.0_WP-wy)*wz*pU(ic,jc,kc+1,uc)+wx*(1.0_WP-wy)*wz*pU(ic+1,jc,kc+1,uc)+(1.0_WP-wx)*wy*wz*pU(ic,jc+1,kc+1,uc)+wx*wy*wz*pU(ic+1,jc+1,kc+1,uc)
            fvel(2)=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz)*pV(ic,jc,kc,vc)+wx*(1.0_WP-wy)*(1.0_WP-wz)*pV(ic+1,jc,kc,vc)+(1.0_WP-wx)*wy*(1.0_WP-wz)*pV(ic,jc+1,kc,vc)+wx*wy*(1.0_WP-wz)*pV(ic+1,jc+1,kc,vc)+(1.0_WP-wx)*(1.0_WP-wy)*wz*pV(ic,jc,kc+1,vc)+wx*(1.0_WP-wy)*wz*pV(ic+1,jc,kc+1,vc)+(1.0_WP-wx)*wy*wz*pV(ic,jc+1,kc+1,vc)+wx*wy*wz*pV(ic+1,jc+1,kc+1,vc)
            fvel(3)=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz)*pW(ic,jc,kc,wc)+wx*(1.0_WP-wy)*(1.0_WP-wz)*pW(ic+1,jc,kc,wc)+(1.0_WP-wx)*wy*(1.0_WP-wz)*pW(ic,jc+1,kc,wc)+wx*wy*(1.0_WP-wz)*pW(ic+1,jc+1,kc,wc)+(1.0_WP-wx)*(1.0_WP-wy)*wz*pW(ic,jc,kc+1,wc)+wx*(1.0_WP-wy)*wz*pW(ic+1,jc,kc+1,wc)+(1.0_WP-wx)*wy*wz*pW(ic,jc+1,kc+1,wc)+wx*wy*wz*pW(ic+1,jc+1,kc+1,wc)
         end if
         ! Density and viscosity
         if (present(cst_rho)) then
            frho=cst_rho
         else
            frho=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz)*pRho(ic,jc,kc,rhoc)+wx*(1.0_WP-wy)*(1.0_WP-wz)*pRho(ic+1,jc,kc,rhoc)+(1.0_WP-wx)*wy*(1.0_WP-wz)*pRho(ic,jc+1,kc,rhoc)+wx*wy*(1.0_WP-wz)*pRho(ic+1,jc+1,kc,rhoc)+(1.0_WP-wx)*(1.0_WP-wy)*wz*pRho(ic,jc,kc+1,rhoc)+wx*(1.0_WP-wy)*wz*pRho(ic+1,jc,kc+1,rhoc)+(1.0_WP-wx)*wy*wz*pRho(ic,jc+1,kc+1,rhoc)+wx*wy*wz*pRho(ic+1,jc+1,kc+1,rhoc)
         end if
         if (present(cst_visc)) then
            fvisc=cst_visc
         else
            fvisc=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz)*pVisc(ic,jc,kc,viscc)+wx*(1.0_WP-wy)*(1.0_WP-wz)*pVisc(ic+1,jc,kc,viscc)+(1.0_WP-wx)*wy*(1.0_WP-wz)*pVisc(ic,jc+1,kc,viscc)+wx*wy*(1.0_WP-wz)*pVisc(ic+1,jc+1,kc,viscc)+(1.0_WP-wx)*(1.0_WP-wy)*wz*pVisc(ic,jc,kc+1,viscc)+wx*(1.0_WP-wy)*wz*pVisc(ic+1,jc,kc+1,viscc)+(1.0_WP-wx)*wy*wz*pVisc(ic,jc+1,kc+1,viscc)+wx*wy*wz*pVisc(ic+1,jc+1,kc+1,viscc)
         end if
         fvisc=fvisc+epsilon(1.0_WP)
         ! Speed of sound (only needed by compressible drag laws)
         fC=0.0_WP
         if (present(cst_snd)) then
            fC=cst_snd
         else if (present(snd)) then
            fC=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz)*pSnd(ic,jc,kc,sndc)+wx*(1.0_WP-wy)*(1.0_WP-wz)*pSnd(ic+1,jc,kc,sndc)+(1.0_WP-wx)*wy*(1.0_WP-wz)*pSnd(ic,jc+1,kc,sndc)+wx*wy*(1.0_WP-wz)*pSnd(ic+1,jc+1,kc,sndc)+(1.0_WP-wx)*(1.0_WP-wy)*wz*pSnd(ic,jc,kc+1,sndc)+wx*(1.0_WP-wy)*wz*pSnd(ic+1,jc,kc+1,sndc)+(1.0_WP-wx)*wy*wz*pSnd(ic,jc+1,kc+1,sndc)+wx*wy*wz*pSnd(ic+1,jc+1,kc+1,sndc)
         end if
         ! Gas properties must not be diluted by zero-valued pure-liquid entries. Keep the
         ! ordinary interpolation unchanged for callers that do not supply a gas fraction.
         do cg=0,1; do b=0,1; do a=0,1
            gas_weights(a,b,cg)=merge(wx,1.0_WP-wx,a.eq.1)*merge(wy,1.0_WP-wy,b.eq.1)*merge(wz,1.0_WP-wz,cg.eq.1)
         end do; end do; end do
         if (present(gas_fraction)) then
            if (.not.all(particle_finite(pGasFraction(ic:ic+1,jc:jc+1,kc:kc+1,1))).or. &
                minval(pGasFraction(ic:ic+1,jc:jc+1,kc:kc+1,1)).lt.0.0_WP) &
               call die('[amrlpt step] invalid gas sampling fraction')
            ! Scale before multiplying: retain tiny positive fractions without an absolute cutoff.
            gas_scale=maxval(pGasFraction(ic:ic+1,jc:jc+1,kc:kc+1,1),mask=gas_weights.gt.0.0_WP)
            if (gas_scale.gt.0.0_WP) then
               do cg=0,1; do b=0,1; do a=0,1
                  if (gas_weights(a,b,cg).gt.0.0_WP) &
                     gas_weights(a,b,cg)=gas_weights(a,b,cg)*(pGasFraction(ic+a,jc+b,kc+cg,1)/gas_scale)
               end do; end do; end do
            else
               gas_weights=0.0_WP
            end if
            gas_sum=sum(gas_weights)
            if (.not.(gas_sum.gt.0.0_WP)) then
               if (this%cull_gas_free) then
                  no_gas=.true.; acc=0.0_WP
                  return
               end if
               ! ENGULFED: the whole stencil is pure liquid, so there is no gas state to build a
               ! drag law from. The particle is made of the same liquid that surrounds it, so it
               ! is neutrally buoyant: use the liquid's own density and viscosity in a Schiller-Naumann law 
               if (this%visc_immersed.le.0.0_WP) &
               &   call die('[amrlpt step] engulfed particle needs visc_immersed (liquid dynamic viscosity)')
               immersed_now=.true.
               frho =this%rho
               fvisc=this%visc_immersed
               vrel =fvel-myp%vel
               Re   =frho*norm2(vrel)*myp%d/fvisc+epsilon(1.0_WP)
               corr =1.0_WP+0.15_WP*Re**(0.687_WP)
               tau  =rho_p*myp%d**2/(18.0_WP*fvisc*corr)
               acc_drag=vrel/tau
               acc  =acc_drag
               if (present(stress_div)) then
                  do c=1,3
                     acc(c)=acc(c)+this%interp(lvl,myp%pos,pStress,c)/rho_p
                  end do
               end if
               myp%dt=tau/real(this%nstep,WP)
               if (.not.immersed_this_particle) then
                  immersed_this_particle=.true.
                  if (this%np_immersed_loc.lt.5) then
                     report_immersed: block
                        use iso_fortran_env, only: output_unit
                        real(WP) :: gasmin,gasmax
                        gasmin=minval(pGasFraction(ic:ic+1,jc:jc+1,kc:kc+1,1))
                        gasmax=maxval(pGasFraction(ic:ic+1,jc:jc+1,kc:kc+1,1))
                        this%immersed_sample=[real(lvl,WP),myp%pos,myp%d,norm2(vrel),gasmin,gasmax]
                        write(output_unit,'(" [",a,"] rank ",i0," engulfed particle: lvl=",i0, &
                        &   " pos=",3es12.5," d=",es12.5," |vrel|=",es12.5," Re_l=",es12.5, &
                        &   " tau=",es12.5," t=",es12.5)') &
                        &   trim(this%name),this%amr%rank,lvl,myp%pos,myp%d,norm2(vrel),Re,tau,this%t
                        flush(output_unit)
                     end block report_immersed
                  end if
               end if
               return
            end if
            gas_weights=gas_weights/gas_sum
            if (.not.present(cst_rho)) frho=sum(gas_weights*pRho(ic:ic+1,jc:jc+1,kc:kc+1,rhoc))
            if (.not.present(cst_visc)) fvisc=sum(gas_weights*pVisc(ic:ic+1,jc:jc+1,kc:kc+1,viscc))+epsilon(1.0_WP)
            if (present(snd).and..not.present(cst_snd)) fC=sum(gas_weights*pSnd(ic:ic+1,jc:jc+1,kc:kc+1,sndc))
         end if
         if (heat_enabled.and.IAND(myp%flag,PART_EXCHANGES).ne.0) then
            Tg=sum(gas_weights*pThermal(ic:ic+1,jc:jc+1,kc:kc+1,1))
            kg=sum(gas_weights*pThermal(ic:ic+1,jc:jc+1,kc:kc+1,2))
            Tp=this%heat_material%get_T_from_rho_e(rho=rho_p,e=myp%eint,y=[1.0_WP])
            cvp=this%heat_material%get_cv_from_rho_e(rho=rho_p,e=myp%eint,y=[1.0_WP])
            if (.not.all(particle_finite([Tg,kg,Tp,cvp,frho,fvisc])).or.min(Tg,Tp,cvp,frho,fvisc).le.0.0_WP.or.kg.lt.0.0_WP) &
               call die('[amrlpt heat] invalid sampled gas properties or particle thermal state')
            ! Gunn (1978), isolated-particle porosity=1, independently of Eulerian gas fraction.
            Reheat=frho*norm2(fvel-myp%vel)*myp%d/fvisc
            Nu=2.0_WP*(1.0_WP+0.7_WP*Reheat**0.2_WP*this%heat_prandtl**(1.0_WP/3.0_WP)) &
               +0.13_WP*Reheat**0.7_WP*this%heat_prandtl**(1.0_WP/3.0_WP)
            H=Pi*myp%d*kg*Nu ! h*A = (Nu*k/d)*(pi*d^2)
            if (.not.particle_finite(H)) call die('[amrlpt heat] non-finite heat conductance')
            if (H.gt.0.0_WP) then
               dedt=H*(Tg-Tp)/myp%mass
               thermal_dt=myp%mass*cvp/H
               if (.not.particle_finite(dedt).or..not.particle_finite(thermal_dt).or.thermal_dt.le.0.0_WP) &
                  call die('[amrlpt heat] invalid heat rate or response time')
            end if
         end if
         fC=fC+epsilon(1.0_WP)
         ! Volume fraction
         pVF=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz)*pVolFrac(ic,jc,kc,1)+wx*(1.0_WP-wy)*(1.0_WP-wz)*pVolFrac(ic+1,jc,kc,1)+(1.0_WP-wx)*wy*(1.0_WP-wz)*pVolFrac(ic,jc+1,kc,1)+wx*wy*(1.0_WP-wz)*pVolFrac(ic+1,jc+1,kc,1)+(1.0_WP-wx)*(1.0_WP-wy)*wz*pVolFrac(ic,jc,kc+1,1)+wx*(1.0_WP-wy)*wz*pVolFrac(ic+1,jc,kc+1,1)+(1.0_WP-wx)*wy*wz*pVolFrac(ic,jc+1,kc+1,1)+wx*wy*wz*pVolFrac(ic+1,jc+1,kc+1,1)
         if (this%drag_vf_override.ge.0.0_WP) pVF=this%drag_vf_override
         pVF=min(max(pVF,0.0_WP),1.0_WP-epsilon(1.0_WP))
         fVF=max(1.0_WP-pVF,epsilon(1.0_WP))
         diag_Re=0.0_WP; diag_Ma=0.0_WP; diag_Kn=0.0_WP
         ! Drag correction factor
         select case(trim(this%drag_model))
         case('None','none')
            corr=epsilon(1.0_WP)
         case('Stokes')
            corr=1.0_WP
         case('Schiller-Naumann','Schiller Naumann','SN')
            Re=frho*norm2(myp%vel-fvel)*myp%d/fvisc+epsilon(1.0_WP)
            corr=1.0_WP+0.15_WP*Re**(0.687_WP)
         case('Tenneti') ! Tenneti and Subramaniam (2011)
            Re=fVF*frho*norm2(myp%vel-fvel)*myp%d/fvisc+epsilon(1.0_WP)
            b1=5.81_WP*pVF/fVF**3+0.48_WP*pVF**(1.0_WP/3.0_WP)/fVF**4
            b2=pVF**3*Re*(0.95_WP+0.61_WP*pVF**3/fVF**2)
            corr=fVF*((1.0_WP+0.15_WP*Re**(0.687_WP))/fVF**3+b1+b2)
         case('Osnes') ! Osnes et al. (2023) IJMF - valid to high Re and transonic Ma
            osnes: block
               real(WP) :: Rep,Rem,Ma,Knp,fKn,CD1,CD2,CD3,sDrag,GM,HM,JM,CM,b3,vp
               Rep=frho*norm2(myp%vel-fvel)*myp%d/fvisc+epsilon(1.0_WP)
               Ma =norm2(myp%vel-fvel)/fC+epsilon(1.0_WP)
               ! Evaluate Ma/Re algebraically from the sampled properties. This remains finite
               ! at zero slip, where forming the ratio from separately epsilon-regularized Ma
               ! and Re would return a spurious O(1) Knudsen number.
               Knp=sqrt(0.5_WP*Pi*this%gas_gamma)*fvisc/(frho*fC*myp%d)
               diag_Re=Rep; diag_Ma=Ma; diag_Kn=Knp
               if (Rep.le.45.0_WP) then
                  ! Rarefaction-dominated regime
                  fKn=1.0_WP/(1.0_WP+Knp*(2.514_WP+0.8_WP*exp(-0.55_WP/Knp)))
                  CD1=(1.0_WP+0.15_WP*Rep**(0.687_WP))*fKn
                  ! The free-molecular expression below contains cancelling O(Ma^-3) terms.
                  ! Its contribution to corr is O(Ma^4), so use the analytic continuum-Mach
                  ! limit when it is negligible rather than evaluating an ill-conditioned form.
                  if (Ma.lt.1.0e-3_WP) then
                     corr=CD1
                  else
                     sDrag=Ma*sqrt(0.5_WP*this%gas_gamma)
                     if (Ma.le.1.0_WP) then
                        JM=2.26_WP-0.1_WP/Ma+0.14_WP/Ma**3
                     else
                        JM=1.6_WP+0.25_WP/Ma+0.11_WP/Ma**2+0.44_WP/Ma**3
                     end if
                     CD2=(1.0_WP+2.0_WP*sDrag**2)*exp(-sDrag**2)/(sDrag**3*sqrt(Pi)+epsilon(1.0_WP)) &
                     &  +(4.0_WP*sDrag**4+4.0_WP*sDrag**2-1.0_WP)*erf(sDrag)/(2.0_WP*sDrag**4+epsilon(1.0_WP)) &
                     &  +2.0_WP/(3.0_WP*sDrag+epsilon(1.0_WP))*sqrt(Pi)
                     CD2=CD2/(1.0_WP+(CD2/JM-1.0_WP)*sqrt(Rep/45.0_WP))
                     corr=CD1/(1.0_WP+Ma**4)+Rep/24.0_WP*Ma**4*CD2/(1.0_WP+Ma**4)
                  end if
                  CD3=corr*24.0_WP/Rep
               else
                  ! Compression-dominated regime
                  if (Ma.lt.1.5_WP) then
                     CM=1.65_WP+0.65_WP*tanh(4.0_WP*Ma-3.4_WP)
                  else
                     CM=2.18_WP-0.13_WP*tanh(0.9_WP*Ma-2.7_WP)
                  end if
                  if (Ma.lt.0.8_WP) then
                     GM=166.0_WP*Ma**3+3.29_WP*Ma**2-10.9_WP*Ma+20.0_WP
                  else
                     GM=5.0_WP+40.0_WP*Ma**(-3.0_WP)
                  end if
                  if (Ma.lt.1.0_WP) then
                     HM=0.0239_WP*Ma**3+0.212_WP*Ma**2-0.074_WP*Ma+1.0_WP
                  else
                     HM=0.93_WP+1.0_WP/(3.5_WP+Ma**5)
                  end if
                  ! CD3 is the drag coefficient, not the Stokes correction factor.
                  ! The Rep/24 conversion is applied once below when corr is formed.
                  CD3=24.0_WP/Rep*(1.0_WP+0.15_WP*Rep**(0.687_WP))*HM+ &
                  &   0.42_WP*CM/(1.0_WP+42500.0_WP/Rep**(1.16_WP*CM)+GM/sqrt(Rep))
               end if
               ! Volume-fraction correction
               vp=min(0.4_WP,pVF)
               CD3=CD3/(1.0_WP-vp)
               Rem=(1.0_WP-vp)*Rep
               b1=(5.81_WP*vp/(1.0_WP-vp)**2+0.48_WP*vp**(1.0_WP/3.0_WP)/(1.0_WP-vp)**3)*(1.0_WP-vp)
               b2=(1.0_WP-vp)*(vp**3*Rem*(0.95_WP+0.61_WP*vp**3/(1.0_WP-vp)**2))*(1.0_WP-vp)
               b3=min(sqrt(20.0_WP*Ma),1.0_WP)*(5.65_WP*vp-22.0_WP*vp**2+23.4_WP*vp**3) &
               &  *(1.0_WP+tanh((Ma-(0.65_WP-0.24_WP*vp))/0.35_WP))
               corr=CD3*Rep/24.0_WP+b1+b2+b3*Rep/24.0_WP
            end block osnes
         case('Beetstra') ! Beetstra et al. (2007)
            Re=fVF*frho*norm2(myp%vel-fvel)*myp%d/fvisc+epsilon(1.0_WP)
            b1=10.0_WP*pVF/fVF**2+fVF**2*(1.0_WP+1.5_WP*sqrt(pVF))
            b2=0.413_WP/24.0_WP*Re/fVF**2*(1.0_WP/fVF+3.0_WP*fVF*pVF+8.4_WP*Re**(-0.343_WP))/(1.0_WP+10.0_WP**(3.0_WP*pVF)*Re**(2.0_WP*fVF-2.5_WP))
            corr=b1+b2
         case default
            corr=1.0_WP
         end select
         ! Particle response time
         tau=rho_p*myp%d**2/(18.0_WP*fvisc*corr)
         ! Return acceleration and update particle timestep size. vrel is handed back to step
         ! by host association so the dissipation it deposits matches this acc exactly.
         vrel=fvel-myp%vel
         acc_drag=vrel/tau
         acc=acc_drag
         if (present(stress_div)) then
            do c=1,3
               acc(c)=acc(c)+this%interp(lvl,myp%pos,pStress,c)/rho_p
            end do
         end if
         myp%dt=tau/real(this%nstep,WP)
         if (heat_enabled) myp%dt=min(myp%dt,0.5_WP*thermal_dt)
         ! Save co-located inputs for the local Re/Ma/Kn extrema. get_info identifies the global
         ! winners and prints only those records in MPI-rank order, avoiding interleaved traces.
         if (trim(this%drag_model).eq.'Osnes') then
            if (this%drag_diagnostics) then
               diag_record=[this%t,dt,dt_done,mydt,real(lvl,WP),myp%pos,myp%vel,fvel,vrel, &
               &            norm2(vrel),myp%d,frho,fvisc,fC,pVF,fVF,diag_Re,diag_Ma,diag_Kn, &
               &            this%gas_gamma,corr,tau]
               if (diag_Re.gt.this%Remax) this%drag_diagnostic_record(:,1)=diag_record
               if (diag_Ma.gt.this%Mamax) this%drag_diagnostic_record(:,2)=diag_record
               if (diag_Kn.gt.this%Knmax) this%drag_diagnostic_record(:,3)=diag_record
            end if
            this%Remax=max(this%Remax,diag_Re)
            this%Mamax=max(this%Mamax,diag_Ma)
            this%Knmax=max(this%Knmax,diag_Kn)
         end if
      end function get_rhs

      !> Deposit a particle source term back to the mesh
      subroutine deposit(val)
         real(WP), dimension(4), intent(in) :: val
         real(WP) :: wx,wy,wz
         real(WP), dimension(2,2,2) :: w
         integer  :: ic,jc,kc,n
         ! Cell-centered stencil indices and weights
         ic=floor((myp%pos(1)-this%amr%xlo)*dxi-0.5_WP); wx=(myp%pos(1)-this%amr%xlo)*dxi-0.5_WP-real(ic,WP)
         jc=floor((myp%pos(2)-this%amr%ylo)*dyi-0.5_WP); wy=(myp%pos(2)-this%amr%ylo)*dyi-0.5_WP-real(jc,WP)
         kc=floor((myp%pos(3)-this%amr%zlo)*dzi-0.5_WP); wz=(myp%pos(3)-this%amr%zlo)*dzi-0.5_WP-real(kc,WP)
         ! Clamp stencil at wall boundaries
         if (this%lo_bc(1).eq.AMRLPT_WALL.and.ic  .lt.this%amr%geom(lvl)%domain%lo(1)) then; ic=this%amr%geom(lvl)%domain%lo(1)  ; wx=0.0_WP; end if
         if (this%hi_bc(1).eq.AMRLPT_WALL.and.ic+1.gt.this%amr%geom(lvl)%domain%hi(1)) then; ic=this%amr%geom(lvl)%domain%hi(1)-1; wx=1.0_WP; end if
         if (this%lo_bc(2).eq.AMRLPT_WALL.and.jc  .lt.this%amr%geom(lvl)%domain%lo(2)) then; jc=this%amr%geom(lvl)%domain%lo(2)  ; wy=0.0_WP; end if
         if (this%hi_bc(2).eq.AMRLPT_WALL.and.jc+1.gt.this%amr%geom(lvl)%domain%hi(2)) then; jc=this%amr%geom(lvl)%domain%hi(2)-1; wy=1.0_WP; end if
         if (this%lo_bc(3).eq.AMRLPT_WALL.and.kc  .lt.this%amr%geom(lvl)%domain%lo(3)) then; kc=this%amr%geom(lvl)%domain%lo(3)  ; wz=0.0_WP; end if
         if (this%hi_bc(3).eq.AMRLPT_WALL.and.kc+1.gt.this%amr%geom(lvl)%domain%hi(3)) then; kc=this%amr%geom(lvl)%domain%hi(3)-1; wz=1.0_WP; end if
         ! Add to source - one trilinear weight set, applied to every component
         w=reshape([(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz),wx*(1.0_WP-wy)*(1.0_WP-wz), &
         &          (1.0_WP-wx)*        wy *(1.0_WP-wz),wx*        wy *(1.0_WP-wz), &
         &          (1.0_WP-wx)*(1.0_WP-wy)*        wz ,wx*(1.0_WP-wy)*        wz , &
         &          (1.0_WP-wx)*        wy *        wz ,wx*        wy *        wz ],[2,2,2])
         do n=1,4
            pSrc(ic:ic+1,jc:jc+1,kc:kc+1,n)=pSrc(ic:ic+1,jc:jc+1,kc:kc+1,n)+w*val(n)
         end do
      end subroutine deposit

   end subroutine step

   !> Update particle volume fraction field based on our current particles
   subroutine update_VF(this)
      use amrex_amr_module, only: amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy,amrex_multifab,amrex_multifab_build,amrex_multifab_destroy,amrex_distromap
      use amrex_distromap_module, only: operator(.eq.)
      use mathtools, only: Pi
      implicit none
      class(amrlpt), intent(inout) :: this
      type(amrex_mfiter) :: mfi
      type(part), dimension(:), pointer :: p
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF
      integer(I8) :: np_
      integer :: lvl,i,ii,jj,kk
      real(WP) :: dxi,dyi,dzi,Vp,wx,wy,wz
      type(amrex_multifab) :: tmpVF
      logical :: dual_dm

      ! Start VF timer
      call this%tmr_vf%start()

      ! Zero VF on all levels
      call this%VF%setval(0.0_WP)
      ! Loop over levels
      do lvl=0,this%amr%clvl()
         ! Get inverse of mesh size
         dxi=1.0_WP/this%amr%dx(lvl)
         dyi=1.0_WP/this%amr%dy(lvl)
         dzi=1.0_WP/this%amr%dz(lvl)
         ! Check if we need dual-grid parallel_copy bridge
         dual_dm=(.not.(this%VF%mf(lvl)%dm.eq.this%get_particle_dm(lvl)))
         if (dual_dm) then
            call amrex_multifab_build(mf=tmpVF,ba=this%amr%ba(lvl),dm=this%get_particle_dm(lvl),nc=this%VF%mf(lvl)%ncomp(),ng=this%VF%mf(lvl)%nghost(),nodal=this%VF%nodal)
            call tmpVF%setval(0.0_WP)
         end if
         ! Loop over tiles
         call this%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointer to data
            if (dual_dm) then
               pVF=>tmpVF%dataptr(mfi)
            else
               pVF=>this%VF%mf(lvl)%dataptr(mfi)
            end if
            ! Loop over particles
            call this%get_particles(lvl=lvl,mfi=mfi,p=p,np=np_)
            do i=1,np_
               ! Skip particles that don't exchange with the fluid
               if (IAND(p(i)%flag,PART_EXCHANGES).eq.0) cycle
               ! Get particle volume
               Vp=Pi/6.0_WP*p(i)%d**3
               ! Get indices and weights
               ii=floor((p(i)%pos(1)-this%amr%xlo)*dxi-0.5_WP); wx=(p(i)%pos(1)-this%amr%xlo)*dxi-0.5_WP-real(ii,WP)
               jj=floor((p(i)%pos(2)-this%amr%ylo)*dyi-0.5_WP); wy=(p(i)%pos(2)-this%amr%ylo)*dyi-0.5_WP-real(jj,WP)
               kk=floor((p(i)%pos(3)-this%amr%zlo)*dzi-0.5_WP); wz=(p(i)%pos(3)-this%amr%zlo)*dzi-0.5_WP-real(kk,WP)
               ! Clamp stencil at wall boundaries
               if (this%lo_bc(1).eq.AMRLPT_WALL.and.ii  .lt.this%amr%geom(lvl)%domain%lo(1)) then; ii=this%amr%geom(lvl)%domain%lo(1)  ; wx=0.0_WP; end if
               if (this%hi_bc(1).eq.AMRLPT_WALL.and.ii+1.gt.this%amr%geom(lvl)%domain%hi(1)) then; ii=this%amr%geom(lvl)%domain%hi(1)-1; wx=1.0_WP; end if
               if (this%lo_bc(2).eq.AMRLPT_WALL.and.jj  .lt.this%amr%geom(lvl)%domain%lo(2)) then; jj=this%amr%geom(lvl)%domain%lo(2)  ; wy=0.0_WP; end if
               if (this%hi_bc(2).eq.AMRLPT_WALL.and.jj+1.gt.this%amr%geom(lvl)%domain%hi(2)) then; jj=this%amr%geom(lvl)%domain%hi(2)-1; wy=1.0_WP; end if
               if (this%lo_bc(3).eq.AMRLPT_WALL.and.kk  .lt.this%amr%geom(lvl)%domain%lo(3)) then; kk=this%amr%geom(lvl)%domain%lo(3)  ; wz=0.0_WP; end if
               if (this%hi_bc(3).eq.AMRLPT_WALL.and.kk+1.gt.this%amr%geom(lvl)%domain%hi(3)) then; kk=this%amr%geom(lvl)%domain%hi(3)-1; wz=1.0_WP; end if
               ! Add to VF
               pVF(ii:ii+1,jj:jj+1,kk:kk+1,1)=pVF(ii:ii+1,jj:jj+1,kk:kk+1,1)+Vp*reshape([(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz),wx*(1.0_WP-wy)*(1.0_WP-wz),(1.0_WP-wx)*wy*(1.0_WP-wz),wx*wy*(1.0_WP-wz),(1.0_WP-wx)*(1.0_WP-wy)*wz,wx*(1.0_WP-wy)*wz,(1.0_WP-wx)*wy*wz,wx*wy*wz],[2,2,2])
            end do
         end do
         call amrex_mfiter_destroy(mfi)
         ! If dual DM, parallel_copy VF back to Euler DM
         if (dual_dm) then
            call this%VF%mf(lvl)%parallel_copy(tmpVF,1,1,this%VF%mf(lvl)%ncomp(),this%VF%mf(lvl)%nghost(),this%VF%mf(lvl)%nghost(),this%amr%geom(lvl))
            call amrex_multifab_destroy(tmpVF)
         end if
      end do
      ! Post-process deposit
      call this%process_deposit(this%VF)
      ! Fill ghost cells
      call this%VF%fill(time=0.0_WP)
      ! Filter
      call this%filter(this%VF)

      ! Stop VF timer
      call this%tmr_vf%stop()

   end subroutine update_VF

   !> CFL based on particle velocities relative to their local cell size
   subroutine get_cfl(this,dt,cflc,cfl)
      use amrex_amr_module, only: amrex_mfiter
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(amrlpt), intent(inout) :: this
      real(WP), intent(in)  :: dt
      real(WP), intent(out) :: cflc
      real(WP), optional    :: cfl
      type(amrex_mfiter) :: mfi
      type(part), dimension(:), pointer :: p
      integer(I8) :: np_
      integer :: lvl,n,ierr
      ! Initialize
      this%CFLp_x=0.0_WP; this%CFLp_y=0.0_WP; this%CFLp_z=0.0_WP; this%CFL_col=0.0_WP
      ! Loop over levels and tiles
      do lvl=0,this%amr%clvl()
         call this%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Loop over particles
            call this%get_particles(lvl,mfi,p,np_)
            do n=1,np_
               ! Skip dead particles
               if (p(n)%flag.eq.PART_IS_DEAD) cycle
               ! Get CFL numbers
               this%CFLp_x =max(this%CFLp_x ,abs(p(n)%vel(1))/this%amr%dx(lvl))
               this%CFLp_y =max(this%CFLp_y ,abs(p(n)%vel(2))/this%amr%dy(lvl))
               this%CFLp_z =max(this%CFLp_z ,abs(p(n)%vel(3))/this%amr%dz(lvl))
               this%CFL_col=max(this%CFL_col,norm2(p(n)%vel)/p(n)%d)
            end do
         end do
         call this%mfiter_destroy(mfi)
      end do
      ! Global MPI max, scale by dt (10 dt for collision CFL to ensure max CFL of 0.1)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLp_x, 1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr); this%CFLp_x=this%CFLp_x*dt
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLp_y, 1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr); this%CFLp_y=this%CFLp_y*dt
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFLp_z, 1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr); this%CFLp_z=this%CFLp_z*dt
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%CFL_col,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr); this%CFL_col=10.0_WP*this%CFL_col*dt
      ! Return convective CFL and optionally overall CFL
      cflc=max(this%CFLp_x,this%CFLp_y,this%CFLp_z)
      if (present(cfl)) cfl=max(cflc,this%CFL_col)
   end subroutine get_cfl

   ! ============================================================================
   ! UTILITIES
   ! ============================================================================

   !> Redistribute particles
   subroutine redistribute(this,minlvl,maxlvl,nover)
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in), optional :: minlvl,maxlvl,nover
      integer :: lmin,lmax,no
      lmin=0;           if (present(minlvl)) lmin=minlvl
      ! Default to the hierarchy as it stands, not to the configured ceiling: AMReX asserts
      ! lev_max <= finestLevel(), which this%maxlvl violates whenever the levels above the
      ! current finest have not been built yet (e.g. right after init_from_scratch)
      lmax=min(this%maxlvl,this%amr%clvl()); if (present(maxlvl)) lmax=maxlvl
      no  =0;           if (present(nover))  no  =nover
      call amrlpt_redistribute(this%pc,lmin,lmax,no)
      call amrlpt_total_np(this%pc,this%np)
   end subroutine redistribute

   !> Enforce the configured particle level cap, including for particles that
   !> already reside above it. This stronger operation is intended for restart
   !> and regrid boundaries; ordinary per-step redistribution only needs to
   !> visit levels [0,maxlvl] once this invariant has been established.
   subroutine enforce_max_level(this)
      implicit none
      class(amrlpt), intent(inout) :: this
      integer :: cap
      cap=min(this%maxlvl,this%amr%clvl())
      call amrlpt_enforce_max_level(this%pc,cap,this%np_demoted)
      call amrlpt_total_np(this%pc,this%np)
   end subroutine enforce_max_level

   !> Return the global number of real particles stored on each AMR level.
   subroutine get_level_counts(this,counts)
      use messager, only: die
      implicit none
      class(amrlpt), intent(in) :: this
      integer(I8), dimension(0:), intent(out) :: counts
      integer :: nlevels
      nlevels=this%amr%clvl()+1
      if (size(counts).ne.nlevels) call die('[amrlpt get_level_counts] output extent must equal clvl()+1')
      call amrlpt_level_counts(this%pc,counts,nlevels)
   end subroutine get_level_counts

   !> Fill ghost particle buffer within no cells
   subroutine fill_ghosts(this,no)
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in) :: no
      call amrlpt_fill_neighbors(this%pc,no)
   end subroutine fill_ghosts

   !> Fill ghost particles, but only communicate particles within radius of tile boundaries
   subroutine fill_ghosts_radius(this,radius)
      use iso_c_binding, only: c_double
      implicit none
      class(amrlpt), intent(inout) :: this
      real(c_double), intent(in) :: radius
      call amrlpt_fill_neighbors_radius(this%pc,radius)
   end subroutine fill_ghosts_radius

   !> Release ghost particle buffer
   subroutine clear_ghosts(this)
      implicit none
      class(amrlpt), intent(inout) :: this
      call amrlpt_clear_neighbors(this%pc)
   end subroutine clear_ghosts

   !> Build explicit pair list within interaction radius rcrit
   subroutine build_neighbor_list(this,rcrit)
      use iso_c_binding, only: c_double
      implicit none
      class(amrlpt), intent(inout) :: this
      real(WP), intent(in) :: rcrit
      call amrlpt_build_neighbor_list(this%pc,real(rcrit,c_double))
   end subroutine build_neighbor_list

   !> Update global particle count
   subroutine get_np(this)
      implicit none
      class(amrlpt), intent(inout) :: this
      integer(c_int64_t) :: np_c
      call amrlpt_total_np(this%pc,np_c)
      this%np=np_c
   end subroutine get_np

   !> Get pointer to particle array for a given tile
   subroutine get_particles(this,lvl,mfi,p,np)
      use amrex_amr_module, only: amrex_mfiter
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(in) :: mfi
      type(part), dimension(:), pointer, intent(out) :: p
      integer(I8), intent(out) :: np
      type(c_ptr) :: dp
      integer(c_int64_t) :: np_c
      call amrlpt_get_particles_mfi(this%pc,lvl,mfi%p,dp,np_c)
      np=np_c
      if (np.gt.0) then
         call c_f_pointer(dp,p,[np])
      else
         nullify(p)
      end if
   end subroutine get_particles

   !> Get pointer to ghost particle array for a given tile (call after fill_ghosts)
   subroutine get_ghosts(this,lvl,mfi,g,ng)
      use amrex_amr_module, only: amrex_mfiter
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(in) :: mfi
      type(part), dimension(:), pointer, intent(out) :: g
      integer(I8), intent(out) :: ng
      type(c_ptr) :: dp
      integer(c_int64_t) :: np_c
      call amrlpt_get_neighbor_particles_mfi(this%pc,lvl,mfi%p,dp,np_c)
      ng=np_c
      if (ng.gt.0) then
         call c_f_pointer(dp,g,[ng])
      else
         nullify(g)
      end if
   end subroutine get_ghosts

   !> Get combined valid+ghost particle array for a tile (call after fill_ghosts)
   subroutine get_all_particles(this,lvl,mfi,p,np_v,np_all)
      use amrex_amr_module, only: amrex_mfiter
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(in) :: mfi
      type(part), dimension(:), pointer, intent(out) :: p
      integer(I8), intent(out) :: np_v,np_all
      type(c_ptr) :: dp
      integer(c_int64_t) :: np_total_c,np_valid_c
      call amrlpt_get_all_particles_mfi(this%pc,lvl,mfi%p,dp,np_total_c,np_valid_c)
      np_v=np_valid_c; np_all=np_total_c
      if (np_total_c.gt.0) then
         call c_f_pointer(dp,p,[np_total_c])
      else
         nullify(p)
      end if
   end subroutine get_all_particles

   !> Get neighbor list pointers for a given tile (call after build_neighbor_list)
   subroutine get_neighbor_list(this,lvl,mfi,off,lst)
      use amrex_amr_module, only: amrex_mfiter
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(in) :: mfi
      integer(c_int32_t), dimension(:), pointer, intent(out) :: off,lst
      type(c_ptr) :: dp_off,dp_lst
      integer(c_int64_t) :: np_c,ntot_c
      call amrlpt_get_neighbor_list_mfi(this%pc,lvl,mfi%p,dp_off,dp_lst,np_c,ntot_c)
      if (ntot_c.gt.0) then
         call c_f_pointer(dp_off,off,[np_c+1_8])
         call c_f_pointer(dp_lst,lst,[ntot_c])
      else
         nullify(off); nullify(lst)
      end if
   end subroutine get_neighbor_list

   !> Interpolate from cell-centered data to a point
   function interp(this,lvl,pos,arr,comp) result(val)
      implicit none
      class(amrlpt), intent(in) :: this
      integer,  intent(in) :: lvl
      real(WP), dimension(3), intent(in) :: pos
      real(WP), dimension(:,:,:,:), contiguous, pointer, intent(in) :: arr
      integer,  intent(in) :: comp
      real(WP) :: val,wx,wy,wz
      integer  :: ii,jj,kk
      ! Get indices and weights
      ii=floor((pos(1)-this%amr%xlo)/this%amr%dx(lvl)-0.5_WP); wx=(pos(1)-this%amr%xlo)/this%amr%dx(lvl)-0.5_WP-real(ii,WP)
      jj=floor((pos(2)-this%amr%ylo)/this%amr%dy(lvl)-0.5_WP); wy=(pos(2)-this%amr%ylo)/this%amr%dy(lvl)-0.5_WP-real(jj,WP)
      kk=floor((pos(3)-this%amr%zlo)/this%amr%dz(lvl)-0.5_WP); wz=(pos(3)-this%amr%zlo)/this%amr%dz(lvl)-0.5_WP-real(kk,WP)
      ! Trilinear interpolation
      val=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-wz)*arr(ii  ,jj  ,kk  ,comp) &
      &  +        wx *(1.0_WP-wy)*(1.0_WP-wz)*arr(ii+1,jj  ,kk  ,comp) &
      &  +(1.0_WP-wx)*        wy *(1.0_WP-wz)*arr(ii  ,jj+1,kk  ,comp) &
      &  +        wx *        wy *(1.0_WP-wz)*arr(ii+1,jj+1,kk  ,comp) &
      &  +(1.0_WP-wx)*(1.0_WP-wy)*        wz *arr(ii  ,jj  ,kk+1,comp) &
      &  +        wx *(1.0_WP-wy)*        wz *arr(ii+1,jj  ,kk+1,comp) &
      &  +(1.0_WP-wx)*        wy *        wz *arr(ii  ,jj+1,kk+1,comp) &
      &  +        wx *        wy *        wz *arr(ii+1,jj+1,kk+1,comp)
   end function interp

   !> Interpolate face velocities to a point
   function interp_face_velocities(this,lvl,pos,pU,uc,pV,vc,pW,wc) result(fvel)
      implicit none
      class(amrlpt), intent(in) :: this
      integer,  intent(in) :: lvl
      real(WP), dimension(3), intent(in) :: pos
      real(WP), dimension(:,:,:,:), contiguous, pointer, intent(in) :: pU,pV,pW
      integer,  intent(in) :: uc,vc,wc
      real(WP), dimension(3) :: fvel
      real(WP) :: dxi,dyi,dzi,wx,wy,wz,fx,fy,fz
      integer  :: ic,jc,kc,ix,iy,iz
      ! Get cell-centered indices and weights
      dxi=1.0_WP/this%amr%dx(lvl); ic=floor((pos(1)-this%amr%xlo)*dxi-0.5_WP); wx=(pos(1)-this%amr%xlo)*dxi-0.5_WP-real(ic,WP)
      dyi=1.0_WP/this%amr%dy(lvl); jc=floor((pos(2)-this%amr%ylo)*dyi-0.5_WP); wy=(pos(2)-this%amr%ylo)*dyi-0.5_WP-real(jc,WP)
      dzi=1.0_WP/this%amr%dz(lvl); kc=floor((pos(3)-this%amr%zlo)*dzi-0.5_WP); wz=(pos(3)-this%amr%zlo)*dzi-0.5_WP-real(kc,WP)
      ! Get face-centered indices and weights
      ix=floor((pos(1)-this%amr%xlo)*dxi); fx=(pos(1)-this%amr%xlo)*dxi-real(ix,WP)
      iy=floor((pos(2)-this%amr%ylo)*dyi); fy=(pos(2)-this%amr%ylo)*dyi-real(iy,WP)
      iz=floor((pos(3)-this%amr%zlo)*dzi); fz=(pos(3)-this%amr%zlo)*dzi-real(iz,WP)
      ! U: staggered in x — cell indices jc/wy, kc/wz; face index ix/fx
      fvel(1)=(1.0_WP-fx)*(1.0_WP-wy)*(1.0_WP-wz)*pU(ix  ,jc  ,kc  ,uc)+fx*(1.0_WP-wy)*(1.0_WP-wz)*pU(ix+1,jc  ,kc  ,uc) &
      &      +(1.0_WP-fx)*        wy *(1.0_WP-wz)*pU(ix  ,jc+1,kc  ,uc)+fx*        wy *(1.0_WP-wz)*pU(ix+1,jc+1,kc  ,uc) &
      &      +(1.0_WP-fx)*(1.0_WP-wy)*        wz *pU(ix  ,jc  ,kc+1,uc)+fx*(1.0_WP-wy)*        wz *pU(ix+1,jc  ,kc+1,uc) &
      &      +(1.0_WP-fx)*        wy *        wz *pU(ix  ,jc+1,kc+1,uc)+fx*        wy *        wz *pU(ix+1,jc+1,kc+1,uc)
      ! V: staggered in y — cell indices ic/wx, kc/wz; face index iy/fy
      fvel(2)=(1.0_WP-wx)*(1.0_WP-fy)*(1.0_WP-wz)*pV(ic  ,iy  ,kc  ,vc)+wx*(1.0_WP-fy)*(1.0_WP-wz)*pV(ic+1,iy  ,kc  ,vc) &
      &      +(1.0_WP-wx)*        fy *(1.0_WP-wz)*pV(ic  ,iy+1,kc  ,vc)+wx*        fy *(1.0_WP-wz)*pV(ic+1,iy+1,kc  ,vc) &
      &      +(1.0_WP-wx)*(1.0_WP-fy)*        wz *pV(ic  ,iy  ,kc+1,vc)+wx*(1.0_WP-fy)*        wz *pV(ic+1,iy  ,kc+1,vc) &
      &      +(1.0_WP-wx)*        fy *        wz *pV(ic  ,iy+1,kc+1,vc)+wx*        fy *        wz *pV(ic+1,iy+1,kc+1,vc)
      ! W: staggered in z — cell indices ic/wx, jc/wy; face index iz/fz
      fvel(3)=(1.0_WP-wx)*(1.0_WP-wy)*(1.0_WP-fz)*pW(ic  ,jc  ,iz  ,wc)+wx*(1.0_WP-wy)*(1.0_WP-fz)*pW(ic+1,jc  ,iz  ,wc) &
      &      +(1.0_WP-wx)*        wy *(1.0_WP-fz)*pW(ic  ,jc+1,iz  ,wc)+wx*        wy *(1.0_WP-fz)*pW(ic+1,jc+1,iz  ,wc) &
      &      +(1.0_WP-wx)*(1.0_WP-wy)*        fz *pW(ic  ,jc  ,iz+1,wc)+wx*(1.0_WP-wy)*        fz *pW(ic+1,jc  ,iz+1,wc) &
      &      +(1.0_WP-wx)*        wy *        fz *pW(ic  ,jc+1,iz+1,wc)+wx*        wy *        fz *pW(ic+1,jc+1,iz+1,wc)
   end function interp_face_velocities

   !> Post-process a deposited extensive field into an intensive composite field
   subroutine process_deposit(this,A)
      use amrex_amr_module, only: amrex_multifab,amrex_multifab_build,amrex_multifab_destroy
      use amrex_interface,  only: amrmfab_sum_downto,amrmfab_interp_from_coarse
      implicit none
      class(amrlpt), intent(inout) :: this
      type(amrdata), intent(inout) :: A
      type(amrex_multifab), dimension(:), allocatable :: tmp
      integer :: lvl,dlvl
      ! Deposit ceiling: particles are capped at this%maxlvl, so levels above it hold NO deposit.
      ! Every pass below must stop there. Walking up to clvl() instead would interpolate from, and
      ! then average back down from, levels that are identically zero - which silently erases the
      ! whole deposit wherever finer levels cover it. Reduces to clvl() when maxlvl is not capped.
      dlvl=min(this%amr%clvl(),this%maxlvl)   ! = clvl unless a caller has capped placement
      ! Convert extensive deposits to intensive
      do lvl=0,this%amr%clvl()
         call A%mf(lvl)%mult(1.0_WP/this%amr%cell_vol(lvl),1,A%ncomp,A%ng)
      end do
      ! Allocate scratch MFs for coarse→fine interpolation
      allocate(tmp(0:this%amr%clvl()))
      do lvl=0,this%amr%clvl()
         call this%amr%mfab_build(lvl,tmp(lvl),ncomp=A%ncomp,nover=0); call tmp(lvl)%setval(0.0_WP)
      end do
      ! Forward pass (coarse to fine)
      do lvl=0,dlvl
         ! SumBoundary: ghost→valid at this level
         call A%syncsum_lvl(lvl)
         ! InterpFromCoarseLevel: propagate coarse deposit to fine tmp
         if (lvl.lt.dlvl) then
            call amrmfab_interp_from_coarse(tmp(lvl+1),A%mf(lvl),[this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],cgeom=this%amr%geom(lvl),fgeom=this%amr%geom(lvl+1),scomp=1,ncomp=A%ncomp)
         end if
         ! Sum_fine_to_coarse: fine→coarse ADD
         if (lvl.gt.0) then
            call amrmfab_sum_downto(A%mf(lvl),A%mf(lvl-1),[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1),fgeom=this%amr%geom(lvl))
         end if
         ! Add interpolated coarse contribution to this level
         call A%mf(lvl)%add(tmp(lvl),1,1,A%ncomp,0)
      end do
      ! Backward pass: average_down fixes double-counting in covered cells. Starts at dlvl-1, not
      ! clvl()-1: levels above dlvl are empty and averaging them down would wipe the deposit.
      do lvl=dlvl-1,0,-1
         call A%average_downto(lvl)
      end do
      ! The AMR fluid advances on the composite hierarchy, so a source deposited on a
      ! particle-coupling level must also exist in every finer valid cell.  This is a
      ! transfer operation, not a smoothing operation: it must happen even when the
      ! optional filter is disabled.  In particular, a capped particle level combined
      ! with filter_width=0 previously left all levels above dlvl at zero; their fine
      ! cells then replaced the covered coarse source in the fluid solve.
      do lvl=dlvl+1,this%amr%clvl()
         call A%fill_from_coarse(lvl,0.0_WP)
      end do
      ! Clean up scratch MFs
      do lvl=0,this%amr%clvl()
         call amrex_multifab_destroy(tmp(lvl))
      end do
      deallocate(tmp)
   end subroutine process_deposit

   !> Explicit diffusion filter for a cell-centered amrdata field
   subroutine filter(this,A)
      use amrex_amr_module, only: amrex_box,amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy,amrex_multifab,amrex_multifab_destroy
      use amrex_interface,  only: amrmfab_average_down_face
      implicit none
      class(amrlpt), intent(inout) :: this
      type(amrdata), intent(inout) :: A
      real(WP) :: alpha,alpha_step
      integer  :: nstep,n,nc,lvl,i,j,k,flvl
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      type(amrex_multifab), dimension(:), allocatable :: Fx,Fy,Fz
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pA,pFx,pFy,pFz
      real(WP) :: dxi,dyi,dzi

      ! Determine filter coefficient from finest level mesh size
      ! Choose the level to diffuse on. The explicit sub-step count is 6*alpha/dx^2 with
      ! alpha ~ filter_width^2, i.e. ~0.54*(filter_width/dx)^2 - so taking dx at the FINEST level
      ! makes the cost grow 4x per AMR level even though filter_width is a fixed physical length.
      ! The deposit is meant to be smooth on that length anyway, so diffuse it on the coarsest
      ! level that still resolves it with filter_ncell cells, then interpolate back up. On a
      ! 6-level hierarchy with filter_width = 17.4 finest cells this goes from 163 sub-steps over
      ! 7 levels to 10 over 5.
      ! Start from the deposit ceiling, not clvl(): there is nothing to filter above this%maxlvl.
      ! The particle placement cap is also the finest level carrying an original
      ! deposit. An optional physical-width filter may choose a coarser level.
      flvl=min(this%amr%clvl(),this%maxlvl)
      do while (flvl.gt.0)
         if (this%filter_width/this%amr%min_meshsize(flvl-1).lt.real(this%filter_ncell,WP)) exit
         flvl=flvl-1
      end do

      alpha=max(this%filter_width**2-this%amr%min_meshsize(flvl)**2,0.0_WP)/(16.0_WP*log(2.0_WP))
      if (alpha.le.0.0_WP) return

      ! Number of explicit sub-steps for stability (dt < dx^2/6 in 3D)
      nstep=ceiling(6.0_WP*alpha/this%amr%min_meshsize(flvl)**2)

      ! Bring the deposit down onto the filtering level (conservative, so nothing is lost). Starts
      ! at the deposit ceiling: above it the field is zero, and averaging that down erases the
      ! deposit rather than preserving it.
      do lvl=min(this%amr%clvl(),this%maxlvl)-1,flvl,-1
         call A%average_downto(lvl)
      end do
      alpha_step=alpha/real(nstep,WP)

      ! Allocate face flux mfabs (no ghost needed for fluxes)
      allocate(Fx(0:this%amr%maxlvl),Fy(0:this%amr%maxlvl),Fz(0:this%amr%maxlvl))
      do lvl=0,flvl
         call this%amr%mfab_build(lvl,Fx(lvl),ncomp=A%ncomp,nover=0,atface=[.true., .false.,.false.]); call Fx(lvl)%setval(0.0_WP)
         call this%amr%mfab_build(lvl,Fy(lvl),ncomp=A%ncomp,nover=0,atface=[.false.,.true., .false.]); call Fy(lvl)%setval(0.0_WP)
         call this%amr%mfab_build(lvl,Fz(lvl),ncomp=A%ncomp,nover=0,atface=[.false.,.false.,.true. ]); call Fz(lvl)%setval(0.0_WP)
      end do

      ! Explicit sub-steps
      do n=1,nstep
         ! Compute diffusive fluxes at each level
         do lvl=0,flvl
            dxi=1.0_WP/this%amr%dx(lvl)
            dyi=1.0_WP/this%amr%dy(lvl)
            dzi=1.0_WP/this%amr%dz(lvl)
            call amrex_mfiter_build(mfi,A%mf(lvl),tiling=.false.)
            do while (mfi%next())
               pA =>A%mf(lvl)%dataptr(mfi)
               pFx=>Fx(lvl)%dataptr(mfi); pFy=>Fy(lvl)%dataptr(mfi); pFz=>Fz(lvl)%dataptr(mfi)
               bx=mfi%nodaltilebox(1)
               do nc=1,A%ncomp; do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pFx(i,j,k,nc)=alpha_step*(pA(i,j,k,nc)-pA(i-1,j,k,nc))*dxi
               end do; end do; end do; end do
               bx=mfi%nodaltilebox(2)
               do nc=1,A%ncomp; do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pFy(i,j,k,nc)=alpha_step*(pA(i,j,k,nc)-pA(i,j-1,k,nc))*dyi
               end do; end do; end do; end do
               bx=mfi%nodaltilebox(3)
               do nc=1,A%ncomp; do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pFz(i,j,k,nc)=alpha_step*(pA(i,j,k,nc)-pA(i,j,k-1,nc))*dzi
               end do; end do; end do; end do
            end do
            call amrex_mfiter_destroy(mfi)
         end do
         ! Average down face fluxes for C/F conservation (finest→coarsest)
         do lvl=flvl,1,-1
            call amrmfab_average_down_face(fmf=Fx(lvl),cmf=Fx(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
            call amrmfab_average_down_face(fmf=Fy(lvl),cmf=Fy(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
            call amrmfab_average_down_face(fmf=Fz(lvl),cmf=Fz(lvl-1),rr=[this%amr%rrefx(lvl-1),this%amr%rrefy(lvl-1),this%amr%rrefz(lvl-1)],cgeom=this%amr%geom(lvl-1))
         end do
         ! Apply divergence to update A
         do lvl=0,flvl
            dxi=1.0_WP/this%amr%dx(lvl)
            dyi=1.0_WP/this%amr%dy(lvl)
            dzi=1.0_WP/this%amr%dz(lvl)
            call amrex_mfiter_build(mfi,A%mf(lvl),tiling=.false.)
            do while (mfi%next())
               pA =>A%mf(lvl)%dataptr(mfi)
               pFx=>Fx(lvl)%dataptr(mfi); pFy=>Fy(lvl)%dataptr(mfi); pFz=>Fz(lvl)%dataptr(mfi)
               bx=mfi%tilebox()
               do nc=1,A%ncomp; do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pA(i,j,k,nc)=pA(i,j,k,nc)+dxi*(pFx(i+1,j,k,nc)-pFx(i,j,k,nc))+dyi*(pFy(i,j+1,k,nc)-pFy(i,j,k,nc))+dzi*(pFz(i,j,k+1,nc)-pFz(i,j,k,nc))
               end do; end do; end do; end do
            end do
            call amrex_mfiter_destroy(mfi)
         end do
         ! Restore coarse valid cells covered by fine, and fill ghosts - within 0..flvl only
         do lvl=flvl-1,0,-1
            call A%average_downto(lvl)
         end do
         do lvl=0,flvl
            call A%fill_lvl(lvl,0.0_WP)
         end do
      end do

      ! Push the filtered field back up. Levels above flvl are replaced by interpolation, which is
      ! the intent - the source is not supposed to carry structure below the filter width.
      do lvl=flvl+1,this%amr%clvl()
         call A%fill_from_coarse(lvl,0.0_WP)
      end do

      ! Destroy flux mfabs
      do lvl=0,flvl
         call amrex_multifab_destroy(Fx(lvl))
         call amrex_multifab_destroy(Fy(lvl))
         call amrex_multifab_destroy(Fz(lvl))
      end do
      deallocate(Fx,Fy,Fz)

   end subroutine filter

   !> Allgather all particles whose position lies inside [lo,hi] across all processors
   !> Returns a newly allocated array pgather(1:ngather) containing copies from every rank
   subroutine gather_region(this,lo,hi,pgather,ngather)
      use mpi_f08, only: MPI_ALLGATHER,MPI_ALLGATHERV,MPI_BYTE,MPI_INTEGER
      use amrex_amr_module, only: amrex_mfiter
      use iso_c_binding, only: c_sizeof
      implicit none
      class(amrlpt), intent(inout) :: this
      real(WP), dimension(3), intent(in) :: lo,hi
      type(part), dimension(:), allocatable, intent(out) :: pgather
      integer(I8), intent(out) :: ngather
      ! Local
      type(amrex_mfiter) :: mfi
      type(part), dimension(:), pointer :: p
      type(part), dimension(:), allocatable :: pbuf
      type(part) :: dummy
      integer(I8) :: np_,m
      integer :: lvl,nlocal,ntotal,part_bytes,ierr,i
      integer, dimension(:), allocatable :: rcounts,rdisps
      ! Size of one particle in bytes
      part_bytes=int(c_sizeof(dummy))
      ! Pass 1: count local particles in box
      nlocal=0
      do lvl=0,this%amr%clvl()
         call this%mfiter_build(lvl,mfi)
         do while(mfi%next())
            call this%get_particles(lvl,mfi,p,np_)
            do m=1,np_
               if (p(m)%pos(1).ge.lo(1).and.p(m)%pos(1).le.hi(1).and. &
               &   p(m)%pos(2).ge.lo(2).and.p(m)%pos(2).le.hi(2).and. &
               &   p(m)%pos(3).ge.lo(3).and.p(m)%pos(3).le.hi(3)) nlocal=nlocal+1
            end do
         end do
         call this%mfiter_destroy(mfi)
      end do
      ! Pack local matches
      allocate(pbuf(max(nlocal,1)))
      nlocal=0
      do lvl=0,this%amr%clvl()
         call this%mfiter_build(lvl,mfi)
         do while(mfi%next())
            call this%get_particles(lvl,mfi,p,np_)
            do m=1,np_
               if (p(m)%pos(1).ge.lo(1).and.p(m)%pos(1).le.hi(1).and. &
               &   p(m)%pos(2).ge.lo(2).and.p(m)%pos(2).le.hi(2).and. &
               &   p(m)%pos(3).ge.lo(3).and.p(m)%pos(3).le.hi(3)) then
                  nlocal=nlocal+1
                  pbuf(nlocal)=p(m)
               end if
            end do
         end do
         call this%mfiter_destroy(mfi)
      end do
      ! Allgatherv: exchange counts
      allocate(rcounts(this%amr%nproc),rdisps(this%amr%nproc))
      call MPI_ALLGATHER(nlocal,1,MPI_INTEGER,rcounts,1,MPI_INTEGER,this%amr%comm,ierr)
      ntotal=sum(rcounts)
      ! Convert counts and displacements to bytes
      rdisps(1)=0
      do i=2,this%amr%nproc
         rdisps(i)=rdisps(i-1)+rcounts(i-1)*part_bytes
      end do
      rcounts=rcounts*part_bytes
      ! Allgatherv particles
      allocate(pgather(max(ntotal,1)))
      call MPI_ALLGATHERV(pbuf,nlocal*part_bytes,MPI_BYTE,pgather,rcounts,rdisps,MPI_BYTE,this%amr%comm,ierr)
      ngather=int(ntotal,I8)
      deallocate(pbuf,rcounts,rdisps)
   end subroutine gather_region

   !> Append a Fortran part array to the particle container (collective)
   subroutine append(this,pnew,n)
      use messager, only: die
      use mpi_f08, only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_MAX,MPI_INTEGER
      implicit none
      class(amrlpt), intent(inout) :: this
      type(part), dimension(:), allocatable, target, intent(in) :: pnew
      integer(I8), intent(in) :: n
      type(c_ptr) :: raw
      integer :: bad,ierr
      integer(I8) :: i
      raw=c_null_ptr; bad=0
      if (n.lt.0_I8) bad=1
      if (n.gt.0_I8) then
         if (.not.allocated(pnew)) then
            bad=1
         else if (int(size(pnew),I8).lt.n) then
            bad=1
         else
            do i=1,n
               if (.not.all(particle_finite([pnew(i)%d,pnew(i)%mass,pnew(i)%eint,pnew(i)%vel]))) bad=1
               if (pnew(i)%d.le.0.0_WP.or.pnew(i)%mass.le.0.0_WP) bad=1
            end do
            raw=c_loc(pnew(1))
         end if
      end if
      call MPI_ALLREDUCE(MPI_IN_PLACE,bad,1,MPI_INTEGER,MPI_MAX,this%amr%comm,ierr)
      if (bad.ne.0) call die('[amrlpt append] invalid array or particle payload: require finite velocity/eint and positive finite diameter/mass')
      call amrlpt_append_particles(this%pc,raw,int(n,c_int64_t))
   end subroutine append

   !> Build MFIter from the particle grid's BA/DM
   subroutine mfiter_build(this,lvl,mfi,tiling)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_mfiter_build
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_mfiter), intent(out) :: mfi
      logical, intent(in), optional :: tiling
      type(amrex_boxarray) :: ba
      type(amrex_distromap) :: dm
      logical :: use_tiling
      use_tiling=.false.; if (present(tiling)) use_tiling=tiling
      call amrlpt_get_particle_boxarray (this%pc,lvl,ba%p)
      call amrlpt_get_particle_distromap(this%pc,lvl,dm%p)
      call amrex_mfiter_build(mfi,ba,dm,tiling=use_tiling)
   end subroutine mfiter_build

   !> Destroy MFIter
   subroutine mfiter_destroy(this,mfi)
      use amrex_amr_module, only: amrex_mfiter,amrex_mfiter_destroy
      implicit none
      class(amrlpt), intent(inout) :: this
      type(amrex_mfiter), intent(inout) :: mfi
      call amrex_mfiter_destroy(mfi)
   end subroutine mfiter_destroy

   !> Set particle-specific BoxArray for a given level
   subroutine set_particle_ba(this,lvl,ba)
      use amrex_amr_module, only: amrex_boxarray
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_boxarray), intent(in) :: ba
      call amrlpt_set_particle_boxarray(this%pc,lvl,ba%p)
   end subroutine set_particle_ba

   !> Set particle-specific DistributionMapping for a given level
   subroutine set_particle_dm(this,lvl,dm)
      use amrex_amr_module, only: amrex_distromap
      implicit none
      class(amrlpt), intent(inout) :: this
      integer, intent(in) :: lvl
      type(amrex_distromap), intent(in) :: dm
      call amrlpt_set_particle_distromap(this%pc,lvl,dm%p)
   end subroutine set_particle_dm

   !> Get particle-specific DistributionMapping for a given level
   function get_particle_dm(this,lvl) result(dm)
      use amrex_amr_module, only: amrex_distromap
      implicit none
      class(amrlpt), intent(in) :: this
      integer, intent(in) :: lvl
      type(amrex_distromap) :: dm
      call amrlpt_get_particle_distromap(this%pc,lvl,dm%p)
      dm%owner=.false.
   end function get_particle_dm

   ! ============================================================================
   ! SOLVER INFO
   ! ============================================================================

   !> Compute particle statistics: np, d/vel min/max/mean/var, Vp_tot, and VF field stats
   subroutine get_info(this)
      use amrex_amr_module, only: amrex_multifab,amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy,amrex_box
      use iso_fortran_env, only: output_unit
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_BARRIER,MPI_SUM,MPI_MIN,MPI_MAX,MPI_IN_PLACE,MPI_INTEGER8,MPI_INTEGER
      use parallel, only: MPI_REAL_WP
      use mathtools, only: Pi
      implicit none
      class(amrlpt), intent(inout) :: this
      type(amrex_mfiter) :: mfi
      type(amrex_box)    :: bx
      type(part), dimension(:), pointer :: p
      real(WP), dimension(:,:,:,:), pointer :: pVF
      integer(I8) :: np_
      integer :: lvl,n,i,j,k,ierr
      real(WP) :: d_sum,d_sq,vx_sum,vx_sq,vy_sum,vy_sq,vz_sum,vz_sq,inv_np,re_loc,ma_loc,kn_loc
      real(WP) :: inventory(6)
      ! Init per-rank accumulators
      inventory=0.0_WP
      this%np=0; this%Vp_tot=0.0_WP
      this%dmin=huge(1.0_WP); this%dmax=-huge(1.0_WP);  d_sum=0.0_WP;  d_sq=0.0_WP
      this%Umin=huge(1.0_WP); this%Umax=-huge(1.0_WP); vx_sum=0.0_WP; vx_sq=0.0_WP
      this%Vmin=huge(1.0_WP); this%Vmax=-huge(1.0_WP); vy_sum=0.0_WP; vy_sq=0.0_WP
      this%Wmin=huge(1.0_WP); this%Wmax=-huge(1.0_WP); vz_sum=0.0_WP; vz_sq=0.0_WP
      ! Loop over all AMR levels and tiles
      do lvl=0,this%amr%clvl()
         call this%mfiter_build(lvl,mfi)
         do while (mfi%next())
            call this%get_particles(lvl,mfi,p,np_)
            do n=1,np_
               if (p(n)%flag.eq.PART_IS_DEAD) cycle
               this%np=this%np+1
               this%Vp_tot=this%Vp_tot+Pi/6.0_WP*p(n)%d**3
               inventory=inventory+[p(n)%mass,p(n)%mass*p(n)%vel, &
                  0.5_WP*p(n)%mass*dot_product(p(n)%vel,p(n)%vel),p(n)%mass*p(n)%eint]
               this%dmin=min(this%dmin,p(n)%d);      this%dmax=max(this%dmax,p(n)%d);       d_sum= d_sum+p(n)%d;       d_sq= d_sq+p(n)%d**2
               this%Umin=min(this%Umin,p(n)%vel(1)); this%Umax=max(this%Umax,p(n)%vel(1)); vx_sum=vx_sum+p(n)%vel(1); vx_sq=vx_sq+p(n)%vel(1)**2
               this%Vmin=min(this%Vmin,p(n)%vel(2)); this%Vmax=max(this%Vmax,p(n)%vel(2)); vy_sum=vy_sum+p(n)%vel(2); vy_sq=vy_sq+p(n)%vel(2)**2
               this%Wmin=min(this%Wmin,p(n)%vel(3)); this%Wmax=max(this%Wmax,p(n)%vel(3)); vz_sum=vz_sum+p(n)%vel(3); vz_sq=vz_sq+p(n)%vel(3)**2
            end do
         end do
         call this%mfiter_destroy(mfi)
      end do
      ! Compute per-rank particle distribution stats
      this%np_loc=this%np; this%np_min=this%np; this%np_max=this%np; this%np_eff=0.0_WP
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%np_min,1,MPI_INTEGER8,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%np_max,1,MPI_INTEGER8,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%np,    1,MPI_INTEGER8,MPI_SUM,this%amr%comm,ierr)
      if (this%np_max.gt.0) this%np_eff=real(this%np,WP)/real(this%np_max,WP)/real(this%amr%nproc,WP)
      call MPI_ALLREDUCE(MPI_IN_PLACE,inventory,6,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
      this%mass_tot=inventory(1); this%momentum_tot=inventory(2:4)
      this%kinetic_tot=inventory(5); this%internal_tot=inventory(6)
      ! Reduce other statistics
      re_loc=this%Remax; ma_loc=this%Mamax; kn_loc=this%Knmax
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Remax, 1,MPI_REAL_WP ,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Mamax, 1,MPI_REAL_WP ,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Knmax, 1,MPI_REAL_WP ,MPI_MAX,this%amr%comm,ierr)
      ! Re/Ma/Kn maxima are independent. When requested, print the full co-located state behind
      ! each global maximum so a monitor value can be reconstructed without ambiguity.
      if (this%drag_diagnostic_pending) then
         if (max(this%Remax,this%Mamax,this%Knmax).gt.0.0_WP) then
            if (this%amr%amRoot) write(output_unit,'(a)') &
            &  'AMRLPT_DRAG_EXTREME columns: kind rank t_lpt step_dt particle_dt_done particle_mydt level x y z vp_x vp_y vp_z ug_x ug_y ug_z slip_x slip_y slip_z slip_mag d rho_g mu c pVF fVF Re Ma Kn gamma corr tau'
            call MPI_BARRIER(this%amr%comm,ierr)
            do i=0,this%amr%nproc-1
               if (this%amr%rank.eq.i) then
                  if (re_loc.eq.this%Remax) write(output_unit,'("AMRLPT_DRAG_EXTREME ReMax ",i0,1x,*(es24.16e3,1x))') this%amr%rank,this%drag_diagnostic_record(:,1)
                  if (ma_loc.eq.this%Mamax) write(output_unit,'("AMRLPT_DRAG_EXTREME MaMax ",i0,1x,*(es24.16e3,1x))') this%amr%rank,this%drag_diagnostic_record(:,2)
                  if (kn_loc.eq.this%Knmax) write(output_unit,'("AMRLPT_DRAG_EXTREME KnMax ",i0,1x,*(es24.16e3,1x))') this%amr%rank,this%drag_diagnostic_record(:,3)
                  flush(output_unit)
               end if
               call MPI_BARRIER(this%amr%comm,ierr)
            end do
         end if
         this%drag_diagnostic_pending=.false.
      end if
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Vp_tot,1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%dmin,  1,MPI_REAL_WP ,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%dmax,  1,MPI_REAL_WP ,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,d_sum,      1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,d_sq,       1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Umin,  1,MPI_REAL_WP ,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Umax,  1,MPI_REAL_WP ,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vx_sum,     1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vx_sq,      1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Vmin,  1,MPI_REAL_WP ,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Vmax,  1,MPI_REAL_WP ,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vy_sum,     1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vy_sq,      1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Wmin,  1,MPI_REAL_WP ,MPI_MIN,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,this%Wmax,  1,MPI_REAL_WP ,MPI_MAX,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vz_sum,     1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vz_sq,      1,MPI_REAL_WP ,MPI_SUM,this%amr%comm,ierr)
      ! Derive mean and variance from single-pass accumulation
      if (this%np.gt.0) then
         inv_np=1.0_WP/real(this%np,WP)
         this%dmean= d_sum*inv_np; this%dvar=max(0.0_WP, d_sq*inv_np-this%dmean**2)
         this%Umean=vx_sum*inv_np; this%Uvar=max(0.0_WP,vx_sq*inv_np-this%Umean**2)
         this%Vmean=vy_sum*inv_np; this%Vvar=max(0.0_WP,vy_sq*inv_np-this%Vmean**2)
         this%Wmean=vz_sum*inv_np; this%Wvar=max(0.0_WP,vz_sq*inv_np-this%Wmean**2)
      end if
      ! Record actual occupied levels and make a violation of maxlvl directly
      ! observable to validation cases and production monitors.
      particle_level_stats: block
         integer(I8), dimension(:), allocatable :: counts
         allocate(counts(0:this%amr%clvl()))
         call this%get_level_counts(counts)
         this%occupied_minlvl=-1; this%occupied_maxlvl=-1
         do lvl=0,this%amr%clvl()
            if (counts(lvl).gt.0_I8) then
               if (this%occupied_minlvl.lt.0) this%occupied_minlvl=lvl
               this%occupied_maxlvl=lvl
            end if
         end do
         this%np_above_maxlvl=0_I8
         if (this%maxlvl.lt.this%amr%clvl()) &
         &   this%np_above_maxlvl=sum(counts(this%maxlvl+1:this%amr%clvl()))
         deallocate(counts)
      end block particle_level_stats
      ! Volume fraction statistics
      vf_stats_block: block
         use amrex_amr_module, only: amrex_imultifab,amrex_imultifab_build,amrex_imultifab_destroy
         use amrex_interface,  only: amrmask_make_fine
         type(amrex_imultifab) :: mask
         integer,  dimension(:,:,:,:), contiguous, pointer :: pMask
         real(WP) :: var_sum
         ! Initialize stats
         this%VFmean=this%VF%get_sum(lvl=0)/real(this%amr%nx*this%amr%ny*this%amr%nz,WP)
         var_sum=0.0_WP; this%VFmin=huge(1.0_WP); this%VFmax=-huge(1.0_WP)
         ! Loop over levels
         do lvl=0,this%amr%clvl()
            ! Build fine mask
            if (lvl.lt.this%amr%clvl()) then
               call amrex_imultifab_build(mask,this%amr%ba(lvl),this%amr%dm(lvl),1,0)
               call amrmask_make_fine(mask,this%amr%ba(lvl+1),[this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],0,1)
            end if
            ! Loop over tiles
            call amrex_mfiter_build(mfi,this%VF%mf(lvl))
            do while (mfi%next())
               ! Get pointer to data
               pVF=>this%VF%mf(lvl)%dataptr(mfi)
               if (lvl.lt.this%amr%clvl()) pMask=>mask%dataptr(mfi)
               ! Loop over tile
               bx=mfi%tilebox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Skip cells covered by finer level
                  if (lvl.lt.this%amr%clvl()) then
                     if (pMask(i,j,k,1).eq.0) cycle
                  end if
                  ! Accumulate statistics
                  this%VFmin=min(this%VFmin,pVF(i,j,k,1)); this%VFmax=max(this%VFmax,pVF(i,j,k,1))
                  var_sum=var_sum+(pVF(i,j,k,1)-this%VFmean)**2*this%amr%cell_vol(lvl)
               end do; end do; end do
            end do
            call amrex_mfiter_destroy(mfi)
            if (lvl.lt.this%amr%clvl()) call amrex_imultifab_destroy(mask)
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%VFmin,1,MPI_REAL_WP,MPI_MIN,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%VFmax,1,MPI_REAL_WP,MPI_MAX,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,var_sum,   1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
         this%VFvar=max(0.0_WP,var_sum/((this%amr%xhi-this%amr%xlo)*(this%amr%yhi-this%amr%ylo)*(this%amr%zhi-this%amr%zlo)))
      end block vf_stats_block

      ! Reduce local accumulators, publish to monitoring fields, then zero for next interval
      reduce_local_stats: block
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%culled_loc,5,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
         this%culled=this%culled+this%culled_loc; this%culled_loc=0.0_WP
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%np_culled_loc,1,MPI_INTEGER,MPI_SUM,this%amr%comm,ierr)
         this%np_culled=this%np_culled_loc; this%np_culled_tot=this%np_culled_tot+this%np_culled_loc
         this%np_culled_loc=0
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%heat_exchange_loc,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
         this%heat_exchange=this%heat_exchange+this%heat_exchange_loc; this%heat_exchange_loc=0.0_WP
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%exchange_loc,5,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
         this%exchange=this%exchange+this%exchange_loc; this%exchange_loc=0.0_WP
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%outflow_loc,5,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
         this%outflow=this%outflow+this%outflow_loc; this%outflow_loc=0.0_WP
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%np_out_loc,1,MPI_INTEGER,MPI_SUM,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%Vp_out_loc,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%np_new_loc,1,MPI_INTEGER,MPI_SUM,this%amr%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%Vp_new_loc,1,MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
         this%np_out=this%np_out_loc; this%Vp_out=this%Vp_out_loc
         this%np_new=this%np_new_loc; this%Vp_new=this%Vp_new_loc
         this%np_out_loc=0; this%Vp_out_loc=0.0_WP
         this%np_new_loc=0; this%Vp_new_loc=0.0_WP
      end block reduce_local_stats

      ! Compute timing statistics across ranks
      call this%tmr_coll%get_stats()
      call this%tmr_coll_%get_stats()
      call this%tmr_fill%get_stats()
      call this%tmr_nbl%get_stats()
      call this%tmr_step%get_stats()
      call this%tmr_step_%get_stats()
      call this%tmr_vf%get_stats()
      call this%tmr_src%get_stats()

      ! Reset timers
      call this%tmr_coll%reset()
      call this%tmr_coll_%reset()
      call this%tmr_fill%reset()
      call this%tmr_nbl%reset()
      call this%tmr_step%reset()
      call this%tmr_step_%reset()
      call this%tmr_vf%reset()
      call this%tmr_src%reset()

   end subroutine get_info

   !> Print solver info
   subroutine print(this)
      use messager, only: log
      use string,   only: str_long
      implicit none
      class(amrlpt), intent(in) :: this
      character(len=str_long) :: message
      call log('AMR Lagrangian particle solver: '//trim(this%name))
      write(message,'("  Particle density  : ",ES12.5)') this%rho
      call log(trim(message))
      write(message,'("  Drag model        : ",a)') trim(this%drag_model)
      call log(trim(message))
      if (this%drag_vf_override.ge.0.0_WP) then
         write(message,'("  Drag VF override  : ",ES12.5)') this%drag_vf_override
         call log(trim(message))
      else
         call log('  Drag VF           : locally deposited particle fraction')
      end if
      if (this%total_energy_coupling) then
         call log('  Energy source     : total energy (mechanical + accepted heat)')
      else
         call log('  Energy source     : legacy drag dissipation')
      end if
      if (this%cull_gas_free) call log('  Gas-free stencil  : cull last accepted state; book mass/momentum/energy loss')
      if (associated(this%heat_material)) call log('  Particle heat     : Gunn, porosity=1, EOS internal energy, dt <= tau_T/2')
      write(message,'("  Gas gamma         : ",ES12.5)') this%gas_gamma
      call log(trim(message))
      write(message,'("  ODE substeps      : ",i0)') this%nstep
      call log(trim(message))
      call log('  Grid: '//trim(this%amr%name))
   end subroutine print

   ! ============================================================================
   ! CHECKPOINT IO
   ! ============================================================================

   !> Write AMReX checkpoint for particles to dirname
   subroutine write(this,dirname)
      implicit none
      class(amrlpt), intent(inout) :: this
      character(len=*), intent(in) :: dirname
      call amrlpt_write(this%pc,trim(dirname)//'/particles'//c_null_char,1_c_int)
   end subroutine write

   !> Read AMReX checkpoint for particles from dirname.
   !> The amrgrid must already have been rebuilt via amr%init_from_checkpoint
   !> before calling this. Syncs the particle container's BA/DM to the restored
   !> grid before AMReX Restart so particles land on the right ranks.
   subroutine read(this,dirname)
      implicit none
      class(amrlpt), intent(inout) :: this
      character(len=*), intent(in) :: dirname
      integer :: lvl
      ! Sync particle container's BA/DM to the restored grid
      do lvl=0,this%amr%clvl()
         call this%set_particle_ba(lvl,this%amr%get_boxarray(lvl))
         call this%set_particle_dm(lvl,this%amr%get_distromap(lvl))
      end do
      call amrlpt_read(this%pc,trim(dirname)//'/particles'//c_null_char)
      call this%enforce_max_level()
      call this%update_VF()
   end subroutine read

end module amrlpt_class
