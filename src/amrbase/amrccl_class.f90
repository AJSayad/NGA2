!> Connected component labeling on an AMR hierarchy
!>
!> build() identifies objects from a Eulerian logical field defined on a SINGLE
!> AMR level - typically amr%maxlvl, where amrvof keeps its PLIC/barycenter data. This is the AMR
!> analogue of cclabel_class, and the union-find machinery is ported from it near-verbatim.
!>
!> Two things make the AMR version simpler than cclabel_class rather than harder:
!>  - Working at one level means there are no covered cells, hence no composite masking anywhere.
!>  - A single fill_boundary() on the label field exchanges labels across every box face at the
!>    level AND wraps periodic ghosts. That collapses cclabel's two separate stitching phases
!>    (interproc_handling for the imin_/jmin_/kmin_ faces, boundary_handling for the periodic
!>    domain edges) into one halo exchange followed by one three-face scan.
!>
!> What is NOT free is that a level generally does not cover the whole structure: amrvof refines
!> only a band around the interface, so a structure wider than that band has an unrefined core and
!> its footprint at this level is a hollow shell. Such structures are flagged struct%open, meaning
!> any quantity integrated over their cells at this level is incomplete.
!> build_composite() instead labels leaf cells across ALL levels and joins coarse/fine faces,
!> so its structures include their coarser interiors without counting covered cells twice.
!>
!> Unlike cclabel_class, no per-structure cell map is stored. Per-structure work is done by the
!> caller as MFIter loops that scatter into id-indexed accumulator arrays, which avoids an
!> O(ncells) allocation and matches how the rest of the AMR code is written.
!>
!> Cross-rank lineage is communicated as a rank-local union forest.  Rank zero resolves the
!> resulting equivalence edges once and broadcasts the final parent map; the iterative full-array
!> allreduces used by cclabel_class are deliberately avoided.
module amrccl_class
   use precision,        only: WP
   use string,           only: str_medium
   use amrgrid_class,    only: amrgrid
   use amrdata_class,    only: amrdata,interp_const
   use amrex_amr_module, only: amrex_mfiter,amrex_imultifab,amrex_box
   implicit none
   private


   ! Expose type/constructor/methods
   public :: amrccl,prepare_fab_ftype,make_label_ftype,same_label_ftype


   ! Some parameters for memory management
   integer , parameter :: min_struct_size=100 !< Default minimum size of structure storage
   real(WP), parameter :: coeff_up=1.5_WP     !< When we run out of structure storage, increase by 50%


   !> Structure object
   type :: struct_type
      integer :: parent                       !< Union-find parent during build; self-referential after
      integer :: n_                           !< Number of rank-local cells contained in struct
      integer, dimension(3) :: per            !< per(dim)=1 if structure wraps the periodic boundary in dim
      logical :: open                         !< True if struct touches the edge of the labeled level
      logical :: has_coarse=.false.            !< Composite structure contains cells below maxlvl
   end type struct_type


   !> amrccl object definition
   type :: amrccl
      ! This is our amrgrid
      class(amrgrid), pointer :: amr=>null()
      ! This is the name of the CCL
      character(len=str_medium) :: name='UNNAMED_AMRCCL'
      ! ID of the structure that contains each cell
      !> Stored as an amrdata so it can be handed to amrviz/amrio. interp_const carries the labels
      !> through a regrid by piecewise-constant interpolation: exact under refinement, one child's
      !> label under coarsening. Correctness does not depend on this - build() resets and rezeroes
      !> the field every time - but interp_none leaves the multifab ALLOCATED AND UNFILLED after a
      !> regrid, so anything reading it before the next build (a visualization write, most
      !> obviously) sees garbage. build() populates one level; build_composite() populates leaf
      !> cells at every level. Covered coarse cells are zero in the composite field.
      type(amrdata) :: id
      ! Reusable integer labeling workspace and coverage mask.  Their layouts change only on
      ! regrid, not every labeling event.  The public real-valued id field is written once after
      ! the final dense map is known, retaining its amrdata/viz interface without using floating
      ! point labels internally.
      type(amrex_imultifab) :: label,cover
      logical :: work_allocated=.false.
      logical :: work_valid=.false.
      logical :: cover_valid=.false.   !< cover holds work_lvl's coverage mask (skipped by composite builds)
      logical :: id_layout_ready=.false.
      integer :: work_lvl=-1
      ! Level the current labeling lives on (-1 if none)
      integer :: lvl=-1
      ! Array of structures, indexed 1:nstruct globally (n_ is the rank-local cell count)
      integer :: nstruct=0
      type(struct_type), dimension(:), allocatable :: struct
      ! Ensure id=1 is always the largest structure in terms of number of cells
      logical :: rename_largest=.true.
   contains
      procedure :: initialize
      procedure :: build
      procedure :: build_composite !< Face-connected leaf cells across the entire AMR hierarchy
      procedure :: empty
      procedure :: finalize
      procedure, private :: allreduce_1d
      procedure, private :: allreduce_2d
      generic   :: allreduce=>allreduce_1d,allreduce_2d
   end type amrccl


   !> User callbacks
   !> prepare_fab is called once per FAB, before that FAB is scanned, so the user can rebind its
   !> array pointers (e.g. pVF=>fs%VF%mf(lvl)%dataptr(mfi)). make_label and same_label then read
   !> those host-associated pointers. NOTE: the user's fields must live on the same BoxArray and
   !> DistributionMapping as amr at that level - true for everything in amrvof/amrmpcomp, but not
   !> for a solver running on its own dual grid.
   abstract interface
      subroutine prepare_fab_ftype(lvl,mfi)
         import :: amrex_mfiter
         integer, intent(in) :: lvl
         type(amrex_mfiter), intent(inout) :: mfi
      end subroutine prepare_fab_ftype
   end interface

   !> Type of the make_label function used to generate a structure
   abstract interface
      logical function make_label_ftype(ind1,ind2,ind3)
         integer, intent(in) :: ind1,ind2,ind3
      end function make_label_ftype
   end interface

   !> Type of the same_label function used to connect two structures
   abstract interface
      logical function same_label_ftype(ind1,ind2,ind3,ind4,ind5,ind6)
         integer, intent(in) :: ind1,ind2,ind3,ind4,ind5,ind6
      end function same_label_ftype
   end interface


contains


   !> Invalidate layout-dependent scratch data after an AMR regrid.
   subroutine amrccl_postregrid(ctx,lbase,time)
      use iso_c_binding, only: c_ptr,c_f_pointer
      implicit none
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      type(amrccl), pointer :: this
      call c_f_pointer(ctx,this)
      this%work_valid=.false.
   end subroutine amrccl_postregrid


   !> Initialization for amrccl class
   subroutine initialize(this,amr,name)
      use amrex_amr_module, only: amrex_bc_foextrap
      use iso_c_binding,    only: c_loc
      implicit none
      class(amrccl), target, intent(inout) :: this
      class(amrgrid), target, intent(in) :: amr
      character(len=*), optional :: name
      ! Set the name for the object
      if (present(name)) this%name=trim(adjustl(name))
      ! Point to amrgrid object
      this%amr=>amr
      ! Create the public label output; one ghost cell supports downstream stencil access.
      call this%id%initialize(amr=amr,name=trim(this%name)//'_id',ncomp=1,ng=1,interp=interp_const)
      ! build() stitches with fill_boundary, which applies no physical BCs at all - ghosts outside
      ! a non-periodic domain face stay zero and never union. These BCs matter only if a caller
      ! ever FillPatches the label field; foextrap is then the self-consistent choice, since a
      ! ghost repeating its interior neighbor's label makes union(x,x), a no-op.
      if (.not.amr%xper) then; this%id%lo_bc(1,1)=amrex_bc_foextrap; this%id%hi_bc(1,1)=amrex_bc_foextrap; end if
      if (.not.amr%yper) then; this%id%lo_bc(2,1)=amrex_bc_foextrap; this%id%hi_bc(2,1)=amrex_bc_foextrap; end if
      if (.not.amr%zper) then; this%id%lo_bc(3,1)=amrex_bc_foextrap; this%id%hi_bc(3,1)=amrex_bc_foextrap; end if
      call this%id%register()
      select type (this)
       type is (amrccl)
         call this%amr%add_postregrid(amrccl_postregrid,c_loc(this))
      end select
      ! Zero structures
      this%nstruct=0
      this%lvl=-1
      this%work_allocated=.false.
      this%work_valid=.false.
      this%cover_valid=.false.
      this%work_lvl=-1
      this%id_layout_ready=.false.
   end subroutine initialize


   !> Build structures at level lvl using the user-set test functions
   subroutine build(this,lvl,prepare_fab,make_label,same_label,preserve_levels)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_ALLGATHER,MPI_GATHER,MPI_GATHERV,MPI_REDUCE,MPI_BCAST, &
         MPI_SUM,MPI_MAX,MPI_INTEGER,MPI_IN_PLACE
      use messager, only: die
      implicit none
      class(amrccl), intent(inout) :: this
      integer, intent(in) :: lvl
      procedure(prepare_fab_ftype) :: prepare_fab
      procedure(make_label_ftype)  :: make_label
      procedure(same_label_ftype)  :: same_label
      logical, optional, intent(in) :: preserve_levels !< Internal composite build: retain other levels
      logical :: clear_all,composite
      integer :: nstruct_,stmin,stmax
      integer, dimension(:), allocatable :: parent             !< Resolving structure id across procs
      integer, dimension(:,:), allocatable :: prov_flags       !< Seven packed boundary flags per local provisional structure
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pId
      integer,  dimension(:,:,:,:), contiguous, pointer :: pLabel,pCov
      integer :: ierr

      ! Sanity checks
      if (.not.associated(this%amr)) call die('[amrccl build] amrccl%initialize must be called first')
      if (lvl.lt.0.or.lvl.gt.this%amr%clvl()) call die('[amrccl build] requested level does not exist')

      ! Start by cleaning up
      call this%empty()
      this%lvl=lvl

      ! A first build may occur when amrccl was initialized after init_from_scratch, so establish
      ! every label MultiFab once.  Thereafter the registered amrdata callbacks keep their layouts
      ! synchronized across regrids; destroying and rebuilding all levels on every CCL event only
      ! adds allocation traffic.  Zero all levels to retain the visualization contract.
      if (.not.this%id_layout_ready) then
         call this%id%reset()
         this%id_layout_ready=.true.
      end if
      ! preserve_levels is passed only by build_composite, which keeps nothing from this build but
      ! the level's id field and nstruct: it recomputes per/open/has_coarse itself and FillPatches
      ! every id ghost afterwards. Composite builds therefore skip the coverage mask, the boundary
      ! flags, their reduction, and the final id halo exchange. Labels are unchanged.
      composite=.false.
      if (present(preserve_levels)) composite=preserve_levels
      clear_all=.not.composite
      if (clear_all) then
         call this%id%setval(0.0_WP)
      else
         call this%id%mf(lvl)%setval(0.0_WP)
      end if

      ! Build reusable integer work fields. Coverage is 1 in valid cells and 0 elsewhere; after
      ! fill_boundary a ghost holds 1 iff it is a valid cell of another box at this level
      ! (periodic images included), and 0 iff it lies outside the level. The label field is
      ! cleared on every build but only reallocated when its AMR layout changes.
      build_work_fields: block
         use amrex_amr_module, only: amrex_imultifab_build,amrex_imultifab_destroy
         use amrex_interface,  only: amrimultifab_fill_boundary
         integer :: i,j,k
         if (.not.this%work_valid.or.this%work_lvl.ne.lvl) then
            if (this%work_allocated) then
               call amrex_imultifab_destroy(this%label)
               call amrex_imultifab_destroy(this%cover)
            end if
            call amrex_imultifab_build(this%label,this%amr%ba(lvl),this%amr%dm(lvl),1,1)
            call amrex_imultifab_build(this%cover,this%amr%ba(lvl),this%amr%dm(lvl),1,1)
            this%work_allocated=.true.
            this%work_lvl=lvl
            this%cover_valid=.false.
            this%work_valid=.true.
         end if
         ! The coverage mask only feeds the open flags, which a composite build discards. It stays
         ! allocated either way; cover_valid makes a later single-level build fill it on demand.
         if (.not.composite.and..not.this%cover_valid) then
            call this%cover%setval(0)
            call this%amr%mfiter_build(lvl,mfi,tiling=.false.)
            do while (mfi%next())
               pCov=>this%cover%dataptr(mfi)
               bx=mfi%validbox()
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pCov(i,j,k,1)=1
               end do; end do; end do
            end do
            call this%amr%mfiter_destroy(mfi)
            call amrimultifab_fill_boundary(this%cover,this%amr%geom(lvl))
            this%cover_valid=.true.
         end if
         call this%label%setval(0)
      end block build_work_fields

      ! Then allocate struct to a default size
      nstruct_=0
      allocate(this%struct(min_struct_size))
      this%struct(:)%parent=0
      this%struct(:)%n_=0
      this%struct(:)%per(1)=0
      this%struct(:)%per(2)=0
      this%struct(:)%per(3)=0
      this%struct(:)%open=.false.

      ! Perform a first pass to build box-local structures and corresponding tree. Neighbors that
      ! fall outside the current valid box are deliberately ignored here - every one of them is
      ! picked up by the stitching pass below, which is what makes this correct for a rank that
      ! owns several boxes.
      first_pass: block
         integer :: i,j,k,ii,jj,kk,dim,nme,nnb
         integer, dimension(3) :: pos
         call this%amr%mfiter_build(lvl,mfi,tiling=.false.)
         do while (mfi%next())
            call prepare_fab(lvl,mfi)
            pLabel=>this%label%dataptr(mfi)
            bx=mfi%validbox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Find next cell in a structure
               if (make_label(i,j,k)) then
                  nme=0
                  ! Loop through one-sided neighbors
                  do dim=1,3
                     pos=0; pos(dim)=-1
                     ii=i+pos(1); jj=j+pos(2); kk=k+pos(3)
                     ! Skip neighbors outside this valid box - the stitching pass handles them
                     if (ii.lt.bx%lo(1).or.jj.lt.bx%lo(2).or.kk.lt.bx%lo(3)) cycle
                     nnb=pLabel(ii,jj,kk,1)
                     if (nnb.gt.0) then
                        if (.not.same_label(i,j,k,ii,jj,kk)) cycle
                        if (nme.eq.0) then
                           nme=nnb
                        else
                           nme=union_struct(nme,nnb)
                        end if
                     end if
                  end do
                  ! If no matching neighbor was labeled, we need a new structure
                  if (nme.eq.0) nme=add()
                  pLabel(i,j,k,1)=nme
               end if
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end block first_pass

      ! Flatten the much smaller local union tree, then compact it and give every surviving root
      ! a globally unique provisional index.  Every label returned by add() is immediately placed
      ! on a cell, so each surviving root is nonempty; a separate full-grid collapse/recount pass
      ! is unnecessary.  Counts are accumulated with the boundary flags in the cell pass below.
      compact_tree: block
         integer :: i,j,k,n,nloc,dim,side,nlabels,gid
         integer, dimension(3) :: pos,nbr,dlo,dhi
         logical, dimension(3) :: per
         integer, dimension(:), allocatable :: all_nstruct,idmap
         type(struct_type), dimension(:), allocatable :: tmp
         ! Flatten local lineage and count exact number of roots.
         nlabels=nstruct_
         do n=1,nlabels
            this%struct(n)%parent=rootify_struct(n)
         end do
         nstruct_=0
         do n=1,nlabels
            if (this%struct(n)%parent.eq.n) nstruct_=nstruct_+1
         end do
         ! Gather one scalar per rank to assign globally unique provisional indices.  The previous
         ! sparse-vector allreduce communicated nproc integers from every rank for this scalar.
         allocate(all_nstruct(0:this%amr%nproc-1))
         call MPI_ALLGATHER(nstruct_,1,MPI_INTEGER,all_nstruct,1,MPI_INTEGER,this%amr%comm,ierr)
         stmin=1
         if (this%amr%rank.gt.0) stmin=stmin+sum(all_nstruct(0:this%amr%rank-1))
         this%nstruct=sum(all_nstruct)
         deallocate(all_nstruct)
         stmax=stmin+nstruct_-1
         ! Generate an index map
         allocate(idmap(max(nlabels,1))); idmap=0
         nstruct_=0
         do n=1,nlabels
            if (this%struct(n)%parent.eq.n) then
               nstruct_=nstruct_+1
               idmap(n)=stmin+nstruct_-1
            end if
         end do
         do n=1,nlabels
            if (this%struct(n)%parent.ne.n) idmap(n)=idmap(this%struct(n)%parent)
         end do
         allocate(prov_flags(7,max(nstruct_,1))); prov_flags=0
         allocate(tmp(stmin:stmax))
         tmp(:)%parent=0; tmp(:)%n_=0; tmp(:)%per(1)=0; tmp(:)%per(2)=0; tmp(:)%per(3)=0; tmp(:)%open=.false.
         dlo=this%amr%geom(lvl)%domain%lo
         dhi=this%amr%geom(lvl)%domain%hi
         per=[this%amr%xper,this%amr%yper,this%amr%zper]
         ! Map box-local labels to globally unique provisional indices.
         call this%amr%mfiter_build(lvl,mfi,tiling=.false.)
         do while (mfi%next())
            pLabel=>this%label%dataptr(mfi)
            if (.not.composite) pCov=>this%cover%dataptr(mfi)
            bx=mfi%validbox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               n=pLabel(i,j,k,1)
               if (n.gt.0) then
                  gid=idmap(n)
                  nloc=gid-stmin+1
                  tmp(gid)%n_=tmp(gid)%n_+1
                  ! Boundary flags feed per/open only, both of which a composite build recomputes
                  if (.not.composite) then
                     if (per(1)) then
                        if (i.eq.dlo(1)) prov_flags(1,nloc)=1
                        if (i.eq.dhi(1)) prov_flags(2,nloc)=1
                     end if
                     if (per(2)) then
                        if (j.eq.dlo(2)) prov_flags(3,nloc)=1
                        if (j.eq.dhi(2)) prov_flags(4,nloc)=1
                     end if
                     if (per(3)) then
                        if (k.eq.dlo(3)) prov_flags(5,nloc)=1
                        if (k.eq.dhi(3)) prov_flags(6,nloc)=1
                     end if
                     do dim=1,3
                        do side=-1,1,2
                           pos=0; pos(dim)=side; nbr=[i,j,k]+pos
                           if (.not.per(dim)) then
                              if (nbr(dim).lt.dlo(dim).or.nbr(dim).gt.dhi(dim)) cycle
                           end if
                           if (pCov(nbr(1),nbr(2),nbr(3),1).eq.0) prov_flags(7,nloc)=1
                        end do
                     end do
                  end if
                  pLabel(i,j,k,1)=gid
               end if
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         deallocate(idmap)
         call move_alloc(tmp,this%struct)
      end block compact_tree

      ! Allocate the global lineage map and fill it with selves.  During stitching it represents
      ! only the equivalences observed by this rank.
      allocate_lineage: block
         integer :: n
         allocate(parent(max(this%nstruct,1))); parent=0
         do n=1,this%nstruct
            parent(n)=n
         end do
      end block allocate_lineage

      ! Stitch box-local structures together. One halo exchange carries every label across every
      ! box face at this level, periodic wraps included.  Scanning only each box's three low faces
      ! visits every adjacency exactly once: an ordinary interface is owned by the high-side box,
      ! and a periodic wrap by the box on the domain-low side.
      stitching: block
         use amrex_interface, only: amrimultifab_fill_boundary
         integer :: i,j,k,dim,nme,nnb
         integer, dimension(3) :: pos,flo,fhi
         call amrimultifab_fill_boundary(this%label,this%amr%geom(lvl))
         call this%amr%mfiter_build(lvl,mfi,tiling=.false.)
         do while (mfi%next())
            call prepare_fab(lvl,mfi)
            pLabel=>this%label%dataptr(mfi)
            bx=mfi%validbox()
            do dim=1,3
               pos=0; pos(dim)=-1
               flo=bx%lo; fhi=bx%hi; fhi(dim)=bx%lo(dim)
               do k=flo(3),fhi(3); do j=flo(2),fhi(2); do i=flo(1),fhi(1)
                  nme=pLabel(i,j,k,1)
                  if (nme.le.0) cycle
                  nnb=pLabel(i+pos(1),j+pos(2),k+pos(3),1)
                  if (nnb.le.0.or.nnb.eq.nme) cycle
                  if (same_label(i,j,k,i+pos(1),j+pos(2),k+pos(3))) call union_parent(nme,nnb)
               end do; end do; end do
            end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end block stitching

      ! Compress each rank's stitching result to a small set of equivalence edges, gather those
      ! edges once, resolve their transitive closure deterministically on rank zero, and broadcast
      ! the final map.  This replaces the old fixed-point loop's 2-4 full-parent allreduces plus a
      ! scalar allreduce per iteration.  At most nstruct-1 forest edges are sent by any rank even
      ! if thousands of face cells repeat the same component adjacency.
      resolve_lineage: block
         integer :: n,nedge,e,total_edges
         integer, dimension(:), allocatable :: edge_counts,edge_displs
         integer, dimension(:,:), allocatable :: local_edges,global_edges

         nedge=0
         do n=1,this%nstruct
            parent(n)=rootify_parent(n)
            if (parent(n).ne.n) nedge=nedge+1
         end do
         allocate(local_edges(2,max(nedge,1)))
         e=0
         do n=1,this%nstruct
            if (parent(n).ne.n) then
               e=e+1
               local_edges(:,e)=[n,parent(n)]
            end if
         end do

         allocate(edge_counts(0:this%amr%nproc-1)); edge_counts=0
         allocate(edge_displs(0:this%amr%nproc-1)); edge_displs=0
         call MPI_GATHER(nedge,1,MPI_INTEGER,edge_counts,1,MPI_INTEGER,0,this%amr%comm,ierr)

         total_edges=0
         if (this%amr%rank.eq.0) then
            total_edges=sum(edge_counts)
            edge_counts=2*edge_counts
            do n=1,this%amr%nproc-1
               edge_displs(n)=edge_displs(n-1)+edge_counts(n-1)
            end do
            allocate(global_edges(2,max(total_edges,1)))
         else
            ! mpi_f08 still requires a valid receive object on non-root ranks even though MPI
            ! ignores it there.
            allocate(global_edges(2,1))
         end if
         call MPI_GATHERV(local_edges,2*nedge,MPI_INTEGER,global_edges,edge_counts,edge_displs, &
            MPI_INTEGER,0,this%amr%comm,ierr)

         if (this%amr%rank.eq.0) then
            do n=1,this%nstruct
               parent(n)=n
            end do
            do e=1,total_edges
               call union_parent(global_edges(1,e),global_edges(2,e))
            end do
            do n=1,this%nstruct
               parent(n)=rootify_parent(n)
            end do
         end if
         call MPI_BCAST(parent,this%nstruct,MPI_INTEGER,0,this%amr%comm,ierr)
         deallocate(local_edges,global_edges,edge_counts,edge_displs)
      end block resolve_lineage

      ! Build the final dense id map directly from the resolved lineage.  Every provisional
      ! structure owns cells, so the surviving ids are exactly the roots in parent; there is no
      ! need to rediscover them by scanning the label field and globally reducing a used-id mask.
      ! Rank-local final cell counts likewise come from the already known provisional counts.
      ! If requested, fold the largest-to-id-1 permutation into this map before touching cells,
      ! avoiding a later full scan and halo exchange.
      compact_struct: block
         integer :: i,j,k,n,nfinal,bigid,tmpcount
         integer, dimension(:), allocatable :: root_dense,dense_map,local_counts,global_counts
         integer, dimension(:,:), allocatable :: dense_flags
         type(struct_type), dimension(:), allocatable :: tmp

         do n=1,this%nstruct
            parent(n)=rootify_parent(n)
         end do
         allocate(root_dense(max(this%nstruct,1))); root_dense=0
         nfinal=0
         do n=1,this%nstruct
            if (parent(n).eq.n) then
               nfinal=nfinal+1
               root_dense(n)=nfinal
            end if
         end do
         allocate(dense_map(max(this%nstruct,1))); dense_map=0
         do n=1,this%nstruct
            dense_map(n)=root_dense(parent(n))
         end do

         allocate(local_counts(max(nfinal,1))); local_counts=0
         allocate(dense_flags(7,max(nfinal,1))); dense_flags=0
         do n=stmin,stmax
            local_counts(dense_map(n))=local_counts(dense_map(n))+this%struct(n)%n_
            dense_flags(:,dense_map(n))=max(dense_flags(:,dense_map(n)),prov_flags(:,n-stmin+1))
         end do

         bigid=1
         if (this%rename_largest.and.nfinal.gt.1) then
            allocate(global_counts(nfinal)); global_counts=0
            call MPI_REDUCE(local_counts,global_counts,nfinal,MPI_INTEGER,MPI_SUM,0,this%amr%comm,ierr)
            if (this%amr%rank.eq.0) bigid=maxloc(global_counts,1)
            call MPI_BCAST(bigid,1,MPI_INTEGER,0,this%amr%comm,ierr)
            deallocate(global_counts)
            if (bigid.ne.1) then
               do n=1,this%nstruct
                  if (dense_map(n).eq.1) then
                     dense_map(n)=bigid
                  else if (dense_map(n).eq.bigid) then
                     dense_map(n)=1
                  end if
               end do
               tmpcount=local_counts(1); local_counts(1)=local_counts(bigid); local_counts(bigid)=tmpcount
               dense_flags(:,[1,bigid])=dense_flags(:,[bigid,1])
            end if
         end if

         ! All zero in a composite build, whose caller recomputes per/open anyway
         if (nfinal.gt.0.and..not.composite) call MPI_ALLREDUCE(MPI_IN_PLACE,dense_flags,7*nfinal,MPI_INTEGER,MPI_MAX,this%amr%comm,ierr)

         allocate(tmp(max(nfinal,1)))
         do n=1,size(tmp)
            tmp(n)%parent=n; tmp(n)%n_=local_counts(n); tmp(n)%per=0; tmp(n)%open=.false.
         end do
         do n=1,nfinal
            if (dense_flags(1,n).eq.1.and.dense_flags(2,n).eq.1) tmp(n)%per(1)=1
            if (dense_flags(3,n).eq.1.and.dense_flags(4,n).eq.1) tmp(n)%per(2)=1
            if (dense_flags(5,n).eq.1.and.dense_flags(6,n).eq.1) tmp(n)%per(3)=1
            tmp(n)%open=(dense_flags(7,n).eq.1)
         end do

         ! This is the only post-lineage cell pass: provisional id -> final dense id.
         call this%amr%mfiter_build(lvl,mfi,tiling=.false.)
         do while (mfi%next())
            pLabel=>this%label%dataptr(mfi)
            pId=>this%id%mf(lvl)%dataptr(mfi)
            bx=mfi%validbox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               n=pLabel(i,j,k,1)
               if (n.gt.0) pId(i,j,k,1)=real(dense_map(n),WP)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         ! build_composite offsets these IDs and then FillPatches every ghost, so skip it there
         if (.not.composite) call this%id%mf(lvl)%fill_boundary(this%amr%geom(lvl))

         this%nstruct=nfinal
         call move_alloc(tmp,this%struct)
         deallocate(root_dense,dense_map,local_counts,dense_flags,prov_flags)
      end block compact_struct
      deallocate(parent)

   contains

      !> Point the local lineage of a structure to its root and return that root.
      function rootify_struct(x) result(y)
         implicit none
         integer, intent(in) :: x
         integer :: y,node,next
         y=x
         do while (y.ne.this%struct(y)%parent)
            y=this%struct(y)%parent
         end do
         node=x
         do while (node.ne.y)
            next=this%struct(node)%parent
            this%struct(node)%parent=y
            node=next
         end do
      end function rootify_struct

      !> This function joins two structures at their roots (the smallest root is chosen and returned)
      function union_struct(x,y) result(rmin)
         implicit none
         integer, intent(in) :: x,y
         integer :: rx,ry,rmin,rmax
         rx=rootify_struct(x); ry=rootify_struct(y)
         rmin=min(rx,ry); rmax=max(rx,ry)
         this%struct(rmax)%parent=rmin
      end function union_struct

      !> This function adds one new root while dynamically handling storage space
      function add() result(x)
         implicit none
         integer :: x
         integer :: size_now,size_new
         type(struct_type), dimension(:), allocatable :: tmp
         ! Check if there is enough room for storing a new structure
         size_now=size(this%struct,dim=1)
         if (nstruct_.eq.size_now) then
            size_new=int(real(size_now,WP)*coeff_up)
            allocate(tmp(size_new))
            tmp(1:nstruct_)=this%struct
            tmp(nstruct_+1:)%parent=0
            tmp(nstruct_+1:)%n_=0
            tmp(nstruct_+1:)%per(1)=0
            tmp(nstruct_+1:)%per(2)=0
            tmp(nstruct_+1:)%per(3)=0
            tmp(nstruct_+1:)%open=.false.
            call move_alloc(tmp,this%struct)
         end if
         ! Add new root
         nstruct_=nstruct_+1
         this%struct(nstruct_)%parent=nstruct_
         this%struct(nstruct_)%n_=0
         this%struct(nstruct_)%per=0
         this%struct(nstruct_)%open=.false.
         x=nstruct_
      end function add

      !> Point global lineage to its root and return that root.
      function rootify_parent(x) result(y)
         implicit none
         integer, intent(in) :: x
         integer :: y,node,next
         y=x
         do while (y.ne.parent(y))
            y=parent(y)
         end do
         node=x
         do while (node.ne.y)
            next=parent(node)
            parent(node)=y
            node=next
         end do
      end function rootify_parent

      !> This function joins two branches at their roots (the smallest root is chosen)
      subroutine union_parent(x,y)
         implicit none
         integer, intent(in) :: x,y
         integer :: rx,ry,rmin,rmax
         rx=rootify_parent(x); ry=rootify_parent(y); rmin=min(rx,ry); rmax=max(rx,ry)
         parent(rmax)=rmin
      end subroutine union_parent

   end subroutine build


   !> Label the non-overlapping (leaf) cells of the hierarchy. Covered coarse cells never
   !> participate, including in connectivity: a coarse average must not weld detached fragments.
   !> First build same-level components with the caller's connectivity test; then FillPatch their
   !> disjoint IDs with piecewise-constant interpolation and union fine/coarse face neighbors.
   !> The callbacks must support every level. same_label is also called at a fine/coarse face
   !> using fine indices (the caller's ghost data describes the coarse neighbor there).
   !> On return id is nonzero ONLY on liquid leaf cells at every level, open is false, and
   !> has_coarse identifies structures which cannot be removed by a finest-level-only transfer.
   !> lvl=-1 denotes a composite build. rename_largest is not applied: select by physical volume.
   subroutine build_composite(this,prepare_fab,make_label,same_label)
      use amrex_amr_module, only: amrex_imultifab_build,amrex_imultifab_destroy
      use amrex_interface, only: amrmask_make_fine
      use mpi_f08, only: MPI_GATHER,MPI_GATHERV,MPI_BCAST,MPI_ALLREDUCE,MPI_IN_PLACE,MPI_INTEGER,MPI_MAX
      implicit none
      class(amrccl), intent(inout) :: this
      procedure(prepare_fab_ftype) :: prepare_fab
      procedure(make_label_ftype) :: make_label
      procedure(same_label_ftype) :: same_label
      type(amrex_imultifab) :: leaf
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pId
      integer, dimension(:,:,:,:), contiguous, pointer :: pLeaf
      integer, allocatable :: parent(:),dense(:),counts(:),displs(:),edges(:,:),all_edges(:,:),flags(:,:)
      integer :: lvl,top,offset,ntotal,i,j,k,n,nb,d,side,nedge,e,total_edges,nfinal,ierr
      integer :: nbr(3),dlo(3),dhi(3)
      logical :: periodic(3)

      top=this%amr%clvl()
      if (.not.this%id_layout_ready) then
         call this%id%reset()
         this%id_layout_ready=.true.
      end if
      call this%id%setval(0.0_WP)
      offset=0
      do lvl=0,top
         call amrex_imultifab_build(leaf,this%amr%ba(lvl),this%amr%dm(lvl),1,0)
         call leaf%setval(1)
         if (lvl.lt.top) call amrmask_make_fine(leaf,this%amr%ba(lvl+1), &
            [this%amr%rrefx(lvl),this%amr%rrefy(lvl),this%amr%rrefz(lvl)],0,1)
         call this%build(lvl,prepare_leaf,make_leaf,same_label,preserve_levels=.true.)
         call this%amr%mfiter_build(lvl,mfi,tiling=.false.)
         do while(mfi%next())
            pId=>this%id%mf(lvl)%dataptr(mfi)
            bx=mfi%validbox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               if (pId(i,j,k,1).gt.0.0_WP) pId(i,j,k,1)=pId(i,j,k,1)+real(offset,WP)
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
         offset=offset+this%nstruct
         call amrex_imultifab_destroy(leaf)
      end do
      ntotal=offset
      allocate(parent(ntotal),dense(ntotal))
      do n=1,ntotal
         parent(n)=n
      end do
      ! Coarse leaf IDs fill only the fine ghosts outside that fine level. Covered coarse IDs
      ! are zero. Proper nesting ensures each fine face is adjacent to level L or level L-1.
      call this%id%fill(0.0_WP)
      do lvl=1,top
         call this%amr%mfiter_build(lvl,mfi,tiling=.false.)
         do while(mfi%next())
            call prepare_fab(lvl,mfi)
            pId=>this%id%mf(lvl)%dataptr(mfi)
            bx=mfi%validbox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               n=nint(pId(i,j,k,1)); if (n.eq.0) cycle
               do d=1,3; do side=-1,1,2
                  nbr=[i,j,k]; nbr(d)=nbr(d)+side
                  nb=nint(pId(nbr(1),nbr(2),nbr(3),1))
                  if (nb.eq.0.or.nb.eq.n) cycle
                  if (same_label(i,j,k,nbr(1),nbr(2),nbr(3))) call join(n,nb)
               end do; end do
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      ! Communicate a forest, rather than one edge for every coarse/fine face cell.
      nedge=0
      do n=1,ntotal
         parent(n)=root(n)
         if (parent(n).ne.n) nedge=nedge+1
      end do
      allocate(edges(2,max(1,nedge))); e=0
      do n=1,ntotal
         if (parent(n).eq.n) cycle
         e=e+1; edges(:,e)=[n,parent(n)]
      end do
      allocate(counts(0:this%amr%nproc-1),displs(0:this%amr%nproc-1))
      counts=0; displs=0
      call MPI_GATHER(2*nedge,1,MPI_INTEGER,counts,1,MPI_INTEGER,0,this%amr%comm,ierr)
      total_edges=0
      if (this%amr%rank.eq.0) then
         total_edges=sum(counts)/2
         do n=1,this%amr%nproc-1
            displs(n)=displs(n-1)+counts(n-1)
         end do
      end if
      allocate(all_edges(2,max(1,total_edges)))
      call MPI_GATHERV(edges,2*nedge,MPI_INTEGER,all_edges,counts,displs,MPI_INTEGER,0,this%amr%comm,ierr)
      if (this%amr%rank.eq.0) then
         do n=1,ntotal
            parent(n)=n
         end do
         do e=1,total_edges
            call join(all_edges(1,e),all_edges(2,e))
         end do
         do n=1,ntotal
            parent(n)=root(n)
         end do
      end if
      call MPI_BCAST(parent,ntotal,MPI_INTEGER,0,this%amr%comm,ierr)
      dense=0; nfinal=0
      do n=1,ntotal
         if (parent(n).ne.n) cycle
         nfinal=nfinal+1; dense(n)=nfinal
      end do
      do n=1,ntotal
         dense(n)=dense(parent(n))
      end do
      call this%empty()
      this%nstruct=nfinal
      allocate(this%struct(nfinal),flags(7,nfinal)); flags=0
      this%struct(:)%n_=0; this%struct(:)%open=.false.
      periodic=[this%amr%xper,this%amr%yper,this%amr%zper]
      do lvl=0,top
         dlo=this%amr%geom(lvl)%domain%lo; dhi=this%amr%geom(lvl)%domain%hi
         call this%amr%mfiter_build(lvl,mfi,tiling=.false.)
         do while(mfi%next())
            pId=>this%id%mf(lvl)%dataptr(mfi)
            bx=mfi%validbox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               n=nint(pId(i,j,k,1)); if (n.eq.0) cycle
               n=dense(n); pId(i,j,k,1)=real(n,WP)
               this%struct(n)%n_=this%struct(n)%n_+1
               if (lvl.lt.this%amr%maxlvl) flags(7,n)=1
               nbr=[i,j,k]
               do d=1,3
                  if (.not.periodic(d)) cycle
                  if (nbr(d).eq.dlo(d)) flags(2*d-1,n)=1
                  if (nbr(d).eq.dhi(d)) flags(2*d,n)=1
               end do
            end do; end do; end do
         end do
         call this%amr%mfiter_destroy(mfi)
      end do
      ! Refresh coarse/fine ghosts as well: provisional IDs must not escape the build.
      call this%id%fill(0.0_WP)
      call MPI_ALLREDUCE(MPI_IN_PLACE,flags,size(flags),MPI_INTEGER,MPI_MAX,this%amr%comm,ierr)
      do n=1,nfinal
         this%struct(n)%parent=n
         this%struct(n)%has_coarse=(flags(7,n).eq.1)
         do d=1,3
            this%struct(n)%per(d)=min(flags(2*d-1,n),flags(2*d,n))
         end do
      end do
   contains
      subroutine prepare_leaf(l,m)
         integer, intent(in) :: l
         type(amrex_mfiter), intent(inout) :: m
         call prepare_fab(l,m)
         pLeaf=>leaf%dataptr(m)
      end subroutine prepare_leaf
      logical function make_leaf(ii,jj,kk)
         integer, intent(in) :: ii,jj,kk
         make_leaf=.false.
         if (pLeaf(ii,jj,kk,1).eq.0) return
         make_leaf=make_label(ii,jj,kk)
      end function make_leaf
      integer function root(x) result(r)
         integer, intent(in) :: x
         integer :: node,next
         r=x
         do while(parent(r).ne.r)
            r=parent(r)
         end do
         node=x
         do while(node.ne.r)
            next=parent(node); parent(node)=r; node=next
         end do
      end function root
      subroutine join(x,y)
         integer, intent(in) :: x,y
         integer :: rx,ry
         rx=root(x); ry=root(y); parent(max(rx,ry))=min(rx,ry)
      end subroutine join
   end subroutine build_composite


   !> Sum-reduce a per-structure array across all ranks - the usual way to turn rank-local
   !> accumulations into global per-structure quantities
   subroutine allreduce_1d(this,vals)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(amrccl), intent(in) :: this
      real(WP), dimension(:), intent(inout) :: vals
      integer :: ierr
      if (size(vals).eq.0) return
      call MPI_ALLREDUCE(MPI_IN_PLACE,vals,size(vals),MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
   end subroutine allreduce_1d

   subroutine allreduce_2d(this,vals)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(amrccl), intent(in) :: this
      real(WP), dimension(:,:), intent(inout) :: vals
      integer :: ierr
      if (size(vals).eq.0) return
      call MPI_ALLREDUCE(MPI_IN_PLACE,vals,size(vals),MPI_REAL_WP,MPI_SUM,this%amr%comm,ierr)
   end subroutine allreduce_2d


   !> Empty structure info
   subroutine empty(this)
      implicit none
      class(amrccl), intent(inout) :: this
      if (allocated(this%struct)) deallocate(this%struct)
      this%nstruct=0
      this%lvl=-1
   end subroutine empty


   !> Finalize CCL object
   subroutine finalize(this)
      use amrex_amr_module, only: amrex_imultifab_destroy
      implicit none
      class(amrccl), intent(inout) :: this
      call this%empty()
      if (this%work_allocated) then
         call amrex_imultifab_destroy(this%label)
         call amrex_imultifab_destroy(this%cover)
      end if
      this%work_allocated=.false.
      this%work_valid=.false.
      this%cover_valid=.false.
      this%work_lvl=-1
      this%id_layout_ready=.false.
      call this%id%finalize()
      nullify(this%amr)
      this%name='UNNAMED_AMRCCL'
   end subroutine finalize


end module amrccl_class
