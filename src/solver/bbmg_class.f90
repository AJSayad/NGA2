!> Black-Box Multi-Grid (BBMG) solver object is defined in this class.
!> It is intended to be driven by an ILS object. It is based on the
!> algorithm by Dendy (1983, 1986), and includes an additional periodic
!> treatment (Dendy, 1988) and Shapira's improvement (2008)
module bbmg_class
   use precision, only: WP
   use string,    only: str_medium
   use mpi_f08
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: bbmg
   
   ! Level data object definition
   type :: lvl_data
      
      ! Global sizes
      integer :: ncell,ncello
      integer :: nx ,imin ,imax ,nxo ,imino ,imaxo
      integer :: ny ,jmin ,jmax ,nyo ,jmino ,jmaxo
      integer :: nz ,kmin ,kmax ,nzo ,kmino ,kmaxo
      
      ! Local sizes
      integer :: ncell_,ncello_
      integer :: nx_,imin_,imax_,nxo_,imino_,imaxo_
      integer :: ny_,jmin_,jmax_,nyo_,jmino_,jmaxo_
      integer :: nz_,kmin_,kmax_,nzo_,kmino_,kmaxo_
      
      ! Operators
      real(WP), dimension(:,:,:,:,:,:), allocatable :: c2f
      real(WP), dimension(:,:,:,:,:,:), allocatable :: f2c
      real(WP), dimension(:,:,:,:,:,:), allocatable :: opr
      real(WP), dimension(:,:,:,:,:,:), allocatable :: oprc2f
      
      ! Vectors
      real(WP), dimension(:,:,:), allocatable :: v
      real(WP), dimension(:,:,:), allocatable :: f
      real(WP), dimension(:,:,:), allocatable :: r
      
      ! Parallel information
      logical, dimension(:), allocatable :: send_xm,send_xp
      logical, dimension(:), allocatable :: send_ym,send_yp
      logical, dimension(:), allocatable :: send_zm,send_zp
      integer :: nsend_xm,nsend_xp
      integer :: nsend_ym,nsend_yp
      integer :: nsend_zm,nsend_zp
      integer :: recv_xm,recv_xp
      integer :: recv_ym,recv_yp
      integer :: recv_zm,recv_zp
      
   end type lvl_data
   
   
   !> BBMG object definition
   type :: bbmg
      
      ! A BBMG has a name
      character(len=str_medium) :: name                               !< Name of solver
      
      ! Parallel information
      type(MPI_Comm) :: comm                                          !< Grid communicator
      logical :: amRoot                                               !< Root process for messaging
      integer :: npx,npy,npz                                          !< Number of processors in each direction
      integer :: nproc                                                !< Number of processors
      integer :: irank                                                !< Rank info
      integer :: iproc,jproc,kproc                                    !< Processor coordinates
      integer, dimension(:,:,:), allocatable :: rank                  !< Cartesian rank information
      
      ! General grid information
      integer :: no                                                   !< Size of the grid overlap
      logical :: xper,yper,zper                                       !< Periodicity info for the problem
      
      ! General operator information
      integer :: nstx                                                 !< Number of diagonal in x
      integer :: nsty                                                 !< Number of diagonal in y
      integer :: nstz                                                 !< Number of diagonal in z
      
      ! Solver parameters
      logical  :: use_pcg                                             !< Use BBMG as a preconditioner to a PCG or directly
      real(WP) :: max_res                                             !< Maximum residual - the solver will attempt to converge below that value
      integer  :: max_ite                                             !< Maximum number of iterations - the solver will not got past that value
      integer  :: cycle=1                                             !< Cycle type: 1=V-cycle (default), 2=W-cycle, ...
      integer  :: relax_pre =1                                        !< Number of pre-sweeps (default is 1)
      integer  :: relax_post=1                                        !< Number of post-sweeps (default is 1)
      integer  :: ncell_direct=0                                      !< Coarse problem size below which a direct solver is employed (default is 0, i.e., no direct solve)
      
      ! Level data management
      integer :: nlvl                                                 !< Number of multigrid levels created - by default all will be used. This can be reduced by the user.
      type(lvl_data), dimension(:), allocatable :: lvl                !< Entire data at each level
      
   contains
      !procedure :: print_cvg                                          !< Print BBMG convergence history
      !procedure :: log_cvg                                            !< Log BBMG convergence history
      procedure :: init_solver                                        !< Solver initialization (at start-up)
      procedure :: update_operator                                    !< Operator update (every time the operator changes)
      procedure :: solve                                              !< Solve the linear system
      
      procedure :: recompute_prolongation                             !< Recompute the prolongation at a given level
      procedure :: recompute_restriction                              !< Recompute the restriction at a given level
      procedure :: recompute_operator                                 !< Recompute the operator at a given level
      
      procedure :: vsync                                              !< Synchronize boundaries for a vector at a given level
      procedure :: msync                                              !< Synchronize boundaries for a matrix at a given level
      
      procedure :: pmodx,pmody,pmodz                                  !< Parity calculation that accounts for periodicity
      
   end type bbmg
   
   
   
   
   
   !> Declare bbmg constructor
   interface bbmg
      procedure bbmg_from_pgrid
   end interface bbmg
   
   
contains
   
   
   !> Division of an integer by 2
   pure integer function div(ind)
     implicit none
     integer, intent(in) :: ind
     div=ind/2+mod(ind,2)
   end function div
   
   
   !> Parity of point i shifted by n to ensure 1<=i<=n
   pure integer function pmodx(this,i,n)
      class(bbmg), intent(in) :: this
      integer, intent(in) :: i,n
      if (this%xper) then
         pmodx=1-mod(mod(i+this%lvl(n)%nx-1,this%lvl(n)%nx)+1,2)
      else
         pmodx=1-mod(i,2)
      end if
   end function pmodx
   pure integer function pmody(this,i,n)
      class(bbmg), intent(in) :: this
      integer, intent(in) :: i,n
      if (this%yper) then
         pmody=1-mod(mod(i+this%lvl(n)%ny-1,this%lvl(n)%ny)+1,2)
      else
         pmody=1-mod(i,2)
      end if
   end function pmody
   pure integer function pmodz(this,i,n)
      class(bbmg), intent(in) :: this
      integer, intent(in) :: i,n
      if (this%zper) then
         pmodz=1-mod(mod(i+this%lvl(n)%nz-1,this%lvl(n)%nz)+1,2)
      else
         pmodz=1-mod(i,2)
      end if
   end function pmodz
   
   
   !> Constructor for a BBMG object - this sets the grid and storage
   function bbmg_from_pgrid(pg,name,nst) result(self)
      use messager,    only: die
      use pgrid_class, only: pgrid
      implicit none
      type(bbmg) :: self
      class(pgrid), intent(in) :: pg
      character(len=*), intent(in) :: name
      integer, dimension(3), intent(in) :: nst
      
      ! Process the operator size
      initialize_overlap: block
         integer :: maxst
         if (nst(1).ne.3.or.nst(2).ne.3.or.nst(3).ne.3) call die('[bbmg constructor] BBMG only supports 3x3x3 operators at this point')
         self%nstx=nst(1); if (mod(self%nstx,2).eq.0) call die('[bbmg constructor] Operator size in x needs to be odd')
         self%nsty=nst(2); if (mod(self%nsty,2).eq.0) call die('[bbmg constructor] Operator size in y needs to be odd')
         self%nstz=nst(3); if (mod(self%nstz,2).eq.0) call die('[bbmg constructor] Operator size in z needs to be odd')
         maxst=max(self%nstx,self%nsty,self%nstz)
         self%no=(maxst-1)/2
      end block initialize_overlap
      
      ! First build the hierarchy of grids
      initialize_grid: block
         integer :: i,n1,n2,n3,nlvl_x,nlvl_y,nlvl_z
         ! Store the name
         self%name  =trim(adjustl(name))
         ! Copy communicator from the pgrid
         self%comm  =pg%comm
         ! Copy root from the pgrid
         self%amRoot=pg%amRoot
         ! Copy number of processors per direction from the pgrid
         self%npx   =pg%npx
         self%npy   =pg%npy
         self%npz   =pg%npz
         ! Copy total number of processors from the pgrid
         self%nproc =pg%nproc
         ! Copy processor rank from the pgrid (we assume irank is starting at 1)
         self%irank =pg%rank+1 !< pgrid starts counting rank at 0
         ! Copy processor coordinates from the pgrid
         self%iproc =pg%iproc
         self%jproc =pg%jproc
         self%kproc =pg%kproc
         ! Copy periodicity from the pgrid
         self%xper=pg%xper
         self%yper=pg%yper
         self%zper=pg%zper
         ! Compute maximum number of levels per direction
         nlvl_x=1; do while (pg%nx.gt.2**(nlvl_x-1)); nlvl_x=nlvl_x+1; end do
         nlvl_y=1; do while (pg%ny.gt.2**(nlvl_y-1)); nlvl_y=nlvl_y+1; end do
         nlvl_z=1; do while (pg%nz.gt.2**(nlvl_z-1)); nlvl_z=nlvl_z+1; end do
         ! Compute maximum number of levels
         n1=nlvl_x; if (nlvl_x.eq.1) n1=huge(nlvl_x)
         n2=nlvl_y; if (nlvl_y.eq.1) n2=huge(nlvl_y)
         n3=nlvl_z; if (nlvl_z.eq.1) n3=huge(nlvl_z)
         self%nlvl=min(n1,n2,n3); if (nlvl_x.eq.1.and.nlvl_y.eq.1.and.nlvl_z.eq.1) self%nlvl=1
         ! Allocate lvl structure
         allocate(self%lvl(self%nlvl))
         ! Generate hierarchy of grids in x from the pgrid
         self%lvl(1)%nx   =pg%nx
         self%lvl(1)%imin =pg%imin -pg%no
         self%lvl(1)%imax =pg%imax -pg%no
         self%lvl(1)%imin_=pg%imin_-pg%no
         self%lvl(1)%imax_=pg%imax_-pg%no
         self%lvl(1)%imino =self%lvl(1)%imin -self%no
         self%lvl(1)%imaxo =self%lvl(1)%imax +self%no
         self%lvl(1)%imino_=self%lvl(1)%imin_-self%no
         self%lvl(1)%imaxo_=self%lvl(1)%imax_+self%no
         self%lvl(1)%nx_ =self%lvl(1)%imax_-self%lvl(1)%imin_+1
         self%lvl(1)%nxo =self%lvl(1)%nx +2*self%no
         self%lvl(1)%nxo_=self%lvl(1)%nx_+2*self%no
         do i=2,self%nlvl
            ! Divide indices and sizes by two
            self%lvl(i)%nx   =div(self%lvl(i-1)%nx     )
            self%lvl(i)%imin =div(self%lvl(i-1)%imin +1)
            self%lvl(i)%imax =div(self%lvl(i-1)%imax   )
            self%lvl(i)%imin_=div(self%lvl(i-1)%imin_+1)
            self%lvl(i)%imax_=div(self%lvl(i-1)%imax_  )
            ! Set secondary info
            self%lvl(i)%imino =self%lvl(i)%imin -self%no
            self%lvl(i)%imaxo =self%lvl(i)%imax +self%no
            self%lvl(i)%imino_=self%lvl(i)%imin_-self%no
            self%lvl(i)%imaxo_=self%lvl(i)%imax_+self%no
            self%lvl(i)%nx_ =self%lvl(i)%imax_-self%lvl(i)%imin_+1
            self%lvl(i)%nxo =self%lvl(i)%nx +2*self%no
            self%lvl(i)%nxo_=self%lvl(i)%nx_+2*self%no
         end do
         ! Generate hierarchy of grids in y from the pgrid
         self%lvl(1)%ny   =pg%ny
         self%lvl(1)%jmin =pg%jmin -pg%no
         self%lvl(1)%jmax =pg%jmax -pg%no
         self%lvl(1)%jmin_=pg%jmin_-pg%no
         self%lvl(1)%jmax_=pg%jmax_-pg%no
         self%lvl(1)%jmino =self%lvl(1)%jmin -self%no
         self%lvl(1)%jmaxo =self%lvl(1)%jmax +self%no
         self%lvl(1)%jmino_=self%lvl(1)%jmin_-self%no
         self%lvl(1)%jmaxo_=self%lvl(1)%jmax_+self%no
         self%lvl(1)%ny_ =self%lvl(1)%jmax_-self%lvl(1)%jmin_+1
         self%lvl(1)%nyo =self%lvl(1)%ny +2*self%no
         self%lvl(1)%nyo_=self%lvl(1)%ny_+2*self%no
         do i=2,self%nlvl
            ! Divide indices and sizes by two
            self%lvl(i)%ny   =div(self%lvl(i-1)%ny     )
            self%lvl(i)%jmin =div(self%lvl(i-1)%jmin +1)
            self%lvl(i)%jmax =div(self%lvl(i-1)%jmax   )
            self%lvl(i)%jmin_=div(self%lvl(i-1)%jmin_+1)
            self%lvl(i)%jmax_=div(self%lvl(i-1)%jmax_  )
            ! Set secondary info
            self%lvl(i)%jmino =self%lvl(i)%jmin -self%no
            self%lvl(i)%jmaxo =self%lvl(i)%jmax +self%no
            self%lvl(i)%jmino_=self%lvl(i)%jmin_-self%no
            self%lvl(i)%jmaxo_=self%lvl(i)%jmax_+self%no
            self%lvl(i)%ny_ =self%lvl(i)%jmax_-self%lvl(i)%jmin_+1
            self%lvl(i)%nyo =self%lvl(i)%ny +2*self%no
            self%lvl(i)%nyo_=self%lvl(i)%ny_+2*self%no
         end do
         ! Generate hierarchy of grids in z from the pgrid
         self%lvl(1)%nz   =pg%nz
         self%lvl(1)%kmin =pg%kmin -pg%no
         self%lvl(1)%kmax =pg%kmax -pg%no
         self%lvl(1)%kmin_=pg%kmin_-pg%no
         self%lvl(1)%kmax_=pg%kmax_-pg%no
         self%lvl(1)%kmino =self%lvl(1)%kmin -self%no
         self%lvl(1)%kmaxo =self%lvl(1)%kmax +self%no
         self%lvl(1)%kmino_=self%lvl(1)%kmin_-self%no
         self%lvl(1)%kmaxo_=self%lvl(1)%kmax_+self%no
         self%lvl(1)%nz_ =self%lvl(1)%kmax_-self%lvl(1)%kmin_+1
         self%lvl(1)%nzo =self%lvl(1)%nz +2*self%no
         self%lvl(1)%nzo_=self%lvl(1)%nz_+2*self%no
         do i=2,self%nlvl
            ! Divide indices and sizes by two
            self%lvl(i)%nz   =div(self%lvl(i-1)%nz     )
            self%lvl(i)%kmin =div(self%lvl(i-1)%kmin +1)
            self%lvl(i)%kmax =div(self%lvl(i-1)%kmax   )
            self%lvl(i)%kmin_=div(self%lvl(i-1)%kmin_+1)
            self%lvl(i)%kmax_=div(self%lvl(i-1)%kmax_  )
            ! Set secondary info
            self%lvl(i)%kmino =self%lvl(i)%kmin -self%no
            self%lvl(i)%kmaxo =self%lvl(i)%kmax +self%no
            self%lvl(i)%kmino_=self%lvl(i)%kmin_-self%no
            self%lvl(i)%kmaxo_=self%lvl(i)%kmax_+self%no
            self%lvl(i)%nz_ =self%lvl(i)%kmax_-self%lvl(i)%kmin_+1
            self%lvl(i)%nzo =self%lvl(i)%nz +2*self%no
            self%lvl(i)%nzo_=self%lvl(i)%nz_+2*self%no
         end do
         ! Count total number of unknowns per level
         do i=1,self%nlvl
            self%lvl(i)%ncell  =max(self%lvl(i)%nx  ,0)*max(self%lvl(i)%ny  ,0)*max(self%lvl(i)%nz  ,0)
            self%lvl(i)%ncello =max(self%lvl(i)%nxo ,0)*max(self%lvl(i)%nyo ,0)*max(self%lvl(i)%nzo ,0)
            self%lvl(i)%ncell_ =max(self%lvl(i)%nx_ ,0)*max(self%lvl(i)%ny_ ,0)*max(self%lvl(i)%nz_ ,0)
            self%lvl(i)%ncello_=max(self%lvl(i)%nxo_,0)*max(self%lvl(i)%nyo_,0)*max(self%lvl(i)%nzo_,0)
         end do
         ! Find direct solve level
         lvl_loop: do i=1,self%nlvl
            if (self%lvl(i)%ncell.le.self%ncell_direct) then
               self%nlvl=i
               exit lvl_loop
            end if
         end do lvl_loop
      end block initialize_grid
      
      ! Now prepare operator storage
      initialize_operator: block
         integer :: n
         do n=1,self%nlvl
            ! Handle empty processors
            if (self%lvl(n)%ncello_.eq.0) cycle
            ! Prolongation operator (from lvl n+1 to lvl n)
            if (n.lt.self%nlvl) then
               allocate(self%lvl(n)%c2f   ( 0:1, 0:1, 0:1,self%lvl(n)%imino_:self%lvl(n)%imaxo_,self%lvl(n)%jmino_:self%lvl(n)%jmaxo_,self%lvl(n)%kmino_:self%lvl(n)%kmaxo_)); self%lvl(n)%c2f=0.0_WP
            end if
            ! Restriction operator (from lvl n-1 to lvl n)
            if (n.gt.1) then
               allocate(self%lvl(n)%f2c   (-1:1,-1:1,-1:1,self%lvl(n)%imino_:self%lvl(n)%imaxo_,self%lvl(n)%jmino_:self%lvl(n)%jmaxo_,self%lvl(n)%kmino_:self%lvl(n)%kmaxo_)); self%lvl(n)%f2c=0.0_WP
            end if
            ! Laplacian operator
            allocate(self%lvl(n)%opr      (-1:1,-1:1,-1:1,self%lvl(n)%imino_:self%lvl(n)%imaxo_,self%lvl(n)%jmino_:self%lvl(n)%jmaxo_,self%lvl(n)%kmino_:self%lvl(n)%kmaxo_)); self%lvl(n)%opr=0.0_WP
            ! Laplacian*Prolongation operator (from lvl n+1 to lvl n)
            if (n.lt.self%nlvl) then
               allocate(self%lvl(n)%oprc2f(-1:1,-1:1,-1:1,self%lvl(n)%imino_:self%lvl(n)%imaxo_,self%lvl(n)%jmino_:self%lvl(n)%jmaxo_,self%lvl(n)%kmino_:self%lvl(n)%kmaxo_)); self%lvl(n)%oprc2f=0.0_WP
            end if
            ! Allocate vectors
            allocate(self%lvl(n)%v(self%lvl(n)%imino_:self%lvl(n)%imaxo_,self%lvl(n)%jmino_:self%lvl(n)%jmaxo_,self%lvl(n)%kmino_:self%lvl(n)%kmaxo_)); self%lvl(n)%v=0.0_WP
            allocate(self%lvl(n)%f(self%lvl(n)%imino_:self%lvl(n)%imaxo_,self%lvl(n)%jmino_:self%lvl(n)%jmaxo_,self%lvl(n)%kmino_:self%lvl(n)%kmaxo_)); self%lvl(n)%f=0.0_WP
            allocate(self%lvl(n)%r(self%lvl(n)%imino_:self%lvl(n)%imaxo_,self%lvl(n)%jmino_:self%lvl(n)%jmaxo_,self%lvl(n)%kmino_:self%lvl(n)%kmaxo_)); self%lvl(n)%r=0.0_WP
         end do
      end block initialize_operator
      
      ! Now prepare communication data
      initialize_comm: block
         integer :: ip,jp,kp,n,ierr
         integer :: my_min,my_max,hisrank
         integer, dimension(3) :: ploc
         integer, dimension(:), allocatable :: ncello_rank
         integer, dimension(:), allocatable :: mino_,min_,max_,maxo_
         integer, dimension(:), allocatable :: buf
         ! Prepare rank array
         allocate(self%rank(self%npx,self%npy,self%npz))
         do kp=1,self%npz
            do jp=1,self%npy
               do ip=1,self%npx
                  ploc=[ip-1,jp-1,kp-1]; call MPI_CART_RANK(self%comm,ploc,self%rank(ip,jp,kp),ierr)
                  self%rank(ip,jp,kp)=self%rank(ip,jp,kp)+1
               end do
            end do
         end do
         ! Allocate index and ncell work arrays
         allocate(buf(self%nproc),mino_(self%nproc),min_(self%nproc),max_(self%nproc),maxo_(self%nproc),ncello_rank(self%nproc))
         ! Work on each level
         do n=1,self%nlvl
            ! Gather mesh sizes at this level
            ncello_rank=0; call MPI_allgather(self%lvl(n)%ncello_,1,MPI_INTEGER,ncello_rank,1,MPI_INTEGER,self%comm,ierr)
            ! X DIRECTION ----------------------------------------------------------------------------------
            ! Allocate send and receive arrays
            allocate(self%lvl(n)%send_xm(self%npx)); self%lvl(n)%send_xm=.false.; self%lvl(n)%recv_xm=MPI_PROC_NULL
            allocate(self%lvl(n)%send_xp(self%npx)); self%lvl(n)%send_xp=.false.; self%lvl(n)%recv_xp=MPI_PROC_NULL
            ! Gather indices
            buf=0; buf(self%irank)=self%lvl(n)%imino_; call MPI_allreduce(buf,mino_,self%nproc,MPI_INTEGER,MPI_SUM,self%comm,ierr)
            buf=0; buf(self%irank)=self%lvl(n)%imin_ ; call MPI_allreduce(buf,min_ ,self%nproc,MPI_INTEGER,MPI_SUM,self%comm,ierr)
            buf=0; buf(self%irank)=self%lvl(n)%imax_ ; call MPI_allreduce(buf,max_ ,self%nproc,MPI_INTEGER,MPI_SUM,self%comm,ierr)
            buf=0; buf(self%irank)=self%lvl(n)%imaxo_; call MPI_allreduce(buf,maxo_,self%nproc,MPI_INTEGER,MPI_SUM,self%comm,ierr)
            ! Check if I'm a sender
            if (self%lvl(n)%ncello_.gt.0.and.self%lvl(n)%imin_.le.self%lvl(n)%imax_) then
               ! Set indices
               my_min=self%lvl(n)%imin_
               my_max=self%lvl(n)%imax_
               ! Find receivers
               do ip=1,self%npx
                  ! Get processor rank
                  hisrank=self%rank(ip,self%jproc,self%kproc)
                  ! Check if he can be a receiver
                  if (ncello_rank(hisrank).le.0) cycle
                  ! Test if he's a receiver - careful about periodicity
                  if (self%xper) then
                     if (mod(maxo_(hisrank),self%lvl(n)%nx).eq.mod(my_min,self%lvl(n)%nx)) self%lvl(n)%send_xm(ip)=.true.
                     if (mod(mino_(hisrank),self%lvl(n)%nx).eq.mod(my_max,self%lvl(n)%nx)) self%lvl(n)%send_xp(ip)=.true.
                  else
                     if (maxo_(hisrank).eq.my_min) self%lvl(n)%send_xm(ip)=.true.
                     if (mino_(hisrank).eq.my_max) self%lvl(n)%send_xp(ip)=.true.
                  end if
               end do
            end if
            ! Count number of send operations
            self%lvl(n)%nsend_xm=0; self%lvl(n)%nsend_xp=0
            do ip=1,self%npx
               if (self%lvl(n)%send_xm(ip)) self%lvl(n)%nsend_xm=self%lvl(n)%nsend_xm+1
               if (self%lvl(n)%send_xp(ip)) self%lvl(n)%nsend_xp=self%lvl(n)%nsend_xp+1
            end do
            ! Check if I'm a receiver
            if (self%lvl(n)%ncello_.gt.0) then
               ! Set indices
               my_max=self%lvl(n)%imaxo_
               my_min=self%lvl(n)%imino_
               ! Find sender
               do ip=1,self%npx
                  ! Get processor rank
                  hisrank=self%rank(ip,self%jproc,self%kproc)
                  ! Check if he can be a sender
                  if (ncello_rank(hisrank).le.0.or.min_(hisrank).gt.max_(hisrank)) cycle
                  ! Test if he's the sender - careful about periodicity
                  if (self%xper) then
                     if (mod(min_(hisrank),self%lvl(n)%nx).eq.mod(my_max,self%lvl(n)%nx)) self%lvl(n)%recv_xm=hisrank-1
                     if (mod(max_(hisrank),self%lvl(n)%nx).eq.mod(my_min,self%lvl(n)%nx)) self%lvl(n)%recv_xp=hisrank-1
                  else
                     if (min_(hisrank).eq.my_max) self%lvl(n)%recv_xm=hisrank-1
                     if (max_(hisrank).eq.my_min) self%lvl(n)%recv_xp=hisrank-1
                  end if
               end do
            end if
            ! Y DIRECTION ----------------------------------------------------------------------------------
            ! Allocate send and receive arrays
            allocate(self%lvl(n)%send_ym(self%npy)); self%lvl(n)%send_ym=.false.; self%lvl(n)%recv_ym=MPI_PROC_NULL
            allocate(self%lvl(n)%send_yp(self%npy)); self%lvl(n)%send_yp=.false.; self%lvl(n)%recv_yp=MPI_PROC_NULL
            ! Gather indices
            buf=0; buf(self%irank)=self%lvl(n)%jmino_; call MPI_allreduce(buf,mino_,self%nproc,MPI_INTEGER,MPI_SUM,self%comm,ierr)
            buf=0; buf(self%irank)=self%lvl(n)%jmin_ ; call MPI_allreduce(buf,min_ ,self%nproc,MPI_INTEGER,MPI_SUM,self%comm,ierr)
            buf=0; buf(self%irank)=self%lvl(n)%jmax_ ; call MPI_allreduce(buf,max_ ,self%nproc,MPI_INTEGER,MPI_SUM,self%comm,ierr)
            buf=0; buf(self%irank)=self%lvl(n)%jmaxo_; call MPI_allreduce(buf,maxo_,self%nproc,MPI_INTEGER,MPI_SUM,self%comm,ierr)
            ! Check if I'm a sender
            if (self%lvl(n)%ncello_.gt.0.and.self%lvl(n)%jmin_.le.self%lvl(n)%jmax_) then
               ! Set indices
               my_min=self%lvl(n)%jmin_
               my_max=self%lvl(n)%jmax_
               ! Find receivers
               do jp=1,self%npy
                  ! Get processor rank
                  hisrank=self%rank(self%iproc,jp,self%kproc)
                  ! Check if he can be a receiver
                  if (ncello_rank(hisrank).le.0) cycle
                  ! Test if he's a receiver - careful about periodicity
                  if (self%yper) then
                     if (mod(maxo_(hisrank),self%lvl(n)%ny).eq.mod(my_min,self%lvl(n)%ny)) self%lvl(n)%send_ym(jp)=.true.
                     if (mod(mino_(hisrank),self%lvl(n)%ny).eq.mod(my_max,self%lvl(n)%ny)) self%lvl(n)%send_yp(jp)=.true.
                  else
                     if (maxo_(hisrank).eq.my_min) self%lvl(n)%send_ym(jp)=.true.
                     if (mino_(hisrank).eq.my_max) self%lvl(n)%send_yp(jp)=.true.
                  end if
               end do
            end if
            ! Count number of send operations
            self%lvl(n)%nsend_ym=0; self%lvl(n)%nsend_yp=0
            do jp=1,self%npy
               if (self%lvl(n)%send_ym(jp)) self%lvl(n)%nsend_ym=self%lvl(n)%nsend_ym+1
               if (self%lvl(n)%send_yp(jp)) self%lvl(n)%nsend_yp=self%lvl(n)%nsend_yp+1
            end do
            ! Check if I'm a receiver
            if (self%lvl(n)%ncello_.gt.0) then
               ! Set indices
               my_max=self%lvl(n)%jmaxo_
               my_min=self%lvl(n)%jmino_
               ! Find sender
               do jp=1,self%npy
                  ! Get processor rank
                  hisrank=self%rank(self%iproc,jp,self%kproc)
                  ! Check if he can be a sender
                  if (ncello_rank(hisrank).le.0.or.min_(hisrank).gt.max_(hisrank)) cycle
                  ! Test if he's the sender - careful about periodicity
                  if (self%yper) then
                     if (mod(min_(hisrank),self%lvl(n)%ny).eq.mod(my_max,self%lvl(n)%ny)) self%lvl(n)%recv_ym=hisrank-1
                     if (mod(max_(hisrank),self%lvl(n)%ny).eq.mod(my_min,self%lvl(n)%ny)) self%lvl(n)%recv_yp=hisrank-1
                  else
                     if (min_(hisrank).eq.my_max) self%lvl(n)%recv_ym=hisrank-1
                     if (max_(hisrank).eq.my_min) self%lvl(n)%recv_yp=hisrank-1
                  end if
               end do
            end if
            ! Z DIRECTION ----------------------------------------------------------------------------------
            ! Allocate send and receive arrays
            allocate(self%lvl(n)%send_zm(self%npz)); self%lvl(n)%send_zm=.false.; self%lvl(n)%recv_zm=MPI_PROC_NULL
            allocate(self%lvl(n)%send_zp(self%npz)); self%lvl(n)%send_zp=.false.; self%lvl(n)%recv_zp=MPI_PROC_NULL
            ! Gather indices
            buf=0; buf(self%irank)=self%lvl(n)%kmino_; call MPI_allreduce(buf,mino_,self%nproc,MPI_INTEGER,MPI_SUM,self%comm,ierr)
            buf=0; buf(self%irank)=self%lvl(n)%kmin_ ; call MPI_allreduce(buf,min_ ,self%nproc,MPI_INTEGER,MPI_SUM,self%comm,ierr)
            buf=0; buf(self%irank)=self%lvl(n)%kmax_ ; call MPI_allreduce(buf,max_ ,self%nproc,MPI_INTEGER,MPI_SUM,self%comm,ierr)
            buf=0; buf(self%irank)=self%lvl(n)%kmaxo_; call MPI_allreduce(buf,maxo_,self%nproc,MPI_INTEGER,MPI_SUM,self%comm,ierr)
            ! Check if I'm a sender
            if (self%lvl(n)%ncello_.gt.0.and.self%lvl(n)%kmin_.le.self%lvl(n)%kmax_) then
               ! Set indices
               my_min=self%lvl(n)%kmin_
               my_max=self%lvl(n)%kmax_
               ! Find receivers
               do kp=1,self%npz
                  ! Get processor rank
                  hisrank=self%rank(self%iproc,self%jproc,kp)
                  ! Check if he can be a receiver
                  if (ncello_rank(hisrank).le.0) cycle
                  ! Test if he's a receiver - careful about periodicity
                  if (self%zper) then
                     if (mod(maxo_(hisrank),self%lvl(n)%nz).eq.mod(my_min,self%lvl(n)%nz)) self%lvl(n)%send_zm(kp)=.true.
                     if (mod(mino_(hisrank),self%lvl(n)%nz).eq.mod(my_max,self%lvl(n)%nz)) self%lvl(n)%send_zp(kp)=.true.
                  else
                     if (maxo_(hisrank).eq.my_min) self%lvl(n)%send_zm(kp)=.true.
                     if (mino_(hisrank).eq.my_max) self%lvl(n)%send_zp(kp)=.true.
                  end if
               end do
            end if
            ! Count number of send operations
            self%lvl(n)%nsend_zm=0; self%lvl(n)%nsend_zp=0
            do kp=1,self%npz
               if (self%lvl(n)%send_zm(kp)) self%lvl(n)%nsend_zm=self%lvl(n)%nsend_zm+1
               if (self%lvl(n)%send_zp(kp)) self%lvl(n)%nsend_zp=self%lvl(n)%nsend_zp+1
            end do
            ! Check if I'm a receiver
            if (self%lvl(n)%ncello_.gt.0) then
               ! Set indices
               my_max=self%lvl(n)%kmaxo_
               my_min=self%lvl(n)%kmino_
               ! Find sender
               do kp=1,self%npz
                  ! Get processor rank
                  hisrank=self%rank(self%iproc,self%jproc,kp)
                  ! Check if he can be a sender
                  if (ncello_rank(hisrank).le.0.or.min_(hisrank).gt.max_(hisrank)) cycle
                  ! Test if he's the sender - careful about periodicity
                  if (self%zper) then
                     if (mod(min_(hisrank),self%lvl(n)%nz).eq.mod(my_max,self%lvl(n)%nz)) self%lvl(n)%recv_zm=hisrank-1
                     if (mod(max_(hisrank),self%lvl(n)%nz).eq.mod(my_min,self%lvl(n)%nz)) self%lvl(n)%recv_zp=hisrank-1
                  else
                     if (min_(hisrank).eq.my_max) self%lvl(n)%recv_zm=hisrank-1
                     if (max_(hisrank).eq.my_min) self%lvl(n)%recv_zp=hisrank-1
                  end if
               end do
            end if
         end do
         ! Deallocate work arrays
         deallocate(buf,mino_,min_,max_,maxo_,ncello_rank)
      end block initialize_comm
      
   end function bbmg_from_pgrid
   
   
   !> Initialization of BBMG object - this sets the solver parameters
   subroutine init_solver(this)
      use messager, only: die
      implicit none
      class(bbmg), intent(inout) :: this
      
      
      
   end subroutine init_solver
   
   
   !> Solve the linear system iteratively
   subroutine solve(this)
      use messager, only: die
      use param,    only: verbose
      implicit none
      class(bbmg), intent(inout) :: this
      
      
      ! If verbose run, log and or print cvg history
      !if (verbose.gt.2) call this%print_cvg
      
   end subroutine solve
   
   
   !> Update of the operators across all levels
   !> This routine assumes that the level(1) opr has been populated in the interior
   subroutine update_operator(this)
      use messager, only: die
      implicit none
      class(bbmg), intent(inout) :: this
      integer :: n,i,j,k
      ! Reset all ghost values to identity and sync boundaries
      do k=this%lvl(1)%kmino_,this%lvl(1)%kmaxo_
         do j=this%lvl(1)%jmino_,this%lvl(1)%jmaxo_
            do i=this%lvl(1)%imino_,this%lvl(1)%imaxo_
               ! Interior cells should have been properly set already
               if (i.ge.this%lvl(1)%imin_.and.i.le.this%lvl(1)%imax_.and. &
               &   j.ge.this%lvl(1)%jmin_.and.j.le.this%lvl(1)%jmax_.and. &
               &   k.ge.this%lvl(1)%kmin_.and.k.le.this%lvl(1)%kmax_) cycle
               ! Set operator to identity
               this%lvl(1)%opr(:,:,:,i,j,k)=0.0_WP
               this%lvl(1)%opr(0,0,0,i,j,k)=1.0_WP
            end do
         end do
      end do
      call this%msync(A=this%lvl(1)%opr,n=1)
      ! Loop over levels
      do n=2,this%nlvl
         ! Recompute prolongation
         call this%recompute_prolongation(n-1)
         ! Recompute restriction
         call this%recompute_restriction(n)
         ! Recompute operator
         call this%recompute_operator(n)
      end do
      ! Prepare direct solver
      !call this%recompute_direct_operator
   end subroutine update_operator
   
   
   !> Recompute prolongation
   subroutine recompute_prolongation(this,n)
      implicit none
      class(bbmg), intent(inout) :: this
      integer, intent(in) :: n
      integer :: n1,n2,n3
      integer :: fi,fj,fk
      integer :: li,lj,lk
      integer :: pi,pj,pk
      integer :: si,sj,sk
      real(WP), dimension(-1:+1)       :: opr1d
      real(WP), dimension(-1:+1,-1:+1) :: opr2d
      real(WP), parameter :: ratio=10.0_WP
      ! Sweep 1 - identity and 1D interpolations
      do fk=this%lvl(n)%kmino_,this%lvl(n)%kmaxo_
         do fj=this%lvl(n)%jmino_,this%lvl(n)%jmaxo_
            do fi=this%lvl(n)%imino_,this%lvl(n)%imaxo_
               ! Check stencil size
               n1=this%pmodx(fi,n); n2=this%pmody(fj,n); n3=this%pmodz(fk,n)
               ! Reset prolongation operator
               this%lvl(n)%c2f(:,:,:,fi,fj,fk)=0.0_WP
               ! Compute prolongation coefficients
               select case (n1+n2+n3)
               case (0) ! Fine and coarse locations coincide => Identity
                  this%lvl(n)%c2f(0,0,0,fi,fj,fk)=1.0_WP
               case (1) ! Fine location lies on coarse x/y/z-line => Which one?
                  if      (n1.eq.1) then ! In x
                     ! Aggregate operator
                     opr1d=0.0_WP
                     do lk=-1,+1
                        do lj=-1,+1
                           do li=-1,+1
                              if (ratio*abs(this%lvl(n)%opr(li,0,0,fi,fj,fk)).ge.abs(this%lvl(n)%opr(li,lj,lk,fi,fj,fk))) then
                                 opr1d(li)=opr1d(li)+this%lvl(n)%opr(li,lj,lk,fi,fj,fk)
                              else
                                 opr1d( 0)=opr1d( 0)+this%lvl(n)%opr(li,lj,lk,fi,fj,fk)
                              end if
                           end do
                        end do
                     end do
                     opr1d(0)=sign(max(abs(opr1d(0)),epsilon(1.0_WP)),opr1d(0))
                     ! Form x-interpolation
                     this%lvl(n)%c2f(0,0,0,fi,fj,fk)=-opr1d(-1)/opr1d(0)
                     this%lvl(n)%c2f(1,0,0,fi,fj,fk)=-opr1d(+1)/opr1d(0)
                  else if (n2.eq.1) then ! In y
                     ! Aggregate operator
                     opr1d=0.0_WP
                     do lk=-1,+1
                        do lj=-1,+1
                           do li=-1,+1
                              if (ratio*abs(this%lvl(n)%opr(0,lj,0,fi,fj,fk)).ge.abs(this%lvl(n)%opr(li,lj,lk,fi,fj,fk))) then
                                 opr1d(lj)=opr1d(lj)+this%lvl(n)%opr(li,lj,lk,fi,fj,fk)
                              else
                                 opr1d( 0)=opr1d( 0)+this%lvl(n)%opr(li,lj,lk,fi,fj,fk)
                              end if
                           end do
                        end do
                     end do
                     opr1d(0)=sign(max(abs(opr1d(0)),epsilon(1.0_WP)),opr1d(0))
                     ! Form y-interpolation
                     this%lvl(n)%c2f(0,0,0,fi,fj,fk)=-opr1d(-1)/opr1d(0)
                     this%lvl(n)%c2f(0,1,0,fi,fj,fk)=-opr1d(+1)/opr1d(0)
                  else if (n3.eq.1) then ! In z
                     ! Aggregate operator
                     opr1d=0.0_WP
                     do lk=-1,+1
                        do lj=-1,+1
                           do li=-1,+1
                              if (ratio*abs(this%lvl(n)%opr(0,0,lk,fi,fj,fk)).ge.abs(this%lvl(n)%opr(li,lj,lk,fi,fj,fk))) then
                                 opr1d(lk)=opr1d(lk)+this%lvl(n)%opr(li,lj,lk,fi,fj,fk)
                              else
                                 opr1d( 0)=opr1d( 0)+this%lvl(n)%opr(li,lj,lk,fi,fj,fk)
                              end if
                           end do
                        end do
                     end do
                     opr1d(0)=sign(max(abs(opr1d(0)),epsilon(1.0_WP)),opr1d(0))
                     ! Form z-interpolation
                     this%lvl(n)%c2f(0,0,0,fi,fj,fk)=-opr1d(-1)/opr1d(0)
                     this%lvl(n)%c2f(0,0,1,fi,fj,fk)=-opr1d(+1)/opr1d(0)
                  end if
               case default ! Do not treat face of cube centers yet
                  cycle
               end select
            end do
         end do
      end do
      ! Sweep 2 - 2D face interpolations
      do fk=this%lvl(n)%kmino_,this%lvl(n)%kmaxo_
         do fj=this%lvl(n)%jmino_,this%lvl(n)%jmaxo_
            do fi=this%lvl(n)%imino_,this%lvl(n)%imaxo_
               ! Check stencil size
               n1=this%pmodx(fi,n); n2=this%pmody(fj,n); n3=this%pmodz(fk,n)
               ! Compute prolongation coefficients
               select case (n1+n2+n3)
               case (2) ! Fine location lies on coarse xy/yz/zx-plane => Which one?
                  if      (n1.eq.0) then ! In y/z
                     ! yz-plane - skip if in y/z ghost cell
                     if (fj.eq.this%lvl(n)%jmino_.or.fj.eq.this%lvl(n)%jmaxo_.or.fk.eq.this%lvl(n)%kmino_.or.fk.eq.this%lvl(n)%kmaxo_) cycle
                     ! Aggregate operator
                     opr2d=0.0_WP
                     do lk=-1,+1
                        do lj=-1,+1
                           do li=-1,+1
                              if (ratio*abs(this%lvl(n)%opr(0,lj,lk,fi,fj,fk)).ge.abs(this%lvl(n)%opr(li,lj,lk,fi,fj,fk))) then
                                 opr2d(lj,lk)=opr2d(lj,lk)+this%lvl(n)%opr(li,lj,lk,fi,fj,fk)
                              else
                                 opr2d( 0, 0)=opr2d( 0, 0)+this%lvl(n)%opr(li,lj,lk,fi,fj,fk)
                              end if
                           end do
                        end do
                     end do
                     opr2d(0,0)=sign(max(abs(opr2d(0,0)),epsilon(1.0_WP)),opr2d(0,0))
                     ! Loop over operator stencil
                     do lk=-1,+1
                        do lj=-1,+1
                           ! Skip middle of stencil
                           if (lj.eq.0.and.lk.eq.0) cycle
                           ! Loop over prolongation stencil
                           do pk=0,1-abs(lk)
                              do pj=0,1-abs(lj)
                                 sj=(lj+1)/2+pj; sk=(lk+1)/2+pk
                                 this%lvl(n)%c2f(0,sj,sk,fi,fj,fk)=this%lvl(n)%c2f(0,sj,sk,fi,fj,fk)-opr2d(lj,lk)/opr2d(0,0)*this%lvl(n)%c2f(0,pj,pk,fi,fj+lj,fk+lk)
                              end do
                           end do
                        end do
                     end do
                  else if (n2.eq.0) then ! In z/x
                     ! zx-plane - skip if in z/x ghost cell
                     if (fk.eq.this%lvl(n)%kmino_.or.fk.eq.this%lvl(n)%kmaxo_.or.fi.eq.this%lvl(n)%imino_.or.fi.eq.this%lvl(n)%imaxo_) cycle
                     ! Aggregate operator
                     opr2d=0.0_WP
                     do lk=-1,+1
                        do lj=-1,+1
                           do li=-1,+1
                              if (ratio*abs(this%lvl(n)%opr(li,0,lk,fi,fj,fk)).ge.abs(this%lvl(n)%opr(li,lj,lk,fi,fj,fk))) then
                                 opr2d(li,lk)=opr2d(li,lk)+this%lvl(n)%opr(li,lj,lk,fi,fj,fk)
                              else
                                 opr2d( 0, 0)=opr2d( 0, 0)+this%lvl(n)%opr(li,lj,lk,fi,fj,fk)
                              end if
                           end do
                        end do
                     end do
                     opr2d(0,0)=sign(max(abs(opr2d(0,0)),epsilon(1.0_WP)),opr2d(0,0))
                     ! Loop over operator stencil
                     do lk=-1,+1
                        do li=-1,+1
                           ! Skip middle of stencil
                           if (li.eq.0.and.lk.eq.0) cycle
                           ! Loop over prolongation stencil
                           do pk=0,1-abs(lk)
                              do pi=0,1-abs(li)
                                 si=(li+1)/2+pi; sk=(lk+1)/2+pk
                                 this%lvl(n)%c2f(si,0,sk,fi,fj,fk)=this%lvl(n)%c2f(si,0,sk,fi,fj,fk)-opr2d(li,lk)/opr2d(0,0)*this%lvl(n)%c2f(pi,0,pk,fi+li,fj,fk+lk)
                              end do
                           end do
                        end do
                     end do
                  else if (n3.eq.0) then ! In x/y
                     ! xy-plane - skip if in x/y ghost cell
                     if (fi.eq.this%lvl(n)%imino_.or.fi.eq.this%lvl(n)%imaxo_.or.fj.eq.this%lvl(n)%jmino_.or.fj.eq.this%lvl(n)%jmaxo_) cycle
                     ! Aggregate operator
                     opr2d=0.0_WP
                     do lk=-1,+1
                        do lj=-1,+1
                           do li=-1,+1
                              if (ratio*abs(this%lvl(n)%opr(li,lj,0,fi,fj,fk)).ge.abs(this%lvl(n)%opr(li,lj,lk,fi,fj,fk))) then
                                 opr2d(li,lj)=opr2d(li,lj)+this%lvl(n)%opr(li,lj,lk,fi,fj,fk)
                              else
                                 opr2d( 0, 0)=opr2d( 0, 0)+this%lvl(n)%opr(li,lj,lk,fi,fj,fk)
                              end if
                           end do
                        end do
                     end do
                     opr2d(0,0)=sign(max(abs(opr2d(0,0)),epsilon(1.0_WP)),opr2d(0,0))
                     ! Loop over operator stencil
                     do lj=-1,+1
                        do li=-1,+1
                           ! Skip middle of stencil
                           if (li.eq.0.and.lj.eq.0) cycle
                           ! Loop over prolongation stencil
                           do pj=0,1-abs(lj)
                              do pi=0,1-abs(li)
                                 si=(li+1)/2+pi; sj=(lj+1)/2+pj
                                 this%lvl(n)%c2f(si,sj,0,fi,fj,fk)=this%lvl(n)%c2f(si,sj,0,fi,fj,fk)-opr2d(li,lj)/opr2d(0,0)*this%lvl(n)%c2f(pi,pj,0,fi+li,fj+lj,fk)
                              end do
                           end do
                        end do
                     end do
                  end if
               case default ! Identity & line already treated, and don't treat cube centers yet
                  cycle
               end select
            end do
         end do
      end do
      ! Sweep 3 - 3D cube interpolation
      do fk=this%lvl(n)%kmin_,this%lvl(n)%kmax_
         do fj=this%lvl(n)%jmin_,this%lvl(n)%jmax_
            do fi=this%lvl(n)%imin_,this%lvl(n)%imax_
               ! Check stencil size
               n1=this%pmodx(fi,n); n2=this%pmody(fj,n); n3=this%pmodz(fk,n)
               ! Compute prolongation coefficients
               select case (n1+n2+n3)
               case (3) ! Fine location lies on xyz-cube
                  ! Loop over operator stencil
                  do lk=-1,+1
                     do lj=-1,+1
                        do li=-1,+1
                           ! Skip middle of stencil
                           if (li.eq.0.and.lj.eq.0.and.lk.eq.0) cycle
                           ! Loop over prolongation stencil
                           do pk=0,1-abs(lk)
                              do pj=0,1-abs(lj)
                                 do pi=0,1-abs(li)
                                    si=(li+1)/2+pi; sj=(lj+1)/2+pj; sk=(lk+1)/2+pk
                                    this%lvl(n)%c2f(si,sj,sk,fi,fj,fk)=this%lvl(n)%c2f(si,sj,sk,fi,fj,fk)-this%lvl(n)%opr(li,lj,lk,fi,fj,fk)&
                                    &/sign(max(abs(this%lvl(n)%opr(0,0,0,fi,fj,fk)),epsilon(1.0_WP)),this%lvl(n)%opr(0,0,0,fi,fj,fk))*this%lvl(n)%c2f(pi,pj,pk,fi+li,fj+lj,fk+lk)
                                 end do
                              end do
                           end do
                        end do
                     end do
                  end do
               case default ! All other cases have been treated already
                  cycle
               end select
            end do
         end do
      end do
      ! Synchronize prolongation
      call this%msync(this%lvl(n)%c2f,n)
   end subroutine recompute_prolongation
   
   
   !> Recompute restriction
   subroutine recompute_restriction(this,n)
      implicit none
      class(bbmg), intent(inout) :: this
      integer, intent(in) :: n
      integer :: fi,fj,fk
      integer :: ci,cj,ck
      integer :: ri,rj,rk
      integer :: pi,pj,pk
      ! Zero out restriction
      this%lvl(n)%f2c=0.0_WP
      ! Loop over coarse cells
      do ck=this%lvl(n)%kmin_,this%lvl(n)%kmax_
         do cj=this%lvl(n)%jmin_,this%lvl(n)%jmax_
            do ci=this%lvl(n)%imin_,this%lvl(n)%imax_
               ! Find corresponding fine cell
               fi=2*ci-1; fj=2*cj-1; fk=2*ck-1
               ! Loop over restriction stencil
               do rk=-this%pmodz(fk-1,n-1),+this%pmodz(fk+1,n-1)
                  do rj=-this%pmody(fj-1,n-1),+this%pmody(fj+1,n-1)
                     do ri=-this%pmodx(fi-1,n-1),+this%pmodx(fi+1,n-1)
                        pi=(1-ri)/2; pj=(1-rj)/2; pk=(1-rk)/2
                        this%lvl(n)%f2c(ri,rj,rk,ci,cj,ck)=this%lvl(n-1)%c2f(pi,pj,pk,fi+ri,fj+rj,fk+rk)
                     end do
                  end do
               end do
            end do
         end do
      end do
      ! Synchronize restriction
      call this%msync(this%lvl(n)%f2c,n)
   end subroutine recompute_restriction
   
   
   !> Recompute operator
   subroutine recompute_operator(this,n)
      implicit none
      class(bbmg), intent(inout) :: this
      integer, intent(in) :: n
      integer :: fi,fj,fk
      integer :: ci,cj,ck
      integer :: ri,rj,rk
      integer :: li,lj,lk
      integer :: pi,pj,pk
      integer :: si,sj,sk
      ! Zero out oprc2f
      this%lvl(n-1)%oprc2f=0.0_WP
      ! Loop over fine cells
      do fk=this%lvl(n-1)%kmin_,this%lvl(n-1)%kmax_
         do fj=this%lvl(n-1)%jmin_,this%lvl(n-1)%jmax_
            do fi=this%lvl(n-1)%imin_,this%lvl(n-1)%imax_
               ! Loop over laplacian stencil
               do lk=-1,+1
                  do lj=-1,+1
                     do li=-1,+1
                        ! Loop over prolongation stencil
                        do pk=0,this%pmodz(fk+lk,n-1)
                           do pj=0,this%pmody(fj+lj,n-1)
                              do pi=0,this%pmodx(fi+li,n-1)
                                 si=pi+(1-this%pmodx(fi,n-1))*min(li,0)+(1-this%pmodx(fi+li,n-1))*max(li,0)
                                 sj=pj+(1-this%pmody(fj,n-1))*min(lj,0)+(1-this%pmody(fj+lj,n-1))*max(lj,0)
                                 sk=pk+(1-this%pmodz(fk,n-1))*min(lk,0)+(1-this%pmodz(fk+lk,n-1))*max(lk,0)
                                 this%lvl(n-1)%oprc2f(si,sj,sk,fi,fj,fk)=this%lvl(n-1)%oprc2f(si,sj,sk,fi,fj,fk)+this%lvl(n-1)%opr(li,lj,lk,fi,fj,fk)*this%lvl(n-1)%c2f(pi,pj,pk,fi+li,fj+lj,fk+lk)
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      ! Synchronize operator
      call this%msync(this%lvl(n-1)%oprc2f,n-1)
      ! Zero out operator
      this%lvl(n)%opr=0.0_WP
      ! Loop over coarse cells
      do ck=this%lvl(n)%kmin_,this%lvl(n)%kmax_
         do cj=this%lvl(n)%jmin_,this%lvl(n)%jmax_
            do ci=this%lvl(n)%imin_,this%lvl(n)%imax_
               ! Find corresponding fine cell
               fi=2*ci-1; fj=2*cj-1; fk=2*ck-1
               ! Loop over restriction stencil
               do rk=-this%pmodz(fk-1,n-1),+this%pmodz(fk+1,n-1)
                  do rj=-this%pmody(fj-1,n-1),+this%pmody(fj+1,n-1)
                     do ri=-this%pmodx(fi-1,n-1),+this%pmodx(fi+1,n-1)
                        ! Loop over lapc2f stencil
                        do lk=-1+this%pmodz(fk+rk,n-1),+1
                           do lj=-1+this%pmody(fj+rj,n-1),+1
                              do li=-1+this%pmodx(fi+ri,n-1),+1
                                 si=(ri-1)/2+li; sj=(rj-1)/2+lj; sk=(rk-1)/2+lk
                                 this%lvl(n)%opr(si,sj,sk,ci,cj,ck)=this%lvl(n)%opr(si,sj,sk,ci,cj,ck)+this%lvl(n)%f2c(ri,rj,rk,ci,cj,ck)*this%lvl(n-1)%oprc2f(li,lj,lk,fi+ri,fj+rj,fk+rk)
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      ! Synchronize operator
      call this%msync(this%lvl(n)%opr,n)
   end subroutine recompute_operator
   
   
   !> Solve the linear system using RBGS
   subroutine smooth_rbgs(this)
      implicit none
      class(bbmg), intent(inout) :: this
      ! integer :: i,j,k,col,st
      ! integer :: colmin,colmax,coldir
      ! real(WP) :: val
      ! ! Reset iteration and error
      ! this%it=0
      ! this%rerr=huge(1.0_WP)
      ! this%aerr=huge(1.0_WP)
      ! ! Loop until done
      ! do while (this%it.lt.this%maxit.and.this%rerr.ge.this%rcvg.and.this%aerr.ge.this%acvg)
      !    ! Increment counter
      !    this%it=this%it+1
      !    ! Alternate sweep direction
      !    if (mod(this%it,2).eq.0) then
      !       colmin=8; colmax=1; coldir=-1 !< Negative sweep
      !    else
      !       colmin=1; colmax=8; coldir=+1 !< Positive sweep
      !    end if
      !    ! Loop over colors
      !    do col=colmin,colmax,coldir
      !       ! Loop over domain
      !       do k=this%cfg%kmin_+mod((col-1)/4,2),this%cfg%kmax_,2
      !          do j=this%cfg%jmin_+mod((col-1)/2,2),this%cfg%jmax_,2
      !             do i=this%cfg%imin_+mod((col-1)/1,2),this%cfg%imax_,2
      !                ! Gauss-Seidel step
      !                if (abs(this%opr(1,i,j,k)).gt.0.0_WP) then
      !                   val=this%rhs(i,j,k)
      !                   do st=2,this%nst
      !                      val=val-this%opr(st,i,j,k)*this%sol(i+this%stc(st,1),j+this%stc(st,2),k+this%stc(st,3))
      !                   end do
      !                   this%sol(i,j,k)=val/this%opr(1,i,j,k)
      !                end if
      !             end do
      !          end do
      !       end do
      !       ! Communicate solution
      !       call this%cfg%sync(this%sol)
      !    end do
      ! end do
   end subroutine smooth_rbgs
   
   
   !> Synchronization of overlap cells for a vector
   subroutine vsync(this,A,n)
      use parallel, only: MPI_REAL_WP
      implicit none
      class(bbmg), intent(inout) :: this
      real(WP), dimension(:,:,:), intent(inout) :: A !< Needs to be (lvl(n)%nxo_,lvl(n)%nyo_,lvl(n)%nzo_)
      integer, intent(in) :: n                       !< Level at which the sync is happening
      real(WP), dimension(:,:,:), allocatable :: buf1,buf2
      integer :: icount,count,dest,ierr
      integer :: no,n1,n2,n3,ip,jp,kp
      type(MPI_Request), dimension(:), allocatable :: request
      ! Shortcut to sizes
      no=this%no
      n1=this%lvl(n)%nxo_; n2=this%lvl(n)%nyo_; n3=this%lvl(n)%nzo_
      ! Communicate in X =======================================================
      ! Initialize buffer and size
      allocate(buf1(no,n2,n3),buf2(no,n2,n3)); icount=no*n2*n3
      ! Allocate request
      allocate(request(max(this%lvl(n)%nsend_xm,this%lvl(n)%nsend_xp)+1))
      ! MINUS DIRECTION ------
      ! Loop over potential destinations
      count=0
      do ip=1,this%npx
         ! Test if I need to send
         if (this%lvl(n)%send_xm(ip)) then
            ! Get processor rank
            dest=this%rank(ip,this%jproc,this%kproc)-1
            ! Copy left buffer
            buf1(:,:,:)=A(1+no:1+no+no-1,:,:)
            ! Send left buffer to destination
            count=count+1
            call MPI_ISEND(buf1,icount,MPI_REAL_WP,dest,0,this%comm,request(count),ierr)
         end if
      end do
      ! Receive from potential sender
      count=count+1
      call MPI_IRECV(buf2,icount,MPI_REAL_WP,this%lvl(n)%recv_xm,0,this%comm,request(count),ierr)
      ! Wait for completion
      call MPI_WAITALL(count,request(1:count),MPI_STATUSES_IGNORE,ierr)
      ! Paste left buffer to right
      if (this%lvl(n)%recv_xm.ne.MPI_PROC_NULL) A(n1-no+1:n1,:,:)=buf2(:,:,:)
      ! PLUS DIRECTION ------
      ! Loop over potential destinations
      count=0
      do ip=1,this%npx
         ! Test if I need to send
         if (this%lvl(n)%send_xp(ip)) then
            ! Get processor rank
            dest=this%rank(ip,this%jproc,this%kproc)-1
            ! Copy right buffer
            buf1(:,:,:)=A(n1-no-no+1:n1-no,:,:)
            ! Send right buffer to destination
            count=count+1
            call MPI_ISEND(buf1,icount,MPI_REAL_WP,dest,0,this%comm,request(count),ierr)
         end if
      end do
      ! Receive from potential sender
      count=count+1
      call MPI_IRECV(buf2,icount,MPI_REAL_WP,this%lvl(n)%recv_xp,0,this%comm,request(count),ierr)
      ! Wait for completion
      call MPI_WAITALL(count,request(1:count),MPI_STATUSES_IGNORE,ierr)
      ! Paste left buffer to right
      if (this%lvl(n)%recv_xp.ne.MPI_PROC_NULL) A(1:no,:,:)=buf2(:,:,:)
      ! Deallocate
      deallocate(buf1,buf2,request)
      ! Communicate in Y =======================================================
      ! Initialize buffer and size
      allocate(buf1(n1,no,n3),buf2(n1,no,n3)); icount=n1*no*n3
      ! Allocate request
      allocate(request(max(this%lvl(n)%nsend_ym,this%lvl(n)%nsend_yp)+1))
      ! MINUS DIRECTION ------
      ! Loop over potential destinations
      count=0
      do jp=1,this%npy
         ! Test if I need to send
         if (this%lvl(n)%send_ym(jp)) then
            ! Get processor rank
            dest=this%rank(this%iproc,jp,this%kproc)-1
            ! Copy left buffer
            buf1(:,:,:)=A(:,1+no:1+no+no-1,:)
            ! Send left buffer to destination
            count=count+1
            call MPI_ISEND(buf1,icount,MPI_REAL_WP,dest,0,this%comm,request(count),ierr)
         end if
      end do
      ! Receive from potential sender
      count=count+1
      call MPI_IRECV(buf2,icount,MPI_REAL_WP,this%lvl(n)%recv_ym,0,this%comm,request(count),ierr)
      ! Wait for completion
      call MPI_WAITALL(count,request(1:count),MPI_STATUSES_IGNORE,ierr)
      ! Paste left buffer to right
      if (this%lvl(n)%recv_ym.ne.MPI_PROC_NULL) A(:,n2-no+1:n2,:)=buf2(:,:,:)
      ! PLUS DIRECTION ------
      ! Loop over potential destinations
      count=0
      do jp=1,this%npy
         ! Test if I need to send
         if (this%lvl(n)%send_yp(jp)) then
            ! Get processor rank
            dest=this%rank(this%iproc,jp,this%kproc)-1
            ! Copy right buffer
            buf1(:,:,:)=A(:,n2-no-no+1:n2-no,:)
            ! Send right buffer to destination
            count=count+1
            call MPI_ISEND(buf1,icount,MPI_REAL_WP,dest,0,this%comm,request(count),ierr)
         end if
      end do
      ! Receive from potential sender
      count=count+1
      call MPI_IRECV(buf2,icount,MPI_REAL_WP,this%lvl(n)%recv_yp,0,this%comm,request(count),ierr)
      ! Wait for completion
      call MPI_WAITALL(count,request(1:count),MPI_STATUSES_IGNORE,ierr)
      ! Paste left buffer to right
      if (this%lvl(n)%recv_yp.ne.MPI_PROC_NULL) A(:,1:no,:)=buf2(:,:,:)
      ! Deallocate
      deallocate(buf1,buf2,request)
      ! Communicate in Z =======================================================
      ! Initialize buffer and size
      allocate(buf1(n1,n2,no),buf2(n1,n2,no)); icount=n1*n2*no
      ! Allocate request
      allocate(request(max(this%lvl(n)%nsend_zm,this%lvl(n)%nsend_zp)+1))
      ! MINUS DIRECTION ------
      ! Loop over potential destinations
      count=0
      do kp=1,this%npz
         ! Test if I need to send
         if (this%lvl(n)%send_zm(kp)) then
            ! Get processor rank
            dest=this%rank(this%iproc,this%jproc,kp)-1
            ! Copy left buffer
            buf1(:,:,:)=A(:,:,1+no:1+no+no-1)
            ! Send left buffer to destination
            count=count+1
            call MPI_ISEND(buf1,icount,MPI_REAL_WP,dest,0,this%comm,request(count),ierr)
         end if
      end do
      ! Receive from potential sender
      count=count+1
      call MPI_IRECV(buf2,icount,MPI_REAL_WP,this%lvl(n)%recv_zm,0,this%comm,request(count),ierr)
      ! Wait for completion
      call MPI_WAITALL(count,request(1:count),MPI_STATUSES_IGNORE,ierr)
      ! Paste left buffer to right
      if (this%lvl(n)%recv_zm.ne.MPI_PROC_NULL) A(:,:,n3-no+1:n3)=buf2(:,:,:)
      ! PLUS DIRECTION ------
      ! Loop over potential destinations
      count=0
      do kp=1,this%npz
         ! Test if I need to send
         if (this%lvl(n)%send_zp(kp)) then
            ! Get processor rank
            dest=this%rank(this%iproc,this%jproc,kp)-1
            ! Copy right buffer
            buf1(:,:,:)=A(:,:,n3-no-no+1:n3-no)
            ! Send right buffer to destination
            count=count+1
            call MPI_ISEND(buf1,icount,MPI_REAL_WP,dest,0,this%comm,request(count),ierr)
         end if
      end do
      ! Receive from potential sender
      count=count+1
      call MPI_IRECV(buf2,icount,MPI_REAL_WP,this%lvl(n)%recv_zp,0,this%comm,request(count),ierr)
      ! Wait for completion
      call MPI_WAITALL(count,request(1:count),MPI_STATUSES_IGNORE,ierr)
      ! Paste left buffer to right
      if (this%lvl(n)%recv_zp.ne.MPI_PROC_NULL) A(:,:,1:no)=buf2(:,:,:)
      ! Deallocate
      deallocate(buf1,buf2,request)
   end subroutine vsync
   
   
   !> Synchronization of overlap cells for a matrix
   subroutine msync(this,A,n)
      use parallel, only: MPI_REAL_WP
      implicit none
      class(bbmg), intent(inout) :: this
      real(WP), dimension(:,:,:,:,:,:), intent(inout) :: A !< Needs to be (:,:,:,lvl(n)%nxo_,lvl(n)%nyo_,lvl(n)%nzo_)
      integer, intent(in) :: n                             !< Level at which the sync is happening
      real(WP), dimension(:,:,:,:,:,:), allocatable :: buf1,buf2
      integer :: icount,count,dest,ierr
      integer :: no,n1,n2,n3,ip,jp,kp,nstx,nsty,nstz
      type(MPI_Request), dimension(:), allocatable :: request
      ! Shortcut to sizes
      no=this%no; nstx=size(A,DIM=1); nsty=size(A,DIM=2); nstz=size(A,DIM=3)
      n1=this%lvl(n)%nxo_; n2=this%lvl(n)%nyo_; n3=this%lvl(n)%nzo_
      ! Communicate in X =======================================================
      ! Initialize buffer and size
      allocate(buf1(nstx,nsty,nstz,no,n2,n3),buf2(nstx,nsty,nstz,no,n2,n3)); icount=nstx*nsty*nstz*no*n2*n3
      ! Allocate request
      allocate(request(max(this%lvl(n)%nsend_xm,this%lvl(n)%nsend_xp)+1))
      ! MINUS DIRECTION ------
      ! Loop over potential destinations
      count=0
      do ip=1,this%npx
         ! Test if I need to send
         if (this%lvl(n)%send_xm(ip)) then
            ! Get processor rank
            dest=this%rank(ip,this%jproc,this%kproc)-1
            ! Copy left buffer
            buf1(:,:,:,:,:,:)=A(:,:,:,1+no:1+no+no-1,:,:)
            ! Send left buffer to destination
            count=count+1
            call MPI_ISEND(buf1,icount,MPI_REAL_WP,dest,0,this%comm,request(count),ierr)
         end if
      end do
      ! Receive from potential sender
      count=count+1
      call MPI_IRECV(buf2,icount,MPI_REAL_WP,this%lvl(n)%recv_xm,0,this%comm,request(count),ierr)
      ! Wait for completion
      call MPI_WAITALL(count,request(1:count),MPI_STATUSES_IGNORE,ierr)
      ! Paste left buffer to right
      if (this%lvl(n)%recv_xm.ne.MPI_PROC_NULL) A(:,:,:,n1-no+1:n1,:,:)=buf2(:,:,:,:,:,:)
      ! PLUS DIRECTION ------
      ! Loop over potential destinations
      count=0
      do ip=1,this%npx
         ! Test if I need to send
         if (this%lvl(n)%send_xp(ip)) then
            ! Get processor rank
            dest=this%rank(ip,this%jproc,this%kproc)-1
            ! Copy right buffer
            buf1(:,:,:,:,:,:)=A(:,:,:,n1-no-no+1:n1-no,:,:)
            ! Send right buffer to destination
            count=count+1
            call MPI_ISEND(buf1,icount,MPI_REAL_WP,dest,0,this%comm,request(count),ierr)
         end if
      end do
      ! Receive from potential sender
      count=count+1
      call MPI_IRECV(buf2,icount,MPI_REAL_WP,this%lvl(n)%recv_xp,0,this%comm,request(count),ierr)
      ! Wait for completion
      call MPI_WAITALL(count,request(1:count),MPI_STATUSES_IGNORE,ierr)
      ! Paste left buffer to right
      if (this%lvl(n)%recv_xp.ne.MPI_PROC_NULL) A(:,:,:,1:no,:,:)=buf2(:,:,:,:,:,:)
      ! Deallocate
      deallocate(buf1,buf2,request)
      ! Communicate in Y =======================================================
      ! Initialize buffer and size
      allocate(buf1(nstx,nsty,nstz,n1,no,n3),buf2(nstx,nsty,nstz,n1,no,n3)); icount=nstx*nsty*nstz*n1*no*n3
      ! Allocate request
      allocate(request(max(this%lvl(n)%nsend_ym,this%lvl(n)%nsend_yp)+1))
      ! MINUS DIRECTION ------
      ! Loop over potential destinations
      count=0
      do jp=1,this%npy
         ! Test if I need to send
         if (this%lvl(n)%send_ym(jp)) then
            ! Get processor rank
            dest=this%rank(this%iproc,jp,this%kproc)-1
            ! Copy left buffer
            buf1(:,:,:,:,:,:)=A(:,:,:,:,1+no:1+no+no-1,:)
            ! Send left buffer to destination
            count=count+1
            call MPI_ISEND(buf1,icount,MPI_REAL_WP,dest,0,this%comm,request(count),ierr)
         end if
      end do
      ! Receive from potential sender
      count=count+1
      call MPI_IRECV(buf2,icount,MPI_REAL_WP,this%lvl(n)%recv_ym,0,this%comm,request(count),ierr)
      ! Wait for completion
      call MPI_WAITALL(count,request(1:count),MPI_STATUSES_IGNORE,ierr)
      ! Paste left buffer to right
      if (this%lvl(n)%recv_ym.ne.MPI_PROC_NULL) A(:,:,:,:,n2-no+1:n2,:)=buf2(:,:,:,:,:,:)
      ! PLUS DIRECTION ------
      ! Loop over potential destinations
      count=0
      do jp=1,this%npy
         ! Test if I need to send
         if (this%lvl(n)%send_yp(jp)) then
            ! Get processor rank
            dest=this%rank(this%iproc,jp,this%kproc)-1
            ! Copy right buffer
            buf1(:,:,:,:,:,:)=A(:,:,:,:,n2-no-no+1:n2-no,:)
            ! Send right buffer to destination
            count=count+1
            call MPI_ISEND(buf1,icount,MPI_REAL_WP,dest,0,this%comm,request(count),ierr)
         end if
      end do
      ! Receive from potential sender
      count=count+1
      call MPI_IRECV(buf2,icount,MPI_REAL_WP,this%lvl(n)%recv_yp,0,this%comm,request(count),ierr)
      ! Wait for completion
      call MPI_WAITALL(count,request(1:count),MPI_STATUSES_IGNORE,ierr)
      ! Paste left buffer to right
      if (this%lvl(n)%recv_yp.ne.MPI_PROC_NULL) A(:,:,:,:,1:no,:)=buf2(:,:,:,:,:,:)
      ! Deallocate
      deallocate(buf1,buf2,request)
      ! Communicate in Z =======================================================
      ! Initialize buffer and size
      allocate(buf1(nstx,nsty,nstz,n1,n2,no),buf2(nstx,nsty,nstz,n1,n2,no)); icount=nstx*nsty*nstz*n1*n2*no
      ! Allocate request
      allocate(request(max(this%lvl(n)%nsend_zm,this%lvl(n)%nsend_zp)+1))
      ! MINUS DIRECTION ------
      ! Loop over potential destinations
      count=0
      do kp=1,this%npz
         ! Test if I need to send
         if (this%lvl(n)%send_zm(kp)) then
            ! Get processor rank
            dest=this%rank(this%iproc,this%jproc,kp)-1
            ! Copy left buffer
            buf1(:,:,:,:,:,:)=A(:,:,:,:,:,1+no:1+no+no-1)
            ! Send left buffer to destination
            count=count+1
            call MPI_ISEND(buf1,icount,MPI_REAL_WP,dest,0,this%comm,request(count),ierr)
         end if
      end do
      ! Receive from potential sender
      count=count+1
      call MPI_IRECV(buf2,icount,MPI_REAL_WP,this%lvl(n)%recv_zm,0,this%comm,request(count),ierr)
      ! Wait for completion
      call MPI_WAITALL(count,request(1:count),MPI_STATUSES_IGNORE,ierr)
      ! Paste left buffer to right
      if (this%lvl(n)%recv_zm.ne.MPI_PROC_NULL) A(:,:,:,:,:,n3-no+1:n3)=buf2(:,:,:,:,:,:)
      ! PLUS DIRECTION ------
      ! Loop over potential destinations
      count=0
      do kp=1,this%npz
         ! Test if I need to send
         if (this%lvl(n)%send_zp(kp)) then
            ! Get processor rank
            dest=this%rank(this%iproc,this%jproc,kp)-1
            ! Copy right buffer
            buf1(:,:,:,:,:,:)=A(:,:,:,:,:,n3-no-no+1:n3-no)
            ! Send right buffer to destination
            count=count+1
            call MPI_ISEND(buf1,icount,MPI_REAL_WP,dest,0,this%comm,request(count),ierr)
         end if
      end do
      ! Receive from potential sender
      count=count+1
      call MPI_IRECV(buf2,icount,MPI_REAL_WP,this%lvl(n)%recv_zp,0,this%comm,request(count),ierr)
      ! Wait for completion
      call MPI_WAITALL(count,request(1:count),MPI_STATUSES_IGNORE,ierr)
      ! Paste left buffer to right
      if (this%lvl(n)%recv_zp.ne.MPI_PROC_NULL) A(:,:,:,:,:,1:no)=buf2(:,:,:,:,:,:)
      ! Deallocate
      deallocate(buf1,buf2,request)
   end subroutine msync
   
   
   ! !> Log ILS info
   ! subroutine ils_log(this)
   !    use string,   only: str_long
   !    use messager, only: log
   !    implicit none
   !    class(ils), intent(in) :: this
   !    character(len=str_long) :: message
   !    if (this%cfg%amRoot) then
   !       write(message,'("Iterative Linear Solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name); call log(message)
   !       write(message,'(" >     method = ",i0)') this%method; call log(message)
   !       write(message,'(" >   it/maxit = ",i0,"/",i0)') this%it,this%maxit; call log(message)
   !       write(message,'(" >  aerr/acvg = ",es12.5,"/",es12.5)') this%aerr,this%acvg; call log(message)
   !       write(message,'(" >  rerr/rcvg = ",es12.5,"/",es12.5)') this%rerr,this%rcvg; call log(message)
   !    end if
   ! end subroutine ils_log
   !
   !
   ! !> Print ILS info to the screen
   ! subroutine ils_print(this)
   !    use, intrinsic :: iso_fortran_env, only: output_unit
   !    implicit none
   !    class(ils), intent(in) :: this
   !    if (this%cfg%amRoot) then
   !       write(output_unit,'("Iterative Linear Solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
   !       write(output_unit,'(" >     method = ",i0)') this%method
   !       write(output_unit,'(" >   it/maxit = ",i0,"/",i0)') this%it,this%maxit
   !       write(output_unit,'(" >  aerr/acvg = ",es12.5,"/",es12.5)') this%aerr,this%acvg
   !       write(output_unit,'(" >  rerr/rcvg = ",es12.5,"/",es12.5)') this%rerr,this%rcvg
   !    end if
   ! end subroutine ils_print
   !
   !
   ! !> Short print of ILS info to the screen
   ! subroutine ils_print_short(this)
   !    use, intrinsic :: iso_fortran_env, only: output_unit
   !    implicit none
   !    class(ils), intent(in) :: this
   !    if (this%cfg%amRoot) write(output_unit,'("Iterative Linear Solver [",a16,"] for config [",a16,"] -> it/maxit = ",i3,"/",i3," and rerr/rcvg = ",es12.5,"/",es12.5)') trim(this%name),trim(this%cfg%name),this%it,this%maxit,this%rerr,this%rcvg
   ! end subroutine ils_print_short
   !
   
end module bbmg_class