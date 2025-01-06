!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private

   !> Single config
   type(config), public :: cfg

   public :: geometry_init

contains

   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read, param_exists
      use parallel,    only: amRoot
      use messager,    only: die
      implicit none
      type(sgrid) :: grid

      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz
         integer :: cpd,cpl,np,nyr,nxr
         real(WP) :: rdx,dx,box_x1,box_x2,box_y1,ddrop
         real(WP), dimension(3) :: dctr
         real(WP) :: Lcalc,r,rold,err,tol
         real(WP) :: Lx,Ly,Lz
         real(WP), dimension(:), allocatable :: x,y,z

         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1))
         dx = Lx/nx
         call param_read('Ly',Ly); call param_read('ny',ny); allocate(y(ny+1))
         call param_read('nz',nz,default=1); allocate(z(nz+1))

         if (nz.eq.1) then
            Lz = dx
            Ly = dx
         else
           call param_read('Lz',Lz)
         end if

         ! Uniform grid
         do i=1,nx+1
            x(i) = real(i-1,WP)*dx
         end do
         do j=1,ny+1
            y(j) = real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do

         ! z is always uniform
         do k=1,nz+1
            z(k) = real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do

         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.true.,zper=.true.,name='ShockDrop')

       end block create_grid

      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition

         ! Read in partition
         call param_read('Partition',partition,short='p')

         ! Create partitioned grid
         cfg=config(grp=group,decomp=partition,grid=grid)

       end block create_cfg

      create_walls: block
         if (cfg%iproc.eq.1) cfg%VF(:cfg%imin-1,:,:)=0.0_WP
         if (cfg%iproc.eq.cfg%npx) cfg%VF(cfg%imax+1:,:,:)=0.0_WP   

       end block create_walls

   end subroutine geometry_init

end module geometry
