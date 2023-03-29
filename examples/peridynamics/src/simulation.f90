!> Various definitions and tools for running an NGA2 simulation
module simulation
   use simsolid_class, only: simsolid
   implicit none
   private
   
   !> Solid simulation
   type(simsolid) :: solid
   
   public :: simulation_init,simulation_run,simulation_final
   
   
contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      implicit none
      
      ! Initialize solid simulation
      call solid%init()
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      
      ! Solid drives overall time integration
      do while (.not.solid%time%done())
         
         ! Advance solid simulation
         call solid%step()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize solid simulation
      call solid%final()
      
   end subroutine simulation_final
   
   
end module simulation