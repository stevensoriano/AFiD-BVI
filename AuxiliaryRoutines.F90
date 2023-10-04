!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: AuxiliaryRoutines.F90                          !
!    CONTAINS: subroutines Allocate*,Destroy*             !
!                                                         ! 
!    PURPOSE: Auxiliary routines used for memory allocs   !
!     and memory freeing                                  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      module AuxiliaryRoutines
      contains 

       subroutine AllocateReal1DArray(var,st1,en1)
       use param, only: ismaster
       implicit none
       integer, intent(in) :: st1,en1
       real,allocatable,dimension(:),intent(inout) :: var
       integer :: alloc_stat
  
       if (.not. allocated(var)) allocate(var(st1:en1), stat=alloc_stat)
      
       if (alloc_stat /= 0) then
        if(ismaster) then 
         write(6,*) 'Memory allocation failed when creating new arrays'
        end if
        call MpiAbort
       end if
  
       var =0.0d0
  
       return
       end subroutine AllocateReal1DArray
      
!===========================================================================

       subroutine AllocateCplx1DArray(var,st1,en1)
       use param, only: ismaster
       implicit none
       integer, intent(in) :: st1,en1
       complex,allocatable,dimension(:),intent(inout) :: var
       integer :: alloc_stat
  
       if (.not. allocated(var)) allocate(var(st1:en1), stat=alloc_stat)
      
       if (alloc_stat /= 0) then
        if(ismaster) then 
         write(6,*) 'Memory allocation failed when creating new arrays'
        end if
        call MpiAbort
       end if
  
       var =0.0d0
  
       return
       end subroutine AllocateCplx1DArray
      
!===========================================================================

       subroutine AllocateInt1DArray(var,st1,en1)
       use param, only: ismaster
       implicit none
       integer, intent(in) :: st1,en1
       integer,allocatable,dimension(:),intent(inout) :: var
       integer :: alloc_stat
  
       if (.not. allocated(var)) allocate(var(st1:en1), stat=alloc_stat)
      
       if (alloc_stat /= 0) then
        if(ismaster) then 
         write(6,*) 'Memory allocation failed when creating new arrays'
        end if
        call MpiAbort
       end if
  
       var =0
  
       return
       end subroutine AllocateInt1DArray
  
!===========================================================================

       subroutine AllocateReal2DArray(var,st1,en1,st2,en2)
       use param, only: ismaster
       implicit none
       integer, intent(in) :: st1,en1,st2,en2
       real,allocatable,dimension(:,:),intent(inout) :: var
       integer :: alloc_stat
  
       if (.not. allocated(var)) then 
         allocate(var(st1:en1,st2:en2), stat=alloc_stat)
       end if
      
       if (alloc_stat /= 0) then
        if(ismaster) then 
         write(6,*) 'Memory allocation failed when creating new arrays'
        end if
        call MpiAbort
       end if
  
       var =0.0
  
       return
       end subroutine AllocateReal2DArray

!===========================================================================

       subroutine AllocateCplx2DArray(var,st1,en1,st2,en2)
       use param, only: ismaster
       implicit none
       integer, intent(in) :: st1,en1,st2,en2
       complex,allocatable,dimension(:,:),intent(inout) :: var
       integer :: alloc_stat
  
       if (.not. allocated(var)) then 
         allocate(var(st1:en1,st2:en2), stat=alloc_stat)
       end if
      
       if (alloc_stat /= 0) then
        if(ismaster) then 
         write(6,*) 'Memory allocation failed when creating new arrays'
        end if
        call MpiAbort
       end if
  
       var =0.0
  
       return
       end subroutine AllocateCplx2DArray

!===========================================================================

       subroutine AllocateInt2DArray(var,st1,en1,st2,en2)
       use param, only: ismaster
       implicit none
       integer, intent(in) :: st1,en1,st2,en2
       integer,allocatable,dimension(:,:),intent(inout) :: var
       integer :: alloc_stat
  
       if (.not. allocated(var)) then 
         allocate(var(st1:en1,st2:en2), stat=alloc_stat)
       end if
      
       if (alloc_stat /= 0) then
        if(ismaster) then 
         write(6,*) 'Memory allocation failed when creating new arrays'
        end if
        call MpiAbort
       end if
  
       var =0
  
       return
       end subroutine AllocateInt2DArray

!===========================================================================

        subroutine AllocateReal3DArray(var,st1,en1,st2,en2,st3,en3)
        use param, only: ismaster
        implicit none
        integer, intent(in) :: st1,en1,st2,en2,st3,en3
        real,allocatable,dimension(:,:,:),intent(inout) :: var
        integer :: alloc_stat
  
        if (.not. allocated(var)) then 
          allocate(var(st1:en1,st2:en2,st3:en3), stat=alloc_stat)
        end if
      
        if (alloc_stat /= 0) then
        if(ismaster) then 
         write(6,*) 'Memory allocation failed when creating new arrays'
        end if
        call MpiAbort
        end if
  
        var =0.0
  
        return
        end subroutine AllocateReal3DArray

!===========================================================================

        subroutine AllocateCplx3DArray(var,st1,en1,st2,en2,st3,en3)
        use param, only: ismaster
        implicit none
        integer, intent(in) :: st1,en1,st2,en2,st3,en3
        complex,allocatable,dimension(:,:,:),intent(inout) :: var
        integer :: alloc_stat
  
        if (.not. allocated(var)) then 
          allocate(var(st1:en1,st2:en2,st3:en3), stat=alloc_stat)
        end if
      
        if (alloc_stat /= 0) then
        if(ismaster) then 
         write(6,*) 'Memory allocation failed when creating new arrays'
        end if
        call MpiAbort
        end if
  
        var =0.0
  
        return
        end subroutine AllocateCplx3DArray

!===========================================================================

        subroutine DestroyCplx1DArray(var)
        implicit none
        complex,allocatable,dimension(:),intent(inout) :: var

        if (allocated(var)) deallocate(var)      

        return
        end subroutine DestroyCplx1DArray
!===========================================================================

        subroutine DestroyReal1DArray(var)
        implicit none
        real,allocatable,dimension(:),intent(inout) :: var

        if (allocated(var)) deallocate(var)      

        return
        end subroutine DestroyReal1DArray

!===========================================================================

        subroutine DestroyInt1DArray(var)
        implicit none
        integer,allocatable,dimension(:),intent(inout) :: var

        if (allocated(var)) deallocate(var)      

        return
        end subroutine DestroyInt1DArray

!===========================================================================

        subroutine DestroyInt2DArray(var)
        implicit none
        integer,allocatable,dimension(:,:),intent(inout) :: var

        if (allocated(var)) deallocate(var)      

        return
        end subroutine DestroyInt2DArray

!===========================================================================

        subroutine DestroyReal2DArray(var)
        implicit none
        real,allocatable,dimension(:,:),intent(inout) :: var

        if (allocated(var)) deallocate(var)      

        return
        end subroutine DestroyReal2DArray

!===========================================================================

        subroutine DestroyCplx2DArray(var)
        implicit none
        complex,allocatable,dimension(:,:),intent(inout) :: var

        if (allocated(var)) deallocate(var)      

        return
        end subroutine DestroyCplx2DArray

!===========================================================================

        subroutine DestroyReal3DArray(var)
        implicit none
        real,allocatable,dimension(:,:,:),intent(inout) :: var

        if (allocated(var)) deallocate(var)      

        return
        end subroutine DestroyReal3DArray

!===========================================================================

        subroutine DestroyCplx3DArray(var)
        implicit none
        complex,allocatable,dimension(:,:,:),intent(inout) :: var

        if (allocated(var)) deallocate(var)      

        return
        end subroutine DestroyCplx3DArray
      
      end module AuxiliaryRoutines
