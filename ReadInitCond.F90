!***********************************************************************
      subroutine ReadInitCond
      use mpih
      use mpi_param, only: kstart,kend
      use local_arrays, only: vx,vy,vz,pr
      use param
      use AuxiliaryRoutines
      IMPLICIT NONE
      character*70 :: filcnw2
      integer :: ihist
      integer :: n1o,n2o,n3o,n3om
      integer :: j,i,kstarto,kendo
      real, allocatable, dimension(:,:,:) :: vyold,vzold,vxold
      integer, parameter :: ghosts = 50 ! min(ghosts) = 1
      integer  :: kstartog,kendog
      integer, allocatable, dimension(:) :: countko

      
!EP   Reading old grid information by rank0
      if (myid .eq. 0) then
      filcnw2 = 'continua_grid.dat'
      open(13,file=filcnw2,status='unknown')
      rewind(13)                                                      
      read(13,*) n1o,n2o,n3o,time
      close(13)
      endif
      
!EP   Bcasting old grid information and time

      call MpiBcastInt(n1o)
      call MpiBcastInt(n2o)
      call MpiBcastInt(n3o)
      call MpiBcastReal(time)
      
!EP   Check whether grid specifications have been updated
      if(n2o.ne.n2.or.n3o.ne.n3.or.n1o.ne.n1) then
      if(ismaster) write(*,*) "Interpolating new grid"
      if(n1.gt.n1o*2.or.n2.gt.n2o*2.or.n3.gt.n3o*2) then
      if(ismaster) write(*,*) "New grid resolution can not be more",&
      " than twice the old resolution"
       call MpiAbort
      endif

      n3om = n3o - 1
      
      call MpiBarrier
      call AllocateInt1DArray(countko,0,numtasks-1)

      call block(n3o-1, numtasks, myid, kstarto, kendo, countko) 

      kstartog = max(kstarto-ghosts,1)
      kendog   = min(kendo+ghosts,n3om)
      

!EP   vx
      allocate(vxold(1:n1o,1:n2o,kstartog-1:kendog+1))

      call mpi_read_continua(n1o,n2o,n3o,kstartog,kendog,1, &
       vxold(1:n1o,1:n2o,kstartog-1:kendog+1))

      if(myid.eq.numtasks-1) then
      do j=1,n2o
      do i=1,n1o
      vxold(i,j,n3o) = 0.0
      enddo
      enddo
      endif

      call interp(vxold,vx(1:n1,1:n2,kstart:kend),n1o,n2o,n3o, &
       1,kstartog,kendog)

      deallocate(vxold)

!EP   vy
      allocate(vyold(1:n1o,1:n2o,kstartog-1:kendog+1))

      call mpi_read_continua(n1o,n2o,n3o,kstartog,kendog,2, &
       vyold(1:n1o,1:n2o,kstartog-1:kendog+1))

      if(myid.eq.numtasks-1) then
      do j=1,n2o
      do i=1,n1o
      vyold(i,j,n3o) = 0.0
      enddo
      enddo
      endif

      call interp(vyold,vy(1:n1,1:n2,kstart:kend),n1o,n2o,n3o, &
       2,kstartog,kendog)

      deallocate(vyold)

!EP   vz
      allocate(vzold(1:n1o,1:n2o,kstartog-1:kendog+1))

      call mpi_read_continua(n1o,n2o,n3o,kstartog,kendog,3, &
       vzold(1:n1o,1:n2o,kstartog-1:kendog+1))

      if(myid.eq.numtasks-1) then
      do j=1,n2o
      do i=1,n1o
      vzold(i,j,n3o) = 0.0
      enddo
      enddo
      endif

      call interp(vzold,vz(1:n1,1:n2,kstart:kend),n1o,n2o,n3o, &
       3,kstartog,kendog)

      deallocate(vzold)

      else

!EP   One to one HDF read
      call mpi_read_continua(n1,n2,n3,kstart,kend,1,vx)
      call mpi_read_continua(n1,n2,n3,kstart,kend,2,vy)
      call mpi_read_continua(n1,n2,n3,kstart,kend,3,vz)
      call mpi_read_continua(n1,n2,n3,kstart,kend,4,pr)

      endif

      if(myid.eq.numtasks-1) then
       vx(1:n1,1:n2,n3) = 0.0
       vy(1:n1,1:n2,n3) = 0.0
       vz(1:n1,1:n2,n3) = 0.0
      endif

      if (ireset.eq.1) then                                             
       ihist=0                                                          
       time=0.
      endif                                                             

      return                                    
      end                                                               
