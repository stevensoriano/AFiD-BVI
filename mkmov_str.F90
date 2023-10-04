      subroutine mkmov_str
      use mpi_param
      use mls_param
      use mpih
      use hdf5
      use param

      IMPLICIT NONE

      integer ic,jc,kc

      integer i,j,k,inp
      integer ndims,itime

      real :: tprfi
      character*70 filnam1,filnamxdm
      character*5 ipfi


      tprfi = 1/tframe
      itime=nint(time*tprfi)
      write(ipfi,82)itime
   82 format(i5.5)

      filnam1='movie/struc'//ipfi//'.dat'

      open(111,file=filnam1)

      do inp=1,Nparticle
       do i=1,maxnv
        write(111,*)xyz(1:3,i,inp)
       end do
      end do

      do inp=1,Nparticle
       do i=1,maxnf
        write(111,*)vert_of_face(1:3,i,inp)
       end do
      end do

      close(111)

      return                                                          
      end                                                             
!     ---------------------------------------------------
      subroutine continua_str1
 
      use param
      use mls_param

      IMPLICIT NONE
 
      integer i,j,k,inp
      integer ndims,itime

      real :: tprfi
      character*70 filnam1,filnam2,filnamxdm
      character*5 ipfi,ipfip

      ipfi = 'resta'

      do inp=1,Nparticle

      write(ipfip,83)inp
   83 format(i5.5)

      filnam1='contstr/str1_'//ipfip//'_'//ipfi//'.dat'

      open(111,file=filnam1)

       do i=1,maxnv
        write(111,*)xyz(1:3,i,inp)
       end do
       do i=1,maxnv
        write(111,*)xyzv(1:3,i,inp)
       end do
       do i=1,maxnv
        write(111,*)xyza(1:3,i,inp)
       end do

      close(111)

      end do
      return
      end


