      subroutine InitRandomForce
      use mpih
      use mpi_param, only: kstart,kend
      use param
      use hdf5
      IMPLICIT NONE
      integer i,j,k,l,m,n
      integer hdf_error, ndims
      character*30 filnambc
      integer(HID_T) :: file_id
      integer(HID_T) :: dset_coef
      integer(HID_T) :: dspace_coef
      integer(HSIZE_T) :: dims(4)
      real kl, km, kn
      complex :: im=(0.,1.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Pre-allocation of Geometrical terms used in the forcing (calcforc.F90)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!! Azimuthal terms
      do i=1,n1m
       do l=1,7
       kl = float(l-4)*2*pi
       term1a(i,l)=exp(im*kl*xm(i))
       term1b(i,l)=exp(im*kl*xc(i))
       enddo
      enddo

!! Radial terms
      do j=1,n2m
       do m=1,7
       km = float(m-4)*2*pi
       term2a(j,m)=exp(im*km*ym(j))
       term2b(j,m)=exp(im*km*yc(j))
       enddo
      enddo

!! Axial terms
      do k=kstart,kend
       do n=1,7
       kn = float(n-4)*2*pi
       term3a(k,n)=exp(im*kn*zm(k))
       term3b(k,n)=exp(im*kn*zc(k))
       end do
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! End pre-allocation of Geometrical terms used in the forcing (calcforc.F90)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      filnambc = 'continua_bcoefs.h5'

      if(ismaster.and.nread.eq.1)then
       ndims = 4

       dims(1)=6
       dims(2)=7
       dims(3)=7
       dims(4)=7

      call h5fopen_f(filnambc, H5F_ACC_RDONLY_F, file_id, hdf_error)

      call h5screate_simple_f(ndims, dims, dspace_coef, hdf_error)

      call h5dopen_f(file_id, 'bcoefs', dset_coef, hdf_error)

      call h5dread_f(dset_coef, H5T_NATIVE_DOUBLE, bcoefs,  &
             dims,hdf_error)

      call h5dclose_f(dset_coef, hdf_error)
      call h5sclose_f(dspace_coef, hdf_error)

      call h5fclose_f(file_id, hdf_error)


      else 

       bcoefs=0.0d0
 
      endif

      return                                    
      end                                                               
