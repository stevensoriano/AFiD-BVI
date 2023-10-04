!----------------------------------------------------
!     Routine for writing plane cuts for movie
!---------------------------------------------------

      subroutine mkmov_flo

      use hdf5
      use param
      use local_arrays, only: vx,vy,vz

      IMPLICIT NONE
      real tprfi
      character*30 :: filename
      character*7 :: ipfi
      integer hdf_error, itime

      integer(HID_T) :: file_xp
      integer(HID_T) :: dset_xp
      integer(HID_T) :: dspace
      integer(HSIZE_T) :: dims(2)

      integer(HID_T) :: dspace_1
      integer(HSIZE_T) :: dims_1(1)
      
      tprfi = 1/tframe
      itime=int(time/tframe)
      write(ipfi,82)itime
   82 format(i5.5)

      filename= 'flowmov/zframe'//trim(ipfi)//'.h5'

!     set up mpi file props

      call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_xp,hdf_error)

      dims(1)=n1m
      dims(2)=n2m
     
      call h5screate_simple_f(2,dims,dspace,hdf_error)
      call h5dcreate_f(file_xp,'vx',H5T_NATIVE_DOUBLE, &
                      dspace,dset_xp,hdf_error)
      call h5dwrite_f(dset_xp,H5T_NATIVE_DOUBLE, &
                        vx(1:n1m,1:n2m,n3m/2.0),&
                        dims,hdf_error)
      call h5dclose_f(dset_xp,hdf_error)
      call h5sclose_f(dspace,hdf_error)

      call h5screate_simple_f(2,dims,dspace,hdf_error)
      call h5dcreate_f(file_xp,'vy',H5T_NATIVE_DOUBLE, &
                      dspace,dset_xp,hdf_error)
      call h5dwrite_f(dset_xp,H5T_NATIVE_DOUBLE, &
                        vy(1:n1m,1:n2m,n3m/2.0), &
                        dims,hdf_error)
      call h5dclose_f(dset_xp,hdf_error)
      call h5sclose_f(dspace,hdf_error)

      call h5screate_simple_f(2,dims,dspace,hdf_error)
      call h5dcreate_f(file_xp,'vz',H5T_NATIVE_DOUBLE, &
                      dspace,dset_xp,hdf_error)
      call h5dwrite_f(dset_xp,H5T_NATIVE_DOUBLE, &
                        vz(1:n1m,1:n2m,n3m/2.0), &
                        dims,hdf_error)
      call h5dclose_f(dset_xp,hdf_error)
      call h5sclose_f(dspace,hdf_error)



      call h5fclose_f(file_xp, hdf_error)

      return
      end
