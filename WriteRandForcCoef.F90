      subroutine WriteRandForcCoef
      use mpih
      use param
      use hdf5
      IMPLICIT NONE
      integer hdf_error, ndims
      character*30 filnambc
      integer(HID_T) :: file_id
      integer(HID_T) :: dset_coef
      integer(HID_T) :: dspace_coef
      integer(HSIZE_T) :: dims_coef(4)


      filnambc = 'continua_bcoefs.h5'

      if(ismaster) then 

      ndims=4
      dims_coef(1)=6
      dims_coef(2)=7
      dims_coef(3)=7
      dims_coef(4)=7

      call h5fcreate_f(filnambc,H5F_ACC_TRUNC_F, file_id, hdf_error)

      call h5screate_simple_f(ndims, dims_coef, dspace_coef, hdf_error)

      call h5dcreate_f(file_id, 'bcoefs', H5T_NATIVE_DOUBLE, &
                      dspace_coef, dset_coef, hdf_error)

      call h5dwrite_f(dset_coef, H5T_NATIVE_DOUBLE, bcoefs, &
             dims_coef,hdf_error)

      call h5dclose_f(dset_coef, hdf_error)
      call h5sclose_f(dspace_coef, hdf_error)

      call h5fclose_f(file_id, hdf_error)

      endif


      return                                    
      end                                                               
