      subroutine write_3dfield 
      use param
      use mpih
      use mpi_param, only: kstart,kend
      use local_arrays, only: vy,vz,vx,pr
      use hdf5
      implicit none

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_pr
      integer(HID_T) :: dset_vx
      integer(HID_T) :: dset_vy
      integer(HID_T) :: dset_vz

      integer(HID_T) :: plist_id

      integer(HSIZE_T), dimension(3) :: dims
      integer(HSIZE_T), dimension(3) :: data_count
      integer(HSSIZE_T), dimension(3) :: data_offset

      integer :: comm, info
      integer :: ndims

      character filnam1*30
      character filnam2*30
      character filnam3*30
      character filnam4*30
      character filnam5*30


!RO   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!RO   Form the name of the file

      write(filnam1,'(a,i6.6,a)')'field/field3d_pres_',int(ntime),'.h5'
      write(filnam2,'(a,i6.6,a)')'field/field3d_velx_',int(ntime),'.h5'
      write(filnam3,'(a,i6.6,a)')'field/field3d_vely_',int(ntime),'.h5'
      write(filnam4,'(a,i6.6,a)')'field/field3d_velz_',int(ntime),'.h5'

!RO   Set offsets and element counts

      ndims = 3

      dims(1)=n1
      dims(2)=n2
      dims(3)=n3m

      data_count(1) = n1
      data_count(2) = n2
      data_count(3) = kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = kstart-1

!EP   dens

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam1, H5F_ACC_TRUNC_F, file_id, &
                       hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'Pr', H5T_NATIVE_DOUBLE, &
                      filespace, dset_pr, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_pr, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dwrite_f(dset_pr, H5T_NATIVE_DOUBLE, &
                      pr(1:n1,1:n2,kstart:kend), dims, &
                      hdf_error, file_space_id = slabspace, &
                      mem_space_id = memspace, &
                      xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_pr, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!EP   vx

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam2, H5F_ACC_TRUNC_F, file_id, &
                       hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'Vx', H5T_NATIVE_DOUBLE, &
                      filespace, dset_vx, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_vx, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dwrite_f(dset_vx, H5T_NATIVE_DOUBLE, &
                      vx(1:n1,1:n2,kstart:kend), dims, &
                      hdf_error, file_space_id = slabspace, &
                      mem_space_id = memspace, &
                      xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_vx, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!EP   vy

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam3, H5F_ACC_TRUNC_F, file_id, &
                       hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'Vy', H5T_NATIVE_DOUBLE, &
                      filespace, dset_vy, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_vy, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dwrite_f(dset_vy, H5T_NATIVE_DOUBLE, &
                       vy(1:n1,1:n2,kstart:kend), dims, &
                       hdf_error, file_space_id = slabspace, &
                       mem_space_id = memspace, &
                       xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_vy, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!EP   vz

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam4, H5F_ACC_TRUNC_F, file_id, &
                       hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'Vz', H5T_NATIVE_DOUBLE, &
                      filespace, dset_vz, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_vz, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
      call h5dwrite_f(dset_vz, H5T_NATIVE_DOUBLE, &
                      vz(1:n1,1:n2,kstart:kend), dims, &
                      hdf_error, file_space_id = slabspace, &
                      mem_space_id = memspace, &
                      xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_vz, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)


      return
      end subroutine write_3dfield


