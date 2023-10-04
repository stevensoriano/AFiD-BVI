      subroutine DumpMovieFiles
      use param
      use mpih
      use mpi_param, only: kstart,kend
      use local_arrays, only: vy,vz,vx,pr
      use local_aux, only: vorx,vory,vorz
      use hdf5
      implicit none

      integer hdf_error
      integer :: i,j,k,ipp,jpp,kpp,kmm,jmm,imm
      integer :: ip,jp,kp,im,jm,km

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_vx
      integer(HID_T) :: dset_vy
      integer(HID_T) :: dset_vz
      integer(HID_T) :: dset_ww

      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info
      integer :: ndims
      integer :: mydata
      integer :: my_down, my_up,tag

      real, allocatable, dimension(:,:,:) :: vort3
      real :: tprfi
      real :: udx1,udx2,udx3,grad_ke
      real :: vortz, vortr, vortth
      integer :: itime
      character*30 filnam1,filnam2,filnam3
      character*30 filnam0,filnam4,filnam5
      character*5 ipfi
      real, allocatable, dimension(:,:,:) :: enstar
 
!RO   Form name

      tprfi = 100.
      itime=nint(time*tprfi)
      write(ipfi,82)itime
   82 format(i5.5)

!RO   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!RO   Form the name of the file

      filnam1 = 'flowmov/vx_'//ipfi//'.h5'
      filnam2 = 'flowmov/vy_'//ipfi//'.h5'
      filnam3 = 'flowmov/vz_'//ipfi//'.h5'
      filnam4 = 'flowmov/ww_'//ipfi//'.h5'
 
      call CalcVorticity

      allocate(enstar(1:n1m,1:n2m,kstart:kend))

      do k=kstart,kend
      do j=1,n2m
      do i=1,n1m
        enstar(i,j,k)=sqrt(vorx(i,j,k)**2.+vory(i,j,k)**2. & 
           +vorz(i,j,k)**2.)
      end do
      end do
      end do


!RO   Set offsets and element counts
   
      ndims = 3

      dims(1)=n1m
      dims(2)=n2m
      dims(3)=n3m

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)

      data_count(1) = n1m
      data_count(2) = n2m
      data_count(3) = kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = kstart-1


!EP   q1
!
!      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id,
!     %    hdf_error)
!
!      call h5pset_fapl_mpio_f(plist_id, comm, info,
!     &  hdf_error)
!
!      call h5fcreate_f(filnam1, H5F_ACC_TRUNC_F, file_id,
!     & hdf_error, access_prp=plist_id)
!
!      call h5pclose_f(plist_id, hdf_error)
!
!      call h5dcreate_f(file_id, 'vx', H5T_NATIVE_DOUBLE,
!     &                filespace,
!     &                dset_vx, hdf_error)
!
!      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 
!
!      call h5dget_space_f(dset_vx, slabspace, hdf_error)
!      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,
!     &                      data_offset, data_count, hdf_error)
!      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
!      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,
!     &                        hdf_error)
!       call h5dwrite_f(dset_vx, H5T_NATIVE_DOUBLE,
!     &   vx(1:n1m,1:n2m,kstart:kend), dims, 
!     &   hdf_error, file_space_id = slabspace, mem_space_id = memspace, 
!     &   xfer_prp = plist_id)
!      call h5pclose_f(plist_id, hdf_error)
!
!      call h5dclose_f(dset_vx, hdf_error)
!
!      call h5sclose_f(memspace, hdf_error)
!      call h5fclose_f(file_id, hdf_error)
!
!
!! =====================================
!
!      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id,
!     %    hdf_error)
!
!      call h5pset_fapl_mpio_f(plist_id, comm, info,
!     &  hdf_error)
!
!      call h5fcreate_f(filnam2, H5F_ACC_TRUNC_F, file_id,
!     & hdf_error, access_prp=plist_id)
!
!      call h5pclose_f(plist_id, hdf_error)
!
!      call h5dcreate_f(file_id, 'vy', H5T_NATIVE_DOUBLE,
!     &                filespace,
!     &                dset_vy, hdf_error)
!
!      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 
!
!      call h5dget_space_f(dset_vy, slabspace, hdf_error)
!      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,
!     &                      data_offset, data_count, hdf_error)
!      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
!      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,
!     &                        hdf_error)
!       call h5dwrite_f(dset_vy, H5T_NATIVE_DOUBLE,
!     &   vy(1:n1m,1:n2m,kstart:kend), dims, 
!     &   hdf_error, file_space_id = slabspace, mem_space_id = memspace, 
!     &   xfer_prp = plist_id)
!      call h5pclose_f(plist_id, hdf_error)
!
!      call h5dclose_f(dset_vy, hdf_error)
!
!      call h5sclose_f(memspace, hdf_error)
!      call h5fclose_f(file_id, hdf_error)
!
!
!! =====================================
!
!      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id,
!     %    hdf_error)
!
!      call h5pset_fapl_mpio_f(plist_id, comm, info,
!     &  hdf_error)
!
!      call h5fcreate_f(filnam3, H5F_ACC_TRUNC_F, file_id,
!     & hdf_error, access_prp=plist_id)
!
!      call h5pclose_f(plist_id, hdf_error)
!
!      call h5dcreate_f(file_id, 'vz', H5T_NATIVE_DOUBLE,
!     &                filespace,
!     &                dset_vz, hdf_error)
!
!      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 
!
!      call h5dget_space_f(dset_vz, slabspace, hdf_error)
!      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,
!     &                      data_offset, data_count, hdf_error)
!      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
!      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,
!     &                        hdf_error)
!       call h5dwrite_f(dset_vz, H5T_NATIVE_DOUBLE,
!     &   vz(1:n1m,1:n2m,kstart:kend), dims, 
!     &   hdf_error, file_space_id = slabspace, mem_space_id = memspace, 
!     &   xfer_prp = plist_id)
!      call h5pclose_f(plist_id, hdf_error)
!
!      call h5dclose_f(dset_vz, hdf_error)
!
!      call h5sclose_f(memspace, hdf_error)
!      call h5fclose_f(file_id, hdf_error)
!
!
! =====================================


      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)

      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

      call h5fcreate_f(filnam4, H5F_ACC_TRUNC_F, file_id, & 
        hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'ww', H5T_NATIVE_DOUBLE, &
                     filespace, dset_vx, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      call h5dget_space_f(dset_vx, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)  
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_vx, H5T_NATIVE_DOUBLE, &
         enstar(1:n1m,1:n2m,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset_vx, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)




      end subroutine 
      
