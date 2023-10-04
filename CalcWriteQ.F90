      subroutine CalcWriteQ
      use param
      use mpih
      use mpi_param, only: kstart,kend
      use local_arrays, only: vy,vz,vx
      use local_aux, only: vorx, vory, vorz
      use hdf5
      implicit none

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_q
      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count
      integer(HSSIZE_T), dimension(3) :: data_offset

      integer :: comm, info
      integer :: ndims

      real, allocatable, dimension(:,:,:) :: qtens
      real, allocatable, dimension(:,:,:) :: enstar
      real :: dvxx1,dvxx2,dvxx3
      real :: dvyx1,dvyx2,dvyx3
      real :: dvzx1,dvzx2,dvzx3
      real :: strn, omeg

      integer :: ic,jc,kc,ip,jp,kp,im,jm,km


      character*30 filnam

      call update_both_ghosts(n1,n2,vx,kstart,kend)
      call update_both_ghosts(n1,n2,vy,kstart,kend)
      call update_both_ghosts(n1,n2,vz,kstart,kend)

      call CalcVorticity

      allocate(qtens(1:n1m,1:n2m,kstart:kend))
      allocate(enstar(1:n1m,1:n2m,kstart:kend))

      do kc=kstart,kend
       kp=kc+1
       km=kc-1
       do jc=1,n2m
        jp=jpv(jc)
        jm=jmv(jc)
         do ic=1,n1m
          ip=ipv(ic)
          im=imv(ic)

          dvxx1=vx(ip,jc,kc)-vx(ic,jc,kc)

          dvxx2=(vx(ic,jp,kc)+vx(ip,jp,kc))- &
                (vx(ic,jm,kc)+vx(ip,jm,kc))*0.25d0

          dvxx3=(vx(ic,jc,kp)+vx(ip,jc,kp))- &
                (vx(ic,jc,km)+vx(ip,jc,km))*0.25d0


          dvyx1=(vy(ip,jc,kc)+vy(ip,jp,kc))- &
                (vy(im,jc,kc)+vy(im,jp,kc))*0.25d0

          dvyx2=vy(ic,jp,kc)-vy(ic,jc,kc)

          dvyx3=(vy(ic,jc,kp)+vy(ic,jp,kp))- &
                (vy(ic,jc,km)+vy(ic,jp,km))*0.25d0

          dvzx1=(vz(ip,jc,kc)+vz(ip,jc,kp))- &
                (vz(im,jc,kc)+vz(im,jc,kp))*0.25d0


          dvzx2=(vz(ic,jp,kc)+vz(ic,jp,kp))- &
                (vz(ic,jm,kc)+vz(ic,jm,kp))*0.25d0

          dvzx3=vz(ic,jc,kp)-vz(ic,jc,kc)

          strn=dvxx1**2.0+dvyx2**2.0+dvzx3**2.0 + &
           0.25*((dvyx1+dvxx2)**2.0+(dvzx1+dvxx3)**2.0 &
             +(dvyx3+dvzx2)**2.0)

          omeg= &
           0.25*((dvyx1-dvxx2)**2.0+(dvzx1-dvxx3)**2.0 &
             +(dvyx3-dvzx2)**2.0)

          qtens(ic,jc,kc)=0.5*(omeg-strn)*dx1*dx1

          enstar(ic,jc,kc)=sqrt(vorx(ic,jc,kc)**2.+vory(ic,jc,kc)**2. & 
               +vorz(ic,jc,kc)**2.)

        end do
       end do
      end do


      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!RO   Form the name of the file

      filnam = 'continua_qvort.h5'


!RO   Set offsets and element counts

      ndims = 3

      dims(1)=n1m
      dims(2)=n2m
      dims(3)=n3m

      call h5screate_simple_f(ndims, dims, &
                              filespace, hdf_error)

      data_count(1) = n1m
      data_count(2) = n2m
      data_count(3) = kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = kstart-1

!RO   Write out Q

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)

      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

      call h5fcreate_f(filnam, H5F_ACC_TRUNC_F, file_id, &
       hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'qtens', H5T_NATIVE_DOUBLE, &
                      filespace, dset_q, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)

      call h5dget_space_f(dset_q, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)

       call h5dwrite_f(dset_q, H5T_NATIVE_DOUBLE, &
         qtens(1:n1m,1:n2m,kstart:kend), dims, &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,&
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset_q, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!RO   Form the name of the file

      filnam = 'continua_enst.h5'

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)

      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

      call h5fcreate_f(filnam, H5F_ACC_TRUNC_F, file_id, &
       hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'enst', H5T_NATIVE_DOUBLE, &
                      filespace, dset_q, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)

      call h5dget_space_f(dset_q, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)

       call h5dwrite_f(dset_q, H5T_NATIVE_DOUBLE, &
         enstar(1:n1m,1:n2m,kstart:kend), dims, &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,&
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset_q, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      if(allocated(qtens)) deallocate(qtens)
      if(allocated(enstar)) deallocate(enstar)


      end subroutine CalcWriteQ

