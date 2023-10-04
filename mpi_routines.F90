!===============================================
      subroutine block(n, p, irank, istart, iend, blcsz)
      implicit none
      integer,intent(in) :: n,p,irank
      integer,intent(out) :: istart,iend
      integer :: i
      integer,dimension(0:p-1),intent(out) :: blcsz
      
      do i=0,p-1
      blcsz(i) = floor(real((n+p-i-1)/p))
      enddo
      istart = sum(blcsz(0:irank))-blcsz(irank)+1
      iend = istart+blcsz(irank)-1

      end subroutine block

!=================================================           
      subroutine mpi_workdistribution
      use param
      use mpih 
      use mpi_param
      implicit none
      integer :: i
      
      if(.not. allocated(countj)) allocate(countj(0:numtasks-1))
      if(.not. allocated(countjp)) allocate(countjp(0:numtasks-1))
      if(.not. allocated(countk)) allocate(countk(0:numtasks-1))

!EP   For PERIODIC pressure solver
      call block(n2+1, numtasks, myid, jstartp, jendp, countjp)
      djp=jendp-jstartp+1

      call block(n2m, numtasks, myid, jstart, jend, countj)
      dj=jend-jstart+1
      
      call block(n3m, numtasks, myid, kstart, kend, countk)
      dk=kend-kstart+1

#ifdef DEBUG
      write(*,*) "jstart: ",jstart
      write(*,*) "jend: ",jend
      write(*,*) "jstartp: ",jstart
      write(*,*) "jendp: ",jend
      write(*,*) "kstart: ",kstart
      write(*,*) "kend: ",kend
#endif

      if( dj .lt. 1 ) then            
       write(6,*)'process ',myid,' has work load <1 cell in j direction'
       write(6,*)"Check grid dimensions and number of processes"
       
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif

      if( dk .lt. 1 ) then            
       write(6,*)'process ',myid,' has work load <1 cell in k direction'
       write(6,*)"Check grid dimensions and number of processes"
       
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif
  
      if(.not. allocated(offsetjp)) allocate(offsetjp(0:numtasks-1))
      if(.not. allocated(offsetj)) allocate(offsetj(0:numtasks-1))
      if(.not. allocated(offsetk)) allocate(offsetk(0:numtasks-1))
      
      offsetjp(:)=0
      offsetj(:)=0
      offsetk(:)=0
      do i=1,numtasks-1
        offsetjp(i)= offsetjp(i-1) + countjp(i-1)
        offsetj(i)= offsetj(i-1) + countj(i-1)
        offsetk(i)= offsetk(i-1) + countk(i-1)
      end do
      
      !-------For MPI-IO--------------------------------
      mydata= n2*dk*n1
      mydatam = n2m*dk*n1m

      if(myid .eq. numtasks-1) mydata = n2*(dk+1)*n1
      
      if(.not. allocated(countf)) allocate(countf(0:numtasks-1))
      if(.not. allocated(offsetf)) allocate(offsetf(0:numtasks-1))
       
      call MPI_ALLGATHER(mydata, 1, MPI_INTEGER, countf, 1, MPI_INTEGER, &
       MPI_COMM_WORLD,ierr)
    
      offsetf(:)=0
      do i=1,numtasks-1
        offsetf(i)= offsetf(i-1) + countf(i-1)
      end do
      
      !------------------------------------------------
      
      end subroutine mpi_workdistribution 
!===============================================
      subroutine update_both_ghosts(n1,n2,q,ks,ke)
      use mpih
      implicit none
      integer, intent(in) :: ks,ke
      real,intent(inout) :: q(n1,n2,ks-lvlhalo:ke+lvlhalo)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag
      
      mydata= n1*n2*lvlhalo
      
      my_down=myid-1
      
      my_up=myid+1

      if(myid .eq. 0) my_down=numtasks-1
      if(myid .eq. numtasks-1) my_up=0

      tag=1
      call MPI_ISEND(q(1,1,ke-lvlhalo+1), mydata, MDP, &
       my_up,tag,MPI_COMM_WORLD,req(1),ierr)
      
      call MPI_ISEND(q(1,1,ks), mydata,  MDP, &
       my_down,tag,MPI_COMM_WORLD,req(2), ierr)
     
      call MPI_IRECV(q(1,1,ks-lvlhalo), mydata,  MDP,  &
       my_down,tag,MPI_COMM_WORLD,req(3),ierr)
     
      call MPI_IRECV(q(1,1,ke+1), mydata,  MDP, &
       my_up, tag,MPI_COMM_WORLD,req(4),ierr)
     
      call MPI_Waitall(4,req,status,ierr)

      end subroutine update_both_ghosts
!=========================================
      subroutine update_upper_ghost(n1,n2,q)
      use mpih
      use mpi_param, only: kstart,kend
      implicit none
      real,intent(inout) :: q(n1,n2,kstart-lvlhalo:kend+lvlhalo)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag
       
      mydata= n1*n2*lvlhalo
      
      my_down= myid-1
      
      my_up= myid+1

      if(myid .eq. 0) my_down=numtasks-1
      if(myid .eq. numtasks-1) my_up= 0
     
      tag=1
      
      call MPI_ISEND(q(1,1,kstart), mydata, MDP, &
       my_down, tag, MPI_COMM_WORLD, req(1), ierr)
      
      call MPI_IRECV(q(1,1,kend+1), mydata, MDP, &
       my_up,tag, MPI_COMM_WORLD, req(2), ierr)
       
      call MPI_Waitall(2,req,status,ierr)

     
      end subroutine update_upper_ghost
!=========================================
      subroutine update_lower_ghost(n1,n2,q)
      use mpih
      use mpi_param, only: kstart,kend
      implicit none
      real,intent(inout) :: q(n1,n2,kstart-lvlhalo:kend+lvlhalo)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag
       
      mydata= n1*n2*lvlhalo
      
      my_down= myid-1
      
      my_up= myid+1

      if(myid .eq. 0) my_down= numtasks-1
      if(myid .eq. numtasks-1) my_up= 0
      
      tag=1
      
      call MPI_ISEND(q(1,1,kend-lvlhalo+1), mydata,  MDP, &
       my_up, tag, MPI_COMM_WORLD, req(1), ierr)
      
      call MPI_IRECV(q(1,1,kstart-lvlhalo), mydata,  MDP, &
       my_down,tag, MPI_COMM_WORLD, req(2), ierr)
       
      call MPI_Waitall(2,req,status,ierr)

      end subroutine update_lower_ghost
!=========================================
      subroutine update_add_upper_ghost(n1,n2,q)
      use mpih
      use mpi_param, only: kstart,kend
      implicit none
      real,intent(inout) :: q(n1,n2,kstart-lvlhalo:kend+lvlhalo-1)
      real :: buf(n1,n2,lvlhalo)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag
      integer :: ic,jc
       
      mydata= n1*n2*lvlhalo
      
      my_down= myid-1
      
      my_up= myid+1

      buf=0.0d0

      if(myid .eq. 0) my_down= numtasks -1
      if(myid .eq. numtasks-1) my_up= 0
     
      tag=1

      call MPI_ISEND(q(1,1,kstart-lvlhalo),mydata,MDP, &
       my_down, tag, MPI_COMM_WORLD, req(1), ierr)
      
      call MPI_IRECV(buf(1,1,1), mydata, MDP, &
       my_up,tag, MPI_COMM_WORLD, req(2), ierr)

      call MPI_Waitall(2,req,status,ierr)

      do ic=1,lvlhalo 
       jc=kend-lvlhalo+ic
       if(jc.eq.kend)then
       q(:,:,jc) = buf(:,:,ic)
       else
       q(:,:,jc) = q(:,:,jc) + buf(:,:,ic)
       end if
      end do

      end subroutine update_add_upper_ghost
!=========================================
      subroutine update_add_lower_ghost(n1,n2,q)
      use mpih
      use mpi_param, only: kstart,kend
      implicit none
      real,intent(inout) :: q(n1,n2,kstart-lvlhalo:kend+lvlhalo-1)
      real :: buf(n1,n2,lvlhalo)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag
      integer :: ic,jc
       
      mydata= n1*n2*(lvlhalo)
      
      my_down= myid-1
      
      my_up= myid+1

      buf=0.0d0

      if(myid .eq. 0) my_down= numtasks-1
      if(myid .eq. numtasks-1) my_up= 0
      
      tag=1
      
      call MPI_ISEND(q(1,1,kend),mydata,MDP, &
       my_up, tag, MPI_COMM_WORLD, req(1), ierr)

      call MPI_IRECV(buf(1,1,1), mydata, MDP, &
       my_down,tag, MPI_COMM_WORLD, req(2), ierr)

      call MPI_Waitall(2,req,status,ierr)

      do ic=1,lvlhalo
       jc=kstart+ic-2
       q(:,:,jc) = q(:,:,jc) + buf(:,:,ic)
      end do

      end subroutine update_add_lower_ghost
!==============================================
        subroutine mpi_globalsum_double_arr(var,nvar)
        use mpih
          implicit none
          real,intent(inout),dimension(nvar) :: var
          real,dimension(nvar) :: var2
          integer,intent(in) :: nvar

          call MPI_ALLREDUCE(var,var2,nvar,MPI_DOUBLE_PRECISION, &
              MPI_SUM,MPI_COMM_WORLD,ierr)          
          var = var2

        end subroutine mpi_globalsum_double_arr
!==============================================
        subroutine mpi_globalsum_double_forc(var)
        use mpih
          implicit none
          real,intent(inout),dimension(6,7,7,7) :: var
          real,dimension(6,7,7,7) :: var2
!          integer :: nvar=1715
          integer :: nvar=2058!6*7*7*7

          call MPI_ALLREDUCE(var,var2,nvar,MPI_DOUBLE_PRECISION, &
              MPI_SUM,MPI_COMM_WORLD,ierr)          
          var = var2

        end subroutine mpi_globalsum_double_forc

!==================================================      
      
      subroutine mem_dealloc
      use local_arrays
      use mpi_param
      use stat_arrays
      implicit none
      
      if(allocated(vx)) deallocate(vx)
      if(allocated(vy)) deallocate(vy)
      if(allocated(vz)) deallocate(vz)
      
      if(allocated(qcap)) deallocate(qcap)
      
      if(allocated(pr)) deallocate(pr)
      
      if(allocated(rhs)) deallocate(rhs)
      
      if(allocated(dph)) deallocate(dph)
      
      if(allocated(ru1)) deallocate(ru1)
      if(allocated(ru2)) deallocate(ru2)
      if(allocated(ru3)) deallocate(ru3)
      
      !---------------------------------------
      if(allocated(countj)) deallocate(countj)
      if(allocated(countk)) deallocate(countk)

      if(allocated(offsetj)) deallocate(offsetj)
      if(allocated(offsetk)) deallocate(offsetk)
      
      if(allocated(countf)) deallocate(countf)
      
      if(allocated(offsetf)) deallocate(offsetf)
      
      if(allocated(vx_me)) deallocate(vx_me)
      if(allocated(vy_me)) deallocate(vy_me)
      if(allocated(vz_me)) deallocate(vz_me)
      if(allocated(vx_rms)) deallocate(vx_rms)
      if(allocated(vy_rms)) deallocate(vy_rms)
      if(allocated(vz_rms)) deallocate(vz_rms)
      if(allocated(pr_rms)) deallocate(pr_rms)
      if(allocated(pr_me)) deallocate(pr_me)
    
      end subroutine mem_dealloc
!================================================
      subroutine mpi_write_continua
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

      integer(HID_T) :: dset_vx
      integer(HID_T) :: dset_vy
      integer(HID_T) :: dset_vz
      integer(HID_T) :: dset_pr
      integer(HID_T) :: dset_enst

      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer(HSIZE_T) :: dims_grid(1)
      integer(HID_T) :: dset_grid
      integer(HID_T) :: dspace_grid

      integer :: comm, info
      integer :: ndims

      character*30 filnam2,filnam3,filnam4
      character*30 filnamgrid,filnam5,filnam6

!RO   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!RO   Form the name of the file

      filnam2 = 'continua_vx.h5'
      filnam3 = 'continua_vy.h5'
      filnam4 = 'continua_vz.h5'
      filnam5 = 'continua_pr.h5'

!RO   Set offsets and element counts
   
      ndims = 3

      dims(1)=n1
      dims(2)=n2
      dims(3)=n3m

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)

      data_count(1) = n1
      data_count(2) = n2
      data_count(3) = kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = kstart-1

!RO   pressure

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, &
          hdf_error)

      call h5pset_fapl_mpio_f(plist_id, comm, info, &
        hdf_error)

      call h5fcreate_f(filnam5, H5F_ACC_TRUNC_F, file_id, &
       hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'pr', H5T_NATIVE_DOUBLE, &
                      filespace, dset_pr, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      call h5dget_space_f(dset_pr, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,&
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_pr, H5T_NATIVE_DOUBLE, &
         pr(1:n1,1:n2,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset_pr, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!EP   vx

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)

      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

      call h5fcreate_f(filnam2, H5F_ACC_TRUNC_F, file_id, &
       hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'Vx', H5T_NATIVE_DOUBLE, &
                      filespace, dset_vx, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      call h5dget_space_f(dset_vx, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_vx, H5T_NATIVE_DOUBLE, &
         vx(1:n1,1:n2,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset_vx, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!EP   vy

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)

      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

      call h5fcreate_f(filnam3, H5F_ACC_TRUNC_F, file_id, &
       hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'Vy', H5T_NATIVE_DOUBLE, &
                      filespace, dset_vy, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      call h5dget_space_f(dset_vy, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_vy, H5T_NATIVE_DOUBLE, &
         vy(1:n1,1:n2,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset_vy, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!EP   vz

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)

      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

      call h5fcreate_f(filnam4, H5F_ACC_TRUNC_F, file_id, &
       hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'Vz', H5T_NATIVE_DOUBLE, &
                      filespace, dset_vz, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      call h5dget_space_f(dset_vz, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dwrite_f(dset_vz, H5T_NATIVE_DOUBLE, &
         vz(1:n1,1:n2,kstart:kend), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset_vz, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)
 
      call CalcWriteQ

      if (myid .eq. 0) then
       open(13,file='continua_grid.dat',status='unknown')
       rewind(13)                                                      
       write(13,*) n1,n2,n3,time
       close(13)
      endif
      
!RO   Write the grid & statistics information
!RO   only if master process

      if (myid.eq.0) then

      ndims=1

      filnamgrid = 'continua_master.h5'
      call h5fcreate_f(filnamgrid,H5F_ACC_TRUNC_F, file_id, hdf_error)

!RO   Write Reynolds number

      dims_grid(1)=1
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

      call h5dcreate_f(file_id, 'Re', H5T_NATIVE_DOUBLE, &
                      dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, ren, &
             dims_grid,hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)
           

!EP   Write Prandtl number

      dims_grid(1)=1
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

      call h5dcreate_f(file_id, 'Pr', H5T_NATIVE_DOUBLE, &
                      dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, pra, &
             dims_grid,hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)
           

!RO   Write the grid information 

      dims_grid(1)=n1
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

      call h5dcreate_f(file_id, 'X_cordin', H5T_NATIVE_DOUBLE, &
                      dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, xc(1:n1), &
             dims_grid,hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=n2
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

      call h5dcreate_f(file_id, 'Y_cordin', H5T_NATIVE_DOUBLE, &
                      dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, yc(1:n2), &
             dims_grid,hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=n3m
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_id, 'Z_cordin', H5T_NATIVE_DOUBLE, &
                      dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zm(1:n3m), &
              dims_grid, hdf_error)


      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

!RO   Close file

      call h5fclose_f(file_id, hdf_error)

      endif
      end subroutine mpi_write_continua
      
!================================================
      
      subroutine mpi_read_continua(n1o,n2o,n3o,ks,ke,intvar,qua)
      use mpih
      use param
      use hdf5
      implicit none
      integer, intent(in) :: ks,ke,n2o,n1o,n3o
      real, dimension(1:n1o,1:n2o,ks-lvlhalo:ke+lvlhalo)::qua

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_qua

      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info
      integer :: ndims

      integer, intent(in) :: intvar
      character*70 :: filnam1
      character*10 :: dsetname

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!EP   Select file and dataset based on intvar

      select case (intvar)
        case (1)
          dsetname = trim('Vx')
          filnam1 = trim('continua_vx.h5')
        case (2)
          dsetname = trim('Vy')
          filnam1 = trim('continua_vy.h5')
        case (3)
          dsetname = trim('Vz')
          filnam1 = trim('continua_vz.h5')
        case (4)
          dsetname = trim('pr')
          filnam1 = trim('continua_pr.h5')
      end select

!RO   Set offsets and element counts
   
      ndims = 3

      dims(1)=n1o
      dims(2)=n2o
      dims(3)=n3o-1


      data_count(1) = n1o
      data_count(2) = n2o
      data_count(3) = ke-ks+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = ks-1


      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)

      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

      call h5fopen_f(filnam1, H5F_ACC_RDONLY_F, file_id, &
       hdf_error, access_prp=plist_id)

      call h5dopen_f(file_id, dsetname, dset_qua, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      call h5dget_space_f(dset_qua, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
                            data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
                              hdf_error)
       call h5dread_f(dset_qua, H5T_NATIVE_DOUBLE, &
         qua(1:n1o,1:n2o,ks:ke), dims,  &
         hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
         xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset_qua, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      end subroutine mpi_read_continua
!================================================

