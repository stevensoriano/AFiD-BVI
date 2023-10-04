      subroutine allocate_mls_local
      use param
      use mls_local
      use mpih
      use mpi_param, only: kstart, kend
      implicit none 
      integer :: merr
      !-------------------------------------------------
      allocate(for_xc(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo-1), &
                                                            stat=merr)
      
      if(merr .ne. 0) then
        write(6,*)"process  ",myid," failed to allocate memory for for"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      !-------------------------------------------------
      allocate(for_yc(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo-1), &
                                                            stat=merr)
      
      if(merr .ne. 0) then
        write(6,*)"process  ",myid," failed to allocate memory for forc"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      !-------------------------------------------------
      allocate(for_zc(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo-1), &
                                                            stat=merr)
      
      if(merr .ne. 0) then
        write(6,*)"process  ",myid," failed to allocate memory for forc"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      !-------------------------------------------------

      return
      end 
