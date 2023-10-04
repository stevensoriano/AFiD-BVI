!***********************************************************************
      subroutine InitArrays
      use param
      use local_arrays
      use local_aux, only: vorx, vory, vorz
      use mpi_param, only: kstart,kend
      use mpih, only: lvlhalo
      use stat_arrays
      use AuxiliaryRoutines
      implicit none
      integer :: j,k,kc,i

      call AllocateReal3DArray(vx,1,n1,1,n2,kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(vy,1,n1,1,n2,kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(vz,1,n1,1,n2,kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(pr,1,n1,1,n2,kstart-lvlhalo,kend+lvlhalo)

      call AllocateReal3DArray(dph,1,n1,1,n2+1, &
          kstart-lvlhalo,kend+lvlhalo)

      call AllocateReal3DArray(vorx,1,n1,1,n2,  &
           kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(vory,1,n1,1,n2, &
           kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(vorz,1,n1,1,n2,  &
           kstart-lvlhalo,kend+lvlhalo)

      call AllocateReal3DArray(forcx,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(forcy,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(forcz,1,n1,1,n2,kstart,kend)

      call AllocateReal3DArray(rhs,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(dq,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(qcap,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(ru1,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(ru2,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(ru3,1,n1,1,n2,kstart,kend)

      call AllocateReal3DArray(vx_me,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(vy_me,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(vz_me,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(pr_me,1,n1,1,n2,kstart,kend)

      call AllocateReal3DArray(vx_rms,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(vy_rms,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(vz_rms,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(pr_rms,1,n1,1,n2,kstart,kend)

      call AllocateCplx2DArray(term1a,1,n1,1,7) 
      call AllocateCplx2DArray(term1b,1,n1,1,7) 
      call AllocateCplx2DArray(term2a,1,n2,1,7) 
      call AllocateCplx2DArray(term2b,1,n2,1,7) 
      call AllocateCplx2DArray(term3a,kstart,kend,1,7) 
      call AllocateCplx2DArray(term3b,kstart,kend,1,7) 

      return 
      end   
