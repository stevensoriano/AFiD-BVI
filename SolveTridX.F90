      subroutine SolveTridX(betadx)
!EP   Solves tridiagonal system in i direction
      use param
      use local_arrays, only : rhs
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc,ic
      real,intent(in) :: betadx
      real, allocatable, dimension(:) :: amil,apil,acil,fil
      real :: ackl_b

      allocate(amil(1:n1))
      allocate(apil(1:n1))
      allocate(acil(1:n1))
      allocate(fil(1:n1))

      do kc=kstart,kend
!$OMP  PARALLEL DO
!$OMP& DEFAULT(NONE)
!$OMP& SHARED(n2m,n1m,kc,betadx,rhs,forclo)
!$OMP& PRIVATE(jc,ic,fil,apil,acil,amil,ackl_b)
          do jc=1,n2m
             do ic=1,n1m
                ackl_b = 1.0/(1.0+2.0*betadx)
                apil(ic)=-betadx*ackl_b
                acil(ic)=1.0d0
                amil(ic)=-betadx*ackl_b
                fil(ic)=rhs(ic,jc,kc)*ackl_b
             enddo
                call PerTridiagSolve(amil,acil,apil,fil,1,n1m,n1)
             do ic=1,n1m
                rhs(ic,jc,kc) = fil(ic)  
             enddo
          end do
!$OMP  END PARALLEL DO
      end do 
 
      if(allocated(amil)) deallocate(amil)
      if(allocated(acil)) deallocate(apil)
      if(allocated(apil)) deallocate(acil)
      if(allocated(fil)) deallocate(fil)

      return
      end
