!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CheckDivergence.F90                            !
!    CONTAINS: subroutine CheckDivergence                 !
!                                                         ! 
!    PURPOSE: Check the maximum divergence of velocity    !
!     in the domain                                       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CheckDivergence(qmax)
      use param
      use local_arrays, only: vy,vz,vx
      use mpih
      use mpi_param, only: kstart,kend
      implicit none
      real,intent(out) :: qmax
      integer :: jc,kc,kp,jp,ic,ip
      real    :: dqcap,dvol,my_qmax
        

      dvol=1.d0/(dx1*dx2*dx3)
      qmax=0.d0                                                     

      do kc=kstart,kend
        kp=kc+1
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc,ic,jp,ip,dqcap)
!$OMP& REDUCTION(max:qmax)
        do jc=1,n2m
          jp=jpv(jc)
            do ic=1,n1m
              ip=ipv(ic)
              dqcap= (vx(ip,jc,kc)-vx(ic,jc,kc))*dx1 &
                    +(vy(ic,jp,kc)-vy(ic,jc,kc))*dx2 &
                    +(vz(ic,jc,kp)-vz(ic,jc,kc))*dx3
              qmax = max(abs(dqcap),qmax)          
      enddo
      enddo
!$OMP  END PARALLEL DO
      enddo
     
      call MpiAllMaxRealScalar(qmax)
      
      return     
      end         
