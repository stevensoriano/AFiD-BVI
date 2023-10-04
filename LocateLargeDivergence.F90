!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: LocateLargeDivergence.F90                      !
!    CONTAINS: subroutine LocateLargeDivergence           !
!                                                         ! 
!    PURPOSE: Debugging routine. Output the location(s)   !
!     of excessive divergence.                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine LocateLargeDivergence
      use param
      use local_arrays, only: vy,vz,vx
      use mpih
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc,kp,jp,ic,ip
      real    :: dqcap
        
!      if(myid.eq.0) write(*,*) "I   J   K   MYID"
      do kc=kstart,kend
        kp=kc+1
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc,ic,jp,ip,dqcap)
        do jc=1,n2m
          jp=jpv(jc)
            do ic=1,n1m
              ip=ipv(ic)
              dqcap= (vx(ip,jc,kc)-vx(ic,jc,kc))*dx1 &
                    +(vy(ic,jp,kc)-vy(ic,jc,kc))*dx2 &
                    +(vz(ic,jc,kp)-vz(ic,jc,kc))*dx3
              if (abs(dqcap).gt.resid) then
!                 write(*,*) ic,jc,kc,myid
              endif
      enddo
      enddo
!$OMP  END PARALLEL DO
      enddo
      
      return     
      end         
