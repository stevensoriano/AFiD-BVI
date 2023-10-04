      subroutine CheckMaxVel
!EP   This routine calculates the maximum velocities 
      use param
      use local_arrays, only: vy,vz,vx
      use mpi_param, only: kstart,kend
      use mpih
      implicit none
      integer :: jc,kc,ic

       vmax=-100.d0

       do kc=kstart,kend
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc,ic)
!$OMP& REDUCTION(max: my_vmax1,my_vmax2,my_vmax3)
        do jc=1,n2m
          do ic=1,n1m
           vmax(1) = max(vmax(1),abs(vx(ic,jc,kc)))
           vmax(2) = max(vmax(2),abs(vy(ic,jc,kc)))
           vmax(3) = max(vmax(3),abs(vz(ic,jc,kc)))
         enddo
        enddo
!$OMP  END PARALLEL DO
       enddo

      call MpiAllMaxRealScalar(vmax(1))
      call MpiAllMaxRealScalar(vmax(2))
      call MpiAllMaxRealScalar(vmax(3))

      return   
      end     
