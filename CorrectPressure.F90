!
!   this subroutine performs the calculation of the corrected pressure.
!   this depends on the fractional step
!
      subroutine CorrectPressure
      use param
      use local_arrays, only: pr,dph
      use mpi_param, only: kstart,kend
      implicit none
      integer :: kp,km,jm,jp,jc,kc,ic,ip,im
      real    :: be
!
!    the pressure is evaluated at the center of the box.
!
!     p^{n+1} = p^{n} + phi^{n+1} - b * Nabla^2 phi^{n+1}
!
      be=al*beta
      do kc=kstart,kend
        kp=kc+1
        km=kc-1
!$OMP  PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(jc,jm,jp,ic,im,ip)
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
          do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)
              pr(ic,jc,kc)=pr(ic,jc,kc)+dph(ic,jc,kc)-be*( &
              (dph(ip,jc,kc)-2.0*dph(ic,jc,kc)+dph(im,jc,kc))*dx1q+ &
              (dph(ic,jp,kc)-2.0*dph(ic,jc,kc)+dph(ic,jm,kc))*dx2q+ &
              (dph(ic,jc,kp)-2.0*dph(ic,jc,kc)+dph(ic,jc,km))*dx3q)
      enddo
      enddo
!$OMP  END PARALLEL DO
      enddo
      return
      end
