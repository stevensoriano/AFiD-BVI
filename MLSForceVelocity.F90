!------------------------------------------------------
!     force the velocity in single phase domain
!------------------------------------------------------
      subroutine MLSForceVelocity
      USE param
      USE mpih
      USE mpi_param, only: kstart, kend
      USE local_arrays, only: vx, vy, vz
      USE mls_local, only: for_xc, for_yc, for_zc
      IMPLICIT NONE


      integer ic,jc,kc
      integer im,jm,km

      call update_add_lower_ghost(n1,n2,for_xc)
      call update_add_lower_ghost(n1,n2,for_yc)
      call update_add_lower_ghost(n1,n2,for_zc)

      call update_add_upper_ghost(n1,n2,for_xc)
      call update_add_upper_ghost(n1,n2,for_yc)
      call update_add_upper_ghost(n1,n2,for_zc)
  
      call MpiBarrier

      do kc=kstart,kend
       km=kc-1
       do jc=1,n2m
       jm=jmv(jc)
        do ic=1,n1m
         im=imv(ic)
         vx(ic,jc,kc)=vx(ic,jc,kc)+(for_xc(ic,jc,kc)+for_xc(im,jc,kc)) &
          *0.5d0
         vy(ic,jc,kc)=vy(ic,jc,kc)+(for_yc(ic,jc,kc)+for_yc(ic,jm,kc)) &
          *0.5d0
         vz(ic,jc,kc)=vz(ic,jc,kc)+(for_zc(ic,jc,kc)+for_zc(ic,jc,km)) &
          *0.5d0
        end do
       end do
      end do
  


      return
      end  
