      subroutine ExplicitTermsVX
      use param
      use local_arrays, only: vx,vy,vz,dq,forcx
      use mpi_param, only: kstart,kend
      implicit none
      integer :: kc,kp,jp,jm,jc,ic,im,ip,km
      real    :: h11,h12,h13,udx1,udx2,udx3

 
      udx1=dx1*0.25
      udx2=dx2*0.25
      udx3=dx3*0.25

      do kc=kstart,kend
      km=kc-1
      kp=kc+1
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc,jm,jp,ic,im,ip)
!$OMP& PRIVATE(h11,h12,h13)
      do jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      do ic=1,n1m
      im=imv(ic)
      ip=ipv(ic)
      
!     vx vx term
!
!
!                 d  q_x q_x 
!                ------------
!                 d   x      
!
      h11=( (vx(ip,jc,kc)+vx(ic,jc,kc)) &
           *(vx(ip,jc,kc)+vx(ic,jc,kc)) &
           -(vx(im,jc,kc)+vx(ic,jc,kc)) &
           *(vx(im,jc,kc)+vx(ic,jc,kc)) &
          )*udx1

!     vx vy term
!
!
!                 d  q_x q_y 
!                ------------
!                 d   y      
!
      h12=( (vy(ic,jp,kc)+vy(im,jp,kc)) &
           *(vx(ic,jp,kc)+vx(ic,jc,kc)) &
           -(vy(ic,jc,kc)+vy(im,jc,kc)) &
           *(vx(ic,jc,kc)+vx(ic,jm,kc)) &
          )*udx2
!
!     vx vz term
!
!
!                 d  q_x q_z 
!                -----------
!                 d   z      
!
      h13=((vz(ic,jc,kp)+vz(im,jc,kp))*(vx(ic,jc,kp)+vx(ic,jc,kc)) &
          -(vz(ic,jc,kc)+vz(im,jc,kc))*(vx(ic,jc,kc)+vx(ic,jc,km)) &
          )*udx3
 
      dq(ic,jc,kc)=-(h11+h12+h13)+forcx(ic,jc,kc)
 
      enddo
      enddo
!$OMP  END PARALLEL DO
      enddo

      return
      end
