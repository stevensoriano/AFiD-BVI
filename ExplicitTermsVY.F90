      subroutine ExplicitTermsVY
      use param
      use local_arrays, only: vx,vy,vz,dph,forcy
      use mpi_param, only: kstart,kend
      implicit none
      integer :: kc,kp,jp,jm,jc,ic,im,ip,km
      real    :: h22,h23,udx1,udx2,h21,udx3

      udx1=dx1*0.25
      udx2=dx2*0.25
      udx3=dx3*0.25
!
!     h term for the vy momentum equation at i+1/2,j,k+1/2
!
      do kc=kstart,kend
      km=kc-1
      kp=kc+1
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc,jm,jp,ic,im,ip)
!$OMP& PRIVATE(h21,h22,h23)
      do jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      do ic=1,n1m
      im=imv(ic)
      ip=ipv(ic)

!     vx vy term
!
!
!                 d  q_t q_r 
!                ------------
!                 d   t      
!
      h21=( (vy(ip,jc,kc)+vy(ic,jc,kc)) &
           *(vx(ip,jc,kc)+vx(ip,jm,kc)) &
           -(vy(ic,jc,kc)+vy(im,jc,kc)) &
           *(vx(ic,jc,kc)+vx(ic,jm,kc)) &
          )*udx1
      
!     vy vy term
!
!
!                 d  q_r q_r 
!                ------------
!                 d   r      
!
      h22=( (vy(ic,jp,kc)+vy(ic,jc,kc)) &
           *(vy(ic,jp,kc)+vy(ic,jc,kc)) &
           -(vy(ic,jm,kc)+vy(ic,jc,kc)) &
           *(vy(ic,jm,kc)+vy(ic,jc,kc)) &
          )*udx2
!
!     vy vz term
!
!
!                 d  q_x q_r 
!                -----------
!                 d   x      
!
      h23=((vz(ic,jc,kp)+vz(ic,jm,kp))*(vy(ic,jc,kp)+vy(ic,jc,kc)) &
          -(vz(ic,jc,kc)+vz(ic,jm,kc))*(vy(ic,jc,kc)+vy(ic,jc,km)) &
          )*udx3

!      dph(ic,jc,kc)=-(h21+h22+h23)
      dph(ic,jc,kc)=-(h21+h22+h23)+forcy(ic,jc,kc)

      enddo
      enddo
!$OMP  END PARALLEL DO
      enddo
      
      return
      end

