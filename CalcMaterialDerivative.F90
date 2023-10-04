      subroutine CalcMaterialDerivative
      use param
      use local_arrays, only: vy,vz,vx,ru3
      use mpi_param, only: kstart,kend
      use local_aux
      implicit none
      integer :: kc,kp,jp,jm,jc,ic,im,ip,km
      real    :: h11,h12,h13,udx1,udx2,udx3
      real    :: h21,h22,h23
      real    :: h31,h32,h33

      udx1=dx1*0.25
      udx2=dx2*0.25
      udx3=dx3*0.25

!     =================================================================
!                 X - material derivative (from hdnl1.F)
!     =================================================================
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
!                 d  q_t q_t 
!                ------------
!                 d   t      
!
      h11=( (vx(ip,jc,kc)+vx(ic,jc,kc)) &
           *(vx(ip,jc,kc)+vx(ic,jc,kc)) &
           -(vx(im,jc,kc)+vx(ic,jc,kc)) &
           *(vx(im,jc,kc)+vx(ic,jc,kc)))*udx1

!     vx vy term
!
!
!                 d  q_t q_r 
!                ------------
!                 d   r      
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
!                 d  q_t q_x 
!                -----------
!                 d   x      
!
      h13=((vz(ic,jc,kp)+vz(im,jc,kp))*(vx(ic,jc,kp)+vx(ic,jc,kc)) &
          -(vz(ic,jc,kc)+vz(im,jc,kc))*(vx(ic,jc,kc)+vx(ic,jc,km)) &
          )*udx3
 
      matderx(ic,jc,kc)=h11+h12+h13 &
                       +(vx(ic,jc,kc)-vxo(ic,jc,kc))/dt
 
      enddo
      enddo
!$OMP  END PARALLEL DO
      enddo
!     =================================================================

!     =================================================================
!                  Y - material derivative (from hdnl2.F)
!     =================================================================

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

      matdery(ic,jc,kc)=h21+h22+h23 &
                       +(vy(ic,jc,kc)-vyo(ic,jc,kc))/dt

      enddo
      enddo
!$OMP  END PARALLEL DO
      enddo

!     =================================================================

!     =================================================================
!                  Z - material derivative (from hdnl3.F)
!     =================================================================

      do kc=kstart,kend
      km=kc-1
      kp=kc+1
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc,jm,jp,im,ip)
!$OMP& PRIVATE(h31,h32,h33,densit)
      do jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      do ic=1,n1m
      im=imv(ic)
      ip=ipv(ic)
!
!
!    vz vx term
!
!
!                d  q_x q_t 
!             -----------
!                d   t      
!
!
      h31=(((vx(ip,jc,kc)+vx(ip,jc,km)) &
           *(vz(ip,jc,kc)+vz(ic,jc,kc))) &
          -((vx(ic,jc,kc)+vx(ic,jc,km)) &
           *(vz(ic,jc,kc)+vz(im,jc,kc))))*udx1
!
!    vz vy term
!
!
!                d  q_x q_r 
!             -----------
!                d   r      
!
      h32=(((vy(ic,jp,kc)+vy(ic,jp,km)) &
           *(vz(ic,jp,kc)+vz(ic,jc,kc))) &
          -((vy(ic,jc,kc)+vy(ic,jc,km)) &
           *(vz(ic,jc,kc)+vz(ic,jm,kc))))*udx2
!
!    vz vz term
!
!
!                 d  q_x q_x 
!                -----------
!                 d   x      
!
      h33=(((vz(ic,jc,kp)+vz(ic,jc,kc)) &
           *(vz(ic,jc,kp)+vz(ic,jc,kc))) &
          -((vz(ic,jc,kc)+vz(ic,jc,km)) &
           *(vz(ic,jc,kc)+vz(ic,jc,km))))*udx3
 
      matderz(ic,jc,kc)=h31+h32+h33  &
                        +(vz(ic,jc,kc)-vzo(ic,jc,kc))/dt

      enddo
      enddo
!$OMP  END PARALLEL DO
      enddo

!     =================================================================

      return
      end
