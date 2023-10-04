      subroutine ExplicitTermsVZ
      use param
      use local_arrays, only: vy,vz,qcap,vx,forcz
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc
      integer :: km,kp,jmm,jpp,ic,im,ip
      real    :: h32,h33,h31
      real    :: udx1,udx2,udx3

      udx1=dx1*0.25
      udx2=dx2*0.25
      udx3=dx3*0.25

      do kc=kstart,kend
      km=kc-1
      kp=kc+1
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc,jmm,jpp,im,ip)
!$OMP& PRIVATE(h31,h32,h33,densit)
      do jc=1,n2m
      jmm=jmv(jc)
      jpp=jpv(jc)
      do ic=1,n1m
      im=imv(ic)
      ip=ipv(ic)
!
!
!    vz vx term
!
!
!                d  q_x q_z 
!             -----------
!                d   x      
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
!                d  q_y q_z 
!             -----------
!                d   y      
!
      h32=(((vy(ic,jpp,kc)+vy(ic,jpp,km)) &
           *(vz(ic,jpp,kc)+vz(ic,jc,kc))) &
          -((vy(ic,jc,kc)+vy(ic,jc,km)) &
           *(vz(ic,jc,kc)+vz(ic,jmm,kc))))*udx2
!
!    vz vz term
!
!
!                 d  q_z q_z 
!                -----------
!                 d   z      
!
      h33=(((vz(ic,jc,kp)+vz(ic,jc,kc)) &
           *(vz(ic,jc,kp)+vz(ic,jc,kc))) &
          -((vz(ic,jc,kc)+vz(ic,jc,km)) &
           *(vz(ic,jc,kc)+vz(ic,jc,km))))*udx3
 
 
      qcap(ic,jc,kc)=-(h31+h32+h33)+forcz(ic,jc,kc)

      enddo
      enddo
!$OMP  END PARALLEL DO
      enddo

      return
      end
