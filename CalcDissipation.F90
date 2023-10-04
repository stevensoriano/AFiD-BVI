!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcDissipation.F90                            !
!    CONTAINS: subroutine CalcDissipation                 !
!                                                         ! 
!    PURPOSE: Calculates the instantaneous kinetic        !
!     energy dissipation and writes it in dissip.out      !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcDissipation
      use mpih
      use param
      use local_arrays,only: vx,vy,vz
      use stat_arrays
      use mpi_param, only: kstart,kend

      implicit none
      integer :: i,j,k
      integer :: imm,ipp,jmm,jpp,kp,km
      real :: eta
      real :: h11,h12,h13,h21,h22,h23,h31,h32,h33
      real :: nute
      real :: udx1,udx2,dissipte,udx3

      
      nute = 0.0d0

      udx1=dx1
      udx2=dx2
      udx3=dx3
      
!================================================================
!      Dissipation rates
!================================================================
!
!                                   1  |         | 2
!                   dissipation:  ---- | nabla  u|
!                                  Re  |         |
!

       do k=kstart,kend
        kp=k+1
        km=k-1

!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(j,i,jmm,jpp,imm,ipp)
!$OMP& PRIVATE(h11,h12,h13)
!$OMP& PRIVATE(h21,h22,h23)
!$OMP& PRIVATE(h31,h32,h33)
!$OMP& PRIVATE(dissipte)
!$OMP& REDUCTION(+: nute)
       do j=1,n2m       
        jmm=jmv(j)
        jpp=jpv(j)

       do i=1,n1m
       imm= imv(i)
       ipp= ipv(i)

       h11=(vx(ipp,j,k)-vx(i,j,k))*udx1
       h12=(vx(i,jpp,k)-vx(i,j,k))*udx2
       h13=(vx(i,j,kp)-vx(i,j,k))*udx3

       h21=(vy(ipp,j,k)-vy(i,j,k))*udx1
       h22=(vy(i,jpp,k)-vy(i,j,k))*udx2
       h23=(vy(i,j,kp)-vy(i,j,k))*udx3

       h31=(vz(ipp,j,k)-vz(i,j,k))*udx1
       h32=(vz(i,jpp,k)-vz(i,j,k))*udx2
       h33=(vz(i,j,kp)-vz(i,j,k))*udx3

       dissipte = 2.0*(h11**2+h22**2+h33**2)+ &
               (h21+h12)**2+ (h31+h13)**2+ (h32+h23)**2


       nute = nute+dissipte
       end do
       end do
!$OMP  END PARALLEL DO
       end do


      call MpiAllSumRealScalar(nute)

      nute = nute/ren/float(n1m*n2m*n3m)
      eta = 1/(ren*ren*ren*nute)**(0.25)
      keta = float(n1m/2)*eta*2*pi
      
      if(ismaster) write(92,*) time,nute,keta

      return   
      end
