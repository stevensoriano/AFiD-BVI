!
!     Creates Initial condition (all zero for now)
!

      subroutine CreateInitialConditions
      use local_arrays, only: vy,vz,vx
      use param
      use mpi_param, only: kstart,kend
      implicit none

      integer :: i,j,k
      real :: aa, circ,xrc0,zrc0,rad,vphi,phi,xx,yy

      vx=0.0d0
      vy=0.0d0
      vz=0.0d0

      aa=vortexrad
      circ = 1.0
      xrc0 = xlen/3.0
      zrc0 = zlen/2.0

      do k=kstart,kend
      do j=1,n2m
       do i=1,n1m
         xx = xc(i)-xrc0
         yy = zm(k)-zrc0
         phi = atan2(yy,xx)
         rad = sqrt((xc(i)-xrc0)**2+(zm(k)-zrc0)**2)
         if(rad.gt.1.e-3) then
          vphi = circ/(pi*2.)*(1.-exp(-(rad**2/aa**2)))/rad
         else
          vphi = 0.
         end if
         vx(i,j,k) =-vphi*sin(phi)
        end do
       end do
      end do

      do k=kstart,kend
      do j=1,n2m
       do i=1,n1m
         xx = xm(i)-xrc0
         yy = zc(k)-zrc0
         phi = atan2(yy,xx)
         rad = sqrt((xm(i)-xrc0)**2+(zc(k)-zrc0)**2)
         if(rad.gt.1.e-3) then
          vphi = circ/(pi*2.)*(1.-exp(-(rad**2/aa**2)))/rad
         else
          vphi = 0.
         end if
         vz(i,j,k) = vphi*cos(phi)
        end do
       end do
      end do

      return                                                            
      end                                                               
