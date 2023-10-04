      subroutine ImplicitAndUpdateVY
      use param
      use local_arrays, only: vy,pr,rhs,dph,ru2
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,im,ip
      real    :: udx2
      real    :: dcvy,dpx22
      real    :: d22vy,d33vy,d11vy
      real    :: alre,udx1q,udx2q,udx3q

      alre=al/ren
      udx2=dx2*al
      udx1q=dx1q
      udx2q=dx2q
      udx3q=dx3q

!
!  compute the rhs of the factored equation
!  everything at i,j+1/2,k+1/2
!
        do kc=kstart,kend
          km=kc-1
          kp=kc+1
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc,jm,jp,ic,ip,im)
!$OMP& PRIVATE(d11vy,d22vy,d33vy,dcvy,dpx22)
          do jc=1,n2m
           jm=jmv(jc)
           jp=jpv(jc)
            do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)

!   viscid terms
!
!   x- second derivative of vy
!
            d11vy=(vy(ip,jc,kc)-2.0*vy(ic,jc,kc)+vy(im,jc,kc))*udx1q

!   y- second derivative of vy

            d22vy=(vy(ic,jp,kc)-2.0*vy(ic,jc,kc)+vy(ic,jm,kc))*udx2q

!   z- second derivative of vy

            d33vy=(vy(ic,jc,kp)-2.0*vy(ic,jc,kc)+vy(ic,jc,km))*udx3q

            dcvy=d11vy+d22vy+d33vy
 
!
!   component of grad(pr) along 2 direction
!
            dpx22=(pr(ic,jc,kc)-pr(ic,jm,kc))*udx2

            rhs(ic,jc,kc)=(ga*dph(ic,jc,kc)+ro*ru2(ic,jc,kc) &
                          +alre*dcvy-dpx22)*dt

            ru2(ic,jc,kc)=dph(ic,jc,kc)
         enddo
       enddo
!$OMP  END PARALLEL DO
      enddo


      call SolveTridX(beta*al*dx1q)
      call SolveTridY(beta*al*dx2q)
      call SolveTridZ(vy(1:n1,1:n2,kstart:kend),beta*al*dx3q)
      
      vy(:,n2,:) = vy(:,1,:)
     
      return
      end
