      subroutine ImplicitAndUpdateVZ
      use param
      use local_arrays, only: vz,qcap, pr,ru3,rhs
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,ip,im
      real    :: udx3
      real    :: dvz2,dvz3,dcvz,dpx33,dvz1
      real    :: alre,udx1q,udx2q,udx3q

      alre=al/ren
      udx1q=dx1q
      udx2q=dx2q
      udx3q=dx3q
      udx3 = al*dx3
!
!  compute the rhs of the factored equation
!  everything at i+1/2,j+1/2,k
!
      do kc=kstart,kend
        km=kc-1
        kp=kc+1
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc,jm,jp,ic,ip,im)
!$OMP& PRIVATE(dvz1,dvz2,dvz3,dcvz,dpx33)
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
            do ic=1,n1m
              im=imv(ic)
              ip=ipv(ic)

!   viscous terms
!
!   x- second derivatives of vz

            dvz1=(vz(im,jc,kc)-2.0*vz(ic,jc,kc)+vz(ip,jc,kc))*udx1q

!   y- second derivatives of vz

            dvz2=(vz(ic,jm,kc)-2.0*vz(ic,jc,kc)+vz(ic,jp,kc))*udx2q

!   z- second derivatives of vz

            dvz3=(vz(ic,jc,kp)-2.0*vz(ic,jc,kc)+vz(ic,jc,km))*udx3q
 
            dcvz=dvz2+dvz3+dvz1
!
!  component of grad(pr) along x3 direction
!
            dpx33=(pr(ic,jc,kc)-pr(ic,jc,km))*udx3

            rhs(ic,jc,kc)=(ga*qcap(ic,jc,kc)+ro*ru3(ic,jc,kc) &
                          +alre*dcvz-dpx33)*dt 

!  updating of the explicit terms

            ru3(ic,jc,kc)=qcap(ic,jc,kc)
         enddo
       enddo
!$OMP  END PARALLEL DO
      enddo

      call SolveTridX(beta*al*dx1q)
      call SolveTridY(beta*al*dx2q)
      call SolveTridZ(vz(1:n1,1:n2,kstart:kend),beta*al*dx3q)
 
      return
      end
