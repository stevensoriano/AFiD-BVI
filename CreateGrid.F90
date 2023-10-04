      subroutine CreateGrid
      use param
      use mls_param,only:Hboxx,celvol
      use mpih
      use AuxiliaryRoutines
      implicit none
      integer :: ic,jc,kc
      real :: x1,x2,x3


!     Create indexing system

      call AllocateInt1DArray(imv,1,n1)
      call AllocateInt1DArray(ipv,1,n1)
      call AllocateInt1DArray(jmv,1,n2)
      call AllocateInt1DArray(jpv,1,n2)
      call AllocateInt1DArray(kmv,1,n3)
      call AllocateInt1DArray(kpv,1,n3)

      call AllocateInt1DArray(jmhv,1,n2+1)

!
!   indices for hor1 direction
!
      do ic=1,n1m
        imv(ic)=ic-1
        ipv(ic)=ic+1
        if(ic.eq.1) imv(ic)=n1m
        if(ic.eq.n1m) ipv(ic)=1
      enddo

!
!   indices for hor2 direction
!

!EP   Threaded for CPU affinity

!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc)
      do jc=1,n2m
        jmv(jc)=jc-1
        jpv(jc)=jc+1
        if(jc.eq.1) jmv(jc)=n2m
        if(jc.eq.n2m) jpv(jc)=1
      enddo
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc)
      do jc = 1,n2+1
       jmhv(jc) = mod(jc,n2m/2+1)
       if(jmhv(jc).eq.0) jmhv(jc) = n2m/2 + 1
      enddo
!$OMP  END PARALLEL DO

!
!     indices for hor3 direction
!
      do kc=1,n3m
        kmv(kc)=kc-1
        kpv(kc)=kc+1
        if(kc.eq.1) kmv(kc)=n3m
        if(kc.eq.n3m) kpv(kc)=1
      end do

!     Create coordinate system

      call AllocateReal1DArray(xc,1,n1)
      call AllocateReal1DArray(yc,1,n2)
      call AllocateReal1DArray(zc,1,n3)

      call AllocateReal1DArray(xm,1,n1)
      call AllocateReal1DArray(ym,1,n2)
      call AllocateReal1DArray(zm,1,n3)

      dx1=real(n1m)/xlen
      dx2=real(n2m)/ylen
      dx3=real(n3m)/zlen
      dx1q=dx1*dx1                                                      
      dx2q=dx2*dx2                                                      
      dx3q=dx3*dx3                                                      

!
!     UNIFORM GRID IN ALL DIRECTIONS
!
      do ic=1,n1
       x1=real(ic-1)/real(n1m)
       xc(ic)= xlen*x1
      end do
      do ic=1,n1m
       xm(ic)=(xc(ic)+xc(ic+1))*0.5d0
      end do

      do jc=1,n2
       x2=real(jc-1)/real(n2m)
       yc(jc)= ylen*x2
      end do
      do jc=1,n2m
       ym(jc)=(yc(jc)+yc(jc+1))*0.5d0
      end do

      do kc=1,n3
       x3=real(kc-1)/real(n3m)
       zc(kc)=zlen*x3
      enddo

      do kc=1,n3m
       zm(kc)=(zc(kc)+zc(kc+1))*0.5d0
      enddo

      Hboxx = (1.0/dx1 + 1.0/dx2 + 1.0/dx3)/3.0
      celvol = (1.0/dx1)*(1.0/dx2)*(1.0/dx3)

      return
      end
