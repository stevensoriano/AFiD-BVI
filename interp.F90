      subroutine interp(arrold,arrnew,n1o,n2o,n3o, &
       intvar,kstarto,kendo)
! Trilinear interpolation to new grid for MPI
! continua.dat files.
! E.P. van der Poel 2012
      use param
      use mpih
      use mpi_param, only: kstart,kend
      implicit none
      integer,intent(in) :: intvar,n2o,n3o,n1o
      integer,intent(in) :: kstarto,kendo
 
      real,intent(out),dimension(1:n1,1:n2,kstart:kend) :: arrnew
      real,dimension(1:n1o,1:n2o,kstarto-1:kendo+1) :: arrold
      real,dimension(0:n1o+1) :: tcold,tmold
      real,dimension(1:n1) :: tcnew,tmnew
      real,dimension(0:n2o+1) :: rcold,rmold
      real,dimension(1:n2) :: rcnew,rmnew
      real,dimension(0:n3o+1) :: zzold,zmold
      real,dimension(1:n3) :: zznew,zmnew

      real, dimension(0:n1o+1) :: xold
      real, dimension(0:n2o+1) :: yold
      real, dimension(0:n3o+1) :: zold
      real, dimension(1:n1) :: xnew
      real, dimension(1:n2) :: ynew
      real, dimension(1:n3) :: znew

      real :: bn(6),an(8)

      real :: bix1,bix2,biy1,biy2,biz1,biz2,bix,biy,biz
      integer :: j,k,i,l
      real :: sfcf,sfff,sccf,scff,sfcc,sffc,sccc,scfc
      integer :: bici,bifi,bicj,bifj,bick,bifk

!EP   Create old grid
      call gridnew(n1o,xlen,0,1.0,tcold(1:n1o),tmold(1:n1o))
      call gridnew(n2o,ylen,0,1.0,rcold(1:n2o),rmold(1:n2o))
      call gridnew(n3o,zlen,0,1.0,zzold(1:n3o),zmold(1:n3o))

!EP   2nd order extrapolation of grid
      tcold(0) = 2*tcold(1)-tcold(2)
      tcold(n1o+1) = 2*tcold(n1o)-tcold(n1o-1)
      tmold(0) = 2*tmold(1)-tmold(2)
      tmold(n1o+1) = 2*tmold(n1o)-tmold(n1o-1)

      rcold(0) = 2*rcold(1)-rcold(2)
      rcold(n2o+1) = 2*rcold(n2o)-rcold(n2o-1)
      rmold(0) = 2*rmold(1)-rmold(2)
      rmold(n2o+1) = 2*rmold(n2o)-rmold(n2o-1)

      zzold(0) = 2*zzold(1)-zzold(2)
      zzold(n3o+1) = 2*zzold(n3o)-zzold(n3o-1)
      zmold(0) = 2*zmold(1)-zmold(2)
      zmold(n3o+1) = 2*zmold(n3o)-zmold(n3o-1)

!EP   Create new grid
      call gridnew(n1,xlen,0,1.0,tcnew,tmnew)
      call gridnew(n2,ylen,0,1.0,rcnew,rmnew)
      call gridnew(n3,zlen,0,1.0,zznew,zmnew)
      
      arrold(n1o,:,:) = arrold(1,:,:)
      arrold(:,n2o,:) = arrold(:,1,:)
      select case (intvar)
          case (1)
             xold=tcold
             yold=rmold
             zold=zmold
             xnew=tcnew
             ynew=rmnew
             znew=zmnew
          case (2)
             xold=tmold
             yold=rcold
             zold=zmold
             xnew=tmnew
             ynew=rcnew
             znew=zmnew
          case (3)
             xold=tmold
             yold=rmold
             zold=zzold
             xnew=tmnew
             ynew=rmnew
             znew=zznew
          case (4)
             xold=tmold
             yold=rmold
             zold=zmold
             xnew=tmnew
             ynew=rmnew
             znew=zmnew
      end select

!EP   Fill top and bottom slab


       if(kstarto.le.1) then
        do i=1,n1o
          do j=1,n2o
            arrold(i,j,0)=2.0*arrold(i,j,1)-arrold(i,j,2)
          enddo
        enddo
       endif

       if(kendo.ge.n3o) then
        do i=1,n1o
          do j=1,n2o
            arrold(i,j,n3o)=2.0*arrold(i,j,n3o-1)-arrold(i,j,n3o-2)
          enddo
        enddo
       endif

!EP   INTERP
      do k=kstart,kend
       do j=1,n2
        do i=1,n1
!    Find nearest grid value
      bix=xnew(i)
      biy=ynew(j)
      biz=znew(k)

      bifi=n1o-1
      bici=n1o
      do l=1,n1o
      if(xold(l).ge.bix) then
      bifi=l-1
      bici=l
      goto 10
      endif
      enddo
10    continue

      bifj=n2o-1
      bicj=n2o
      do l=1,n2o
      if(yold(l).ge.biy) then
      bifj=l-1
      bicj=l
      goto 20
      endif
      enddo
20    continue

      bifk=n3o-1
      bick=n3o
      do l=1,n3o
      if(zold(l).ge.biz) then
      bifk=l-1
      bick=l
      goto 30
      endif
      enddo
30    continue

      if(myid.eq.numtasks-1) then
      if(bifk.lt.kstarto-1) then
      write(*,*) "ERROR: INSUFFICIENT GHOSTS in inirea.F ", &
            "for myid: ",myid
      write(*,*) "INCREASE BY: ",kstarto-bifk
      call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
      endif
      elseif(myid.eq.0) then
      if(bick.gt.kendo+1) then
      write(*,*) "ERROR: INSUFFICIENT GHOSTS in inirea.F ", &
            "for myid: ",myid
      write(*,*) "INCREASE BY: ",bick-kendo
      call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
      endif
      else
      if(bifk.lt.kstarto-1.or.bick.gt.kendo+1) then
      write(*,*) "ERROR: INSUFFICIENT GHOSTS in inirea.F ", &
            "for myid: ",myid
      write(*,*) "INCREASE BY: ",max(kstarto-bifk,bick-kendo)
      call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
      endif
      endif

!     Define
      bix1 = xold(bifi) 
      bix2 = xold(bici)
      biy1 = yold(bifj) 
      biy2 = yold(bicj)
      biz1 = zold(bifk) 
      biz2 = zold(bick)

      if(bifi.eq.0) bifi = n1o-1
      if(bifj.eq.0) bifj = n2o-1

!EP   Send data
       sfcf = arrold(bifi,bicj,bifk)
       sfff = arrold(bifi,bifj,bifk)
       sccf = arrold(bici,bicj,bifk)
       scff = arrold(bici,bifj,bifk)

       sfcc = arrold(bifi,bicj,bick)
       sffc = arrold(bifi,bifj,bick)
       sccc = arrold(bici,bicj,bick)
       scfc = arrold(bici,bifj,bick)

          an(1) = sfcf
          an(2) = sfff
          an(3) = sccf
          an(4) = scff
          an(5) = sfcc
          an(6) = sffc
          an(7) = sccc
          an(8) = scfc
          bn(1) = bix1
          bn(2) = bix2
          bn(3) = biy1
          bn(4) = biy2
          bn(5) = biz1
          bn(6) = biz2
          call interptrilin(bix,biy,biz,an,bn,arrnew(i,j,k))

        enddo
       enddo
      enddo

      end

      subroutine interptrilin(bix,biy,biz,an,bn,ans)
      implicit none
      
      real, intent(in) :: bix,biy,biz
      real,intent(in) :: bn(6),an(8)
      real, intent(out) :: ans

      real :: bix1,bix2,biy1,biy2,biz1,biz2
      real :: afifjck,acifjck,aficjck,acicjck
      real :: afifjfk,acifjfk,aficjfk,acicjfk,dxdydz
      

      aficjfk = an(1)
      afifjfk = an(2)
      acicjfk = an(3)
      acifjfk = an(4)
      aficjck = an(5)
      afifjck = an(6)
      acicjck = an(7)
      acifjck = an(8)
      
      bix1 = bn(1)
      bix2 = bn(2)
      biy1 = bn(3)
      biy2 = bn(4)
      biz1 = bn(5)
      biz2 = bn(6)

      dxdydz = 1.0/((bix2-bix1)*(biy2-biy1)*(biz2-biz1))

      ans = (aficjfk*(bix2-bix)*(biy-biy1)*(biz2-biz) &
            +afifjfk*(bix2-bix)*(biy2-biy)*(biz2-biz) &
            +acicjfk*(bix-bix1)*(biy-biy1)*(biz2-biz) &
            +acifjfk*(bix-bix1)*(biy2-biy)*(biz2-biz) &
            +aficjck*(bix2-bix)*(biy-biy1)*(biz-biz1) &
            +afifjck*(bix2-bix)*(biy2-biy)*(biz-biz1) &
            +acicjck*(bix-bix1)*(biy-biy1)*(biz-biz1) &
            +acifjck*(bix-bix1)*(biy2-biy)*(biz-biz1) &
            )*dxdydz

      end

      subroutine gridnew(n,rext,istr,str,rc,rm)
      implicit none
      double precision,intent(in) :: rext,str
      double precision,dimension(1:n),intent(out) :: rc,rm
      double precision,dimension(1:n) :: etaz
      double precision,dimension(1:n+400) :: etazm
      double precision :: x3,etain,delet,pi
      integer,intent(in) :: n,istr
      integer :: i,n3mo,nclip

      if (istr.eq.0) then
        do i=1,n
          x3=real(i-1)/real(n-1)
          rc(i)=rext*x3
        enddo
       endif

      if(istr.eq.6) then
      pi=2.0d0*asin(1.0d0)
      nclip = int(str)
      n3mo = n+nclip+nclip
      do i=1,n3mo
        etazm(i)=+cos(pi*(float(i)-0.5)/float(n3mo))
      end do
      do i=1,n
        etaz(i)=etazm(i+nclip)
      end do
      delet = etaz(1)-etaz(n)
      etain = etaz(1)
      do i=1,n
        etaz(i)=etaz(i)/(0.5*delet)
      end do
      rc(1) = 0.0
      do i=2,n-1
        rc(i) = rext*(1.-etaz(i))*0.5
      end do
      rc(n) = rext
      endif


      do i=1,n-1
        rm(i)=(rc(i)+rc(i+1))*0.5d0
      enddo
      rm(n) = 2*rc(n)-rm(n-1)

      return
      end

