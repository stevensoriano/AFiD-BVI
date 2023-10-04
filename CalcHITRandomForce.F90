      subroutine CalcHITRandomForce
      use param
      use local_arrays, only: forcx,forcy,forcz
      use mpih
      use mpi_param
      use stat_arrays
      implicit none
      integer :: j,k,i,l,n,m
      integer :: kc,jc,ic
      real :: kl,kn,km
      real :: wlow,whigh
      real :: waven
      real :: scalpr, scalpc
      real :: u1, u2
      real, dimension(6,7,7,7) :: fcoefs
      complex :: fcoefsx,fcoefsy,fcoefsz
      complex :: dummy1,dummy3,dummy1b,dummy2b,dummy3b

      wlow=0.01
      whigh=kfmax*2*pi
      fcoefs=0.0d0

      if(myid.eq.0) then
      do n=1,7
       kn = float(n-4)*2*pi
       do m=1,7
        km = float(m-4)*2*pi
        do l=1,4
         kl = float(l-4)*2*pi
         waven=sqrt(kl**2+km**2+kn**2)
         if((waven.gt.wlow).and.(waven.lt.(whigh))) then
         do i=1,3
          call random_number(u1)
          call random_number(u2)
          bcoefs(2*i-1,l,m,n) = bcoefs(2*i-1,l,m,n)*(1-dt/tl) +  & 
     &     sqrt(2*epsstar*dt/(tl*tl))*sqrt(-2.0*log(u1))*cos(2*pi*u2)
          bcoefs(2*i,l,m,n) = bcoefs(2*i,l,m,n)*(1-dt/tl) +  & 
     &     sqrt(2*epsstar*dt/(tl*tl))*sqrt(-2.0*log(u1))*cos(2*pi*u2)
         end do
         scalpr=bcoefs(1,l,m,n)*kl + bcoefs(2,l,m,n)*km + &
     &          bcoefs(3,l,m,n)*kn
         scalpc=bcoefs(4,l,m,n)*kl + bcoefs(5,l,m,n)*km + &
     &          bcoefs(6,l,m,n)*kn
         fcoefs(1,l,m,n) =bcoefs(1,l,m,n)-scalpr*kl/(waven**2)
         fcoefs(2,l,m,n) =bcoefs(2,l,m,n)-scalpr*km/(waven**2)
         fcoefs(3,l,m,n) =bcoefs(3,l,m,n)-scalpr*kn/(waven**2)
         fcoefs(4,l,m,n) =bcoefs(4,l,m,n)-scalpc*kl/(waven**2)
         fcoefs(5,l,m,n) =bcoefs(5,l,m,n)-scalpc*km/(waven**2)
         fcoefs(6,l,m,n) =bcoefs(6,l,m,n)-scalpc*kn/(waven**2)
         end if
        end do
       end do
      end do

!RO     Complex conjugate
        do k=1,7
        kc = 8-k
        do j=1,7
        jc = 8-j
        do i=5,7
         ic = 8-i
         fcoefs(1,i,j,k) =fcoefs(1,ic,jc,kc)
         fcoefs(2,i,j,k) =fcoefs(2,ic,jc,kc)
         fcoefs(3,i,j,k) =fcoefs(3,ic,jc,kc)
         fcoefs(4,i,j,k)=-fcoefs(4,ic,jc,kc)
         fcoefs(5,i,j,k)=-fcoefs(5,ic,jc,kc)
         fcoefs(6,i,j,k)=-fcoefs(6,ic,jc,kc)
        end do
        end do
        end do

!RO     Special case i=4 (k=0)
        i=4 
        do k=1,7
         do j=1,7
         fcoefs(4,i,j,k)=0.
         fcoefs(5,i,j,k)=0.
         fcoefs(6,i,j,k)=0.
         end do
        end do


        do k=1,7
         kc=8-k
         do j=5,7
         jc=8-j
         fcoefs(1,i,j,k)=fcoefs(1,i,jc,kc)
         fcoefs(2,i,j,k)=fcoefs(2,i,jc,kc)
         fcoefs(3,i,j,k)=fcoefs(3,i,jc,kc)
         end do
        end do

        j=4
        do k=5,7
         kc=8-k
         fcoefs(1,i,j,k)=fcoefs(1,i,j,kc)
         fcoefs(2,i,j,k)=fcoefs(2,i,j,kc)
         fcoefs(3,i,j,k)=fcoefs(3,i,j,kc)
        end do

!RO    0 constant force

       fcoefs(:,4,4,4) = 0.0d0

       end if

!RS    Check the data range that actually needs to be send in this routine
       call mpi_globalsum_double_forc(fcoefs)

  
       forcx=0.0d0 ! Initial force as we are using these matrix elements to sum
       do n=2,6 ! The ring values (n=1,n=7) are zero, so exclude from summation
       do m=2,6 ! The ring values (m=1,m=7) are zero, so exclude from summation
       do l=2,6 ! The ring values (l=1,l=7) are zero, so exclude from summation
        fcoefsx=cmplx(fcoefs(1,l,m,n),fcoefs(4,l,m,n))
         do k=kstart,kend
         dummy1=term3a(k,n) ! use pre-allocated axial terms
         do j=1,n2m
          dummy1b=dummy1*term2a(j,m)!*term3a(k,n) ! use pre-allocated radial terms
          do i=1,n1m
           forcx(i,j,k) = forcx(i,j,k)+real(fcoefsx*dummy1b*term1b(i,l)) ! sum force
          end do 
         end do
         end do
        end do
        end do
        end do

!      Seperate loops per component, seem to benefit data locality, could be system dependent
       forcy=0.0d0 ! Initial force as we are using these matrix elements to sum
       do n=2,6 ! The ring values (n=1,n=7) are zero, so exclude from summation
       do m=2,6 ! The ring values (m=1,m=7) are zero, so exclude from summation
       do l=2,6 ! The ring values (l=1,l=7) are zero, so exclude from summation
        fcoefsy=cmplx(fcoefs(2,l,m,n),fcoefs(5,l,m,n))
         do k=kstart,kend
         dummy1=term3a(k,n) ! use pre-allocated axial terms
         do j=1,n2m 
         dummy2b=dummy1*term2b(j,m) ! use pre-allocated radial terms
          do i=1,n1m
           forcy(i,j,k)=forcy(i,j,k)+real(fcoefsy*dummy2b*term1a(i,l)) ! sum force
          end do
         end do
         end do
       end do
       end do
       end do

!      Seperate loops per component, seem to benefit data locality, could be system dependent
       forcz=0.0d0 ! Initial force as we are using these matrix elements to sum
       do n=2,6 ! The ring values (n=1,n=7) are zero, so exclude from summation
       do m=2,6 ! The ring values (m=1,m=7) are zero, so exclude from summation
       do l=2,6 ! The ring values (l=1,l=7) are zero, so exclude from summation
        fcoefsz=cmplx(fcoefs(3,l,m,n),fcoefs(6,l,m,n))
         do k=kstart,kend
         dummy3=term3b(k,n) ! use pre-allocated axial terms
         do j=1,n2m
          dummy3b=dummy3*term2a(j,m) ! use pre-allocated radial terms
          do i=1,n1m
           forcz(i,j,k)=forcz(i,j,k)+real(fcoefsz*dummy3b*term1a(i,l)) ! sum force
          end do
         end do
         end do
      end do
      end do
      end do

      return
      end
