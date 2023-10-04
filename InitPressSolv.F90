!***********************************************************************
!  this subroutine perform the calculation of trigz for temperton fft
!
      subroutine InitFourierMetric
      use param
      implicit none
      integer :: n2mp,j,i,n1mp
      integer :: k, n3mh, n3mp


!      n1mh=n1m/2+1  Already defined in papero.F
      n1mp=n1mh+1
!      n2mh=n2m/2+1  Already ddefined in papero.F
      n2mp=n2mh+1
      n3mh=n3m/2+1
      n3mp=n3mh+1

      allocate(ao(1:n1))
      allocate(ap(1:n2))
      allocate(af(1:n3))

      allocate(ak1(1:n1))
      allocate(ak2(1:n2))
      allocate(ak3(1:n3))
!
!     wave number definition
!
      do i=1,n1mh
        ao(i)=(i-1)*2.d0*pi
      enddo
      do i=n1mp,n1m
        ao(i)=-(n1m-i+1)*2.d0*pi
      enddo
      do i=1,n1m
        ak1(i)=2.d0*(1.d0-dcos(ao(i)/n1m))*(float(n1m)/xlen)**2
      enddo

      do j=1,n2mh
        ap(j)=(j-1)*2.d0*pi
      enddo
      do j=n2mp,n2m
        ap(j)=-(n2m-j+1)*2.d0*pi
      enddo
      do j=1,n2m
        ak2(j)=2.d0*(1.d0-dcos(ap(j)/n2m))*(float(n2m)/ylen)**2
      enddo

      do k=1,n3mh
        af(k)=(k-1)*2.d0*pi
      enddo
      do k=n3mp,n3m
        af(k)=-(n3m-k+1)*2.d0*pi
      enddo
      do k=1,n3m
        ak3(k)=2.d0*(1.d0-dcos(af(k)/n3m))*(float(n3m)/zlen)**2
      enddo
 
      return
      end
!=====================================================
      subroutine InitPressSolv
      use param
      implicit none
      integer FFTW_EXHAUSTIVE
      parameter(FFTW_EXHAUSTIVE=64)
      real, allocatable, dimension(:,:) :: xr
      complex, allocatable, dimension(:,:) :: xa
      complex, allocatable, dimension(:) :: xb
#ifdef _OPENMP
      integer :: nt,fftw_info
#endif

      allocate(xr(n1m,n2m))
      allocate(xa(n1mh,n2m))
      allocate(xb(n3m))
    
      call InitFourierMetric

#ifdef _OPENMP
!$OMP  PARALLEL
!$OMP$ REDUCTION(+:nt)
      nt = nt + 1
!$OMP  END PARALLEL
      call dfftw_init_threads(fftw_info)
      if(fftw_info.eq.0) write(*,*) "ERROR: FFTW THREAD INIT FAIL"
      call dfftw_plan_with_nthreads(nt)
#endif

      call dfftw_plan_dft_r2c_2d(fwd_plan,n1m,n2m,xr,xa,FFTW_EXHAUSTIVE)
      call dfftw_plan_dft_c2r_2d(bck_plan,n1m,n2m,xa,xr,FFTW_EXHAUSTIVE)

      call dfftw_plan_dft_1d(fwdplan_1d,n3m,xb,xb,-1,64)
      call dfftw_plan_dft_1d(bckplan_1d,n3m,xb,xb,+1,64)

      if(allocated(xr)) deallocate(xr)
      if(allocated(xa)) deallocate(xa)
      if(allocated(xb)) deallocate(xb)

      return
      end
      
      
