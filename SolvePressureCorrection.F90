!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolvePressureCorrection.F90                    !
!    CONTAINS: subroutine SolvePressureCorrection         !
!                                                         ! 
!    PURPOSE: Compute the pressure correction by solving  !
!     a Poisson equation                                  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SolvePressureCorrection
      use param
      use local_arrays, only: dph
      use mpi_param
      use AuxiliaryRoutines
      implicit none
      integer :: j,k,i
      complex,allocatable,dimension(:) :: xb
      real,allocatable,dimension(:,:) :: xr
      complex,allocatable,dimension(:,:) :: xa
      complex,allocatable,dimension(:,:,:) :: buf
      complex,allocatable,dimension(:,:,:) :: buft

      call AllocateReal2DArray(xr,1,n1m,1,n2m)
      call AllocateCplx1DArray(xb,1,n3m)
      call AllocateCplx2DArray(xa,1,n1mh,1,n2m)
      call AllocateCplx3DArray(buf,1,n1,1,n2,kstart,kend)
      call AllocateCplx3DArray(buft,1,n3,1,n1,jstart,jend)

!     Fourier transform in all three directions

      do k=kstart,kend
       xr(1:n1m,1:n2m) = dph(1:n1m,1:n2m,k)
       call dfftw_execute_dft_r2c(fwd_plan,xr,xa)
       buf(1:n1mh,1:n2m,k) = xa(1:n1mh,1:n2m)
      end do

      call PackZ_UnpackR_C(buf(:,:,kstart:kend),buft(:,:,jstart:jend))

      do j=jstart,jend
       do i=1,n1m
        xb(1:n3m) = buft(1:n3m,i,j)
        call dfftw_execute_dft(fwdplan_1d,xb,xb)
        buft(1:n3m,i,j) = xb(1:n3m)
       end do
      end do

      buft=buft/float(n1m*n2m*n3m)

!     Solve the Poisson equation in frequency space
!     d^2/dx^2 = -k_x^2, etc...

      do j=jstart,jend
       do i=1,n1m
        do k=1,n3m
          buft(k,i,j) = buft(k,i,j)/(-ak3(k)-ak1(i)-ak2(j))
        end do
       end do
      end do
 
!     1,1,1 term would be divided by zero as all ks are zero
!     so manually set it to zero to avoid NaNs

      if(jstart.eq.1) buft(1,1,1) = 0.0d0

!     Inverse Fourier transform in all three directions

      do j=jstart,jend
       do i=1,n1m
        xb(1:n3m) = buft(1:n3m,i,j)
        call dfftw_execute_dft(bckplan_1d,xb,xb)
        buft(1:n3m,i,j) = xb(1:n3m)
       end do
      end do

      call PackR_UnpackZ_C(buft(:,:,jstart:jend),buf(:,:,kstart:kend))

       
      do k=kstart,kend
       xa(1:n1mh,1:n2m) = buf(1:n1mh,1:n2m,k)
       call dfftw_execute_dft_c2r(bck_plan,xa,xr)
       dph(1:n1m,1:n2m,k) = xr(1:n1m,1:n2m)
      end do


! =====================================

      call DestroyReal2DArray(xr)
      call DestroyCplx1DArray(xb)
      call DestroyCplx2DArray(xa)
      call DestroyCplx3DArray(buf)
      call DestroyCplx3DArray(buft)
 
      return
      end
