      subroutine SolveTridZ(q,betadx)
      use param
      use local_arrays, only : rhs
      use mpi_param
      use mpih
      implicit none
      real, intent(inout) :: q(1:n1,1:n2,kstart:kend)
      real, allocatable, dimension(:) :: amkl,apkl,ackl,fkl
      integer :: jc,kc,ic
      real :: betadx,ackl_b
      real,allocatable :: rhst(:,:,:)

      allocate(rhst(1:n3,1:n1,jstart:jend))

      allocate(amkl(1:n3))
      allocate(apkl(1:n3))
      allocate(ackl(1:n3))
      allocate(fkl(1:n3))

!$OMP  PARALLEL
!$OMP  SINGLE
      call PackZ_UnpackR(rhs(:,:,kstart:kend),rhst(:,:,jstart:jend))
!$OMP  END SINGLE NOWAIT
!$OMP  END PARALLEL 

      do jc=jstart,jend
!$OMP  PARALLEL DO
!$OMP& DEFAULT(NONE)
!$OMP& SHARED(n1m,n3m,betadx,jc)
!$OMP& SHARED(rhst)
!$OMP& PRIVATE(ic,kc,ackl_b,fkl,info,appk,ipkv)
!$OMP& PRIVATE(amkl,ackl,apkl)
         do ic=1,n1m
          do kc=1,n3m
            ackl_b=1.0d0/(1.+2.0*betadx)
            amkl(kc)=-betadx*ackl_b
            ackl(kc)=1.0d0
            apkl(kc)=-betadx*ackl_b
            fkl(kc)=rhst(kc,ic,jc)*ackl_b
          end do

          call PerTridiagSolve(amkl,ackl,apkl,fkl,1,n3m,n3)

          do kc=1,n3m
            rhst(kc,ic,jc)=fkl(kc)
          end do
         enddo
!$OMP  END PARALLEL DO
      end do

      call PackR_UnpackZ(rhst(:,:,jstart:jend),rhs(:,:,kstart:kend))

!     Update velocities

      do kc=kstart,kend
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(ic,jc)
      do jc=1,n2m
      do ic=1,n1m
      q(ic,jc,kc) = q(ic,jc,kc) + rhs(ic,jc,kc)
      enddo
      enddo
!$OMP  END PARALLEL DO
      enddo

      if(allocated(rhst)) deallocate(rhst)

      if(allocated(amkl)) deallocate(amkl)
      if(allocated(ackl)) deallocate(apkl)
      if(allocated(apkl)) deallocate(ackl)
      if(allocated(fkl)) deallocate(fkl)

      return
      end
