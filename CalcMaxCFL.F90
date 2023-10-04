!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcMaxCFL.F90                                 !
!    CONTAINS: subroutine CalcMaxCFL                      !
!                                                         ! 
!    PURPOSE: Calculates the maximum value of U/dx to     !
!     determine the new value of the time-step            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcMaxCFL(cflm)
      use param
      use local_arrays, only: vx,vy,vz
      use mpih
      use mpi_param, only: kstart,kend
      implicit none
      real,intent(inout)    :: cflm
      integer :: j,k,jp,kp,i,ip
      real :: qcf
      
      cflm=0.00000001d0
                                                                       
      do k=kstart,kend
        kp=k+1
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(j,i,jp,ip,qcf)
!$OMP& REDUCTION(max: cflm)
        do j=1,n2m
          jp=j+1
          do i=1,n1m
            ip=i+1
            qcf=( abs((vx(i,j,k)+vx(ip,j,k))*0.5d0*dx1)  &
                 +abs((vy(i,j,k)+vy(i,jp,k))*0.5d0*dx2)  &
                 +abs((vz(i,j,k)+vz(i,j,kp))*0.5d0*dx3))

            cflm = max(cflm,qcf)
      enddo
      enddo
!$OMP  END PARALLEL DO
      enddo
            
      call MpiAllMaxRealScalar(cflm)


      return  
      end                                                               
