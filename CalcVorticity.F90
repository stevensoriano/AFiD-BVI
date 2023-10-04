!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcVorticity.F90                              !
!    CONTAINS: subroutine CalcMaxCFL                      !
!                                                         ! 
!    PURPOSE: Calculates the vorticity modulus and        !
!     stores it in enstro array                           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcVorticity
      use mpih
      use param
      use local_arrays, only: vx,vy,vz
      use local_aux, only: vorx, vory, vorz
      use mpi_param, only: kstart,kend
      implicit none
      integer :: kc,kp,jp,jm,jc,ic,im,ip,km

      do kc=kstart,kend
      kp=kc+1
      km=kc-1
       do jc=1,n2m
        jp=jpv(jc)
        jm=jmv(jc)
        do ic=1,n1m
         ip=ipv(ic)
         im=imv(ic)

         vorz(ic,jc,kc) =   & 
          ((vy(ip,jc,kc)+vy(ip,jp,kc)-vy(im,jc,kc)-vy(im,jp,kc))*dx1- &
           (vx(ip,jp,kc)+vx(ic,jp,kc)-vx(ip,jm,kc)-vx(ic,jm,kc))*dx2)*0.25

         vory(ic,jc,kc) =   &
           ((vx(ip,jc,kp)+vx(ic,jc,kp)-vx(ip,jc,km)-vx(ic,jc,km))*dx3- &
            (vz(ip,jc,kc)+vz(ip,jc,kp)-vz(im,jc,kc)-vz(im,jc,kp))*dx1)*0.25

          vorx(ic,jc,kc) =  & 
           ((vz(ic,jp,kc)+vz(ic,jp,kp)-vz(ic,jm,kc)-vz(ic,jm,kp))*dx2- &
            (vy(ic,jp,kp)+vy(ic,jc,kp)-vy(ic,jp,km)-vy(ic,jc,km))*dx3)*0.25

          end do
        end do
      end do

      return
      end


