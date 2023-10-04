!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ReadPpartInput.F90                             !
!    CONTAINS: subroutine ReadPpartInput                  !
!                                                         ! 
!    PURPOSE: Read parameters from ppart.in file          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ReadPpartInput
      use param
      use mpih
      use pointparticle
      implicit none
      character(len=4) :: dummy

      open(unit=15,file='ppart.in',status='old')
        read(15,301) dummy
        read(15,*) Npointpart,timeONp,iresetp,toutpp,p1p2
        read(15,301) dummy
        read(15,*) usfroude,dbd,dbd2,rhohat1,rhohat2,stokes1,stokes2
        read(15,301) dummy
        read(15,*) cpi(1),cpi(2),cpi(3)
        read(15,301) dummy
        read(15,*) cpf(1),cpf(2),cpf(3)
301     format(a4)                
      close(15)

      if(Npointpart.gt.0) withppart = .true.

      return
      end
