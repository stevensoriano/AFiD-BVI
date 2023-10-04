!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ReadMLSInput.F90                               !
!    CONTAINS: subroutine ReadMLSInput                    !
!                                                         ! 
!    PURPOSE: Read parameters from pert.in file           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ReadMLSInput
      use param
      use mpih
      use mls_param
      implicit none
      logical :: fexist
      character(len=4) :: dummy

      open(unit=15,file='mlspart.in',status='old')
        read(15,301) dummy       
        read(15,*) imlsfor, imlsstr, imlsref,pread
        read(15,301) dummy       
        read(15,*) wcheck, fcheck, wcub, wexp, nel
        read(15,301) dummy       
        read(15,*) wcon,wscl,sclf,rhop
        read(15,301) dummy       
        read(15,*) ke,kb,kv,kat,kal,kvis,Cv
        read(15,301) dummy       
        read(15,*) gtsfx,datfx
      close(15)

301     format(a4)                

      if(imlsfor.gt.0) mlsforcing = .true.
      if(imlsstr.gt.0) solvestructure = .true.


      return
      end

