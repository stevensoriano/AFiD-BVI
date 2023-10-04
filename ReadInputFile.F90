!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ReadInputFile.F90                              !
!    CONTAINS: subroutine ReadInputFile                   !
!                                                         ! 
!    PURPOSE: Read parameters from bou.in file            !
!     and pert.in file                                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ReadInputFile
      use param
      use mpih
      use mls_param
      implicit none
      logical :: fexist
      character(len=4) :: dummy

      open(unit=15,file='bou.in',status='old')
        read(15,301) dummy
        read(15,*) n1m,n2m,n3m,nsst,nwrit,nread
        read(15,301) dummy
        read(15,*) ntst,tframe,tpin,tmax,walltimemax,ireset
        read(15,301) dummy
        read(15,*) xlen,ylen,zlen
        read(15,301) dummy
        read(15,*) ren,dt,resid,cflmax
        read(15,301) dummy
        read(15,*) tsta,starea
        read(15,301) dummy
        read(15,*) tl, epsstar, kfmax
        read(15,301) dummy       
        read(15,*) idtv,dtmax,cfllim  
        read(15,301) dummy       
        read(15,*) imovie,tmovie
        read(15,301) dummy       
        read(15,*) vortexrad, wirev
301     format(a4)                
      close(15)

      n1=n1m+1                                                          
      n2=n2m+1  
      n3=n3m+1
      n1mh = n1m/2 + 1
      n2mh = n2m/2 + 1

      return
      end
