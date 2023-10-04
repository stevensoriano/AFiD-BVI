      subroutine OpenLogs
      use param
      implicit none

       open(113,file='fforce.out',status='unknown',position='append' &
               ,access='sequential')
     
       open(114,file='fforce_nor.out',status='unknown',position='append' &
                ,access='sequential')
    
       open(115,file='viscf.out',status='unknown',position='append' &
               ,access='sequential')


!EP   dissipation.out in balance.F
      open(92,file='dissipation.out',status='unknown', &
       access='sequential',position='append')

!EP   rms_vel.out in stst.F
       open(94,file='rms_vel.out',status='unknown',position='append', &
       access='sequential')

      if(ireset.eq.1) then    
       rewind(92)
       rewind(94)
       rewind(113)
       rewind(114)
       rewind(115)
      endif

      return
      end   
      


      subroutine CloseLogs
      implicit none
      close(94)
      close(92)
      close(113)
      close(114)
      close(115)
      return 
      end
