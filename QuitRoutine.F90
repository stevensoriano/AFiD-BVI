!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: QuitRoutine.F90                                !
!    CONTAINS: subroutine QuitRoutine, NotifyError        !
!                                                         ! 
!    PURPOSE: Routines to exit the program and write the  !
!     data if necessary                                   !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine QuitRoutine(tin,normalexit,errorcode)
      use hdf5
      use mpih
      use param
      use pointparticle, only: withppart
      implicit none
      logical, intent(in) :: normalexit
      integer :: errorcode
      real :: tin(3)

      if(errorcode.ne.100) then !EP skip if already finalized

      tin(3) = MPI_WTIME()
      if(ismaster) write(6,'(a,f10.2,a)') 'Total Iteration Time = ',tin(3) -tin(2),' sec.'

      if(ismaster) then
       call NotifyError(errorcode)
      endif

      if(normalexit) then
       call WriteStats
       call CalcVorticity
       call mpi_write_continua
       call WriteRandForcCoef
      else
       call MpiAbort
      endif

      call dfftw_destroy_plan(fwd_plan)
      call dfftw_destroy_plan(bck_plan)

      if(withppart) call part_write_continua(.false.)

      call HdfClose

      if(ismaster) then
        open(27,file="Total_time.out")
        write(27,*)"Total simulation time in sec.: ", tin(3)-tin(1)
        close(27)
      endif
 

      if(ismaster) call CloseLogs

      call mem_dealloc
      call FinalizeMPI

      endif

      end subroutine QuitRoutine

!====================================================================================

      subroutine NotifyError(errorcode)
      use param
      implicit none
      integer, intent(in) :: errorcode

      if(errorcode.eq.166) then
        write(6,168) dt
168     format(10x,'dt too small, DT= ',e14.7)
      else if(errorcode.eq.165) then
        write(6,164)
164     format(10x,'cfl too large  ')
      else if(errorcode.eq.266) then
        write(6,268)
268     format(10x,'velocities diverged')
      else if(errorcode.eq.169) then
        write(6,178)
        write(6,179)
        write(6,180)
178     format(10x,'too large local residue for mass conservation at:')
179     format(10x,'Probably the matrix in SolvePressureCorrection')
180     format(10x,'is singular. Try changing nxm or str3')
        call LocateLargeDivergence
      else if(errorcode.eq.333) then
         write(*,*) "time greater than tmax"
         write(*,*) "statistics and continuation updated"
      else if(errorcode.eq.334) then
         write(*,*) "walltime greater than walltimemax"
         write(*,*) "statistics and continuation updated"
      else if(errorcode.eq.444) then
         write(*,*) "FFT size in ny or nz is not efficient"
      else
         write(*,*) "Maximum number of timesteps reached"
      end if

      return
      end

