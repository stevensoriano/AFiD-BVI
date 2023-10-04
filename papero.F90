      program papero
      use mpih
      use mpi_param
      use param
      use mls_param, only: ke,kb,kv,kat,kal,kvis,Cv
      use local_arrays, only: vy,vz,pr,vx
      use pointparticle, only: iresetp,toutpp,withppart
      implicit none
      character(len=4) :: dummy
      integer :: n,ns,nt,errorcode

      integer :: ntstf
      real    :: cflm,dmax,minwtdt
      real    :: ti(2), tin(3)

      call InitializeMPI

!RO   Calculate number of OpenMP Threads

      nt = 0
!$OMP  PARALLEL
!$OMP& REDUCTION(+:nt)
      nt = nt + 1
!$OMP  END PARALLEL

      if (ismaster) then 
        write(6,*) 'MPI tasks=', numtasks
        write(6,*) 'OMP threads per task=', nt
        write(6,*) 'No. of processors=', nt*numtasks
      end if

      call ReadInputFile
      call ReadMLSInput
      call ReadPpartInput
      
      call MpiBarrier

!RO   Initialize output files

      if(ismaster) call OpenLogs

!RO   Definition of PI

      pi=2.d0*dasin(1.d0)                          

      if(ismaster) then
!m====================================================                                                                             
      write(6,112)xlen/zlen,ylen/zlen
  112 format(//,20x,'H I T',//,10x, &
       '3D Cube with aspect-ratio:  L1/L3 = ',f5.2,'   L2/L3 = ',f5.2)
      write(6,202) ren
  202 format(/,5x,'Parameters: ',' Re=',e10.3)
      if(idtv.eq.1) then
         write(6,204) cflmax
  204 format(/,5x,'Variable dt and fixed cfl= ', &
       e11.4,/ )            
      else 
         write(6,205) dtmax,cfllim
  205 format(/,5x,'Fixed dt= ',e11.4,' and maximum cfl=', &
        e11.4,/ )            
      endif
!m====================================================    
      endif
!m======================================================
      
      call InitTimeMarchScheme

      call mpi_workdistribution

      call InitArrays

!======================================================
!
!     the solution of the problem starts
!
!======================================================

      tin(1) = MPI_WTIME()

!     Initialize MLS routines
      call tri_geo
      call allocate_mls_local
      call write_to_vtk

      call HdfStart

      call InitStats
      call CreateGrid 
      call InitPressSolv

      time=0.d0
      vmax=0.0d0

      if(myid.eq.0) then
      write(6,754)n1,n2,n3                                              
  754 format(/,5x,'grid resolution: ',' n1= ',i5,' n2= ',i5, &
       ' n3= ',i5/)                       
      write(6,755) 1.d0/dx1,1.d0/dx2,1.d0/dx3,dt,ntst                  
  755 format(/,2x,' dx1=',e10.3,' dx2=',e10.3,' dx3=',e10.3,' dt=' &
       ,e10.3,' ntst=',i7,/)
      endif

!===================================                                                                       
!   Read or create initial fields                                       
!===================================  


!     create the initial conditions

      if(nread.eq.0) then

       if(ismaster) write(6,*)' nread=0 ---> create initial conditions'

       ntime=0                                                         
       time=0.d0
       cflm=0.d0
         
       call CreateInitialConditions

       call update_both_ghosts(n1,n2,pr,kstart,kend)
       call update_both_ghosts(n1,n2,vx,kstart,kend)
       call update_both_ghosts(n1,n2,vy,kstart,kend)
       call update_both_ghosts(n1,n2,vz,kstart,kend)
        

!=================================
!     Or read initial condition
!=================================

      else

       if(ismaster) write(6,*)' nread=1 ---> Read initial conditions'

       call ReadInitCond

       call update_both_ghosts(n1,n2,pr,kstart,kend)
       call update_both_ghosts(n1,n2,vx,kstart,kend)
       call update_both_ghosts(n1,n2,vy,kstart,kend)
       call update_both_ghosts(n1,n2,vz,kstart,kend)
       tmax = tmax + time


      endif                                                             

!=================================
! Check velocity divergence of initial conditions
!=================================

      call CheckDivergence(dmax)

      if(ismaster) write(6,*)' initial divg dmax  ',dmax

!=================================
      ntstf=ntst                                                   
!m================================
      if(ismaster) then
       write(6,711) tframe,ntstf,tpin
711    format(3x,'check in cond : tframe =',f10.1,  &
             '  ntstf =',i8,2x,'tpin =',f10.1//)
      endif
!m================================ 
!       If variable timestep

      if(idtv.eq.1) then

       if(ismaster) then
        write(6,*)ntime,time,vmax(1),vmax(2),vmax(3),&
         dt,dmax
        endif

      else

       if(ismaster) then
        write(6,*)ntime,time,vmax(1),vmax(2),vmax(3),&
        cflm,dmax
       endif
     
       cflm=cflm*dt

       endif

      call InitRandomForce

      tin(2) = MPI_WTIME()

      if(ismaster) then
       write(6,*) 'Initialization Time = ', tin(2) -tin(1), ' sec.'
      endif

!
!  ********* starts the time dependent calculation ***

      errorcode = 0 ! set errocode to 0 (OK)
      minwtdt = huge(0.0d0) ! initialize minimum time step walltime

      do ntime=0,ntstf                                           

       ti(1) = MPI_WTIME()

       call CalcMaxCFL(cflm)

!      Adjust timestep or exit if too small 
            
       if(idtv.eq.1) then
         if(ntime.ne.1) then
            dt=cflmax/cflm
            if(dt.gt.dtmax) dt=dtmax
         endif
         if(dt.lt. 0.00000001d0) errorcode=166 
       else
         cflm=cflm*dt
         if(cflm.gt.cfllim) errorcode=165 
       endif


       call CalcHITRandomForce

!RO    March in time one step

       call TimeMarcher

       time=time+dt

!RO    Here we come if on first time-step or every tpin and do some
!RO    routine checks and output

       if((ntime.eq.1).or.(mod(time,tpin).lt.dt)) then 

        call write_to_vtk

        call CheckMaxVel !Make sure velocities are not exceeding maximum
        if(vmax(1).gt.1000.d0.and.vmax(2).gt.1000.d0) errorcode=266 

        call CalcMaxCFL(cflm) !Recalculate CFL 
        if(idtv.ne.1) cflm=cflm*dt 

        call CheckDivergence(dmax) !Make sure velocity is solenoidal
        if(dmax.gt.resid) errorcode=169 !If velocity field not solenoidal, stop

!RO    Calculate statistics if Time > TSTA
!RO     (Initialization/Transient)

        if(time.gt.tsta) then
          call UpdateStats
          call CalcDissipation
        endif

       end if

       if((imovie.eq.1).and.(mod(time,tmovie).lt.dt)) call DumpMovieFiles

       if(withppart.and.(mod(time,toutpp).lt.dt)) then 
          call part_write_continua(.true.)
       end if

       if(time.gt.tmax) errorcode=333

!RO    Output wall-time of time-step

       ti(2) = MPI_WTIME()
       minwtdt = min(minwtdt,ti(2) - ti(1))

        if(mod(time,tpin).lt.dt) then
          if(ismaster) then
          write(6,*) 'Maximum divergence = ', dmax
          write(6,*) 'ntime - time - vmax(1) - vmax(2) - vmax(3)'
          write(6,*) ntime,time,vmax(1),vmax(2),vmax(3)
          write(6,'(a,f8.3,a)') 'Minimum Iteration Time = ', minwtdt,  &
                  ' sec.'
          endif
          minwtdt = huge(0.0d0)
        endif

      if( (ti(2) -tin(1)) .gt. walltimemax) errorcode=334

      if( ntime .eq. ntst ) errorcode = 555

      call MpiBcastInt(errorcode)

      if(errorcode.ne.0) then

!EP    dt too small
        if(errorcode.eq.166) call QuitRoutine(tin,.false.,errorcode)

!EP   cfl too high    
        if(errorcode.eq.165) call QuitRoutine(tin,.false.,errorcode)

!EP   velocities diverged
        if(errorcode.eq.266) call QuitRoutine(tin,.false.,errorcode)

!EP   mass not conserved
        if(errorcode.eq.169) call QuitRoutine(tin,.false.,errorcode)

!EP   Physical time exceeded tmax, no error; normal quit
        if(errorcode.eq.333) call QuitRoutine(tin,.true.,errorcode)

!EP   walltime exceeded walltimemax, no error; normal quit
        if(errorcode.eq.334) call QuitRoutine(tin,.true.,errorcode)

!RS   FFT input not correct
        if(errorcode.eq.444) call QuitRoutine(tin,.false.,errorcode)

!RS   maximum number of timesteps reached, no error; normal quit
        if(errorcode.eq.555) call QuitRoutine(tin,.true.,errorcode)

        errorcode = 100 !EP already finalized

        exit

      endif
 
      end do


!     Store structure output
      call continua_str1

      call QuitRoutine(tin,.true.,errorcode)

      stop
      end
