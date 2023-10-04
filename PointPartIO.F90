!----------------------------------------------------
!     Routine for reading in particle position 
!     and velocity data.
!----------------------------------------------------

      subroutine part_read_continua

      USE hdf5
      USE param
      USE pointparticle

      IMPLICIT NONE
      character*30 :: filename,dsetnam

      filename = trim('outputdir/parthist.h5')
 
      if(ismaster) then

       dsetnam=trim('xp')
       call HdfSerialReadReal2D(dsetnam,filename,xp,Npointpart,3)
 
       dsetnam=trim('q1p')
       call HdfSerialReadReal1D(dsetnam,filename,q1p,Npointpart)
       dsetnam=trim('q2p')
       call HdfSerialReadReal1D(dsetnam,filename,q2p,Npointpart)
       dsetnam=trim('q3p')
       call HdfSerialReadReal1D(dsetnam,filename,q3p,Npointpart)

       dsetnam=trim('vort1')
       call HdfSerialReadReal1D(dsetnam,filename,vort1,Npointpart)
       dsetnam=trim('vort2')
       call HdfSerialReadReal1D(dsetnam,filename,vort2,Npointpart)
       dsetnam=trim('vort3')
       call HdfSerialReadReal1D(dsetnam,filename,vort3,Npointpart)

       dsetnam=trim('kalb1')
       call HdfSerialReadReal1D(dsetnam,filename,kalb1,Npointpart)
       dsetnam=trim('kalb2')
       call HdfSerialReadReal1D(dsetnam,filename,kalb2,Npointpart)
       dsetnam=trim('kalb3')
       call HdfSerialReadReal1D(dsetnam,filename,kalb3,Npointpart)

      end if

      call MpiBcastReal2DArray(xp,Npointpart,3)
      call MpiBcastReal1DArray(q1p,Npointpart)
      call MpiBcastReal1DArray(q2p,Npointpart)
      call MpiBcastReal1DArray(q3p,Npointpart)
      call MpiBcastReal1DArray(vort1,1,Npointpart)
      call MpiBcastReal1DArray(vort2,1,Npointpart)
      call MpiBcastReal1DArray(vort3,1,Npointpart)
      call MpiBcastReal1DArray(kalb1,1,Npointpart)
      call MpiBcastReal1DArray(kalb2,1,Npointpart)
      call MpiBcastReal1DArray(kalb3,1,Npointpart)

      call MpiBarrier

      return
      end subroutine part_read_continua

!----------------------------------------------------
!     Routine for writing particle position 
!     and velocity data.
!----------------------------------------------------

      subroutine part_write_continua(ismidrun)

      USE hdf5
      USE param
      USE pointparticle

      IMPLICIT NONE
      character*50 :: filename,dsetnam
      character*5  :: ipfi
      integer hdf_error, itime
      real tprfi
      logical :: ismidrun

      if(ismaster) then 

      if(ismidrun) then 
        tprfi = 1/tframe
        itime=nint(time*tprfi)
        write(ipfi,82)itime
   82   format(i5.5)
        filename = trim('outputdir/flowmov/parthist_'//ipfi//'.h5')
      else 
        filename = trim('outputdir/parthist.h5')
      end if

      call HdfCreateBlankFile(filename)

      dsetnam=trim('xp')
      call HdfSerialWriteReal2D(dsetnam,filename,xp,Npointpart,3)
 
      dsetnam=trim('renp')
      call HdfSerialWriteReal1D(dsetnam,filename,renp,Npointpart)

      dsetnam=trim('q1p')
      call HdfSerialWriteReal1D(dsetnam,filename,q1p,Npointpart)
      dsetnam=trim('q2p')
      call HdfSerialWriteReal1D(dsetnam,filename,q2p,Npointpart)
      dsetnam=trim('q3p')
      call HdfSerialWriteReal1D(dsetnam,filename,q3p,Npointpart)

      dsetnam=trim('vort1')
      call HdfSerialWriteReal1D(dsetnam,filename,vort1,Npointpart)
      dsetnam=trim('vort2')
      call HdfSerialWriteReal1D(dsetnam,filename,vort2,Npointpart)
      dsetnam=trim('vort3')
      call HdfSerialWriteReal1D(dsetnam,filename,vort3,Npointpart)

      dsetnam=trim('kalb1')
      call HdfSerialWriteReal1D(dsetnam,filename,kalb1,Npointpart)
      dsetnam=trim('kalb2')
      call HdfSerialWriteReal1D(dsetnam,filename,kalb2,Npointpart)
      dsetnam=trim('kalb3')
      call HdfSerialWriteReal1D(dsetnam,filename,kalb3,Npointpart)


      dsetnam=trim('aap')
      call HdfSerialWriteReal2D(dsetnam,filename,aap,Npointpart,3)
      dsetnam=trim('facc_for')
      call HdfSerialWriteReal2D(dsetnam,filename,facc_for,Npointpart,3)
      dsetnam=trim('drag_for')
      call HdfSerialWriteReal2D(dsetnam,filename,drag_for,Npointpart,3)
      dsetnam=trim('lift_for')
      call HdfSerialWriteReal2D(dsetnam,filename,lift_for,Npointpart,3)

      end if
      call MpiBarrier

      return
      end subroutine part_write_continua
