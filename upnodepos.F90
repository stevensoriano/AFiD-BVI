!----------------------------------------------------
!     Update position of nodes from forces
!----------------------------------------------------

      subroutine upnodepos
      USE param
      USE mls_param
      USE mpi_param, only: kstart, kend
      USE mls_local, only: q1c, q2c, q3c, for_xc, for_yc, for_zc
      IMPLICIT NONE

      integer i,inp

      do inp=1,Nparticle
       do i=1,maxnv

      xyza(inp,1:3,i)=fxyz(inp,1:3,i)*uspm(inp)
!
      xyzv(inp,1:3,i)=xyzv(inp,1:3,i)+0.50*(dt/float(nsstep))
     %                            *(xyza(inp,1:3,i)+xyza0(inp,1:3,i))
!
      xyz(inp,1:3,i)=xyz(inp,1:3,i)+0.50*(dt/float(nsstep))
     %                     *(xyzv(inp,1:3,i)+xyzv0(inp,1:3,i))
!


       end do
      end do

#ifdef MLSDEBUG      
      print*,'Positions',maxval(xyz(:,:,:)) 
      print*,'Vels',maxval(xyzv(:,:,:)) 
#endif

      xyza0(:,:,:) = xyza(:,:,:)
      xyzv0(:,:,:) = xyzv(:,:,:)


      return
      end 
