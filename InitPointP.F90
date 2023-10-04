!-----------------------------------------------------------------------
!     Initialization of particles by master node.
!     - Find grid positions of each particle and communicate to all procs.
!------------------------------------------------------------------------

      subroutine InitPointP

      USE param
      USE mls_param
      USE mpih
      USE mpi_param
      USE pointparticle

      IMPLICIT NONE

      real      :: gp1, gp2, d1, d2
      real      :: len_beg(3),len_end(3)
      real      :: pos(3),qInt(7)
      real      :: dbub1,dbub2,Volbub,Volcell,volfr
      integer   :: bubind(6)
      integer   :: inp,seed,i,j,k,numpart1

      call AllocatePointPArrays


      gp1 = 3.d0/(1.d0+2.d0*rhohat1)
      gp2 = 3.d0/(1.d0+2.d0*rhohat2)
      d1 = dbd
      d2 = dbd2 
!      d2 = sqrt(gp2/gp1)*d1

      numpart1=floor(dble(Npointpart*p1p2))

      if(ismaster)then
       do i=1,Npointpart
         if(i.le.numpart1) then
         gammap(i) = gp1
         dbd12(i) = d1
         stokes(i) = stokes1
        else 
         gammap(i) = gp2
         dbd12(i) = d2
         stokes(i) = stokes2
        end if
       end do
      end if

      call MpiBcastReal1DArray(gammap,Npointpart)
      call MpiBcastReal1DArray(dbd12,Npointpart)
      call MpiBarrier


      ! CS  Specificy Stokes number based on ratio of particle to gap width
      ! St = d_p^2/(12*gamma*\nu*\tau_s)
      !    where, \tau_s \equiv \tau_\eta = 1 (say)
!      do i=1,Npointpart
!       stokes(i)=dbd12(i)**2*ren/(gammap(i)*12.d0)
!      end do

      ! CS  If loop to reset particles
      IF (iresetp.eq.1) then
!     -----------------------------------------------------
!     Master process initializing random particle positions
!     -----------------------------------------------------

      if(ismaster)then
       print*,'  --------------------------------------------'
       write(*,753)Npointpart 
       write(*,758)minval(gammap(:)),maxval(gammap(:))
       write(*,755)minval(stokes(:)),maxval(stokes(:))
       write(*,756)usfroude
       print*,'  --------------------------------------------'
753   format(4x,'Particle initiated in code --- No. of Bubbles: ',i5)

      seed = 361
      dbub1=dbd12(1)*zlen
      dbub2=dbd12(Npointpart)*zlen

      Volbub = float(numpart1)*4.d0/3.d0*pi*dbub1**3/8.d0 &
              +float(Npointpart-numpart1)*4.d0/3.d0*pi*dbub2**3/8.d0
      Volcell = xlen*ylen*zlen
      volfr = Volbub/Volcell
      volfr = 100.0*volfr
      write(*,757)volfr

!     ------------------------------------------
      ! CS These 6 lines are for the range of seeding, nothing special
      len_end(1)=cpf(1)*xlen
      len_end(2)=cpf(2)*ylen
      len_end(3)=cpf(3)*zlen

      len_beg(1)=cpi(1)*xlen
      len_beg(2)=cpi(2)*ylen
      len_beg(3)=cpi(3)*zlen

!     -------------------Random Distribution------------------

      do inp=1,Npointpart

       xp(inp,1) = len_beg(1) + ran(seed)*(len_end(1)-len_beg(1))
       xp(inp,2) = len_beg(2) + ran(seed)*(len_end(2)-len_beg(2))
       xp(inp,3) = len_beg(3) + ran(seed)*(len_end(3)-len_beg(3))

       seed=seed+(seed*inp)

       xpo(inp,1)=xp(inp,1)
       xpo(inp,2)=xp(inp,2)
       xpo(inp,3)=xp(inp,3) 

      end do

      print*,'   Bubbles initiated at random rositions by rank 0'
      end if

!     End of master process initializing

      call MpiBarrier

      call MpiBcastReal1DArray(xp, 3*Npointpart)
      call MpiBcastReal1DArray(xpo, 3*Npointpart)


!     -----------------------------------------------------------------
!     Initiate particle velocity
!     -----------------------------------------------------------------
      do inp=1,Npointpart
      
        pos(1)=xp(inp,1);pos(2)=xp(inp,2);pos(3)=xp(inp,3)

        ! CS Important routine to get particle index
        call indpos(pos,inp,bubind) ! CS now you have the index of particle inp

        ! CS Calculate only if particle is in slab of processor
        if(bubind(6).ge.kstart-1.and.bubind(3).le.kend-1)then
!CS       TODO bubind(3) or bubind(6)?

          pos(1)=xp(inp,1);pos(2)=xp(inp,2);pos(3)=xp(inp,3)

          ! CS Important routine to get interpolates particle velocities
          call interpol(pos,bubind,inp,qInt)

          q1p(inp)=qInt(1);q2p(inp)=qInt(2);q3p(inp)=qInt(3)
          qVal1(inp)=qInt(1);qVal2(inp)=qInt(2);qVal3(inp)=qInt(3)

        else

          q1p  (inp)=0.d0; q2p  (inp)=0.d0; q3p  (inp)=0.d0
          qVal1(inp)=0.d0; qVal2(inp)=0.d0; qVal3(inp)=0.d0

        end if
     
      end do

      vort1(:)=0.d0 ; vort2(:)=0.d0 ; vort3(:)=0.d0 
      kalb1(:)=0.d0 ; kalb2(:)=0.d0 ; kalb3(:)=0.d0

      call MpiBarrier
      call MpiAllSumReal1D(q1p,Npointpart)
      call MpiAllSumReal1D(q2p,Npointpart)
      call MpiAllSumReal1D(q3p,Npointpart)

!     --------------------------------------------------------------------
      ELSE
!     -----------------------------------------------------
!     Read previous particle positions
!     -----------------------------------------------------
      
      if(ismaster) then
     
       print*,'  --------------------------------------------'
       write(*,754)Npointpart 
       write(*,755)minval(stokes(:)),maxval(stokes(:))
       write(*,756)usfroude
       print*,'  --------------------------------------------'

754    format(4x,'Particle data from file --- No. of Bubbles: ',i5)
758    format(4x,'Particle mass ratios, min:',e10.3,' max:',e10.3)
755    format(4x,'Particle time scale, min:',e10.3,' max:',e10.3)
756    format(4x,'Inverse Froude number: ',e10.3)

       seed = 361
       dbub1=dbd12(1)*zlen
       dbub2=dbd12(Npointpart)*zlen

       Volbub = float(numpart1)*4.d0/3.d0*pi*dbub1**3/8.d0 &
               +float(Npointpart-numpart1)*4.d0/3.d0*pi*dbub2**3/8.d0
       Volcell = xlen*ylen*zlen
       volfr = Volbub/Volcell
       volfr = 100.0*volfr
       write(*,757)volfr
757    format(4x,'Volume fraction of bubbles (%): ',e10.3)

      end if

      END IF

      return
      end
