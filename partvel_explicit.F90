!----------------------------------------------------------------------
!     explicit  solver for bubble velocity (drag gravity,lift,addded mass)
!     solved by each node from the position info received from master
!-----------------------------------------------------------------------

      subroutine CalcPointPVel
      USE param
      USE mpih
      USE mls_param
      USE mpi_param
      use pointparticle

      IMPLICIT NONE

      integer inp
      integer bubind(6)
      real pos(3),qInt(7),dmatInt(3)
      real uxn,uyn,uzn,vxn,vyn,vzn,omexn,omeyn,omezn
      real uxn1,uyn1,uzn1,vxn1,vyn1,vzn1,omexn1,omeyn1,omezn1
      real kalbxn,kalbyn,kalbzn,kalbxn1,kalbyn1,kalbzn1
      real tau_b,fdrag,renp_temp
      real Forxn,Foryn,Forzn,Forxn1,Foryn1,Forzn1
      real Adn1                 ! Diagonal matrix entries
      real Aoxn1,Aoyn1,Aozn1    ! Off-diagonal matrix entries
      real RHSx,RHSy,RHSz

      integer i,j,k,ist,jst,kst
      integer dum,inb,li,lj,lk,istwoway
      real fBox(8,3),matderpart(3),posp(3)
      real Volbub,dbub,Volcell
      real dummy

!     print*,'Beginning bubblesolve'      
      aap(1:Npointpart,1:3) = 0.d0

!      buoy_for(1:Npointpart,1:3) = 0.d0
      facc_for(1:Npointpart,1:3) = 0.d0
      drag_for(1:Npointpart,1:3) = 0.d0
      lift_for(1:Npointpart,1:3) = 0.d0

      renp(:) = 0.d0

      istwoway = 0

      do inp=1,Npointpart

!     ------------------------------------------------
      Forxn=0.0 ; Foryn=0.0 ; Forzn=0.0
      Forxn1=0.0; Foryn1=0.0; Forzn1=0.0
      Aoxn1=0.0 ; Aoyn1=0.0 ; Aozn1=0.0
!     ------------------------------------------------
      ! Get positions
      pos(1)=xp(inp,1)
      pos(2)=xp(inp,2)
      pos(3)=xp(inp,3)

      ! Get indices from indpos
      call indpos(pos,inp,bubind)

      ! Only specific processors will calculate particles
      if(bubind(6).ge.kstart-1.and.bubind(6).le.kend-1)then

      ! Interpolate velocity and vorticity onto bubble position
      call interpol(pos,bubind,inp,qInt)

      ! Interpolate material derivatives onto bubble position
      call interdmat(pos,bubind,inp,dmatInt)

      ! Store values for simplicity
      qVal1(inp)=qInt(1)    !u_{f,1,interp}
      qVal2(inp)=qInt(2)    !u_{f,2,interp}
      qVal3(inp)=qInt(3)    !u_{f,3,interp}

      vort1(inp)=qInt(5)    !vort_q1_interp
      vort2(inp)=qInt(6)
      vort3(inp)=qInt(7)

      kalb1(inp)=dmatInt(1) !D/dt
      kalb2(inp)=dmatInt(2)
      kalb3(inp)=dmatInt(3)

!     -----------------------------------------------------------------
      uxn=qValo1(inp) ; uxn1=qVal1(inp) !Set old and new flow velocities
      uyn=qValo2(inp) ; uyn1=qVal2(inp)
      uzn=qValo3(inp) ; uzn1=qVal3(inp)

      omexn=vorto1(inp) ; omexn1=vort1(inp)
      omeyn=vorto2(inp) ; omeyn1=vort2(inp)
      omezn=vorto3(inp) ; omezn1=vort3(inp)

      kalbxn=kalbo1(inp) ; kalbxn1=kalb1(inp)
      kalbyn=kalbo2(inp) ; kalbyn1=kalb2(inp)
      kalbzn=kalbo3(inp) ; kalbzn1=kalb3(inp)
      
      vxn=q1po(inp) ; vyn=q2po(inp) ; vzn=q3po(inp)
!     -----------------------------------------------------------------
!      fdrag = gammap(inp)/stokes(inp)

      renp_temp=sqrt(((uxn-vxn)**2)+((uyn-vyn)**2)+((uzn-vzn)**2))* &
                  (dbd12(inp)*ren)
      renp(inp)=renp_temp

!     -------------------calculate [A]--------------------------------

      Adn1 = 1.d0
      Aoxn1 = 0.d0
      Aoyn1 = 0.d0
      Aozn1 = 0.d0

!     ------------ Check each component ----------------------------
      facc_for(inp,1) = kalbxn*gammap(inp)
      facc_for(inp,2) = kalbyn*gammap(inp)
      facc_for(inp,3) = kalbzn*gammap(inp)

!      drag_for(inp,1) = (uxn-vxn)*gammap(inp)/stokes(inp)
!      drag_for(inp,2) = (uyn-vyn)*gammap(inp)/stokes(inp)
!      drag_for(inp,3) = (uzn-vzn)*gammap(inp)/stokes(inp)

      drag_for(inp,1) = (uxn-vxn)/stokes(inp)
      drag_for(inp,2) = (uyn-vyn)/stokes(inp)
      drag_for(inp,3) = (uzn-vzn)/stokes(inp)

      lift_for(inp,1) = 0.d0
      lift_for(inp,2) = 0.d0
      lift_for(inp,3) = 0.d0
!!      lift_for(inp,1) = -((vyn-uyn)*omezn   &
!!                        - (vzn-uzn)*omeyn)*gammap(inp)/3.d0
!!      lift_for(inp,2) = -((vzn-uzn)*omexn   &
!!                        - (vxn-uxn)*omezn)*gammap(inp)/3.d0 
!!      lift_for(inp,3) = -((vxn-uxn)*omeyn   &
!!                        - (vyn-uyn)*omexn)*gammap(inp)/3.d0

!     --------------------calculate [Fn,Fn+1]----------------------------

      Forxn  =  facc_for(inp,1) + drag_for(inp,1) + lift_for(inp,1)
      Foryn  =  facc_for(inp,2) + drag_for(inp,2) + lift_for(inp,2)
      Forzn  =  facc_for(inp,3) + drag_for(inp,3) + lift_for(inp,3)

!     ------------ Calculate [RHS] implicitly ----------------------

      RHSx = vxn + dt*(Forxn)
      RHSy = vyn + dt*(Foryn)
      RHSz = vzn + dt*(Forzn)

      aap(inp,1) = Forxn
      aap(inp,2) = Foryn
      aap(inp,3) = Forzn

!     --------------------calculate [V]-----------------------------

      vxn1 = (RHSx*(Adn1**2+Aoxn1**2)           &
            + RHSy*(Aoxn1*Aoyn1-Adn1*Aozn1)     &
            + RHSz*(Adn1*Aoyn1+Aoxn1*Aozn1))    &
            /(Adn1*(Adn1**2+Aoxn1**2+Aoyn1**2+Aozn1**2))

      vyn1 = (RHSx*(Aoxn1*Aoyn1+Adn1*Aozn1)     &
            + RHSy*(Adn1**2+Aoyn1**2)           &
            + RHSz*(-Adn1*Aoxn1+Aoyn1*Aozn1))   &
            /(Adn1*(Adn1**2+Aoxn1**2+Aoyn1**2+Aozn1**2))

      vzn1 = (RHSx*(-Adn1*Aoyn1+Aoxn1*Aozn1)    &
            + RHSy*(Adn1*Aoxn1+Aoyn1*Aozn1)     &
            + RHSz*(Adn1**2+Aozn1**2))          &
            /(Adn1*(Adn1**2+Aoxn1**2+Aoyn1**2+Aozn1**2))

!      print*,'Bubble Vel',vxn1,vyn1,vzn1

      q1p(inp)=vxn1
      q2p(inp)=vyn1
      q3p(inp)=vzn1

      if(istwoway.eq.1) then
!     -----------------------------------------------------------------
!     ================================================================
!     two way coupling
!     ================================================================

!     ------------------ Bubble indices -------------------------------
      i=bubind(1)   ; j=bubind(2)   ; k=bubind(3)
      ist=bubind(4) ; jst=bubind(5) ; kst=bubind(6)

!     ----- Percentage calculations for trilinear interpolation -------
      
      posp(1)=(pos(1)-xm(ist))/(xm(2)-xm(1))
      posp(2)=(pos(2)-ym(jst))/(ym(2)-ym(1))
      posp(3)=(pos(3)-zm(kst))/(zm(2)-zm(1))

      if(jst.eq.0.or.jst.eq.n2m)then
      print*,'Interpol Error - parvel.F',jst
      stop
      end if

      posp(1)=1.d0-posp(1); posp(2)=1.d0-posp(2); posp(3)=1.d0-posp(3)

!     ---------------- Cell Volume -------------------------
      Volcell = 1.d0/(dx1*dx2*dx3)

      dbub=dbd12(inp)
      Volbub = 1.d0/6.d0*pi*dbub*dbub*dbub

!      print*,'Volumes',Volbub,Volcell

!     ------ Forcing term at particle position [units = velocity] ------

!CS   Based on no-slip velocity (Stokes drag)

      matderpart(1) = (Volbub/Volcell)*((kalb1(inp)+kalbo1(inp))*0.5d0  &
                            + (gammap(inp)-1.d0)*usfroude)
      matderpart(2) = (Volbub/Volcell)*(kalb2(inp)+kalbo2(inp))*0.5d0
      matderpart(3) = (Volbub/Volcell)*(kalb3(inp)+kalbo3(inp))*0.5d0

!     ------------allocating all 3 forces to fBox------------
      do dum=1,3
      
       fBox(1,dum)=posp(1)*posp(2)*posp(3)*matderpart(dum)
       fBox(2,dum)=(1-posp(1))*posp(2)*posp(3)*matderpart(dum)
       fBox(3,dum)=posp(1)*(1-posp(2))*posp(3)*matderpart(dum)
       fBox(4,dum)=(1-posp(1))*(1-posp(2))*posp(3)*matderpart(dum)
       fBox(5,dum)=posp(1)*posp(2)*(1-posp(3))*matderpart(dum)
       fBox(6,dum)=(1-posp(1))*posp(2)*(1-posp(3))*matderpart(dum)
       fBox(7,dum)=posp(1)*(1-posp(2))*(1-posp(3))*matderpart(dum)
       fBox(8,dum)=(1-posp(1))*(1-posp(2))*(1-posp(3))*matderpart(dum)

      end do
!     -------------------------------------------------------
      inb=1
      do lk=0,1
        do lj=0,1
          do li=0,1

!         X-forcing
          for_xc_part(ist+li,jst+lj,kst+lk)=    &
          for_xc_part(ist+li,jst+lj,kst+lk) + fBox(inb,1)

!         Y-forcing
          for_yc_part(ist+li,jst+lj,kst+lk)=    &
          for_yc_part(ist+li,jst+lj,kst+lk) + fBox(inb,2)

!         Z-forcing
          for_zc_part(ist+li,jst+lj,kst+lk)=    &
          for_zc_part(ist+li,jst+lj,kst+lk) + fBox(inb,3)
    

          inb=inb+1
          end do
        end do
      end do
    
      end if    ! End istwoway
!     ----------------------------------------------------------------------

      else

      q1p(inp)=0.d0 ; q2p(inp)=0.d0 ; q3p(inp)=0.d0 
      
      qVal1(inp)=0.d0 ; qVal2(inp)=0.d0 ; qVal3(inp)=0.d0
      vort1(inp)=0.d0 ; vort2(inp)=0.d0 ; vort3(inp)=0.d0
      kalb1(inp)=0.d0 ; kalb2(inp)=0.d0 ; kalb3(inp)=0.d0
       
      end if

      end do

      call MpiBarrier

      call MpiAllSumReal1D(q1p,Npointpart)
      call MpiAllSumReal1D(q2p,Npointpart)
      call MpiAllSumReal1D(q3p,Npointpart)

      call MpiAllSumReal1D(qVal1,Npointpart)
      call MpiAllSumReal1D(qVal2,Npointpart)
      call MpiAllSumReal1D(qVal3,Npointpart)

      call MpiAllSumReal1D(vort1,Npointpart)
      call MpiAllSumReal1D(vort2,Npointpart)
      call MpiAllSumReal1D(vort3,Npointpart)

      call MpiAllSumReal1D(kalb1,Npointpart)
      call MpiAllSumReal1D(kalb2,Npointpart)
      call MpiAllSumReal1D(kalb3,Npointpart)
    
      call MpiAllSumReal2D(aap,Npointpart,3)
      call MpiAllSumReal2D(facc_for,Npointpart,3)
      call MpiAllSumReal2D(drag_for,Npointpart,3)
      call MpiAllSumReal2D(lift_for,Npointpart,3)
      call MpiAllSumReal1D(renp,Npointpart)

      return
      end
!----------------------------------------------------------------------------------
