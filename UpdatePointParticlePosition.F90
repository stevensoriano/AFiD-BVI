!------------------------------------------------------------
!     update particle position based on velocity
!     find out indices and distribute to slaves
!------------------------------------------------------------
      subroutine UpdatePointParticlePosition

      USE param
      USE pointparticle

      IMPLICIT NONE

      real      :: pos(3),qInt(7)
      real      :: vx_new,vy_new,vz_new
      real      :: vx_old,vy_old,vz_old,rbub
      integer   :: inp,merr

      rbub = 0.5*dbd*zlen

      do inp=1,Npointpart

       xpo(inp,1) = xp(inp,1)
       xpo(inp,2) = xp(inp,2)
       xpo(inp,3) = xp(inp,3) 

       vx_new=q1p (inp); vy_new=q2p (inp); vz_new=q3p (inp)
       vx_old=q1po(inp); vy_old=q2po(inp); vz_old=q3po(inp)      

       if(ONparticle.eq.0)then
        xp(inp,1) = xp(inp,1) + dt*vx_new
        xp(inp,2) = xp(inp,2) + dt*vy_new
        xp(inp,3) = xp(inp,3) + dt*vz_new
       else
        xp(inp,1) = xp(inp,1) + 0.5d0*dt*(3.d0*vx_new - vx_old)
        xp(inp,2) = xp(inp,2) + 0.5d0*dt*(3.d0*vy_new - vy_old)
        xp(inp,3) = xp(inp,3) + 0.5d0*dt*(3.d0*vz_new - vz_old)
       end if

!     =================================================================
!               Periodic boundary corrections for particles
!     =================================================================

!RO    Is this correct? doesn't seem like it to me. 

!      X-boundary check
       if(xp(inp,1).lt.xm(1))   xp(inp,1)=xm(n1m)-(xm(1)-xp(inp,1))
       if(xp(inp,1).gt.xm(n1m)) xp(inp,1)=xm(1)+(xp(inp,1)-xm(n1m))

!      Y-boundary check          
       if(xp(inp,2).lt.ym(1))   xp(inp,2)=ym(n2m)-(ym(1)-xp(inp,2))
       if(xp(inp,2).gt.ym(n2m)) xp(inp,2)=ym(1)+(xp(inp,2)-ym(n2m))

!      Z-boundary check
       if(xp(inp,3).lt.zm(1))   xp(inp,3)=zm(n3m)-(zm(1)-xp(inp,3))
       if(xp(inp,3).gt.zm(n3m)) xp(inp,3)=zm(1)+(xp(inp,3)-zm(n3m))

!     =================================================================
!                Update old variables by all processors
!     =================================================================

       q1po(inp)=q1p(inp)
       q2po(inp)=q2p(inp)
       q3po(inp)=q3p(inp)

       kalbo1(inp)=kalb1(inp)
       kalbo2(inp)=kalb2(inp)
       kalbo3(inp)=kalb3(inp)

       qValo1(inp)=qVal1(inp)
       qValo2(inp)=qVal2(inp)
       qValo3(inp)=qVal3(inp)

       vorto1(inp)=vort1(inp)
       vorto2(inp)=vort2(inp)
       vorto3(inp)=vort3(inp) 

      end do
     
      return
      end
