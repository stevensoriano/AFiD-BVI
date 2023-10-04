!------------------------------------------------------------------
!	calculate bubble indices
!------------------------------------------------------------------

      subroutine indpos(pos,inp,bubind)

      USE param
      IMPLICIT NONE
      
      real pos(3)
      integer bubind(6)
      integer i1,j1,k1,ist,jst,kst,inp

! CS  Just get i-j-k index (bubind) of particle (inp) using its position (pos) 

!     --------------------------------------------------------
!     X - indices
      i1=FLOOR(pos(1)*dx1) + 1

!     staggered
      if(pos(1).gt.xm(i1))ist=i1
      if(pos(1).le.xm(i1))ist=i1-1

      if(ist.eq.0)ist=n1m 

!     --------------------------------------------------------
!     Y - indices
      j1=FLOOR(pos(2)*dx2) + 1

!     staggered
      if(pos(2).gt.ym(j1))jst=j1
      if(pos(2).le.ym(j1))jst=j1-1

      if(jst.eq.0)jst=n2m

!     --------------------------------------------------------
!     Z - indices

      k1=FLOOR(pos(3)*dx3) + 1

!     staggered
      if(pos(3).gt.zm(k1))kst=k1
      if(pos(3).le.zm(k1))kst=k1-1

      if(kst.eq.0)kst=n3m

!     -----------------------------------------------------
      bubind(1)=i1  ; bubind(2)=j1  ; bubind(3)=k1
      bubind(4)=ist ; bubind(5)=jst ; bubind(6)=kst
!     -------------------------------------------------------
      if(i1.gt.n1m.or.j1.gt.n2m.or.k1.gt.n3m)then
      print*,'error in indic',i1,j1,k1
      print*,'posn indic err',pos(1),pos(2),pos(3)
      stop        
      end if

!     -------------------------------------------------------

      return
      end

!-----------------------------------------------------------------------
!     Tri-linear interpolation of fluid velocity and vorticity
!-----------------------------------------------------------------------

      subroutine interpol(pos,bubind,inp,qInt)

      USE param
      use local_aux
      USE pointparticle
      IMPLICIT NONE
      
      real      :: qBox(8,7),pos(3),posp(3),qInt(7)
      integer   :: inb,li,lj,lk,inp 
      integer   :: i,j,k,ist,jst,kst
      integer   :: bubind(6)

!     ------------------ Bubble indices -------------------------------
      i=bubind(1)   ; j=bubind(2)   ; k=bubind(3)
      ist=bubind(4) ; jst=bubind(5) ; kst=bubind(6)

!     ----- Percentage calculations for trilinear interpolation -------

      posp(1)=(pos(1)-xm(ist))/(xm(2)-xm(1)) !CS TODO: check definition
      posp(2)=(pos(2)-ym(jst))/(ym(2)-ym(1)) !CS TODO: check definition
      posp(3)=(pos(3)-zm(kst))/(zm(2)-zm(1)) !CS TODO: check definition

      posp(1)=1.d0-posp(1); posp(2)=1.d0-posp(2); posp(3)=1.d0-posp(3)
!     -----------------------------------------------------------------

      inb=1

      do lk=0,1
       do lj=0,1
        do li=0,1

         qBox(inb,1) = vxc(ist+li,jst+lj,kst+lk)
         qBox(inb,2) = vyc(ist+li,jst+lj,kst+lk)
         qBox(inb,3) = vzc(ist+li,jst+lj,kst+lk)

         qBox(inb,5) = vorx(ist+li,jst+lj,kst+lk)
         qBox(inb,6) = vory(ist+li,jst+lj,kst+lk)
         qBox(inb,7) = vorz(ist+li,jst+lj,kst+lk)

         inb=inb+1

         end do
        end do
      end do

!     X-velocity 
      qInt(1)=posp(1)*posp(2)*posp(3)*qBox(1,1)     &
                  + (1-posp(1))*posp(2)*posp(3)*qBox(2,1)     &
                  + posp(1)*(1-posp(2))*posp(3)*qBox(3,1)     &
                  + (1-posp(1))*(1-posp(2))*posp(3)*qBox(4,1)     &
                  + posp(1)*posp(2)*(1-posp(3))*qBox(5,1)     &
                  + (1-posp(1))*posp(2)*(1-posp(3))*qBox(6,1)     &
                  + posp(1)*(1-posp(2))*(1-posp(3))*qBox(7,1)     &
                  + (1-posp(1))*(1-posp(2))*(1-posp(3))*qBox(8,1)    

!     Y-velocity 
      qInt(2)=posp(1)*posp(2)*posp(3)*qBox(1,2)     &
                  + (1-posp(1))*posp(2)*posp(3)*qBox(2,2)     &
                  + posp(1)*(1-posp(2))*posp(3)*qBox(3,2)     &
                  + (1-posp(1))*(1-posp(2))*posp(3)*qBox(4,2)     &
                  + posp(1)*posp(2)*(1-posp(3))*qBox(5,2)     &
                  + (1-posp(1))*posp(2)*(1-posp(3))*qBox(6,2)     &
                  + posp(1)*(1-posp(2))*(1-posp(3))*qBox(7,2)     &
                  + (1-posp(1))*(1-posp(2))*(1-posp(3))*qBox(8,2)
    

!     Z-velocity 
      qInt(3)=posp(1)*posp(2)*posp(3)*qBox(1,3)  &
                  + (1-posp(1))*posp(2)*posp(3)*qBox(2,3)  &
                  + posp(1)*(1-posp(2))*posp(3)*qBox(3,3)  &
                  + (1-posp(1))*(1-posp(2))*posp(3)*qBox(4,3)  &
                  + posp(1)*posp(2)*(1-posp(3))*qBox(5,3)  &
                  + (1-posp(1))*posp(2)*(1-posp(3))*qBox(6,3)  &
                  + posp(1)*(1-posp(2))*(1-posp(3))*qBox(7,3)  &
                  + (1-posp(1))*(1-posp(2))*(1-posp(3))*qBox(8,3) 
 

!     X-vorticity
      qInt(5)=posp(1)*posp(2)*posp(3)*qBox(1,5)  &
                  + (1-posp(1))*posp(2)*posp(3)*qBox(2,5)  &
                  + posp(1)*(1-posp(2))*posp(3)*qBox(3,5)  &
                  + (1-posp(1))*(1-posp(2))*posp(3)*qBox(4,5)  &
                  + posp(1)*posp(2)*(1-posp(3))*qBox(5,5)  &
                  + (1-posp(1))*posp(2)*(1-posp(3))*qBox(6,5)  &
                  + posp(1)*(1-posp(2))*(1-posp(3))*qBox(7,5)  &
                  + (1-posp(1))*(1-posp(2))*(1-posp(3))*qBox(8,5)

!     Y-vorticity
      qInt(6)=posp(1)*posp(2)*posp(3)*qBox(1,6)  &
                  + (1-posp(1))*posp(2)*posp(3)*qBox(2,6)  &
                  + posp(1)*(1-posp(2))*posp(3)*qBox(3,6)  &
                  + (1-posp(1))*(1-posp(2))*posp(3)*qBox(4,6)  &
                  + posp(1)*posp(2)*(1-posp(3))*qBox(5,6)  &
                  + (1-posp(1))*posp(2)*(1-posp(3))*qBox(6,6)  &
                  + posp(1)*(1-posp(2))*(1-posp(3))*qBox(7,6)  &
                  + (1-posp(1))*(1-posp(2))*(1-posp(3))*qBox(8,6)

!     Z-vorticity
      qInt(7)=posp(1)*posp(2)*posp(3)*qBox(1,7)  &
                  + (1-posp(1))*posp(2)*posp(3)*qBox(2,7)  &
                  + posp(1)*(1-posp(2))*posp(3)*qBox(3,7)  &
                  + (1-posp(1))*(1-posp(2))*posp(3)*qBox(4,7)  &
                  + posp(1)*posp(2)*(1-posp(3))*qBox(5,7)  &
                  + (1-posp(1))*posp(2)*(1-posp(3))*qBox(6,7)  &
                  + posp(1)*(1-posp(2))*(1-posp(3))*qBox(7,7)  &
                  + (1-posp(1))*(1-posp(2))*(1-posp(3))*qBox(8,7)
     

      return
      end
        
!------------------------------------------------------------------
!     calculate values at the centre of cell volumes
!     used for interpolation
!------------------------------------------------------------------

      subroutine CentreVariables
      USE mpi_param, only: kstart, kend
      USE mpih
      USE local_aux
      USE local_arrays, only: vx, vy, vz
      USE param
      IMPLICIT NONE

      integer li,lj,lk
      integer lbk,ubk

      do lk=kstart-lvlhalo,kend+lvlhalo-1
       do lj=1,n2m
        do li=1,n1m

         vxc(li,lj,lk)=(vx(li,lj,lk)+vx(li+1,lj,lk))*0.5d0
         vyc(li,lj,lk)=(vy(li,lj,lk)+vy(li,lj+1,lk))*0.5d0
         vzc(li,lj,lk)=(vz(li,lj,lk)+vz(li,lj,lk+1))*0.5d0

        matderxc(li,lj,lk) = (matderx(li  ,lj  ,lk  ) &
                            + matderx(li+1,lj  ,lk  ))*0.5d0
        matderyc(li,lj,lk) = (matdery(li  ,lj  ,lk  ) &
                            + matdery(li  ,lj+1,lk  ))*0.5d0
        matderzc(li,lj,lk) = (matderz(li  ,lj  ,lk  ) &
                            + matderz(li  ,lj  ,lk+1))*0.5d0


        end do
       end do      
      end do


!-----------------------------------------------------------------

      return
      end

!----------------------------------------------------------------------------
!     Tri-linear interpolation of material derivative at particle position
!----------------------------------------------------------------------------

      subroutine interdmat(pos,bubind,inp,dmatInt)

      USE param
      USE local_aux, only: matderxc,matderyc,matderzc
      IMPLICIT NONE

      real      :: mBox(8,3),pos(3),posp(3),dmatInt(3)
      integer   :: inb,li,lj,lk,inp 
      integer   :: i,j,k,ist,jst,kst
      integer   :: bubind(6)

!     ------------------ Bubble indices -------------------------------
      i=bubind(1)   ; j=bubind(2)   ; k=bubind(3)
      ist=bubind(4) ; jst=bubind(5) ; kst=bubind(6)

!     ----- Percentage calculations for trilinear interpolation -------
      posp(1)=(pos(1)-xm(ist))/(xm(2)-xm(1)) !CS TODO: check definition
      posp(2)=(pos(2)-ym(jst))/(ym(2)-ym(1)) !CS TODO: check definition
      posp(3)=(pos(3)-zm(kst))/(zm(2)-zm(1)) !CS TODO: check definition

      posp(1)=1.d0-posp(1); posp(2)=1.d0-posp(2); posp(3)=1.d0-posp(3)
!     -----------------------------------------------------------------

      inb=1

      do lk=0,1
        do lj=0,1
         do li =0,1

         mBox(inb,1) = matderxc(ist+li,jst+lj,kst+lk)
         mBox(inb,2) = matderyc(ist+li,jst+lj,kst+lk)
         mBox(inb,3) = matderzc(ist+li,jst+lj,kst+lk)

         inb=inb+1

        end do
       end do
      end do

!     X-material derivative
      dmatInt(1)=posp(1)*posp(2)*posp(3)*mBox(1,1)  &
                  +(1-posp(1))*posp(2)*posp(3)*mBox(2,1)   &
                  + posp(1)*(1-posp(2))*posp(3)*mBox(3,1)   &
                  + (1-posp(1))*(1-posp(2))*posp(3)*mBox(4,1)   &
                  + posp(1)*posp(2)*(1-posp(3))*mBox(5,1)   &
                  + (1-posp(1))*posp(2)*(1-posp(3))*mBox(6,1)   &
                  + posp(1)*(1-posp(2))*(1-posp(3))*mBox(7,1)   &
                  + (1-posp(1))*(1-posp(2))*(1-posp(3))*mBox(8,1)

!     Y-material derivative
      dmatInt(2)=posp(1)*posp(2)*posp(3)*mBox(1,2)   &
                  + (1-posp(1))*posp(2)*posp(3)*mBox(2,2)   &
                  + posp(1)*(1-posp(2))*posp(3)*mBox(3,2)   &
                  + (1-posp(1))*(1-posp(2))*posp(3)*mBox(4,2)   &
                  + posp(1)*posp(2)*(1-posp(3))*mBox(5,2)   &
                  + (1-posp(1))*posp(2)*(1-posp(3))*mBox(6,2)   &
                  + posp(1)*(1-posp(2))*(1-posp(3))*mBox(7,2)   &
                  + (1-posp(1))*(1-posp(2))*(1-posp(3))*mBox(8,2)

!     Z-material derivative
      dmatInt(3)=posp(1)*posp(2)*posp(3)*mBox(1,3)   &
                  + (1-posp(1))*posp(2)*posp(3)*mBox(2,3)   &
                  + posp(1)*(1-posp(2))*posp(3)*mBox(3,3)  &
                  + (1-posp(1))*(1-posp(2))*posp(3)*mBox(4,3)   &
                  + posp(1)*posp(2)*(1-posp(3))*mBox(5,3)   &
                  + (1-posp(1))*posp(2)*(1-posp(3))*mBox(6,3)   &
                  + posp(1)*(1-posp(2))*(1-posp(3))*mBox(7,3)   &
                  + (1-posp(1))*(1-posp(2))*(1-posp(3))*mBox(8,3)

      
      return
      end

!	---------------------------------------------------------------

      subroutine storeold
      USE param
      USE mpih
      USE mpi_param, only: kstart, kend
      USE local_arrays, only: vx, vy, vz
      USE local_aux
      USE mls_local
      use pointparticle
      IMPLICIT NONE


      vxo(1:n1,1:n2,kstart:kend) = &
                vx(1:n1,1:n2,kstart:kend)
      vyo(1:n1,1:n2,kstart:kend) =  &
                vy(1:n1,1:n2,kstart:kend)
      vzo(1:n1,1:n2,kstart:kend) = &
                vz(1:n1,1:n2,kstart:kend)

      for_xc_part(1:n1m,1:n2m,kstart-lvlhalo:kend)=0.d0
      for_yc_part(1:n1m,1:n2m,kstart-lvlhalo:kend)=0.d0
      for_zc_part(1:n1m,1:n2m,kstart-lvlhalo:kend)=0.d0

      for_x_part(1:n1m,1:n2m,kstart:kend)=0.d0
      for_y_part(1:n1m,1:n2m,kstart:kend)=0.d0
      for_z_part(1:n1m,1:n2m,kstart:kend)=0.d0

      matderx(1:n1,1:n2,kstart:kend)=0.d0 
      matdery(1:n1,1:n2,kstart:kend)=0.d0 
      matderz(1:n1,1:n2,kstart:kend)=0.d0 

      matderxc(1:n1m,1:n2m,kstart-lvlhalo:kend)=0.d0 
      matderyc(1:n1m,1:n2m,kstart-lvlhalo:kend)=0.d0 
      matderzc(1:n1m,1:n2m,kstart-lvlhalo:kend)=0.d0 



      return
      end

!	---------------------------------------------------------------

      subroutine AllocatePointPArrays
      USE param
      USE mpih
      USE mpi_param, only: kstart, kend
      USE local_arrays, only: vx, vy, vz
      USE local_aux
      USE mls_local
      use pointparticle
      use AuxiliaryRoutines
      IMPLICIT NONE

      call AllocateReal3DArray(vxo,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(vyo,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(vzo,1,n1,1,n2,kstart,kend)

      call AllocateReal3DArray(vxc,1,n1,1,n2,  &
               kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(vyc,1,n1,1,n2,  &
               kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(vzc,1,n1,1,n2,  &
               kstart-lvlhalo,kend+lvlhalo)

      call AllocateReal3DArray(matderx,1,n1,1,n2,  &
               kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(matdery,1,n1,1,n2,  &
               kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(matderz,1,n1,1,n2,  &
               kstart-lvlhalo,kend+lvlhalo)

      call AllocateReal3DArray(matderx,1,n1,1,n2,  &
               kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(matdery,1,n1,1,n2,  &
               kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(matderz,1,n1,1,n2,  &
               kstart-lvlhalo,kend+lvlhalo)

      call AllocateReal3DArray(matderxc,1,n1,1,n2,  &
               kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(matderyc,1,n1,1,n2,  &
               kstart-lvlhalo,kend+lvlhalo)
      call AllocateReal3DArray(matderzc,1,n1,1,n2,  &
               kstart-lvlhalo,kend+lvlhalo)

      call AllocateReal3DArray(for_xc_part,1,n1,1,n2,kstart-lvlhalo,kend)
      call AllocateReal3DArray(for_yc_part,1,n1,1,n2,kstart-lvlhalo,kend)
      call AllocateReal3DArray(for_zc_part,1,n1,1,n2,kstart-lvlhalo,kend)

      call AllocateReal3DArray(for_x_part,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(for_y_part,1,n1,1,n2,kstart,kend)
      call AllocateReal3DArray(for_z_part,1,n1,1,n2,kstart,kend)

!     Allocate particle parameters

      !-------- Now allocate particle parameters ---------
      call AllocateReal2DArray(xp,1,Npointpart,1,3)
      call AllocateReal2DArray(xpo,1,Npointpart,1,3)

      call AllocateReal1DArray(q1p,1,Npointpart)
      call AllocateReal1DArray(q2p,1,Npointpart)
      call AllocateReal1DArray(q3p,1,Npointpart)

      call AllocateReal1DArray(q1po,1,Npointpart)
      call AllocateReal1DArray(q2po,1,Npointpart)
      call AllocateReal1DArray(q3po,1,Npointpart)

      call AllocateReal2DArray(aap,1,Npointpart,1,3)
      call AllocateReal2DArray(facc_for,1,Npointpart,1,3)
      call AllocateReal2DArray(lift_for,1,Npointpart,1,3)
      call AllocateReal2DArray(drag_for,1,Npointpart,1,3)

      call AllocateReal1DArray(renp,1,Npointpart)

      call AllocateReal1DArray(qVal1,1,Npointpart)
      call AllocateReal1DArray(qVal2,1,Npointpart)
      call AllocateReal1DArray(qVal3,1,Npointpart)

      call AllocateReal1DArray(qValo1,1,Npointpart)
      call AllocateReal1DArray(qValo2,1,Npointpart)
      call AllocateReal1DArray(qValo3,1,Npointpart)

      call AllocateReal1DArray(kalb1,1,Npointpart)
      call AllocateReal1DArray(kalb2,1,Npointpart)
      call AllocateReal1DArray(kalb3,1,Npointpart)

      call AllocateReal1DArray(kalbo1,1,Npointpart)
      call AllocateReal1DArray(kalbo2,1,Npointpart)
      call AllocateReal1DArray(kalbo3,1,Npointpart)

      call AllocateReal1DArray(vort1,1,Npointpart)
      call AllocateReal1DArray(vort2,1,Npointpart)
      call AllocateReal1DArray(vort3,1,Npointpart)

      call AllocateReal1DArray(vorto1,1,Npointpart)
      call AllocateReal1DArray(vorto2,1,Npointpart)
      call AllocateReal1DArray(vorto3,1,Npointpart)

      call AllocateReal1DArray(dbd12,1,Npointpart)
      call AllocateReal1DArray(stokes,1,Npointpart)
      call AllocateReal1DArray(gammap,1,Npointpart)

      return
      end
