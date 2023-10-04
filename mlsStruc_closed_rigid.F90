!-----------------------------------------------------------------------
!     read in position of a individual trial marker for MLS and compute 
!     compute support domain, shape function and interpolate
!     
!     compute internal forces on the nodes 
!------------------------------------------------------------------------

      subroutine mlsStruc_closed_rigid
      USE mpih
      USE param
      USE mls_param
      USE local_arrays, only: vx,vy,vz
      USE mpi_param, only: kstart, kend
      USE mls_local, only: for_xc, for_yc, for_zc
      USE local_arrays, only: pr
      IMPLICIT NONE

      integer inp,seed,merr,ntr
      real pos_MLS(3)
      integer pind_i(6),pind_o(6)

      integer inw,i,j,k,cmp
      real norp(3),el_c(3),el_cx(3),Hbox(3)
      real wbet,dismaxpro,elmag,norpd,normd,epsw
      real pinvA(4,4),invA(4,4),B(4,nel)
      real Bx(4,nel),By(4,nel),Bz(4,nel)
      real pxk(4,1),pxkx(4,1),pxky(4,1),pxkz(4,1)
      real ui(nel,3),prp(nel),Wt(3,nel),Wtx(nel)
      real ptx(1,4),Bu(4,1),ABu(4,1),pABu(1,1)
      real ptxA(1,4),ptxAB(1,nel),ptxABu(1,1)
      real pxx(4,1),pxy(4,1),pxz(4,1)
      real dWt(3,nel),dWtx(nel),dWty(nel),dWtz(nel)     
      real piAx1(4,4),piAx2(4,4),pinvAx(4,4)
      real piAy1(4,4),piAy2(4,4),pinvAy(4,4)
      real piAz1(4,4),piAz2(4,4),pinvAz(4,4)
      real Gmat(4,1),Gmatx(4,1),Gmaty(4,1),Gmatz(4,1)
      real Gmatx1(4,1),Gmatx2(4,1),Gmatx3(4,1),Gmaty1(4,1),Gmaty2(4,1)
      real Gmaty3(4,1),Gmatz1(4,1),Gmatz2(4,1),Gmatz3(4,1)
      real Gtxb(1,nel),Gtbx(1,nel),PhiTx(1,nel)
      real Gtyb(1,nel),Gtby(1,nel),PhiTy(1,nel)
      real Gtzb(1,nel),Gtbz(1,nel),PhiTz(1,nel)
      
      real pr_p,press,dotp

      real sr_xx,sr_yy,sr_zz,sr_xy,sr_yz,sr_zx,sr_yx,sr_zy,sr_xz
      real,dimension(1,1)::dUxx,dUxy,dUxz,dUyx,dUyy,dUyz,dUzx,dUzy,dUzz
      real tau_f(3),pos_pro(3),fncy(3),fnca(3)
      real xyzp1,xyzp2,xyzp3,rtzp1,rtzp2,rtzp3

      integer::ind_pal
      integer::f1,f2,v1,v2,v3,v4,iv,ie
      real,dimension(3)::tvec1,tvec2,tvec3,tvec4,csi,zet,Vij
      real,dimension(3)::a32,a13,a34,a21,a42,a23,a31,a24,csi1,zet1,tcv
      real::modcsi,modzet,betab,betav,alphat,alphal,b11,b12,b22,tdum,VV
      real,dimension(3)::gravc,frinv

      integer li,lj,lk,wrongn
!     --------------------------------------------------------
      fpxyz(:,:,:)=0.0
      press_face=0.0d0
      vforc_face=0.0d0
      sca_nod1=0.0d0
      sca_nod2=0.0d0
      sca_nod3=0.0d0


      do inp=1,Nparticle

       wrongn=0

!
       do ntr=1,maxnf

!     ------checking for normal direction (not needed)------
!       dotp = (tri_bar(1,ntr,inp)-0.5)*tri_nor(1,ntr,inp) +  &
!              (tri_bar(2,ntr,inp)-0.5)*tri_nor(2,ntr,inp) +  &
!              (tri_bar(3,ntr,inp)-0.5)*tri_nor(3,ntr,inp)  
!     
!      if (dotp.lt.0.) then
!!       wrongn=wrongn+1
!       tri_nor(1:3,ntr,inp)=tri_nor(1:3,ntr,inp)*(-1.)
!      end if
!     ------------------------------------------------------

       if(dum_for(ntr,inp).eq.0)then

!     Initialise positions
        pos_MLS(1)=tri_bar(1,ntr,inp)
        pos_MLS(2)=tri_bar(2,ntr,inp)
        pos_MLS(3)=tri_bar(3,ntr,inp)
 

!     Collision force - allocate and compute

!     Compute positions of probes from centroids
      xyzp1=tri_bar(1,ntr,inp)+tri_nor(1,ntr,inp)*Hboxx
      xyzp2=tri_bar(2,ntr,inp)+tri_nor(2,ntr,inp)*Hboxx
      xyzp3=tri_bar(3,ntr,inp)+tri_nor(3,ntr,inp)*Hboxx


!     Support domain for probe
      pos_pro(1)=xyzp1;pos_pro(2)=xyzp2;pos_pro(3)=xyzp3

      call partindicesMLS(pos_pro,pind_i,pind_o,dismaxpro)

!     initialise pre-factor matrix

      ptx(1,1)=1.d0;ptx(1,2)=pos_pro(1);ptx(1,3)=pos_pro(2)
      ptx(1,4)=pos_pro(3)


      pxx(1,1)=0.0d0;pxx(2,1)=1.0d0;pxx(3,1)=0.0d0;pxx(4,1)=0.0d0
      pxy(1,1)=0.0d0;pxy(2,1)=0.0d0;pxy(3,1)=1.0d0;pxy(4,1)=0.0d0
      pxz(1,1)=0.0d0;pxz(2,1)=0.0d0;pxz(3,1)=0.0d0;pxz(4,1)=1.0d0

!     --------------------------------------------------------

      Hbox(1)=dx1
      Hbox(2)=dx2
      Hbox(3)=dx3

!     --------------------------------------------------------

      inw=1
      epsw = 1e-10
      dismaxpro = dismaxpro+epsw

!     spline weights for all three components and initialise ui matrix
      ui(:,:) = 0.0d0

      do k=pind_i(3),pind_o(3)
       do j=pind_i(2),pind_o(2)
        do i=pind_i(1),pind_o(1)
 
          el_cx(1)=xm(i)-pos_pro(1)
          el_cx(2)=ym(j)-pos_pro(2)
          el_cx(3)=zm(k)-pos_pro(3)

          elmag=sqrt(el_cx(1)**2+el_cx(2)**2+el_cx(3)**2)

          do cmp=1,3
            norp(cmp) = abs(el_cx(cmp))*Hbox(cmp)/wscl
        
!     Computing weight functions

!     ----------------EXPONENTIAL SPLINES------------------
            if(wexp.eq.1)then 
              if(norp(cmp).le.1.0d0)then
                Wt(cmp,inw)=exp(-(norp(cmp)/wcon)**2)
                dWt(cmp,inw)=(-2.d0/(wcon**2))*norp(cmp)* &
                          exp(-(norp(cmp)/wcon)**2) &
                         *(el_cx(cmp)/abs(el_cx(cmp)))*Hbox(cmp)/wscl
              else
                Wt(cmp,inw)=0.d0
                dWt(cmp,inw)=0.d0
              end if
            end if
!     ----------------CUBIC SPLINES---------------------
            if(wcub.eq.1)then
              if(norp(cmp).le.0.5d0)then
                Wt(cmp,inw)=(2.d0/3.d0)-4.d0*(norp(cmp)**2)+ &
                                    4.d0*(norp(cmp)**3)
                dWt(cmp,inw)=(-8.d0*(norp(cmp))+12.d0*(norp(cmp)**2)) &
                        *(el_cx(cmp)/abs(el_cx(cmp)))*Hbox(cmp)/wscl
                elseif(norp(cmp).le.1.d0)then
                Wt(cmp,inw)=(4.d0/3.d0)*(1.d0-norp(cmp)**3)- &
                              4.d0*(norp(cmp)-norp(cmp)**2)
                dWt(cmp,inw)=(-4.d0+8.d0*norp(cmp)-4.d0*(norp(cmp)**2)) &
                        *(el_cx(cmp)/abs(el_cx(cmp)))*Hbox(cmp)/wscl
               else
                Wt(cmp,inw)=0.d0
                dWt(cmp,inw)=0.d0
             end if
            end if
           
          end do

          Wtx(inw) = Wt(1,inw)*Wt(2,inw)*Wt(3,inw)
          dWtx(inw) = dWt(1,inw)*Wt(2,inw)*Wt(3,inw)
          dWty(inw) = Wt(1,inw)*dWt(2,inw)*Wt(3,inw)
          dWtz(inw) = Wt(1,inw)*Wt(2,inw)*dWt(3,inw)
!     ------------------------------------------------

!     Allocating variables for support domain
          prp(inw) = pr(i,j,k) 

          ui(inw,1) = 0.5*(vx(i,j,k)+vx(i+1,j,k))
          ui(inw,2) = 0.5*(vy(i,j,k)+vy(i,j+1,k))
          ui(inw,3) = 0.5*(vz(i,j,k)+vz(i,j,k+1))
         
          inw = inw + 1
         end do
        end do
       end do

!     --------------------------------------------------------
      
      B(:,:)=0.0d0;Bx(:,:)=0.0d0;By(:,:)=0.0d0;Bz(:,:)=0.0d0
      pinvA(:,:)=0.0d0; invA(:,:)=0.0d0
      pinvAx(:,:)=0.0d0;pinvAy(:,:)=0.0d0;pinvAz(:,:)=0.0d0

      inw = 1

!     pre-inverse matrix A(x) and B(x)
      do k=pind_i(3),pind_o(3)      
       do j=pind_i(2),pind_o(2)
        do i=pind_i(1),pind_o(1)

        pxk(1,1)=1.0d0;pxk(2,1)=xm(i)
        pxk(3,1)=ym(j);pxk(4,1)=zm(k)

        pxkx(1,1)=0.0d0;pxkx(2,1)=1.0d0;pxkx(3,1)=0.0d0
        pxkx(4,1)=0.0d0
        pxky(1,1)=0.0d0;pxky(2,1)=0.0d0;pxky(3,1)=1.0d0
        pxky(4,1)=0.0d0
        pxkz(1,1)=0.0d0;pxkz(2,1)=0.0d0;pxkz(3,1)=0.0d0
        pxkz(4,1)=1.0d0

        call DGEMM('N','T',4,4,1,Wtx(inw),pxk,4,pxk,4, &
                                 1.0d0,pinvA,4)

!     ------ derivative of A matrix - req for Gam derivative ------------
        call DGEMM('N','T',4,4,1,dWtx(inw),pxk,4,pxk,4,1.0d0,pinvAx,4)

        call DGEMM('N','T',4,4,1,dWty(inw),pxk,4,pxk,4,1.0d0,pinvAy,4)

        call DGEMM('N','T',4,4,1,dWtz(inw),pxk,4,pxk,4,1.0d0,pinvAz,4)
!     -------------------------------------------------------------------
        B(:,inw)=Wtx(inw)*pxk(:,1)
!     -------derivative of B matrix - req for shape derivative---------
        Bx(:,inw) = dWtx(inw)*pxk(:,1) 
        By(:,inw) = dWty(inw)*pxk(:,1) 
        Bz(:,inw) = dWtz(inw)*pxk(:,1) 
!     -----------------------------------------------------------------
        inw = inw + 1
  
        end do
       end do
      end do

!     calling routine to compute inverse
      call inverseLU(pinvA,invA)
            
!     Compute Gamma=inva*p for derivatives of shape functions
      call DGEMM('N','T',4,1,4,1.0d0,invA,4,ptx,1,0.0d0,Gmat,4)
 
      call DGEMM('N','N',4,1,4,1.0d0,invA,4,pxx,4,0.0d0,Gmatx1,4)  
      call DGEMM('N','N',4,1,4,1.0d0,pinvAx,4,Gmat,4,0.0d0,Gmatx2,4)  
      call DGEMM('N','N',4,1,4,1.0d0,invA,4,Gmatx2,4,0.0d0,Gmatx3,4)  
      Gmatx = Gmatx1-Gmatx3  
      call DGEMM('N','N',4,1,4,1.0d0,invA,4,pxy,4,0.0d0,Gmaty1,4)  
      call DGEMM('N','N',4,1,4,1.0d0,pinvAy,4,Gmat,4,0.0d0,Gmaty2,4)  
      call DGEMM('N','N',4,1,4,1.0d0,invA,4,Gmaty2,4,0.0d0,Gmaty3,4)  
      Gmaty = Gmaty1-Gmaty3 
      call DGEMM('N','N',4,1,4,1.0d0,invA,4,pxz,4,0.0d0,Gmatz1,4)  
      call DGEMM('N','N',4,1,4,1.0d0,pinvAz,4,Gmat,4,0.0d0,Gmatz2,4)  
      call DGEMM('N','N',4,1,4,1.0d0,invA,4,Gmatz2,4,0.0d0,Gmatz3,4)  
      Gmatz = Gmatz1-Gmatz3 

!     Compute shape function derivatives
      call DGEMM('T','N',1,nel,4,1.0d0,Gmatx,4,B,4,0.0d0,Gtxb,1)
      call DGEMM('T','N',1,nel,4,1.0d0,Gmat,4,Bx,4,0.0d0,Gtbx,1)
      PhiTx = Gtxb + Gtbx
      call DGEMM('T','N',1,nel,4,1.0d0,Gmaty,4,B,4,0.0d0,Gtyb,1)
      call DGEMM('T','N',1,nel,4,1.0d0,Gmat,4,By,4,0.0d0,Gtby,1)
      PhiTy = Gtyb + Gtby
      call DGEMM('T','N',1,nel,4,1.0d0,Gmatz,4,B,4,0.0d0,Gtzb,1)
      call DGEMM('T','N',1,nel,4,1.0d0,Gmat,4,Bz,4,0.0d0,Gtbz,1)
      PhiTz = Gtzb + Gtbz

      if(wcheck.eq.1)then
        if(abs(sum(PhiTx)).ge.0.10.or.abs(sum(PhiTy)).ge.0.10.or. &
           abs(sum(PhiTz)).ge.0.10)then
        print*,'Shape deriv. sum',sum(PhiTx),sum(PhiTy),sum(PhiTz)
        call MPI_ABORT(MPI_COMM_WORLD,ierr)
        end if
      end if

!     ------Velocity gradient interpolation-------------------

      call DGEMM('N','N',1,1,nel,1.0d0,PhiTx,1,ui(:,1),nel,0.0d0, &
                                                      dUxx,1) 
      call DGEMM('N','N',1,1,nel,1.0d0,PhiTy,1,ui(:,1),nel,0.0d0, &
                                                      dUxy,1) 
      call DGEMM('N','N',1,1,nel,1.0d0,PhiTz,1,ui(:,1),nel,0.0d0, &
                                                      dUxz,1) 
!
      call DGEMM('N','N',1,1,nel,1.0d0,PhiTx,1,ui(:,2),nel,0.0d0, &
                                                      dUyx,1) 
      call DGEMM('N','N',1,1,nel,1.0d0,PhiTy,1,ui(:,2),nel,0.0d0, &
                                                      dUyy,1) 
      call DGEMM('N','N',1,1,nel,1.0d0,PhiTz,1,ui(:,2),nel,0.0d0, &
                                                      dUyz,1) 
!
      call DGEMM('N','N',1,1,nel,1.0d0,PhiTx,1,ui(:,3),nel,0.0d0, &
                                                      dUzx,1) 
      call DGEMM('N','N',1,1,nel,1.0d0,PhiTy,1,ui(:,3),nel,0.0d0, &
                                                      dUzy,1) 
      call DGEMM('N','N',1,1,nel,1.0d0,PhiTz,1,ui(:,3),nel,0.0d0, &
                                                      dUzz,1) 
      
!     ------------------------------------------------------
!     matrix multiplications for final interpolation
!     DGEMM(transA,transB,m,n,k,alpha,A,LDA,b,LDB,beta,c,LDC)
!     C = alpha * A * B + beta * C

!     ---------------Shape function calculation---------------
      call DGEMM('N','N',1,4,4,1.0d0,ptx,1,invA,4,0.0d0,ptxA,1) 
      call DGEMM('N','N',1,nel,4,1.0d0,ptxA,1,B,4,0.0d0,ptxAB,1) 

!     --------------Pressure Interpolation at probe from \phi---------
      call DGEMM('N','N',1,1,nel,1.0d0,ptxAB,1,prp(:),nel,0.0d0, &
                                                      ptxABu,1) 
      
      if(wcheck.eq.1)then
        if(sum(ptxAB).le.0.95.or.sum(ptxAB).gt.1.05)then
        print*,'Shape function sum mlsStruc',sum(ptxAB)
        end if
        if(isnan(sum(ptxAB))) then
        write(*,*)'NaN Shape function mlsStruc',pos_pro(1:3)
        end if
      end if

      pr_p=ptxABu(1,1)

!     -----------------------------------------------------------------
!     Face normals 
      fnca(1:3) = tri_nor(1:3,ntr,inp)
            
!     Compute pressure at marker from probe
      press=pr_p+Hboxx*(acc_tri(1,ntr,inp)*fnca(1)+ &
              acc_tri(2,ntr,inp)*fnca(2)+acc_tri(3,ntr,inp)*fnca(3))

!     Compute shear stress from vel.gradients

      sr_xx = dUxx(1,1); sr_yy = dUyy(1,1); sr_zz = dUzz(1,1)
      sr_xy = 0.50*(dUxy(1,1)+dUyx(1,1))
      sr_yz = 0.50*(dUyz(1,1)+dUzy(1,1))
      sr_zx = 0.50*(dUzx(1,1)+dUxz(1,1))

      sr_yx=sr_xy;sr_zy=sr_yz;sr_xz=sr_zx

      tau_f(1)=(2.0/ren)*(sr_xx*fnca(1)+sr_xy*fnca(2)+sr_xz*fnca(3))
      tau_f(2)=(2.0/ren)*(sr_yx*fnca(1)+sr_yy*fnca(2)+sr_yz*fnca(3))
      tau_f(3)=(2.0/ren)*(sr_zx*fnca(1)+sr_zy*fnca(2)+sr_zz*fnca(3))

!     -----------------------------------------------------------------

      press_face(ntr,inp) = -press*fnca(1)*sur(ntr,inp)
      vforc_face(ntr,inp) = tau_f(1)*sur(ntr,inp)

      sca_nod1(ntr,inp) = -press*fnca(1)*sur(ntr,inp)+tau_f(1)*sur(ntr,inp)
      sca_nod2(ntr,inp) = -press*fnca(2)*sur(ntr,inp)+tau_f(2)*sur(ntr,inp) 
      sca_nod3(ntr,inp) = -press*fnca(3)*sur(ntr,inp)+tau_f(3)*sur(ntr,inp) 

       else
      
      press_face(ntr,inp) = 0.0d0
      vforc_face(ntr,inp) = 0.0d0

      sca_nod1(ntr,inp) = 0.0d0
      sca_nod2(ntr,inp) = 0.0d0
      sca_nod3(ntr,inp) = 0.0d0
      
       end if 


       end do

      end do
 
      return
      end
