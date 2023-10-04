!-----------------------------------------------------------------------
!     read in position of a individual trial marker for MLS and compute 
!     compute support domain, shape function and interpolate
!     
!     compute internal forces on the nodes 
!------------------------------------------------------------------------

      subroutine mlsStruc_open_rigid
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
      integer pind_i_i(6),pind_o_i(6),pind_i_o(6),pind_o_o(6)

      integer inw,i,j,k,cmp
      real norp(3),el_c(3),el_cx(3),Hbox(3)
      real wbet,dismaxpro_i,dismaxpro_o,elmag,norpd,normd,epsw
      real pinvA(4,4),invA(4,4),B(4,nel)
      real Bx(4,nel),By(4,nel),Bz(4,nel)
      real pxk(4,1),pxkx(4,1),pxky(4,1),pxkz(4,1)
      real ui_o(nel,3),ui_i(nel,3),prp_o(nel),prp_i(nel),Wt(3,nel)
      real Wtx_i(nel),Wtx_o(nel)
      real ptx_o(1,4),ptx_i(1,4),Bu(4,1),ABu(4,1),pABu(1,1)
      real ptxA(1,4),ptxAB(1,nel),ptxABu(1,1)
      real pxx(4,1),pxy(4,1),pxz(4,1)
      real dWt(3,nel)
      real dWtx_o(nel),dWty_o(nel),dWtz_o(nel)     
      real dWtx_i(nel),dWty_i(nel),dWtz_i(nel)     
      real piAx1(4,4),piAx2(4,4),pinvAx(4,4)
      real piAy1(4,4),piAy2(4,4),pinvAy(4,4)
      real piAz1(4,4),piAz2(4,4),pinvAz(4,4)
      real Gmat(4,1),Gmatx(4,1),Gmaty(4,1),Gmatz(4,1)
      real Gmatx1(4,1),Gmatx2(4,1),Gmatx3(4,1),Gmaty1(4,1),Gmaty2(4,1)
      real Gmaty3(4,1),Gmatz1(4,1),Gmatz2(4,1),Gmatz3(4,1)
      real Gtxb(1,nel),Gtbx(1,nel),PhiTx(1,nel)
      real Gtyb(1,nel),Gtby(1,nel),PhiTy(1,nel)
      real Gtzb(1,nel),Gtbz(1,nel),PhiTz(1,nel)
      
      real pr_p_i,pr_p_o,press_i,press_o

      real sr_xx,sr_yy,sr_zz,sr_xy,sr_yz,sr_zx,sr_yx,sr_zy,sr_xz
      real,dimension(1,1)::dUxx,dUxy,dUxz,dUyx,dUyy,dUyz,dUzx,dUzy,dUzz
      real tau_f_i(3),tau_f_o(3),pos_pro_i(3),pos_pro_o(3),fnca(3)
      real xyzp1_o,xyzp2_o,xyzp3_o,xyzp1_i,xyzp2_i,xyzp3_i

      integer::ind_pal
      integer::f1,f2,v1,v2,v3,v4,iv,ie
      real,dimension(3)::tvec1,tvec2,tvec3,tvec4,csi,zet,Vij
      real,dimension(3)::a32,a13,a34,a21,a42,a23,a31,a24,csi1,zet1,tcv
      real::modcsi,modzet,betab,betav,alphat,alphal,b11,b12,b22,tdum,VV
      real,dimension(3)::gravc
!     --------------------------------------------------------
      fpxyz(:,:,:)=0.0
      press_face = 0.0
      gravc(:)=0.0
      gravc(1)=-9.8

      do inp=1,Nparticle

       do ntr=1,maxnf
     
       if(dum_for(ntr,inp).eq.0)then

        pos_MLS(1)=tri_bar(1,ntr,inp)
        pos_MLS(2)=tri_bar(2,ntr,inp)
        pos_MLS(3)=tri_bar(3,ntr,inp)

!     Compute positions of probes from centroids
      xyzp1_o=tri_bar(1,ntr,inp)+tri_nor(1,ntr,inp)*Hboxx
      xyzp2_o=tri_bar(2,ntr,inp)+tri_nor(2,ntr,inp)*Hboxx
      xyzp3_o=tri_bar(3,ntr,inp)+tri_nor(3,ntr,inp)*Hboxx

      xyzp1_i=tri_bar(1,ntr,inp)-tri_nor(1,ntr,inp)*Hboxx
      xyzp2_i=tri_bar(2,ntr,inp)-tri_nor(2,ntr,inp)*Hboxx
      xyzp3_i=tri_bar(3,ntr,inp)-tri_nor(3,ntr,inp)*Hboxx

!     Support domain for probe
      pos_pro_o(1)=xyzp1_o;pos_pro_o(2)=xyzp2_o;pos_pro_o(3)=xyzp3_o
      pos_pro_i(1)=xyzp1_i;pos_pro_i(2)=xyzp2_i;pos_pro_i(3)=xyzp3_i

      call partindicesMLS(pos_pro_o,pind_i_o,pind_o_o,dismaxpro_o)
      call partindicesMLS(pos_pro_i,pind_i_i,pind_o_i,dismaxpro_i)

!     initialise pre-factor matrix

      ptx_o(1,1)=1.d0;ptx_o(1,2)=pos_pro_o(1);ptx_o(1,3)=pos_pro_o(2)
      ptx_o(1,4)=pos_pro_o(3)
      ptx_i(1,1)=1.d0;ptx_i(1,2)=pos_pro_i(1);ptx_i(1,3)=pos_pro_i(2)
      ptx_i(1,4)=pos_pro_i(3)

      pxx(1,1)=0.0d0;pxx(2,1)=1.0d0;pxx(3,1)=0.0d0;pxx(4,1)=0.0d0
      pxy(1,1)=0.0d0;pxy(2,1)=0.0d0;pxy(3,1)=1.0d0;pxy(4,1)=0.0d0
      pxz(1,1)=0.0d0;pxz(2,1)=0.0d0;pxz(3,1)=0.0d0;pxz(4,1)=1.0d0

!     --------------------------------------------------------

      Hbox(1)=dx1
      Hbox(2)=dx2
      Hbox(3)=dx3

!     --------------------------------------------------------
      
      epsw = 1e-10
      dismaxpro_o = dismaxpro_o+epsw
      dismaxpro_i = dismaxpro_i+epsw

!     Outer probe calculation
!     spline weights for all three components and initialise ui matrix
      inw=1
      ui_o(:,:) = 0.0d0

      do k=pind_i_o(3),pind_o_o(3)
       do j=pind_i_o(2),pind_o_o(2)
        do i=pind_i_o(1),pind_o_o(1)
 
          el_cx(1)=xm(i)-pos_pro_o(1)
          el_cx(2)=ym(j)-pos_pro_o(2)
          el_cx(3)=zm(k)-pos_pro_o(3)

          elmag=sqrt(el_cx(1)**2+el_cx(2)**2+el_cx(3)**2)

          do cmp=1,3
            norp(cmp) = abs(el_cx(cmp))*Hbox(cmp)/wscl
        
!     Computing weight functions

!     ----------------EXPONENTIAL SPLINES------------------
            if(wexp.eq.1)then 
              if(norp(cmp).le.1.0d0)then
                Wt(cmp,inw)=exp(-(norp(cmp)/wcon)**2)
                dWt(cmp,inw)=(-2.d0/(wcon**2))*norp(cmp)*
     %                    exp(-(norp(cmp)/wcon)**2)
     %                   *(el_cx(cmp)/abs(el_cx(cmp)))*Hbox(cmp)/wscl
              else
                Wt(cmp,inw)=0.d0
                dWt(cmp,inw)=0.d0
              end if
            end if
!     ----------------CUBIC SPLINES---------------------
            if(wcub.eq.1)then
              if(norp(cmp).le.0.5d0)then
                Wt(cmp,inw)=(2.d0/3.d0)-4.d0*(norp(cmp)**2)+
     %                              4.d0*(norp(cmp)**3)
                dWt(cmp,inw)=(-8.d0*(norp(cmp))+12.d0*(norp(cmp)**2))
     %                  *(el_cx(cmp)/abs(el_cx(cmp)))*Hbox(cmp)/wscl
                elseif(norp(cmp).le.1.d0)then
                Wt(cmp,inw)=(4.d0/3.d0)*(1.d0-norp(cmp)**3)-
     %                        4.d0*(norp(cmp)-norp(cmp)**2)
                dWt(cmp,inw)=(-4.d0+8.d0*norp(cmp)-4.d0*(norp(cmp)**2))
     %                  *(el_cx(cmp)/abs(el_cx(cmp)))*Hbox(cmp)/wscl
               else
                Wt(cmp,inw)=0.d0
                dWt(cmp,inw)=0.d0
             end if
            end if
           
          end do

          Wtx_o(inw) = Wt(1,inw)*Wt(2,inw)*Wt(3,inw)
          dWtx_o(inw) = dWt(1,inw)*Wt(2,inw)*Wt(3,inw)
          dWty_o(inw) = Wt(1,inw)*dWt(2,inw)*Wt(3,inw)
          dWtz_o(inw) = Wt(1,inw)*Wt(2,inw)*dWt(3,inw)
!     ------------------------------------------------

!     Allocating variables for support domain
          prp_o(inw) = pr(i,j,k) 

          ui_o(inw,1) = 0.5*(vx(i,j,k)+vx(i+1,j,k))
          ui_o(inw,2) = 0.5*(vy(i,j,k)+vy(i,j+1,k))
          ui_o(inw,3) = 0.5*(vz(i,j,k)+vz(i,j,k+1))
         
          inw = inw + 1
         end do
        end do
       end do

!     Inner probe calculation
!     spline weights for all three components and initialise ui matrix
      inw=1
      ui_i(:,:) = 0.0d0

      do k=pind_i_i(3),pind_o_i(3)
       do j=pind_i_i(2),pind_o_i(2)
        do i=pind_i_i(1),pind_o_i(1)
 
          el_cx(1)=xm(i)-pos_pro_i(1)
          el_cx(2)=ym(j)-pos_pro_i(2)
          el_cx(3)=zm(k)-pos_pro_i(3)

          elmag=sqrt(el_cx(1)**2+el_cx(2)**2+el_cx(3)**2)

          do cmp=1,3
            norp(cmp) = abs(el_cx(cmp))*Hbox(cmp)/wscl
        
!     Computing weight functions

!     ----------------EXPONENTIAL SPLINES------------------
            if(wexp.eq.1)then 
              if(norp(cmp).le.1.0d0)then
                Wt(cmp,inw)=exp(-(norp(cmp)/wcon)**2)
                dWt(cmp,inw)=(-2.d0/(wcon**2))*norp(cmp)*
     %                    exp(-(norp(cmp)/wcon)**2)
     %                   *(el_cx(cmp)/abs(el_cx(cmp)))*Hbox(cmp)/wscl
              else
                Wt(cmp,inw)=0.d0
                dWt(cmp,inw)=0.d0
              end if
            end if
!     ----------------CUBIC SPLINES---------------------
            if(wcub.eq.1)then
              if(norp(cmp).le.0.5d0)then
                Wt(cmp,inw)=(2.d0/3.d0)-4.d0*(norp(cmp)**2)+
     %                              4.d0*(norp(cmp)**3)
                dWt(cmp,inw)=(-8.d0*(norp(cmp))+12.d0*(norp(cmp)**2))
     %                  *(el_cx(cmp)/abs(el_cx(cmp)))*Hbox(cmp)/wscl
                elseif(norp(cmp).le.1.d0)then
                Wt(cmp,inw)=(4.d0/3.d0)*(1.d0-norp(cmp)**3)-
     %                        4.d0*(norp(cmp)-norp(cmp)**2)
                dWt(cmp,inw)=(-4.d0+8.d0*norp(cmp)-4.d0*(norp(cmp)**2))
     %                  *(el_cx(cmp)/abs(el_cx(cmp)))*Hbox(cmp)/wscl
               else
                Wt(cmp,inw)=0.d0
                dWt(cmp,inw)=0.d0
             end if
            end if
           
          end do

          Wtx_i(inw) = Wt(1,inw)*Wt(2,inw)*Wt(3,inw)
          dWtx_i(inw) = dWt(1,inw)*Wt(2,inw)*Wt(3,inw)
          dWty_i(inw) = Wt(1,inw)*dWt(2,inw)*Wt(3,inw)
          dWtz_i(inw) = Wt(1,inw)*Wt(2,inw)*dWt(3,inw)
!     ------------------------------------------------

!     Allocating variables for support domain
          prp_i(inw) = pr(i,j,k) 

          ui_i(inw,1) = 0.5*(vx(i,j,k)+vx(i+1,j,k))
          ui_i(inw,2) = 0.5*(vy(i,j,k)+vy(i,j+1,k))
          ui_i(inw,3) = 0.5*(vz(i,j,k)+vz(i,j,k+1))
         
          inw = inw + 1
         end do
        end do
       end do


!     --------------------------------------------------------
!     --------------------------------------------------------

!     Second part - outer probe calculations      
      B(:,:)=0.0d0;Bx(:,:)=0.0d0;By(:,:)=0.0d0;Bz(:,:)=0.0d0
      pinvA(:,:)=0.0d0;invA(:,:)=0.0d0
      pinvAx(:,:)=0.0d0;pinvAy(:,:)=0.0d0;pinvAz(:,:)=0.0d0

      inw = 1

!     pre-inverse matrix A(x) and B(x)
      do k=pind_i_o(3),pind_o_o(3)      
       do j=pind_i_o(2),pind_o_o(2)
        do i=pind_i_o(1),pind_o_o(1)

        pxk(1,1)=1.0d0;pxk(2,1)=xm(i)
        pxk(3,1)=ym(j);pxk(4,1)=zm(k)

        pxkx(1,1)=0.0d0;pxkx(2,1)=1.0d0;pxkx(3,1)=0.0d0
        pxkx(4,1)=0.0d0
        pxky(1,1)=0.0d0;pxky(2,1)=0.0d0;pxky(3,1)=1.0d0
        pxky(4,1)=0.0d0
        pxkz(1,1)=0.0d0;pxkz(2,1)=0.0d0;pxkz(3,1)=0.0d0
        pxkz(4,1)=1.0d0

        call DGEMM('N','T',4,4,1,Wtx_o(inw),pxk,4,pxk,4,
     %                           1.0d0,pinvA,4)

!     -------------------------------------------------------------------
        B(:,inw)=Wtx_o(inw)*pxk(:,1)

        inw = inw + 1
  
        end do
       end do
      end do

!     calling routine to compute inverse
      call inverseLU(pinvA,invA)
            
!     ------------------------------------------------------
!     matrix multiplications for final interpolation
!     DGEMM(transA,transB,m,n,k,alpha,A,LDA,b,LDB,beta,c,LDC)
!     C = alpha * A * B + beta * C

!     ---------------Shape function calculation---------------
      call DGEMM('N','N',1,4,4,1.0d0,ptx_o,1,invA,4,0.0d0,ptxA,1) 
      call DGEMM('N','N',1,nel,4,1.0d0,ptxA,1,B,4,0.0d0,ptxAB,1) 

!     --------------Pressure Interpolation at probe from \phi---------
      call DGEMM('N','N',1,1,nel,1.0d0,ptxAB,1,prp_o(:),nel,0.0d0,
     %                                                ptxABu,1) 
      
      if(wcheck.eq.1)then
        if(sum(ptxAB).le.0.95.or.sum(ptxAB).gt.1.05)then
        print*,'Shape function sum mlsStruc',sum(ptxAB),ptxA(:,:),
     %            prp_o(:),Wtx_o(:),pinvA(:,:),invA(:,:)
        stop
        end if
        if(isnan(sum(ptxAB))) then
        write(*,*)'NaN Shape function mlsStruc',pos_pro_o(1:3)
        end if
      end if

      pr_p_o=ptxABu(1,1)

!     Face normals 
      fnca(1:3) = tri_nor(1:3,ntr,inp)
            
!     Compute pressure at marker from probe
      press_o=pr_p_o+Hboxx*(acc_tri(1,ntr,inp)*fnca(1)+
     %        acc_tri(2,ntr,inp)*fnca(2)+acc_tri(3,ntr,inp)*fnca(3))


!     -----------------------------------------------------------------
!     -----------------------------------------------------------------
!     Second part - inner probe calculations      
      B(:,:)=0.0d0;Bx(:,:)=0.0d0;By(:,:)=0.0d0;Bz(:,:)=0.0d0
      pinvA(:,:)=0.0d0;invA(:,:)=0.0d0
      pinvAx(:,:)=0.0d0;pinvAy(:,:)=0.0d0;pinvAz(:,:)=0.0d0

      inw = 1

!     pre-inverse matrix A(x) and B(x)
      do k=pind_i_i(3),pind_o_i(3)      
       do j=pind_i_i(2),pind_o_i(2)
        do i=pind_i_i(1),pind_o_i(1)

        pxk(1,1)=1.0d0;pxk(2,1)=xm(i)
        pxk(3,1)=ym(j);pxk(4,1)=zm(k)

        pxkx(1,1)=0.0d0;pxkx(2,1)=1.0d0;pxkx(3,1)=0.0d0
        pxkx(4,1)=0.0d0
        pxky(1,1)=0.0d0;pxky(2,1)=0.0d0;pxky(3,1)=1.0d0
        pxky(4,1)=0.0d0
        pxkz(1,1)=0.0d0;pxkz(2,1)=0.0d0;pxkz(3,1)=0.0d0
        pxkz(4,1)=1.0d0

        call DGEMM('N','T',4,4,1,Wtx_i(inw),pxk,4,pxk,4,
     %                           1.0d0,pinvA,4)

        B(:,inw)=Wtx_i(inw)*pxk(:,1)

        inw = inw + 1
  
        end do
       end do
      end do

!     calling routine to compute inverse
      call inverseLU(pinvA,invA)
            
!     ------------------------------------------------------
!     matrix multiplications for final interpolation
!     DGEMM(transA,transB,m,n,k,alpha,A,LDA,b,LDB,beta,c,LDC)
!     C = alpha * A * B + beta * C

!     ---------------Shape function calculation---------------
      call DGEMM('N','N',1,4,4,1.0d0,ptx_i,1,invA,4,0.0d0,ptxA,1) 
      call DGEMM('N','N',1,nel,4,1.0d0,ptxA,1,B,4,0.0d0,ptxAB,1) 

!     --------------Pressure Interpolation at probe from \phi---------
      call DGEMM('N','N',1,1,nel,1.0d0,ptxAB,1,prp_i(:),nel,0.0d0,
     %                                                ptxABu,1) 
      
      if(wcheck.eq.1)then
        if(sum(ptxAB).le.0.95.or.sum(ptxAB).gt.1.05)then
        print*,'Shapeifunction sum mlsStruc',sum(ptxAB),ptxA(:,:),
     %            prp_o(:),Wtx_o(:),pinvA(:,:),invA(:,:)
        stop
        end if
        if(isnan(sum(ptxAB))) then
        write(*,*)'NaN Shape function mlsStruc',pos_pro_i(1:3)
        end if
      end if

      pr_p_i=ptxABu(1,1)

!     -----------------------------------------------------------------
!     Face normals 
      fnca(1:3) = tri_nor(1:3,ntr,inp)
            
!     Compute pressure at marker from probe
      press_i=pr_p_i-Hboxx*(acc_tri(1,ntr,inp)*fnca(1)+
     %        acc_tri(2,ntr,inp)*fnca(2)+acc_tri(3,ntr,inp)*fnca(3))

!     -----------------------------------------------------------------
      press_face(ntr,inp) = (press_o-press_i)*sur(ntr,inp)


       else
      
      press_face(ntr,inp) = 0.0d0      
       end if 


       end do
      end do

c     --------------------------------------------------------
!     --------------------------------------------------------
 
      return
      end
