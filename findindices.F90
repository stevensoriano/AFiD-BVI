!------------------------------------------------------------------
!     Routine to compute indices of all centroids on all particles
!     Only starting and ending indices of the support domain are 
!     reqiured. Then loop over it.	
!------------------------------------------------------------------

      subroutine findindices
      USE param
      USE mpih
      USE mls_param, only : pind,bboxind,dismax,Nparticle,maxnf,tri_bar
      IMPLICIT NONE
      
      real pos(3)
      integer i1,j1,k1,ist,jst,kst,inp,ntr
      integer ii,jj,kk
      real tstr,rcdp,x2dp
      real a1,a2,alp,alp1
!     --------INNER AND OUTER VARIABLES----------
      integer i1i,j1i,k1i,isti,jsti,ksti
      integer i1o,j1o,k1o,isto,jsto,ksto
      
      real tstri,rcdpi,x2dpi
      real tstro,rcdpo,x2dpo

      real a1i,a2i,alpi,alp1i
      real a1o,a2o,alpo,alp1o
!     --------------------------------------------
      real disdum,distsdum
      integer ind_pal

      do inp=1,Nparticle

       do ntr=1,maxnf

       pos(1)=tri_bar(1,ntr,inp)
       pos(2)=tri_bar(2,ntr,inp)
       pos(3)=tri_bar(3,ntr,inp)
!     ++++++++Indices of the marker++++++++++++++++++++++++	
!     X - indices
      i1=FLOOR(pos(1)*dx1) + 1

      if(i1.gt.n1m)i1=1
      if(i1.le.1)i1=n1m
!     staggered
      if(pos(1).gt.xm(i1))ist=i1
      if(pos(1).le.xm(i1))ist=i1-1

      if(ist.eq.0)ist=n1m 

!     Y - indices
      j1=FLOOR(pos(2)*dx2) + 1

      if(j1.gt.n2m) j1=1
      if(j1.le.1) j1=n2m
!     staggered
      if(pos(2).gt.ym(j1))jst=j1
      if(pos(2).le.ym(j1))jst=j1-1

      if(jst.eq.0)jst=n2m 

!     Z - indices
!     
     
      k1=FLOOR(pos(3)*dx3) + 1

      if(k1.gt.n3m) k1=1
      if(k1.le.1) k1=n3m
!     staggered
      if(pos(3).gt.zm(k1))kst=k1
      if(pos(3).le.zm(k1))kst=k1-1

      if(kst.eq.0)kst=n3m 



!     OUTER and INNER INDICES
      isto = ist+1 ; isti = ist-1
      jsto = jst+1 ; jsti = jst-1
      ksto = kst+1 ; ksti = kst-1

!     dummies
      i1i = 1; i1o = 1
      j1i = 1; j1o = 1
      k1i = 1; k1o = 1

!     -------------------------------------------------------------
      pind(1,ntr,inp)=i1 ; pind(2,ntr,inp)=j1 ; pind(3,ntr,inp)=k1
      pind(4,ntr,inp)=ist; pind(5,ntr,inp)=jst; pind(6,ntr,inp)=kst
      if(myid.eq.0)then
!            write(*,*) pind(1:6,ntr,inp)
      end if
!     -------------------------------------------------------------

!      disdum = -100000.0
!      do kk=ksti,ksto
!       do jj=jsti,jsto
!        do ii=isti,isto
! 
!      distsdum=(pos(1)-xm(ii))**2+(pos(2)-rm(jj))**2+(pos(3)-zm(kk))**2
!
!      disdum = max(disdum,distsdum)
!
!        end do
!       end do
!      end do
!
!      dismax(ntr,inp) = disdum

       end do

!     Compute average axial index for allocating particle to processor
      bboxind(1,1,inp) = MINVAL(pind(4,:,inp)) 
      bboxind(1,2,inp) = MAXVAL(pind(4,:,inp))
      bboxind(2,1,inp) = MINVAL(pind(5,:,inp))
      bboxind(2,2,inp) = MAXVAL(pind(5,:,inp))
      bboxind(3,1,inp) = MINVAL(pind(6,:,inp))
      bboxind(3,2,inp) = MAXVAL(pind(6,:,inp))
      
      ind_pal = 0.5*(bboxind(3,1,inp)+bboxind(3,2,inp)) - 0
      
      if(ind_pal.lt.1.or.ind_pal.gt.n3m)then
      write(*,*),'Cannot allocate processor',ind_pal,bboxind(3,1,inp), &
                        bboxind(3,2,inp)
      write(*,*),bboxind(1,1,inp),bboxind(1,2,inp)
      write(*,*),bboxind(2,1,inp),bboxind(2,2,inp)
      write(*,*),bboxind(3,1,inp),bboxind(3,2,inp)
      stop
      end if      

      end do

      return
      end


