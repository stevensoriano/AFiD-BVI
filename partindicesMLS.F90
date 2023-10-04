!------------------------------------------------------------------
!     Given position of lagrangian marker - this routine computes
!     the indices and the indices of the support domain
!     Only starting and ending indices of the support domain are 
!     reqiured. Then loop over it.	
!------------------------------------------------------------------

      subroutine partindicesMLS(pos,pind_i,pind_o,dismaxpro)
      USE param
      IMPLICIT NONE
      
      real, intent(in) ::  pos(3)
      integer, dimension(6), intent(inout) :: pind_i, pind_o
      integer i1,j1,k1,ist,jst,kst
      integer ii,jj,kk
      real tstr,rcdp,x2dp
      real a1,a2,alp,alp1
!     --------INNER AND OUTER VARIABLES----------
      real pos_i(3),pos_o(3),dists(8),dismaxpro
      real disdum,distsdum
      integer i1i,j1i,k1i,isti,jsti,ksti
      integer i1o,j1o,k1o,isto,jsto,ksto
      
      real tstri,rcdpi,x2dpi
      real tstro,rcdpo,x2dpo

      real a1i,a2i,alpi,alp1i
      real a1o,a2o,alpo,alp1o
!     --------------------------------------------

!     ++++++++Indices of the marker++++++++++++++++++++++++	
!     X indices
      i1=FLOOR(pos(1)*dx1) + 1

      if(i1.gt.n1m)i1=1
      if(i1.le.1)i1=n1m
!     staggered
      if(pos(1).gt.xm(i1))ist=i1
      if(pos(1).le.xm(i1))ist=i1-1

      if(ist.eq.0)ist=n1m 

!     Y indices
      j1=FLOOR(pos(2)*dx2) + 1

      if(j1.gt.n2m) j1=1
      if(j1.le.1) j1=n2m
!     staggered
      if(pos(2).gt.ym(j1))jst=j1
      if(pos(2).le.ym(j1))jst=j1-1

      if(jst.eq.0)jst=n2m 

!     Z indices
     
!     uniform grid : wall-normal direction
      k1=pos(3)*dx3 + 1
      if(k1.gt.n3m)write(6,*)'Wall normal index error in findindices.F'
      if(k1.lt.1)write(6,*)  'Wall normal index error in findindices.F'
!     staggered
      if(pos(3).gt.zm(k1))kst=k1
      if(pos(3).le.zm(k1))kst=k1-1

!     OUTER and INNER INDICES
      if(nel.eq.27)then

      i1o = i1+1 ; i1i = i1-1
      j1o = j1+1 ; j1i = j1-1
      k1o = k1+1 ; k1i = k1-1

      isto = ist+1 ; isti = ist-1
      jsto = jst+1 ; jsti = jst-1
      ksto = kst+1 ; ksti = kst-1

      elseif(nel.eq.125)then
      i1o = ist+2 ; i1i = ist-2
      j1o = jst+2 ; j1i = jst-2
      k1o = kst+2 ; k1i = kst-2

      isto = ist+1 ; isti = ist-1
      jsto = jst+1 ; jsti = jst-1
      ksto = kst+1 ; ksti = kst-1


      end if     

!     -------------------------------------------------------

      pind_i(1)=i1i ; pind_i(2)=j1i ; pind_i(3)=k1i
      pind_i(4)=isti; pind_i(5)=jsti; pind_i(6)=ksti

      pind_o(1)=i1o ; pind_o(2)=j1o ; pind_o(3)=k1o
      pind_o(4)=isto; pind_o(5)=jsto; pind_o(6)=ksto
      
!     ----------------------------------------------------------

      disdum = -100000.0
      do kk=ksti,ksto
       do jj=jsti,jsto
        do ii=isti,isto
 
        distsdum=(pos(1)-xm(ii))**2+(pos(2)-ym(jj))**2 &
            +(pos(3)-zm(kk))**2

        disdum = max(disdum,distsdum)

        end do
       end do
      end do


      dismaxpro = disdum


      return
      end


