!***********************************************************************************
          subroutine ibgrid3d (istl,ni,nj,nk,th,rr,zz &
                               ,mpugeo,nb,xyzb,nxyz,bar,area  &
                               ,mpun,mpuni,nins,nintf,indgeoi &
                               ,indgeoe,indgeoee,distb,omega,vel,cen &
                               ,vb,vbi,ntrib &
                               ,nii,nif,nji,njf,nki,nkf)

          implicit none
          integer ni,nj,nk,nb,mpun,mpuni,mpugeo,indgeoi(4,mpuni,3)     
          integer nins(4),indgeoe(4,mpun,3),indgeoee(4,mpun,3),nintf(4)
          integer ntrib(4,mpun),nii,nif,nji,njf,nki,nkf,niff,njff,nkff,ica
          real th(ni),rr(nj),zz(nk),xyzb(mpugeo,9),nxyz(mpugeo,3)
          real bar(mpugeo,3),area(mpugeo),distb(4,mpun),vb(4,mpun),vbi(4,mpuni)
          real omega(3),vel(3),cen(3)
          real,    dimension(:,:  ), allocatable :: xp,yp!,zp
          real,    dimension(:    ), allocatable :: zp
          real,    dimension(:    ), allocatable :: coefang
          real,    dimension(:,:  ), allocatable :: probe,deltav
          integer, dimension(:,:  ), allocatable :: iposbar,iposprobe
          integer, dimension(:,:,:), allocatable :: ibpoint,ibcell1,ibcell2,ibcell3,ibniente
          character(50) string
          integer nx,ny,nz,isclive3d,isclive3dp,istl,mins,minter,idummy
          integer i,j,k,l,loop,iloop,itmp,ie,je,ke,ntrp
          real degr,fx,fy,fz,anp,tras,delta,wp
          real vx,vy,vz
          real xm,ym,zm,vtrz
          real csi,xni,eta,fn1,fn2,fn3,fn4,fn5,fn6,fn7,fn8
          real pig,erre,theta
          integer ip,jp,kp,ibsum,ibv1,ibv2,ibv3,ibv4,ibv5,ibv6,nzmin,nzmax

!...tag cell centers and vertices
!
          pig=2.*asin(1.)
          nins=0
          nintf=0
!
          do loop=1,4
             write(*,*)'Loop: ',loop
             if (loop.eq.1) then
                nz=nk-1
                ny=nj-1
                nx=ni
                niff = ni-1 !marcoooo
                njff = nj-1
                nkff = nk-1
                allocate(ibcell1(nx,ny,nz),ibpoint(nx,ny,nz),xp(nx,ny),yp(nx,ny),zp(nz))
                do k=1,nz
                   zp(k)=0.5*(zz(k)+zz(k+1))
                enddo
                do j=1,ny
                   do i=1,nx
                      xp(i,j)=(th(i))
                      yp(i,j)=0.5*(rr(j)+rr(j+1))
                   enddo
                enddo
             elseif (loop.eq.2) then
                nz=nk-1
                ny=nj
                nx=ni-1
                niff = ni-1
                njff = nj
                nkff = nk-1
                allocate(ibcell2(nx,ny,nz),ibpoint(nx,ny,nz),xp(nx,ny),yp(nx,ny),zp(nz))
                do k=1,nz
                   zp(k)=0.5*(zz(k)+zz(k+1))
                enddo
                do j=1,ny
                   do i=1,nx
                      xp(i,j)=(0.5*(th(i)+th(i+1)))
                      yp(i,j)=(rr(j))
                   enddo
                enddo
             elseif (loop.eq.3) then
                nz=nk
                ny=nj-1
                nx=ni-1
                niff = ni-1
                njff = nj-1
                nkff = nk
                allocate(ibcell3(nx,ny,nz),ibpoint(nx,ny,nz),xp(nx,ny),yp(nx,ny),zp(nz))
                zp=zz
                do j=1,ny
                   do i=1,nx
                      xp(i,j)= (0.5*(th(i)+th(i+1)))
                      yp(i,j)=(0.5*(rr(j)+rr(j+1)))
                   enddo
                enddo
             elseif (loop.eq.4) then
                nz=nk-1
                ny=nj-1
                nx=ni-1
                niff = ni-1
                njff = nj-1
                nkff = nk-1
                allocate(ibpoint(nx,ny,nz),xp(nx,ny),yp(nx,ny),zp(nz))
                do k=1,nz
                   zp(k)=0.5*(zz(k)+zz(k+1))
                enddo
                do j=1,ny
                   do i=1,nx
                      xp(i,j)= (0.5*(th(i)+th(i+1)))
                      yp(i,j)=(0.5*(rr(j)+rr(j+1)))
                   enddo
                enddo
          endif
!
          nzmin=nki
          nzmax=nkff

          ibpoint=1

          call tag3d(xp,yp,zp,nx,ny,nz,mpugeo,nb,xyzb,ibpoint,nii,niff,nji,njff,nki,nkff)

          if(loop.eq.1) ibcell1=ibpoint
          if(loop.eq.2) ibcell2=ibpoint
          if(loop.eq.3) ibcell3=ibpoint
 
          do k=nki,nkff
             do j=nji,njff
                do i=nii,niff
                   itmp = isclive3d(i,j,k,nx,ny,nz,ibpoint)
                   if (itmp.eq. 0) then
                      nintf(loop) = nintf(loop) + 1
                      call int3d(i,j,k,nx,ny,nz,xp,yp,zp,mpugeo,nb,xyzb,nxyz,ibpoint,ie,je,ke,wp,xm,ym,zm,ntrp)
 
                      if (loop.eq.4)then
                         ip = i+1
                         jp = j+1
                         kp = k+1

                         if (i.eq.nx) ip = 1

                         ibv1  = ibcell1 (i ,j ,k )
                         ibv2  = ibcell1 (ip,j ,k )
                         ibv3  = ibcell2 (i ,j ,k )
                         ibv4  = ibcell2 (i ,jp,k )
                         ibv5  = ibcell3 (i ,j ,k )
                         ibv6  = ibcell3 (i ,j ,kp)
                         ibsum = ibv1+ibv2+ibv3+ibv4+ibv5+ibv6
!                         if(ibsum.eq.6)then
!                            nintf(loop) = nintf(loop) - 1   !
!                            ntrib(loop,nintf(loop))=0
!                            goto 729
!                         endif
                      endif
 
                      indgeoe (loop,nintf(loop),1)=i
                      indgeoe (loop,nintf(loop),2)=j
                      indgeoe (loop,nintf(loop),3)=k
                      indgeoee(loop,nintf(loop),1)=ie
                      indgeoee(loop,nintf(loop),2)=je
                      indgeoee(loop,nintf(loop),3)=ke
                      distb(loop,nintf(loop))=wp
                      ntrib(loop,nintf(loop))=ntrp
                      vb(loop,nintf(loop))=0.0
729                end if
 
                   if (itmp.eq.-1) then
                      nins(loop) = nins(loop) + 1
                      indgeoi(loop,nins(loop),1)=i
                      indgeoi(loop,nins(loop),2)=j
                      indgeoi(loop,nins(loop),3)=k
                      vbi(loop,nins(loop))=0.0
                   end if
                enddo
             enddo
          enddo

945       deallocate(ibpoint,xp,yp,zp)
          enddo

          deallocate(ibcell1,ibcell2,ibcell3)
          return
          end subroutine ibgrid3d
!...................................................................................
!***********************************************************************************
          subroutine intersect_triangle(orig,dir,vert0,vert1,vert2,t,u,v,intersect)

          implicit none
          integer intersect
          real u,v,t,det,inv_det
          real, dimension (3) ::  edge1,edge2,tvec,pvec,qvec,orig,dir,vert0,vert1,vert2
          real epsilon

          epsilon = 1.0e-12

            t = 0.
            u = 0.
            v = 0.

!.....find vectors for two edges sharing vert0 
            call sub(edge1, vert1, vert0)
            call sub(edge2, vert2, vert0)

!.....begin calculating determinant - also used to calculate U parameter 
            call cross(pvec, dir, edge2)

!.....if determinant is near zero, ray lies in plane of triangle 
            call dot(det, edge1, pvec)

            if (det.gt.-epsilon.and.det.lt.epsilon) then
               intersect = 2
               return
            endif
            inv_det = 1.0 / det

!......calculate distance from vert0 to ray origin 
            call sub(tvec, orig, vert0)

!.....calculate U parameter and test bounds 
            call dot(u, tvec, pvec)
            u = u * inv_det

            if (u.lt.0.0.or.u.gt.1.0) then
               intersect = 3
               return
            endif

!.....prepare to test V parameter 
            call cross(qvec, tvec, edge1)

!.....calculate V parameter and test bounds 
            call dot(v, dir, qvec)
            v = v * inv_det

            if (v.lt.0.0.or.(u+v).gt.1.0) then
               intersect = 3
               return
            endif

!.....calculate t, ray intersects triangle 
            call dot(t, edge2, qvec)
            t = t * inv_det

            intersect = 1

            return
          end subroutine intersect_triangle
!...................................................................................
!***********************************************************************************
          subroutine segtriint(x1,y1,z1,x2,y2,z2,xa,ya,za,xb,yb,zb,xc,yc,zc,xm,ym,zm,iflag)
          implicit none
          real x1,y1,z1,x2,y2,z2,xa,ya,za,xb,yb,zb,xc,yc,zc,xm,ym,zm
          real q(3),r(3),va(3),vb(3),vc(3),t,u,v
          integer iflag,code
!
!...redefine ray as origin plus lenght vector
!
          q(1) = x1
          q(2) = y1
          q(3) = z1
          r(1) = x2-x1
          r(2) = y2-y1
          r(3) = z2-z1
          va(1) = xa
          va(2) = ya
          va(3) = za
          vb(1) = xb
          vb(2) = yb
          vb(3) = zb
          vc(1) = xc
          vc(2) = yc
          vc(3) = zc

          iflag = 0.
          call intersect_triangle(q,r,va,vb,vc,t,u,v,code)
          if (code.eq.2) call intersect_triangle(q,r,va,vc,vb,t,u,v,code)
           
          xm = q(1)+t*r(1)
          ym = q(2)+t*r(2)
          zm = q(3)+t*r(3)

          if (code.eq.1) iflag = 1
          return
          end subroutine segtriint
!...................................................................................
!***********************************************************************************
          integer function isclive3d(i,j,k,ni,nj,nk,ibcc)

          implicit none
          integer i,j,k,ni,nj,nk,ibcc(ni,nj,nk)
          integer ibv0,ibv1,ibv2,ibv3,ibv4,ibv5,ibv6,ibsum
          integer im,ip,jm,jp,km,kp

          im = max(1 ,i-1)
          ip = min(ni,i+1)
          jm = max(1 ,j-1)
          jp = min(nj,j+1)
          km = max(1 ,k-1)
          kp = min(nk,k+1)
          if (i.eq.1 ) im = ni !---periodicità su X
          if (i.eq.ni) ip = 1    !---periodicità su X
          if (j.eq.1 ) jm = nj !---periodicità su Y
          if (j.eq.nj) jp = 1    !---periodicità su Y


          ibv0  = ibcc (i ,j ,k )
          ibv1  = ibcc (im,j ,k )
          ibv2  = ibcc (ip,j ,k )
          ibv3  = ibcc (i ,jm,k )
          ibv4  = ibcc (i ,jp,k )
          ibv5  = ibcc (i ,j ,km)
          ibv6  = ibcc (i ,j ,kp)
          ibsum = 0
          ibsum = ibv1+ibv2+ibv3+ibv4+ibv5+ibv6

          if(abs(ibsum).eq.0) goto 123
          if(abs(ibsum).eq.2) goto 123
          if(abs(ibsum).eq.4) goto 123
          if(abs(ibsum).eq.6) goto 123
          write(*,*) 'RAELLY WARNING',ibsum
!          pause
123       continue

          if ((ibv0.eq.-1).and.(ibsum.ge.-6)) isclive3d = -1  ! ins
          if ((ibv0.eq.1).and.(ibsum.eq. 6)) isclive3d =  1  ! ext                
          if ((ibv0.eq.1).and.(ibsum.lt. 6)) isclive3d =  0  ! int            

          if (isclive3d.eq.0) then
!             if (i.eq.1.or.i.eq.ni) isclive3d = -1
!             if (j.eq.1.or.j.eq.nj) isclive3d = -1
!             if (k.eq.1.or.k.eq.nk) isclive3d = -1
          endif

          return
          end function isclive3d
!...................................................................................
!***********************************************************************************
          subroutine int3d(i,j,k,nx,ny,nz,xp,yp,zp,mpugeo,nb,xyzb,nxyz,ibcc,ie,je,ke,wp,xm,ym,zm,ntrp)

          implicit none
          integer i,j,k,nx,ny,nz,nb,id,n,mpugeo,ntr1,ntr2,ntr3,ntrp
          integer im,ip,jm,jp,km,kp,ie,je,ke
          real xp(nx,ny),yp(nx,ny),zp(nz),w(7),wi(7),wp,wpi
          real xyzb(mpugeo,9),nxyz(mpugeo,3),norm(6),normp
          real xcc,xccm,xccp,ycc,yccm,yccp,zcc,zccp,zccm
          real di,dj,dk,dm,wmax,wmin
          real d,wfip,wfjp,wfim,wfjm,wsum,wwal
          real xm,ym,zm,xmt,ymt,zmt,xm1,ym1,zm1,xm2,ym2,zm2,xm3,ym3,zm3
          real dist1,dist2p,dist2m,distmp,distmt,distpt
          integer iflag,jflag,kflag,iflagt
          integer ibcc(nx,ny,nz),ibv0,ibv1,ibv2,ibv3,ibv4,ibv5,ibv6,ibsum
          real xa,ya,za,xb,yb,zb,xc,yc,zc
          real minxaxbxc,minyaybyc,minzazbzc
          real maxxaxbxc,maxyaybyc,maxzazbzc
          real ycen,zcen,omegat,erre,theta,pig,vty
          real vyi,vzi,vri,vthi
          real vyj,vzj,vrj,vthj
          real vyk,vzk,vrk,vthk

          w=0.
          wmax=0.
          wmin=0.

          im = max(1 ,i-1)
          ip = min(nx,i+1)
          jm = max(1 ,j-1)
          jp = min(ny,j+1)
          km = max(1 ,k-1)
          kp = min(nz,k+1)
          if (i.eq.1 ) im = nx   !---periodicità su X    
          if (i.eq.nx) ip = 1
          if (j.eq.1 ) jm = ny   !---periodicità su Y
          if (j.eq.ny) jp = 1


          ibv0  = ibcc (i ,j ,k )
          ibv1  = ibcc (im,j ,k )
          ibv2  = ibcc (ip,j ,k )
          ibv3  = ibcc (i ,jm,k )
          ibv4  = ibcc (i ,jp,k )
          ibv5  = ibcc (i ,j ,km)
          ibv6  = ibcc (i ,j ,kp)
          ibsum = ibv1+ibv2+ibv3+ibv4+ibv5+ibv6

          if (ibv0 .ne. 1) write(*,*)'WARNING IBV0 ',ibv0
          if (abs(ibsum) .eq. 6) write(*,*)'WARNING IBVSUM ',ibsum

          zcc  = zp(k)
          ycc  = yp(i,j)
          xcc  = xp(i,j)

          if (ibv1.ne.ibv2) then
             yccm = yp(im,j)
             yccp = yp(ip,j)
             xccm = xp(im,j)
             xccp = xp(ip,j)
             do n=1,nb
                xa = xyzb(n,1)
                ya = xyzb(n,2)
                za = xyzb(n,3)
                xb = xyzb(n,4)
                yb = xyzb(n,5)
                zb = xyzb(n,6)
                xc = xyzb(n,7)
                yc = xyzb(n,8)
                zc = xyzb(n,9)

                minxaxbxc = min(xa,xb,xc)
                minyaybyc = min(ya,yb,yc)
                minzazbzc = min(za,zb,zc)
                maxxaxbxc = max(xa,xb,xc)
                maxyaybyc = max(ya,yb,yc)
                maxzazbzc = max(za,zb,zc)

                if (zcc.ge.minzazbzc.and.zcc.le.maxzazbzc) then
                   call segtriint(xccm,yccm,zcc,xccp,yccp,zcc,xa,ya,za,xb,yb,zb,xc,yc,zc,xmt,ymt,zmt,iflagt)
                   if (iflagt.eq.1) then
                      distmp=sqrt((yccm-yccp)**2+(xccm-xccp)**2)
                      distmt=sqrt((yccm-ymt )**2+(xccm-xmt )**2)
                      distpt=sqrt((yccp-ymt )**2+(xccp-xmt )**2)
                      if (1.001*distmp.ge.(distmt+distpt))then
                         xm1 = xmt
                         ym1 = ymt
                         zm1 = zmt
                         dist1 =sqrt((ycc -ym1)**2+(xcc -xm1)**2)
                         dist2p=sqrt((yccp-ym1)**2+(xccp-xm1)**2)
                         dist2m=sqrt((yccm-ym1)**2+(xccm-xm1)**2)
                         if (ibv1.eq.1) w(1) = dist1/dist2m
                         if (ibv2.eq.1) w(2) = dist1/dist2p
                         if (w(1).ne.0.0.and.w(2).ne.0.0) write(*,*) 'REALLY WARNING:iflag'
                         ntr1=n
                      endif
                   endif
                endif
               enddo
            endif
!'
            if (ibv3.ne.ibv4) then
               yccm = yp(i,jm)
               yccp = yp(i,jp)
               xccm = xp(i,jm)
               xccp = xp(i,jp)
               do n=1,nb
                  xa = xyzb(n,1)
                  ya = xyzb(n,2)
                  za = xyzb(n,3)
                  xb = xyzb(n,4)
                  yb = xyzb(n,5)
                  zb = xyzb(n,6)
                  xc = xyzb(n,7)
                  yc = xyzb(n,8)
                  zc = xyzb(n,9)
!                  
                  minxaxbxc = min(xa,xb,xc)
                  minyaybyc = min(ya,yb,yc)
                  minzazbzc = min(za,zb,zc)
                  maxxaxbxc = max(xa,xb,xc)
                  maxyaybyc = max(ya,yb,yc)
                  maxzazbzc = max(za,zb,zc)

                  if (zcc.ge.minzazbzc.and.zcc.le.maxzazbzc) then
                     call segtriint(xccm,yccm,zcc,xccp,yccp,zcc,xa,ya,za,xb,yb,zb,xc,yc,zc,xmt,ymt,zmt,iflagt)

                     if (iflagt.eq.1) then
                        distmp=sqrt((yccm-yccp)**2+(xccm-xccp)**2)
                        distmt=sqrt((yccm-ymt )**2+(xccm-xmt )**2)
                        distpt=sqrt((yccp-ymt )**2+(xccp-xmt )**2)

                        if (1.001*distmp.ge.(distmt+distpt))then
                           xm2 = xmt
                           ym2 = ymt
                           zm2 = zmt
                           dist1 =sqrt((ycc -ym2)**2+(xcc -xm2)**2)
                           dist2p=sqrt((yccp-ym2)**2+(xccp-xm2)**2)
                           dist2m=sqrt((yccm-ym2)**2+(xccm-xm2)**2)
                           if (ibv3.eq.1) w(3) = dist1/dist2m
                           if (ibv4.eq.1) w(4) = dist1/dist2p
                           if (w(3).ne.0.0.and.w(4).ne.0.0) write(*,*) 'REALLY WARNING:jflag'
                           ntr2=n
                        endif
                     endif
                  endif
               enddo
            endif
!'
            if (ibv5.ne.ibv6) then
               zccm = zp(km)
               zccp = zp(kp)
               do n=1,nb
                  xa = xyzb(n,1)
                  ya = xyzb(n,2)
                  za = xyzb(n,3)
                  xb = xyzb(n,4)
                  yb = xyzb(n,5)
                  zb = xyzb(n,6)
                  xc = xyzb(n,7)
                  yc = xyzb(n,8)
                  zc = xyzb(n,9)
                  
                  minxaxbxc = min(xa,xb,xc)
                  minyaybyc = min(ya,yb,yc)
                  minzazbzc = min(za,zb,zc)
                  maxxaxbxc = max(xa,xb,xc)
                  maxyaybyc = max(ya,yb,yc)
                  maxzazbzc = max(za,zb,zc)

                  if ((ycc.ge.minyaybyc.and.ycc.le.maxyaybyc).and.(xcc.ge.minxaxbxc.and.xcc.le.maxxaxbxc)) then
                     call segtriint(xcc,ycc,zccm,xcc,ycc,zccp,xa,ya,za,xb,yb,zb,xc,yc,zc,xmt,ymt,zmt,iflagt)

                     if (iflagt.eq.1) then
                        if (zmt.ge.zccm.and.zmt.le.zccp) then
                           xm3 = xmt
                           ym3 = ymt
                           zm3 = zmt
                           if (ibv5.eq.1) w (5) = (zcc -zm3)/(zccm-zm3)
                           if (ibv6.eq.1) w (6) = (zcc -zm3)/(zccp-zm3)
                           if (w(5).ne.0.0.and.w(6).ne.0.0) write(*,*) 'REALLY WARNING:kflag'
                           ntr3=n
                        endif
                     endif
                  endif
               enddo
            endif
!'
            wmax = maxval(w)
            wmin = minval(w)
            if (wmax.gt.1.0) write(*,*) 'REALLY WARNING:w>1.0',wmax
            if ((xcc.ne.xmt).and.(ycc.ne.ymt).and.(zcc.ne.zmt)) write(*,*) 'REALLY WARNING:w=0.0'

            if (wmin.lt.0.0) write(*,*) 'REALLY WARNING:w<0.0',wmin

            ie=i
            je=j
            ke=k
            id=0
            wp=1.0

            do n=1,6
               if (w(n).gt.0.000001.and.w(n).lt.wp) then
                  wp = w(n)
                  id=n
               endif
7382        enddo

!          if (id.eq.0) write(*,*)'Error' 

          if (id.eq.1)ie=im
          if (id.eq.2)ie=ip
          if (id.eq.3)je=jm
          if (id.eq.4)je=jp
          if (id.eq.5)ke=km
          if (id.eq.6)ke=kp
          if (id.eq.1.or.id.eq.2)then
             xm = xm1
             ym = ym1
             zm = zm1
             ntrp=ntr1
          endif
          if (id.eq.3.or.id.eq.4)then
             xm = xm2
             ym = ym2
             zm = zm2
             ntrp=ntr2
          endif
          if (id.eq.5.or.id.eq.6)then
             xm = xm3
             ym = ym3
             zm = zm3
             ntrp=ntr3
          endif

          return
          end subroutine int3d
!...................................................................................
!***********************************************************************************
          subroutine tag3d(x,y,z,ni,nj,nk,mpugeo,nb,xyzb,ibs,nii,nif,nji,njf,nki,nkf)

          implicit none
          integer nb,ni,nj,nk,iflag,icross,ncross,i1,i2,isave,n,nins,next,mpugeo,found
          integer ibs(ni,nj,nk),i,j,k,nii,nif,nji,njf,nki,nkf,ica
          real xyzb(mpugeo,9),x(ni,nj),y(ni,nj),z(nk)
          real xcc,ycc,zcc,xm,ym,zm,zcross(100),zcrosso(100),ztmin
          real xa,ya,za,xb,yb,zb,xc,yc,zc

          nins=0
          next=0

          do i=nii,nif
             do j=nji,njf
                xcc = x(i,j)
                ycc = y(i,j)
                icross = 0
                do n=1,nb
                   xa = xyzb(n,1)
                   ya = xyzb(n,2)
                   za = xyzb(n,3)
                   xb = xyzb(n,4)
                   yb = xyzb(n,5)
                   zb = xyzb(n,6)
                   xc = xyzb(n,7)
                   yc = xyzb(n,8)
                   zc = xyzb(n,9)

                   if ((ycc.ge.min(ya,yb,yc).and.ycc.le.max(ya,yb,yc)).and.       &
                       (xcc.ge.min(xa,xb,xc).and.xcc.le.max(xa,xb,xc))) then
                      call segtriint(xcc,ycc,-100.,xcc,ycc,100.,xa,ya,za,xb,yb,zb,xc,yc,zc,xm,ym,zm,iflag)
                      if (iflag.eq.1) then
                         icross = icross+1
                         zcross(icross) = zm
                      endif
                   endif
                enddo

                ncross = icross

!   eliminate duplicate intersection

                if (ncross.gt.2) then
                   icross  = 1
                   zcrosso(1)=zcross(1)
                   do i1=2,ncross
                      found = 0
                      do i2=1,icross
                         if (abs(zcrosso(i2)-zcross(i1)).lt.1.0E-12) found = 1
                      enddo
                      if (found.eq.0) then
                         icross = icross+1
                         zcrosso(icross) = zcross(i1)
                      endif
                   enddo
                   ncross = icross
                   zcross(1:ncross) = zcrosso(1:ncross)
                endif

                do i1=1,ncross
                   ztmin = 1.0e10
                   do i2=1,ncross
                      if (zcross(i2).le.ztmin) then
                         ztmin = zcross(i2)
                         isave = i2;
                      endif
                   enddo
                   zcrosso(i1) = ztmin
                   zcross (isave) = 2.0e10
                enddo

                do k=nki,nkf
                   zcc = z(k)
                   do icross=1,ncross
                      if (zcrosso(icross) .le. zcc) ibs(i,j,k) = -ibs(i,j,k)
                   enddo
                   if(ibs(i,j,k).eq.-1) nins=nins+1
                   if(ibs(i,j,k).eq.1)  next=next+1
                enddo

               enddo
            enddo

          return
          end subroutine tag3d
!...................................................................................
!***********************************************************************************
          subroutine readgeo3ddim(nb,namefile)

          implicit none
          character(50)  namefile,string
          integer nb

          open(11,file=namefile)
          nb = 0
10        read(11,*,err=98)string

          if (string.eq.'solid') goto 10
          if (string.eq.'endsolid') goto 99
          if (string.eq.'facet') goto 10
          if (string.eq.'endfacet') goto 10

          if (string.eq.'outer') then
             nb = nb+1
             read(11,*)string
             read(11,*)string
             read(11,*)string
             goto 10
          endif

          if (string.eq.'endloop') goto 10
98        write(*,*)' Error in STL file '
99        close(11)
          write(*,*)' STL with ',nb,' surface triangles '

          return
          end subroutine readgeo3ddim
!...................................................................................
!***********************************************************************************
          subroutine readgeo3d(mpugeo,xyzb,bar,area,nb,namefile)

          implicit none
          integer j,nb,i,l,mpugeo
          real xyzb(mpugeo,9),bar(mpugeo,3),area(mpugeo)
          character(50)  namefile,string
          real d12,d23,d31,sp

          open(11,file=namefile)
          j = 0
10        read(11,*,err=98)string

          if (string.eq.'solid') goto 10
          if (string.eq.'endsolid') goto 99
          if (string.eq.'facet') goto 10
          if (string.eq.'endfacet') goto 10

          if (string.eq.'outer') then
             j = j+1
             read(11,*)string,(xyzb(j,i),i=1,3)
             read(11,*)string,(xyzb(j,i),i=4,6)
             read(11,*)string,(xyzb(j,i),i=7,9)
             goto 10
          endif

          if (string.eq.'endloop') goto 10
98        write(*,*)' Error in STL file '
99        close(11)

          do l=1,nb
             bar(l,1)=(xyzb(l,1)+xyzb(l,4)+xyzb(l,7))/3.
             bar(l,2)=(xyzb(l,2)+xyzb(l,5)+xyzb(l,8))/3.
             bar(l,3)=(xyzb(l,3)+xyzb(l,6)+xyzb(l,9))/3.
     
             d12=sqrt( (xyzb(l,1)-xyzb(l,4))**2. &
                      +(xyzb(l,2)-xyzb(l,5))**2. &
                      +(xyzb(l,3)-xyzb(l,6))**2. )
             d23=sqrt( (xyzb(l,4)-xyzb(l,7))**2. &
                      +(xyzb(l,5)-xyzb(l,8))**2. &
                      +(xyzb(l,6)-xyzb(l,9))**2. )
             d31=sqrt( (xyzb(l,7)-xyzb(l,1))**2. &
                      +(xyzb(l,8)-xyzb(l,2))**2. &
                      +(xyzb(l,9)-xyzb(l,3))**2. )
             sp=(d12+d23+d31)/2.
             area(l)=sqrt(sp*(sp-d12)*(sp-d23)*(sp-d31))

!             write(23,823) l,(bar(l,i),i=1,3)
823          format(i5,3(2x,e14.7))
          enddo
          close(23)

          return
          end subroutine readgeo3d
!...................................................................................
!***********************************************************************************
          subroutine readnormals(mpugeo,nxyz,nb,namefile)

          implicit none
          integer j,nb,i,mpugeo
          real nxyz(mpugeo,3)
          character(50)  namefile,string

          open(11,file=namefile)
          read(11,*)string
          do j=1,nb
             read(11,*)string,string,(nxyz(j,i),i=1,3)
             read(11,*)string
             read(11,*)string
             read(11,*)string
             read(11,*)string
             read(11,*)string
             read(11,*)string
824          format(i5,3(2x,e14.7))
          enddo
          read(11,*)string
          close(11)
          close(24)

          return
          end subroutine readnormals
!...................................................................................
