!----------------------------------------------------------------------
!     Auxillary routines for triangulated geometry
!----------------------------------------------------------------------
        subroutine calculate_normal(nv,nf,xyz, &
                                   vert_of_face, tri_nor)

        implicit none
        integer::nv,nf,v1,v2,v3,i
        integer,dimension(3,nf)::vert_of_face
        real,dimension(3,nf):: face_normal
        real,dimension(3,nv)::xyz
        real,dimension(3)::ve1,ve2
        real,dimension(3,nf)::tri_nor

      do i=1,nf
           v1=vert_of_face(1,i)
           v2=vert_of_face(2,i)
           v3=vert_of_face(3,i)

           ve1(1:3)=xyz(1:3,v2)-xyz(1:3,v1)
           ve2(1:3)=xyz(1:3,v3)-xyz(1:3,v1)

           face_normal(1,i)=ve1(2)*ve2(3) - ve1(3)*ve2(2)
           face_normal(2,i)=ve1(3)*ve2(1) - ve1(1)*ve2(3)
           face_normal(3,i)=ve1(1)*ve2(2) - ve1(2)*ve2(1)

         tri_nor(1:3,i)=face_normal(1:3,i)/ &
                                    sqrt(sum(face_normal(1:3,i)**2))
        enddo

        return
        end subroutine calculate_normal
!------------------------------------------------------
        subroutine calculate_volume (Volume,nv,nf,xyz,vert_of_face)

        implicit none
        integer :: nv,nf,v1,v2,v3,i
        integer, dimension (3,nf) :: vert_of_face
        real, dimension (3,nv) ::xyz
        real :: Volume

        Volume=0.0
        do i=1,nf
           v1=vert_of_face(1,i)
           v2=vert_of_face(2,i)
           v3=vert_of_face(3,i)

        Volume = Volume + (xyz(1,v1)*(xyz(2,v2)*xyz(3,v3)- &
                           xyz(3,v2)*xyz(2,v3)) +   &
                          xyz(1,v2)*(xyz(2,v3)*xyz(3,v1)- &
                           xyz(3,v3)*xyz(2,v1)) +      &
                          xyz(1,v3)*(xyz(2,v1)*xyz(3,v2)- &
                           xyz(3,v1)*xyz(2,v2)))
        enddo
        Volume=Volume/6.

        return
        end subroutine calculate_volume
!------------------------------------------------------
        subroutine calculate_area (Surface,nv,nf,xyz,vert_of_face,sur)

        implicit none
        integer :: nv,nf,v1,v2,v3,i
        integer, dimension (3,nf) :: vert_of_face
        real, dimension (3,nv) ::xyz
        real, dimension (nf) :: sur
        real :: Surface,d12,d23,d31,sp

        Surface=0.0
        do i=1,nf
           v1=vert_of_face(1,i)
           v2=vert_of_face(2,i)
           v3=vert_of_face(3,i)

           d12=sqrt( (xyz(1,v1)-xyz(1,v2))**2. &
                    +(xyz(2,v1)-xyz(2,v2))**2. &
                    +(xyz(3,v1)-xyz(3,v2))**2. )
           d23=sqrt( (xyz(1,v2)-xyz(1,v3))**2. &
                    +(xyz(2,v2)-xyz(2,v3))**2. &
                    +(xyz(3,v2)-xyz(3,v3))**2. )
           d31=sqrt( (xyz(1,v3)-xyz(1,v1))**2. &
                    +(xyz(2,v3)-xyz(2,v1))**2. &
                    +(xyz(3,v3)-xyz(3,v1))**2. )
           sp=(d12+d23+d31)/2.
           sur(i)=sqrt(sp*(sp-d12)*(sp-d23)*(sp-d31))
           Surface=Surface+sur(i)
        enddo

        return
        end subroutine calculate_area
!------------------------------------------------------
        subroutine calculate_angle(theta,ne,nf,face_of_edge,face_normal)

        implicit none
        integer :: ne,nf,f1,f2,i
        integer, dimension (2,ne) :: face_of_edge
        real, dimension (3,nf) :: face_normal
        real, dimension (ne) :: theta
        real, dimension (3) :: n1,n2
        real :: pvx,pvy,pvz,tmp,PI,th,tmp1,tmp2

        PI=acos(-1.0)

        theta(1:ne)=0.0
        do i=1,ne
           f1=face_of_edge(1,i)
           f2=face_of_edge(2,i)
           if (f1.ne.0.and.f2.ne.0) then
              n1(1:3)=face_normal(1:3,f1)
              n2(1:3)=face_normal(1:3,f2)
              pvx = n1(2)*n2(3) - n1(3)*n2(2)
              pvy = n1(3)*n2(1) - n1(1)*n2(3)
              pvz = n1(1)*n2(2) - n1(2)*n2(1)

              tmp1 = sqrt(pvx*pvx+pvy*pvy+pvz*pvz)
              tmp2 = (n1(1)*n2(1)+n1(2)*n2(2)+n1(3)*n2(3))
              th = atan2(tmp1,tmp2)
              theta(i)=th
           endif
        enddo

        return
        end subroutine calculate_angle
!------------------------------------------------------
        subroutine calculate_distance (dist,nv,ne,xyz,vert_of_edge)

        implicit none
        integer :: nv,ne,v1,v2,i
        integer, dimension (2,ne) :: vert_of_edge
        real, dimension (3,nv) ::xyz
        real, dimension (ne) :: dist

        do i=1,ne
           v1=vert_of_edge(1,i)
           v2=vert_of_edge(2,i)
        dist(i)=sqrt((xyz(1,v1)-xyz(1,v2))**2+(xyz(2,v1)-xyz(2,v2))**2+ &
                     (xyz(3,v1)-xyz(3,v2))**2) 
        enddo

        return
        end subroutine calculate_distance
!------------------------------------------------------
        subroutine find_quartet (v1234,nv,ne,nf,xyz,face_of_edge, &
                              edge_of_face,vert_of_edge,face_normal)

        implicit none
        integer :: nv,ne,nf,i,f1,f2,e1,e2,e3,v,v1,v2,v3,v4,vtmp
        integer, dimension (2,ne) :: face_of_edge,vert_of_edge
        integer, dimension (3,nf) :: edge_of_face
        real, dimension (3,nf) :: face_normal
        real, dimension (3,nv) ::xyz
        integer, dimension (4,ne) :: v1234
        real, dimension (3) :: csi,zet,csi1,zet1,a21,a31,a24,a34
        real :: res

        do i=1,ne
           f1=face_of_edge(1,i)
           f2=face_of_edge(2,i)
           if (f1.ne.0.and.f2.ne.0) then
              ! get opposite vertex of edge
              e1=edge_of_face(1,f1)
              e2=edge_of_face(2,f1)
              e3=edge_of_face(3,f1)
              if (e1.eq.i) then
                 v=vert_of_edge(1,e2)
              if(v.ne.vert_of_edge(1,i).and.v.ne.vert_of_edge(2,i))then
                    v1=v
                 else
                    v1=vert_of_edge(2,e2)
              endif
              endif
              if (e2.eq.i) then
                 v=vert_of_edge(1,e1)
              if(v.ne.vert_of_edge(1,i).and.v.ne.vert_of_edge(2,i))then
                    v1=v
                 else
                    v1=vert_of_edge(2,e1)
              endif
              endif
              if (e3.eq.i) then
                 v=vert_of_edge(1,e2)
              if(v.ne.vert_of_edge(1,i).and.v.ne.vert_of_edge(2,i))then
                    v1=v
                 else
                    v1=vert_of_edge(2,e2)
              endif
              endif
              ! get opposite vertex of edge
              e1=edge_of_face(1,f2)
              e2=edge_of_face(2,f2)
              e3=edge_of_face(3,f2)
              if (e1.eq.i) then
                 v=vert_of_edge(1,e2)
              if(v.ne.vert_of_edge(1,i).and.v.ne.vert_of_edge(2,i))then
                    v4=v
                 else
                    v4=vert_of_edge(2,e2)
              endif
              endif
              if (e2.eq.i) then
                 v=vert_of_edge(1,e1)
              if(v.ne.vert_of_edge(1,i).and.v.ne.vert_of_edge(2,i))then
                    v4=v
                 else
                    v4=vert_of_edge(2,e1)
              endif
              endif
              if (e3.eq.i) then
                 v=vert_of_edge(1,e2)
              if(v.ne.vert_of_edge(1,i).and.v.ne.vert_of_edge(2,i))then
                    v4=v
                 else
                 v4=vert_of_edge(2,e2)
              endif
              endif
              ! Other two vertices, on edge
              v2=vert_of_edge(1,i)
              v3=vert_of_edge(2,i)

              csi(1:3)=face_normal(1:3,f1)
              zet(1:3)=face_normal(1:3,f2)

              a21(1:3)=xyz(1:3,v2)-xyz(1:3,v1)
              a31(1:3)=xyz(1:3,v3)-xyz(1:3,v1)
              a34(1:3)=xyz(1:3,v3)-xyz(1:3,v4)
              a24(1:3)=xyz(1:3,v2)-xyz(1:3,v4)

              call cross(csi1,a21,a31)
              call cross(zet1,a34,a24)

              call dot(res,csi,csi1)
              if (res.lt.0.0) then
                 vtmp=v2
                 v2=v3
                 v3=vtmp
              endif

              v1234(1,i)=v1
              v1234(2,i)=v2
              v1234(3,i)=v3
              v1234(4,i)=v4
           endif
        enddo

        return
        end subroutine find_quartet     
!     ---------------------------------------------------------------------
        subroutine convert_geo(nv,ne,nf,xyz,xyzv,xyza, &
                   vert_of_face,tri_ver,tri_bar, &
                   vel_tri,acc_tri)

        implicit none
        integer :: nv,ne,nf,v1,v2,v3,i
        real, dimension (3,nv) :: xyz,xyzv,xyza
        integer, dimension (3,nf) :: vert_of_face
        real, dimension(9,nf) :: tri_ver
        real, dimension(3,nf) :: tri_bar,vel_tri,acc_tri


        do i=1,nf
           v1=vert_of_face(1,i)
           v2=vert_of_face(2,i)
           v3=vert_of_face(3,i)

           tri_ver(1:3,i) = xyz(1:3,v1)
           tri_ver(4:6,i) = xyz(1:3,v2)
           tri_ver(7:9,i) = xyz(1:3,v3)

          ! Find triangle's velocity and acceleration 

           vel_tri(1,i)=(xyzv(1,v1)+xyzv(1,v2)+xyzv(1,v3))/3.
           vel_tri(2,i)=(xyzv(2,v1)+xyzv(2,v2)+xyzv(2,v3))/3.
           vel_tri(3,i)=(xyzv(3,v1)+xyzv(3,v2)+xyzv(3,v3))/3.
                      
           acc_tri(1,i)=(xyza(1,v1)+xyza(1,v2)+xyza(1,v3))/3.
           acc_tri(2,i)=(xyza(2,v1)+xyza(2,v2)+xyza(2,v3))/3.
           acc_tri(3,i)=(xyza(3,v1)+xyza(3,v2)+xyza(3,v3))/3.
        enddo

        ! Find triangles' baricentre

        do i=1,nf
          tri_bar(1,i)=(tri_ver(1,i)+tri_ver(4,i)+tri_ver(7,i))/3.
          tri_bar(2,i)=(tri_ver(2,i)+tri_ver(5,i)+tri_ver(8,i))/3.
          tri_bar(3,i)=(tri_ver(3,i)+tri_ver(6,i)+tri_ver(9,i))/3.
        enddo
        end subroutine convert_geo
!     ----------------------------------------------------------------
          subroutine dot(xy,x,y)

          implicit none
          real x(3),y(3),xy

          xy = x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
          return
          end subroutine dot
!     ----------------------------------------------------------------
          subroutine  sub(xmy,x,y)

          implicit none
          real x(3),y(3),xmy(3)

          xmy = x-y
          return
          end subroutine sub
!     ----------------------------------------------------------------
          subroutine  cross(xcy,x,y)
          implicit none
          real x(3),y(3),xcy(3)

          xcy(1) = x(2)*y(3)-y(2)*x(3)
          xcy(2) = x(3)*y(1)-y(3)*x(1)
          xcy(3) = x(1)*y(2)-y(1)*x(2)
          return
          end subroutine cross
!     ----------------------------------------------------------------

