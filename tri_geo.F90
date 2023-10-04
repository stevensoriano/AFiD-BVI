!-------------------------------------------------------
!     reading in gts file and storing geometry data
!-------------------------------------------------------
      subroutine tri_geo
      
      USE param
      USE mpih
      USE mpi_param, only: kstart, kend
      USE mls_param
      IMPLICIT NONE

      character*50 filename, strucfilename
      character*5 ipfi, ipfip

      integer :: i,j,k,v1,v2,v3,inp,e1,e2,e3,count
      real :: totmass
      real :: xyz_tr(3,Nparticle)
      integer,dimension(:),allocatable :: nv,ne,nf
      real :: xmax,xmin,ymax,ymin,zmax,zmin
      real :: usVolbub

      filename = gtsfx

      maxnv=-1 ; maxne=-1 ; maxnf=-1
      xmax = -1 ; ymax = -1 ; zmax = -1
      xmin = 1e3 ; ymin = 1e3 ; zmin = 1e3


      allocate(nv(Nparticle),ne(Nparticle),nf(Nparticle))


!VS     Reading in number of vertices, edges and faces from first line of the gts file
!VS     Use different filenames if there are different objects - for now monodispersed

      do inp=1,Nparticle

        open(109,file=filename)
        read(109,*)nv(inp),ne(inp),nf(inp)
        maxnv=max(maxnv,nv(inp))
        maxne=max(maxne,ne(inp))
        maxnf=max(maxnf,nf(inp))

        if(myid.eq.0.and.inp.eq.1)then
          write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++++'
          write(6,*)'reading mesh',filename
          write(6,*)'N vertices,N edges,N faces',maxnv,maxne,maxnf
          write(6,*)'Total number of particles ',Nparticle
          write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++++'
        end if

        close(109)
      end do

      usVolbub = 0.75/(pi*(sclf**3))
      usVolele = float(maxnf)*usVolbub

!     ------------------------------------------------------------------

!     Allocations using vertices,edges and faces as parameters     

      call allocate_trigeo


!     Read in centres of particles from a input file
!     Trivial method for now - has to be changed
      do inp=1,Nparticle
       open(108,file='spos.in')
       read(108,*)xyz_tr(1,inp),xyz_tr(2,inp),xyz_tr(3,inp)
      end do


    

!     elastic constants for the springs
!     Read in from input file
      
       thickness = 1.0
       rhos = 1.0
!     Read in the vertices, edges and faces from the file
!     Also make connections between the three

      do inp=1,Nparticle

        open(109,file=filename)
        read(109,*)nv(inp),ne(inp),nf(inp)
        
        do i=1,maxnv
          read(109,*)xyz0(1,i,inp),xyz0(2,i,inp),xyz0(3,i,inp)
          xmin=min(xmin,xyz0(1,i,inp))
          ymin=min(ymin,xyz0(2,i,inp))
          zmin=min(zmin,xyz0(3,i,inp))
          xmax=max(xmax,xyz0(1,i,inp))
          ymax=max(ymax,xyz0(2,i,inp))
          zmax=max(zmax,xyz0(3,i,inp))
        end do
      
        if(myid.eq.0.and.inp.eq.1) then
        write(6,*) 'Bounding box of object before scaling',inp
        write(6,*) 'xmin = ',xmin, ' xmax =', xmax
        write(6,*) 'ymin = ',ymin, ' ymax =', ymax
        write(6,*) 'zmin = ',zmin, ' zmax =', zmax
        endif

!      Translate the sphere from 0,0,0 and scale radius

!     Scaling factor

       do i=1,maxnv
         xyz0(1,i,inp) = xyz0(1,i,inp)*sclf + xyz_tr(1,inp)
         xyz0(2,i,inp) = xyz0(2,i,inp)*sclf + xyz_tr(2,inp)
         xyz0(3,i,inp) = xyz0(3,i,inp)*sclf + xyz_tr(3,inp)
       end do

       xmax = -1 ; ymax = -1 ; zmax = -1
       xmin = 1e3 ; ymin = 1e3 ; zmin = 1e3

        do i=1,maxnv
          xmin=min(xmin,xyz0(1,i,inp))
          ymin=min(ymin,xyz0(2,i,inp))
          zmin=min(zmin,xyz0(3,i,inp))
          xmax=max(xmax,xyz0(1,i,inp))
          ymax=max(ymax,xyz0(2,i,inp))
          zmax=max(zmax,xyz0(3,i,inp))
        end do

        if(myid.eq.0.and.inp.eq.1) then
        write(6,*) 'Bounding box of object after scaling',inp
        write(6,*) 'xmin = ',xmin, ' xmax =', xmax
        write(6,*) 'ymin = ',ymin, ' ymax =', ymax
        write(6,*) 'zmin = ',zmin, ' zmax =', zmax

        write(6,*) 'Average grid length ',1.0/float(n1m)
        endif

!     ---------------------------------------------------------------
!     Read from structural file to load in previous particle positions
      if(pread.eq.1)then
      write(ipfip,43)inp
   43 format(i5.5)
      ipfi = trim(datfx)

      strucfilename = 'contstr/str1_'//ipfip//'_'//ipfi//'.dat' 
           open(909,file=strucfilename)
           do i=1,maxnv
                 read(909,*)xyz(1,i,inp),xyz(2,i,inp),xyz(3,i,inp)
           end do
           do i=1,maxnv
                 read(909,*)xyzv(1,i,inp),xyzv(2,i,inp),xyzv(3,i,inp)
           end do
           do i=1,maxnv
                 read(909,*)xyza(1,i,inp),xyza(2,i,inp),xyza(3,i,inp)
           end do
           close(909)

      else

         xyz(1:3,1:maxnv,inp) = xyz0(1:3,1:maxnv,inp)
         xyzv(1:3,1:maxnv,inp) = 0.0
         xyza(1:3,1:maxnv,inp) = 0.0
      end if

      xyzv1(1:3,1:maxnv,inp) = 0.0
      xyza0(1:3,1:maxnv,inp) = 0.0
      fxyz(1:3,1:maxnv,inp) = 0.0

!      Finished translation and scaling        

        
        n_edge_of_vert(1:maxnv,inp)=0

        do i=1,maxne
          read(109,*)v1,v2
          vert_of_edge(1,i,inp)=v1
          vert_of_edge(2,i,inp)=v2
          n_edge_of_vert(v1,inp)=n_edge_of_vert(v1,inp)+1
          n_edge_of_vert(v2,inp)=n_edge_of_vert(v2,inp)+1
          vert_of_vert(n_edge_of_vert(v1,inp),v1,inp)=v2           
          vert_of_vert(n_edge_of_vert(v2,inp),v2,inp)=v1           
          edge_of_vert(n_edge_of_vert(v1,inp),v1,inp)=i
          edge_of_vert(n_edge_of_vert(v2,inp),v2,inp)=i
        enddo

        do i=1,maxnf
          read(109,*)edge_of_face(1,i,inp),edge_of_face(2,i,inp), &
                     edge_of_face(3,i,inp)
        end do
 
        do i=1,maxnf
           e1=edge_of_face(1,i,inp)
           e2=edge_of_face(2,i,inp)

           if (vert_of_edge(2,e1,inp).eq.vert_of_edge(1,e2,inp)) then
              v1=vert_of_edge(1,e1,inp)
              v2=vert_of_edge(2,e1,inp)
              v3=vert_of_edge(2,e2,inp)
           elseif(vert_of_edge(2,e1,inp).eq.vert_of_edge(2,e2,inp)) then
              v1=vert_of_edge(1,e1,inp)
              v2=vert_of_edge(2,e1,inp)
              v3=vert_of_edge(1,e2,inp)
           elseif(vert_of_edge(1,e1,inp).eq.vert_of_edge(1,e2,inp)) then
              v1=vert_of_edge(2,e1,inp)
              v2=vert_of_edge(1,e1,inp)
              v3=vert_of_edge(2,e2,inp)
           else 
              v1=vert_of_edge(2,e1,inp)
              v2=vert_of_edge(1,e1,inp)
              v3=vert_of_edge(1,e2,inp)
           endif 

              vert_of_face(1,i,inp)=v1
              vert_of_face(2,i,inp)=v2
              vert_of_face(3,i,inp)=v3
           enddo

          close(109)
!VS   Completed reading the gts file

!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Check - vertex cannot be connected to itself
           do i=1,maxnv
              do j=1,n_edge_of_vert(i,inp)
                 if (vert_of_vert(j,i,inp).eq.i)  &
                        write(*,*)'Error ',vert_of_vert(j,i,inp),i
              enddo
           enddo
!     Check - for faces and edges
           face_of_edge(1:2,1:maxne,inp)=0
           do i=1,maxnf
              e1=edge_of_face(1,i,inp)
              e2=edge_of_face(2,i,inp)
              e3=edge_of_face(3,i,inp)
              if (face_of_edge(1,e1,inp).eq.0) then
                 face_of_edge(1,e1,inp)=i
              elseif (face_of_edge(2,e1,inp).eq.0) then
                 face_of_edge(2,e1,inp)=i
              else
                 write(*,*)'Edge error', i,e1,e2,e3
              endif
              if (face_of_edge(1,e2,inp).eq.0) then
                 face_of_edge(1,e2,inp)=i
              elseif (face_of_edge(2,e2,inp).eq.0) then
                 face_of_edge(2,e2,inp)=i
              else
                 write(*,*)'Edge error'
              endif
              if (face_of_edge(1,e3,inp).eq.0) then
                 face_of_edge(1,e3,inp)=i
              elseif (face_of_edge(2,e3,inp).eq.0) then
                 face_of_edge(2,e3,inp)=i
              else
                 write(*,*)'Edge error'
              endif
           enddo 

           !Check
           count=0
           do i=1,maxne
              if (face_of_edge(1,i,inp).eq.face_of_edge(2,i,inp)) then
                 write(*,*)'Error on edges '
              endif
              if (face_of_edge(1,i,inp).eq.0.or.face_of_edge(2,i,inp).eq.0) then
                 count=count+1
              endif
           enddo
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        call calculate_normal(maxnv,maxnf, &
                       xyz0(:,:,inp),vert_of_face(:,:,inp), &
                             tri_nor(:,:,inp))

        call calculate_angle(theta0(:,inp),maxne,maxnf, &
                        face_of_edge(:,:,inp),tri_nor(:,:,inp))

        call calculate_distance(dist0(:,inp),maxnv,maxne,xyz0(:,:,inp) &
                              ,vert_of_edge(:,:,inp))

        if(myid.eq.0.and.inp.eq.1)write(*,*)'Min edge length', &
                                                 minval(dist0(:,inp))
        if(myid.eq.0.and.inp.eq.1)write(*,*)'Max edge length', &
                                                 maxval(dist0(:,inp))
        if(myid.eq.0.and.inp.eq.1)write(*,*)'Average edge length', &
                                                 sum(dist0(:,inp))/size(dist0(:,inp))

        call calculate_area(Surface0(inp),maxnv,maxnf,xyz0(:,:,inp), &
                        vert_of_face(:,:,inp),sur0(:,inp))
        if(myid.eq.0.and.inp.eq.1)write(*,*)'Total surface: ', &
                                                          Surface0(inp)
        sur(:,inp)=sur0(:,inp)

        call calculate_volume(Volume0(inp),maxnv,maxnf,xyz0(:,:,inp), &
                              vert_of_face(:,:,inp))
        if(myid.eq.0.and.inp.eq.1)write(*,*)'Volume: ', &
                                                            Volume0(inp)

        call find_quartet(v1234(:,:,inp),maxnv,maxne,maxnf,xyz0(:,:,inp) &
                          ,face_of_edge(:,:,inp),edge_of_face(:,:,inp), &
                          vert_of_edge(:,:,inp),tri_nor(:,:,inp))

        if(myid.eq.0.and.inp.eq.1)write(*,*)'Finished reading geo.'


      totmass=Surface0(inp)*thickness*rhos
      uspm = 1.0

      if(myid.eq.0.and.inp.eq.1)print*,'Point mass ',uspm
 
      call convert_geo(maxnv,maxne,maxnf,xyz(:,:,inp),xyzv(:,:,inp), &
           xyza(:,:,inp),vert_of_face(:,:,inp), &
           tri_ver(:,:,inp),tri_bar(:,:,inp), &
           vel_tri(:,:,inp),acc_tri(:,:,inp))

      end do
      
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 

!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        return
        end subroutine tri_geo

