!-------------------------------------------------
!     write structure data to vtk file
!-------------------------------------------------

      subroutine write_to_vtk
      use hdf5
      use param
      use mls_param
      use mpih
      implicit none
     
      integer :: i,inp,itime,ntr,v1,v2,v3
      integer, dimension(maxnf) :: cell_type      
!      real, dimension(maxnf,Nparticle) :: sca_nod1,sca_nod2,sca_nod3
      real :: tprfi
      character*70 filname
      character*5 ipfi,ipfip

      cell_type(:) = 5 !5 is for triangular elements
      

!      do inp=1,Nparticle
!       do ntr=1,maxnf

!       v1 = vert_of_face(1,ntr,inp)
!       v2 = vert_of_face(2,ntr,inp)
!       v3 = vert_of_face(3,ntr,inp)

!       sca_nod1(ntr,inp)=(fpxyz(1,v1,inp)+fpxyz(1,v2,inp)+ &
!                        fpxyz(1,v3,inp))/(3.0*sur(ntr,inp))
!       sca_nod2(ntr,inp)=(fpxyz(2,v1,inp)+fpxyz(2,v2,inp)+ &
!                        fpxyz(2,v3,inp))/(3.0*sur(ntr,inp))
!       sca_nod3(ntr,inp)=(fpxyz(3,v1,inp)+fpxyz(3,v2,inp)+ &
!                        fpxyz(3,v3,inp))/(3.0*sur(ntr,inp))

!        sca_nod1(ntr,inp) = press_face(ntr,inp)+vforce_face(ntr,inp)

!       end do
!      end do

      tprfi = 1000.
      itime = nint(time*tprfi)
      write(ipfi,93) itime
   93 format(i5.5)

      do inp=1,Nparticle

      write(ipfip,94)inp
   94 format(i5.5)

      filname = 'vtkfiles/struc_'//ipfip//'_'//ipfi//'.vtk'

      open(121,file = filname)
!Header     
        write(121,'(a)')adjustl('# vtk DataFile Version 3.1')
        write(121,'(a)')adjustl('Stores triangular gts mesh')
        write(121,'(a)')adjustl('ASCII')
        write(121,'(a)')adjustl('DATASET UNSTRUCTURED_GRID')
        write(121,*)''
        write(121,*)'POINTS ',maxnv,' FLOAT'
      do i=1,maxnv
        write(121,*)xyz(1:3,i,inp)
      end do
        write(121,*)''
        write(121,*)'CELLS ',maxnf, 4*maxnf
      do i=1,maxnf
        write(121,*)'3 ',vert_of_face(1:3,i,inp)-1
      end do
        write(121,*)''
        write(121,*)'CELL_TYPES ',maxnf
        write(121,*)cell_type(:)
        
        write(121,*)''
        write(121,*)'CELL_DATA ',maxnf
        write(121,*)'Scalars Ex-For1 FLOAT'
        write(121,*)'LOOKUP_TABLE default'
        write(121,*)sca_nod1(:,inp)
        write(121,*)'Scalars Ex-For2 FLOAT'
        write(121,*)'LOOKUP_TABLE default'
        write(121,*)sca_nod2(:,inp)
         write(121,*)'Scalars Ex-For3 FLOAT'
        write(121,*)'LOOKUP_TABLE default'
        write(121,*)sca_nod3(:,inp)
     
      close(121)

      end do


      return
      end
