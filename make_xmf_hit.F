      program make_xmf
      use hdf5
      implicit none
      character*30 filnamk, filnamx
      character*5 ipfi
      character*30 filnamgrid

      integer :: n1m, n1, n2m, n2, n3m, n3
      integer :: i
      real    :: dx1, dx2, dx3
      real, allocatable, dimension(:) :: tc, tm
      logical :: fexist


      integer hdf_error
      integer(HSIZE_T) :: dims_grid(1)
      integer(HID_T) :: file_grid
      integer(HID_T) :: dset_grid
      integer(HID_T) :: dspace_grid

      n1m = 624
      n1= n1m+1

      n2m = 624
      n2 = n2m+1

      n3m = 1872
      n3 = n3m+1

      dx1 = 1./real(n1m)
      dx2 = 1./real(n2m)
      dx3 = 1./real(n3m)

      allocate(tc(1:n1),tm(1:n1))

!     MAKE COORDINATE SYSTEM

      do i=1,n1                                                    
        tc(i) = (real(i-1))*dx1 
      end do

      do i=1,n1m                                                     
        tm(i) = (real(i-1)+0.5d0)*dx1 
      end do

!RO   Write reference for opening in visit

      call h5open_f(hdf_error)

!RO   Write the grid information 

      filnamgrid = 'continua_master.h5'
      call h5fcreate_f(filnamgrid,H5F_ACC_TRUNC_F, file_grid, hdf_error)

      dims_grid(1)=n1m
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'xm', H5T_NATIVE_DOUBLE,
     &                dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, tc(1:n1),
     &        dims_grid, hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5fclose_f(file_grid, hdf_error)



!RO   Write xmf file



      do i=1,100000

       write(ipfi,82)i
   82  format(i5.5)
       filnamk = './ww_'//ipfi//'.h5'
       inquire(file=filnamk, exist=fexist)
       if(fexist) then 
 
        filnamx = './ww_'//ipfi//'.xmf'
        write(*,*) filnamx 
        open(45,file=filnamx,status='unknown')
        rewind(45)
      write(45,'("<?xml version=""1.0"" ?>")')
      write(45,'("<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>")')
      write(45,'("<Xdmf Version=""2.0"">")')
      write(45,'("<Domain>")')
      write(45,'("<Grid Name=""VortexR"" GridType=""Uniform"">")')
      write(45,'("<Topology TopologyType=""3DRectMesh"" 
     &NumberOfElements=""",i4," ",i4," ",i4,"""/>")') n3m,n2m,n1m
      write(45,'("<Geometry GeometryType=""VXVYVZ"">")')
      write(45,'("<DataItem Dimensions=""",i4,"""
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1m
      write(45,'("continua_master.h5:/xm")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,"""
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2m
      write(45,'("continua_master.h5:/ym")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,"""
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m
      write(45,'("continua_master.h5:/zm")')
      write(45,'("</DataItem>")')
      write(45,'("</Geometry>")')
      write(45,'("<Attribute Name=""Vorticity modulus""
     & AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,"""
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')
     &n3m,n2m,n1m
      write(45,'("ww_",i5.5,".h5:/Ww")') i
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
       write(45,'("</Grid>")')
       write(45,'("</Domain>")')
       write(45,'("</Xdmf>")')
       close(45)

       end if

      end do
    

      stop
      end


