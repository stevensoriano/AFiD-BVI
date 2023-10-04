      program make_xmf
      use hdf5
      implicit none
      character*30 filnamk, filnamx
      character*5 ipfi
      character*30 filnamgrid

      character(len=4) :: dummy
      integer :: nxm, nym, nzm, na
      integer :: i, j, k
      real    :: xlen,ylen,zlen
      real    :: x2, x3, x1
      real, allocatable, dimension(:) :: xm,ym,zm
      logical :: fexist


      integer hdf_error
      integer(HSIZE_T) :: dims_grid(1)
      integer(HID_T) :: file_grid
      integer(HID_T) :: dset_grid
      integer(HID_T) :: dspace_grid



      open(unit=15,file='bou.in',status='old')
        read(15,301) dummy
        read(15,*) nxm,nym,nzm,na,na,na
        read(15,301) dummy
        read(15,*) na,na,na,na,na,na
        read(15,301) dummy
        read(15,*) xlen,ylen,zlen
       close(15)
301     format(a4)                

      allocate(xm(1:nxm))
      allocate(ym(1:nym))
      allocate(zm(1:nzm))


!     MAKE COORDINATE SYSTEM

      do i=1,nxm
       x1=real(i-1)/real(nxm)
       xm(i)= xlen*x1
      end do

      do j=1,nym
       x2=real(j-1)/real(nym)
       ym(j)= ylen*x2
      end do

      do k=1,nzm
       x3=real(k-1)/real(nzm)
       zm(k)= zlen*x3
      end do


!     WRITE COORDINATE SYSTEM

!RO   Write reference for opening in visit

      call h5open_f(hdf_error)

!RO   Write the grid information 

      filnamgrid = 'continua_master.h5'
      call h5fcreate_f(filnamgrid,H5F_ACC_TRUNC_F, file_grid, hdf_error)

      dims_grid(1)=nxm
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'xm', H5T_NATIVE_DOUBLE,
     &                dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, xm(1:nxm),
     &        dims_grid, hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=nym
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)

      call h5dcreate_f(file_grid, 'ym', H5T_NATIVE_DOUBLE,
     &                dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, ym(1:nym),
     &       dims_grid,hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=nzm
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'zm', H5T_NATIVE_DOUBLE,
     &                dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zm(1:nzm),
     &        dims_grid, hdf_error)


      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5fclose_f(file_grid, hdf_error)



!RO   Write xmf file



      do i=1,99999

       write(ipfi,82)i
   82  format(i5.5)
       filnamk = './flowmov/ww_'//ipfi//'.h5'
       inquire(file=filnamk, exist=fexist)
       if(fexist) then 
 
        filnamx = './flowmov/ww_'//ipfi//'.xmf'
        write(*,*) filnamx 
        open(45,file=filnamx,status='unknown')
        rewind(45)
      write(45,'("<?xml version=""1.0"" ?>")')
      write(45,'("<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>")')
      write(45,'("<Xdmf Version=""2.0"">")')
      write(45,'("<Domain>")')
      write(45,'("<Grid Name=""VortexR"" GridType=""Uniform"">")')
      write(45,'("<Topology TopologyType=""3DRectMesh"" 
     &NumberOfElements=""",i4," ",i4," ",i4,"""/>")') nzm,nym,nxm
      write(45,'("<Geometry GeometryType=""VXVYVZ"">")')
      write(45,'("<DataItem Dimensions=""",i4,"""
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')nxm
      write(45,'("continua_master.h5:/xm")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,"""
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')nym
      write(45,'("continua_master.h5:/ym")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,"""
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')nzm
      write(45,'("continua_master.h5:/zm")')
      write(45,'("</DataItem>")')
      write(45,'("</Geometry>")')
      write(45,'("<Attribute Name=""Vorticity modulus""
     & AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,"""
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')
     &nzm,nym,nxm
      write(45,'("ww_",i5.5,".h5:/ww")') i
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


