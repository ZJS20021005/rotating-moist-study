      subroutine mkmov_qvap

      use local_arrays, only: qvap
      use mpi_param
      use mpih
      use hdf5
      use param

      IMPLICIT NONE

      integer hdf_error
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: file_qvapv

      integer(HID_T) :: dset_qvapv

      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: file_plist
      integer(HID_T) :: slab_plist
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info

      integer ndims,itime

      real :: tprfi
      character(70) filnam1,filnamxdm
      character(5) ipfi
      character(len=256) :: line
      integer :: ios
      character(len=256) :: iomsg

      !RO   Sort out MPI definitions and file names

      tprfi = 1/tframe
      itime=nint(time*tprfi)
      write(ipfi,'(i5.5)')itime

      filnam1='movie/qvap'//ipfi//'.h5'
      filnamxdm = 'movie/qvap'//ipfi//'.xmf' 

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

      !RO   Set offsets and element counts

      ndims=3

      dims(1)=n1mr
      dims(2)=n2mr
      dims(3)=n3mr

      data_count(1) = n1mr
      data_count(2) = n2mr
      data_count(3) = kendr-kstartr+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = kstartr-1 


      call h5open_f(hdf_error)

      !RO   Set up MPI file properties
      call h5pcreate_f(H5P_FILE_ACCESS_F, file_plist, hdf_error)
      call h5pset_fapl_mpio_f(file_plist, comm, info, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, slab_plist, hdf_error) 
      call h5pset_dxpl_mpio_f(slab_plist, H5FD_MPIO_COLLECTIVE_F, hdf_error)

      !RO   Create dataspace
      call h5screate_simple_f(ndims, dims, filespace, hdf_error)

      !RO   Create dataspace in memory
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      !RO   Open first continua file for qvap
      call h5fcreate_f(filnam1, H5F_ACC_TRUNC_F, file_qvapv, hdf_error, access_prp=file_plist)

      !RO   Create dataset on file

      call h5dcreate_f(file_qvapv, 'S', H5T_NATIVE_DOUBLE,filespace, dset_qvapv, hdf_error)

      !RO   Set hyperslab
      call h5dget_space_f(dset_qvapv, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)

      call h5dwrite_f(dset_qvapv, H5T_NATIVE_DOUBLE,&
     &   qvap(1:n1mr,1:n2mr,kstartr:kendr), dims, &
     &   hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
     &   xfer_prp = slab_plist)

      !RO   Close dataset and file for qvap

      call h5dclose_f(dset_qvapv, hdf_error)
      call h5fclose_f(file_qvapv, hdf_error)

      !RO   Close all other stuff

      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5pclose_f(file_plist, hdf_error)
      call h5pclose_f(slab_plist, hdf_error)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      !EP   Write the xdm

      if (myid.eq.0) then

      open(45,file=filnamxdm,status='replace',action='write',iostat=ios,iomsg=iomsg)
      if (ios .ne. 0) then
        write(*,'(a,1x,a,1x,i0,1x,a)') 'Error opening XMF file:', trim(filnamxdm), ios, trim(iomsg)
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      endif
      rewind(45)
      write(45,'("<?xml version=""1.0"" ?>")')
      write(45,'("<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>")')
      write(45,'("<Xdmf Version=""2.0"">")')
      write(45,'("<Domain>")')
      write(45,'("<Grid Name=""RB Cartesian"" GridType=""Uniform"">")')
      ! Build and write the topology line
      write(line,'(a,i4,a,i4,a,i4,a)') '<Topology TopologyType="3DRectMesh" NumberOfElements="', n3mr, ' ', n2mr, ' ', n1mr, '"/>'
      write(45,'(a)') trim(line)
      ! Geometry
      write(45,'(a)') '<Geometry GeometryType="VXVYVZ">'
      ! DataItem for x
      write(line,'(a,i4,a)') '<DataItem Dimensions="', n1mr, '" NumberType="Float" Precision="4" Format="HDF">'
      write(45,'(a)') trim(line)
      write(45,'(a)') 'cordin_info.h5:/x'
      write(45,'(a)') '</DataItem>'
      ! DataItem for y
      write(line,'(a,i4,a)') '<DataItem Dimensions="', n2mr, '" NumberType="Float" Precision="4" Format="HDF">'
      write(45,'(a)') trim(line)
      write(45,'(a)') 'cordin_info.h5:/y'
      write(45,'(a)') '</DataItem>'
      ! DataItem for z
      write(line,'(a,i4,a)') '<DataItem Dimensions="', n3mr, '" NumberType="Float" Precision="4" Format="HDF">'
      write(45,'(a)') trim(line)
      write(45,'(a)') 'cordin_info.h5:/z'
      write(45,'(a)') '</DataItem>'
      write(45,'(a)') '</Geometry>'
      write(45,'(a)') '<Attribute Name="Vapor" AttributeType="Scalar" Center="Node">'
      ! Attribute DataItem dimensions (n3mr n2mr n1mr)
      write(line,'(a,i4,a,i4,a,i4,a)') '<DataItem Dimensions="', n3mr, ' ', n2mr, ' ', n1mr, '" NumberType="Float" Precision="4" Format="HDF">'
      write(45,'(a)') trim(line)
      write(45,'(a)') trim('qvap'//adjustl(ipfi)//'.h5:/S')
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Time Value=""",e12.5,""" />")') time
      write(45,'("</Grid>")')
      write(45,'("</Domain>")')
      write(45,'("</Xdmf>")')
      close(45)

      endif

      call h5close_f(hdf_error)

      return                                                          
      end                                                             
