      subroutine mkmov_hdf_xcut

      use param
      use hdf5
      use mpih
      use mpi_param, only: kstart,kend
      use local_arrays, only: q2,q3,q1,pr
      use mgrd_arrays, only: dsalc
      use local_arrays, only: hsal,rusal,rhsr,dsal

      IMPLICIT none

      integer im,jm,km
      integer ic,jc,kc
      integer ip,jp,kp

      integer hdf_error
      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_q1v
      integer(HID_T) :: dset_q2v
      integer(HID_T) :: dset_q3v
      integer(HID_T) :: dset_prv
      integer(HID_T) :: dset_sal

      integer(HSIZE_T) :: dims(2)


      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(2) :: data_count  
      integer(HSSIZE_T), dimension(2) :: data_offset 

      integer :: comm, info
      integer :: ndims

      real tprfi
      real prx(n2m,n3m),v1(n2m,n3m),v2(n2m,n3m),v3(n2m,n3m)
      real dsalx(n2m,n3m)
      integer itime

      character(70) namfile,xdmnam
      character(5) ipfi


      call h5open_f(hdf_error)
      !RO   first of all, fill up velocity arrays
      !RO   so that we have node-centred velocities

      ndims=2
      ic=n1m/2
      ip=ic+1
      
      do kc=kstart,kend
       kp=kc+1 

       do jc=1,n2m                                                    
        jp=jc+1

        v1(jc,kc) = (q1(ic,jc,kc)+q1(ip,jc,kc))*0.5d0
        v2(jc,kc) = (q2(ic,jc,kc)+q2(ic,jp,kc))*0.5d0
        v3(jc,kc) = (q3(ic,jc,kc)+q3(ic,jc,kp))*0.5d0
        !prx(jc,kc)= pr(ic,jc,kc)
        dsalx(jc,kc) = dsalc(ic,jc,kc)
       end do
      end do

      !RO   File writing part

      !RO   Form the name of the file

      tprfi = 1.d0/tframe
      itime=nint(time*tprfi)
      write(ipfi,'(i5.5)')itime

      namfile='flowmov/frame_'//ipfi//'_xcut.h5'
      xdmnam='flowmov/frame_'//ipfi//'_xcut.xmf'

      !RO   Sort out MPI definitions and open file

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

      call h5fcreate_f(namfile, H5F_ACC_TRUNC_F, file_id, hdf_error,access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      !RO   Create dataspace

      dims(1)=n2m
      dims(2)=n3m
      call h5screate_simple_f(ndims, dims,filespace, hdf_error)
         

      !RO   Create the dataset with default properties.
      call h5dcreate_f(file_id, 'Vx', H5T_NATIVE_REAL,filespace,dset_q1v, hdf_error)
      call h5dcreate_f(file_id, 'Vy', H5T_NATIVE_REAL,filespace,dset_q2v, hdf_error)
      call h5dcreate_f(file_id, 'Vz', H5T_NATIVE_REAL,filespace,dset_q3v, hdf_error)
      !call h5dcreate_f(file_id, 'Pr', H5T_NATIVE_REAL,filespace,dset_prv, hdf_error)
      call h5dcreate_f(file_id, 'Sal', H5T_NATIVE_REAL,filespace,dset_sal, hdf_error)

      !RO   Set offsets and element counts
      data_count(1)=n2m
      data_count(2)=kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = kstart-1

      !RO   Create dataspace in memory

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      !RO   Select hyperslab  and then write it

      !v1
      call h5dget_space_f(dset_q1v, filespace, hdf_error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_q1v, H5T_NATIVE_DOUBLE,v1(1:n2m,kstart:kend), dims, hdf_error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      !v2
      call h5dget_space_f(dset_q2v, filespace, hdf_error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_q2v, H5T_NATIVE_DOUBLE,v2(1:n2m,kstart:kend), dims, hdf_error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      !v3
      call h5dget_space_f(dset_q3v, filespace, hdf_error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
       call h5dwrite_f(dset_q3v, H5T_NATIVE_DOUBLE, v3(1:n2m,kstart:kend), dims, hdf_error, file_space_id = filespace, mem_space_id = memspace,xfer_prp = plist_id)

!      !Pressure
!      call h5dget_space_f(dset_prv, filespace, hdf_error)
!      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
!      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
!      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
!      call h5dwrite_f(dset_prv, H5T_NATIVE_DOUBLE, prx(1:n2m,kstart:kend), dims, hdf_error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      !Dsal
      call h5dget_space_f(dset_sal, filespace, hdf_error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
      call h5dwrite_f(dset_sal, H5T_NATIVE_DOUBLE, dsalx(1:n2m,kstart:kend), dims,  hdf_error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)


      !RO   Close properties and file

      call h5dclose_f(dset_q1v, hdf_error)
      call h5dclose_f(dset_q2v, hdf_error)
      call h5dclose_f(dset_q3v, hdf_error)
      !call h5dclose_f(dset_prv, hdf_error)
      call h5dclose_f(dset_sal, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5pclose_f(plist_id, hdf_error)
      call h5fclose_f(file_id, hdf_error)
       
      if (myid.eq.0) then

      open(45,file=xdmnam,status='unknown')
      rewind(45)
      write(45,'("<?xml version=""1.0"" ?>")')
      write(45,'("<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>")')
      write(45,'("<Xdmf Version=""2.0"">")')
      write(45,'("<Domain>")')
      write(45,'("<Grid Name=""xcut"" GridType=""Uniform"">")')
      write(45,'("<Topology TopologyType=""2DSMesh"" NumberOfElements=&
     &""",i4," ",i4,"""/>")') n3m,n2m
      write(45,'("<Geometry GeometryType=""X_Y"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n2m
      write(45,'("cordin_xcut_info.h5:/y")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n2m
      write(45,'("cordin_xcut_info.h5:/z")')
      write(45,'("</DataItem>")')
      write(45,'("</Geometry>")')
      write(45,'("<Attribute Name=""X-velcoity""&
     & AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n2m
      write(45,'("frame_",i5.5,"_xcut.h5:/Vx")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""Y-Velocity""&
     & AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n2m
      write(45,'("frame_",i5.5,"_xcut.h5:/Vy")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""Z-Velocity""&
     & AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n2m
      write(45,'("frame_",i5.5,"_xcut.h5:/Vz")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
!      write(45,'("<Attribute Name=""Pressure""&
!     & AttributeType=""Scalar"" Center=""Node"">")')
!      write(45,'("<DataItem Dimensions=""",i4," ",i4,"""&
!     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n2m
!      write(45,'("frame_",i5.5,"_xcut.h5:/Pr")') itime
!      write(45,'("</DataItem>")')
!      write(45,'("</Attribute>")')
      write(45,'("<Attribute Name=""Salinity""&
     & AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m,n2m
      write(45,'("frame_",i5.5,"_xcut.h5:/Sal")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')

      write(45,'("<Time Value=""",e12.5""" />")')time
      write(45,'("</Grid>")')
      write(45,'("</Domain>")')
      write(45,'("</Xdmf>")')
      close(45)

      end if


      call h5close_f(hdf_error)

      return
      end
