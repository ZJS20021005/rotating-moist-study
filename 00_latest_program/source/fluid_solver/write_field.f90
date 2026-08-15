      subroutine write_grid_info
      use mpih
      use param
      use hdf5

      IMPLICIT none

      integer hdf_error

      integer(HID_T) :: file_grid
      integer(HID_T) :: dset_grid
      integer(HID_T) :: dspace_grid

      integer(HSIZE_T) :: dims_grid(1)

      character namfile*70

      call h5open_f(hdf_error)

      IF(myid.eq.0)THEN      

      namfile='field_gridc.h5'

      call h5fcreate_f(namfile,H5F_ACC_TRUNC_F,file_grid,hdf_error)

      dims_grid(1)=n1

      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'xf', H5T_NATIVE_DOUBLE,dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, xc(1:n1),dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'xc', H5T_NATIVE_DOUBLE,dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, xm(1:n1),dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=n2
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'yf', H5T_NATIVE_DOUBLE,dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, yc(1:n2),dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'yc', H5T_NATIVE_DOUBLE,dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, ym(1:n2),dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=n3m
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'zf', H5T_NATIVE_DOUBLE,dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zc(1:n3m),dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'zc', H5T_NATIVE_DOUBLE,dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zm(1:n3m),dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5fclose_f(file_grid, hdf_error)

      
      namfile='field_gridd.h5'

      call h5fcreate_f(namfile,H5F_ACC_TRUNC_F,file_grid,hdf_error)

      dims_grid(1)=n1r
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'xf', H5T_NATIVE_DOUBLE,dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, xcr(1:n1r),dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'xc', H5T_NATIVE_DOUBLE, dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, xmr(1:n1r),dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=n2r
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'yf', H5T_NATIVE_DOUBLE,dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, ycr(1:n2r),dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'yc', H5T_NATIVE_DOUBLE, dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, ymr(1:n2r), dims_grid,hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=n3mr
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'zf', H5T_NATIVE_DOUBLE,dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zcr(1:n3mr),dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'zc', H5T_NATIVE_DOUBLE,dspace_grid, dset_grid, hdf_error)
      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zmr(1:n3mr), dims_grid, hdf_error)
      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5fclose_f(file_grid, hdf_error)


      !  xmf file for continua data files.
      write(namfile,'(a)')'field_ms.xmf'
      open(45,file=namfile,status='unknown')
      rewind(45)
      write(45,'("<?xml version=""1.0"" ?>")')
      write(45,'("<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>")')
      write(45,'("<Xdmf Version=""2.0"">")')
      write(45,'("<Domain>")')

      write(45,'("<Grid Name=""DDC_u"" GridType=""Uniform"">")')
      write(45,'("<Topology TopologyType=""3DRectMesh"" &
     & NumberOfElements=""",i4," ",i4," ",i4,"""/>")') n3m,n2,n1
      write(45,'("<Geometry GeometryType=""VXVYVZ"">")')
      write(45,'("<DataItem Dimensions=""",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1
      write(45,'("field_gridc.h5:/xf")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2
      write(45,'("field_gridc.h5:/yc")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m
      write(45,'("field_gridc.h5:/zc")')
      write(45,'("</DataItem>")')
      write(45,'("</Geometry>")')
      write(45,'("<Attribute Name=""u""&
     & AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')&
     &n3m,n2,n1
      write(45,'("continua_q1.h5:/Vth")')
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Time Value=""",e12.5,""" />")') time
      write(45,'("</Grid>")')

      write(45,'("<Grid Name=""DDC_v"" GridType=""Uniform"">")')
      write(45,'("<Topology TopologyType=""3DRectMesh"" &
     & NumberOfElements=""",i4," ",i4," ",i4,"""/>")') n3m,n2,n1
      write(45,'("<Geometry GeometryType=""VXVYVZ"">")')
      write(45,'("<DataItem Dimensions=""",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1
      write(45,'("field_gridc.h5:/xc")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2
      write(45,'("field_gridc.h5:/yf")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m
      write(45,'("field_gridc.h5:/zc")')
      write(45,'("</DataItem>")')
      write(45,'("</Geometry>")')
      write(45,'("<Attribute Name=""v""&
     & AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')&
     &n3m,n2,n1
      write(45,'("continua_q2.h5:/Vr")')
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Time Value=""",e12.5,""" />")') time
      write(45,'("</Grid>")')

      write(45,'("<Grid Name=""DDC_w"" GridType=""Uniform"">")')
      write(45,'("<Topology TopologyType=""3DRectMesh"" &
     & NumberOfElements=""",i4," ",i4," ",i4,"""/>")') n3m,n2,n1
      write(45,'("<Geometry GeometryType=""VXVYVZ"">")')
      write(45,'("<DataItem Dimensions=""",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1
      write(45,'("field_gridc.h5:/xc")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2
      write(45,'("field_gridc.h5:/yc")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m
      write(45,'("field_gridc.h5:/zf")')
      write(45,'("</DataItem>")')
      write(45,'("</Geometry>")')
      write(45,'("<Attribute Name=""w""&
     & AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')&
     &n3m,n2,n1
      write(45,'("continua_q3.h5:/Vz")')
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Time Value=""",e12.5,""" />")') time
      write(45,'("</Grid>")')

      write(45,'("<Grid Name=""DDC_S"" GridType=""Uniform"">")')
      write(45,'("<Topology TopologyType=""3DRectMesh"" &
     & NumberOfElements=""",i4," ",i4," ",i4,"""/>")') n3mr,n2r,n1r
      write(45,'("<Geometry GeometryType=""VXVYVZ"">")')
      write(45,'("<DataItem Dimensions=""",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1r
      write(45,'("field_gridd.h5:/xc")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2r
      write(45,'("field_gridd.h5:/yc")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3mr
      write(45,'("field_gridd.h5:/zc")')
      write(45,'("</DataItem>")')
      write(45,'("</Geometry>")')
      write(45,'("<Attribute Name=""S""&
     & AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')&
     & n3mr,n2r,n1r
      write(45,'("continua_dsal.h5:/dsal")')
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Time Value=""",e12.5,""" />")') time
      write(45,'("</Grid>")')
      write(45,'("</Domain>")')
      write(45,'("</Xdmf>")')
      close(45) 

      ENDIF

      call h5close_f(hdf_error)

      return
      end subroutine write_grid_info

!================================================
      subroutine write_base
      use param
      use mpih
      use mpi_param, only: kstart,kend
      use local_arrays, only:q1, q2, q3
      use hdf5
      implicit none

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_phi

      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info
      integer :: ndims

      character filnam1*30
      character filnam2*30
      character filnam3*30

      call h5open_f(hdf_error)

      !RO   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

      !RO   Form the name of the file

      write(filnam1,'(a,i5.5,a)')'field_q1c_',int(ntime),'.h5'
      write(filnam2,'(a,i5.5,a)')'field_q2c_',int(ntime),'.h5'
      write(filnam3,'(a,i5.5,a)')'field_q3c_',int(ntime),'.h5'

      !RO   Set offsets and element counts
   
      ndims = 3

      dims(1)=n1
      dims(2)=n2
      dims(3)=n3m

      call h5screate_simple_f(ndims, dims,filespace, hdf_error)

      data_count(1) = n1
      data_count(2) = n2
      data_count(3) = kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = kstart-1

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

      call h5fcreate_f(filnam1, H5F_ACC_TRUNC_F, file_id, hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'uc', H5T_NATIVE_DOUBLE,filespace, dset_phi, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      call h5dget_space_f(dset_phi, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_phi, H5T_NATIVE_DOUBLE,&
     &                q1(:,:,kstart:kend), dims, &
     &                hdf_error, file_space_id = slabspace,&
     &                mem_space_id = memspace, &
     &                xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_phi, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

      call h5fcreate_f(filnam2, H5F_ACC_TRUNC_F, file_id,  hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'vc', H5T_NATIVE_DOUBLE, filespace, dset_phi, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)

      call h5dget_space_f(dset_phi, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
      call h5dwrite_f(dset_phi, H5T_NATIVE_DOUBLE,&
     &                q2(:,:,kstart:kend), dims,&
     &                hdf_error, file_space_id = slabspace,&
     &                mem_space_id = memspace,&
     &                xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_phi, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)

      call h5fcreate_f(filnam3, H5F_ACC_TRUNC_F, file_id,   hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'wc', H5T_NATIVE_DOUBLE, filespace, dset_phi, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)

      call h5dget_space_f(dset_phi, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_phi, H5T_NATIVE_DOUBLE,&
     &                q3(:,:,kstart:kend), dims,&
     &                hdf_error, file_space_id = slabspace,&
     &                mem_space_id = memspace,&
     &                xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_phi, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      call h5sclose_f(filespace, hdf_error)
      call h5close_f(hdf_error)

      return
      end subroutine write_base


!================================================
      subroutine write_fine
      use param
      use mpih
      use mpi_param, only: kstartr,kendr
      use mgrd_arrays, only:q1lr, q2lr, q3lr
      use hdf5
      implicit none

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_phi

      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info
      integer :: ndims

      character filnam1*30
      character filnam2*30
      character filnam3*30

      call h5open_f(hdf_error)

      !RO   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

      !RO   Form the name of the file

      write(filnam1,'(a,i5.5,a)')'field_q1r_',int(ntime),'.h5'
      write(filnam2,'(a,i5.5,a)')'field_q2r_',int(ntime),'.h5'
      write(filnam3,'(a,i5.5,a)')'field_q3r_',int(ntime),'.h5'


      !RO   Set offsets and element counts
   
      ndims = 3

      dims(1)=n1r
      dims(2)=n2r
      dims(3)=n3mr

      call h5screate_simple_f(ndims, dims,  filespace, hdf_error)

      data_count(1) = n1r
      data_count(2) = n2r
      data_count(3) = kendr-kstartr+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = kstartr-1

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam1, H5F_ACC_TRUNC_F, file_id, hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'ur', H5T_NATIVE_DOUBLE,filespace, dset_phi, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)

      call h5dget_space_f(dset_phi, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_phi, H5T_NATIVE_DOUBLE,&
     &                q1lr(:,:,kstartr:kendr), dims,&
     &                hdf_error, file_space_id = slabspace,&
     &                mem_space_id = memspace,&
     &                xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_phi, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)


      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam2, H5F_ACC_TRUNC_F, file_id, hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'vr', H5T_NATIVE_DOUBLE, filespace, dset_phi, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)

      call h5dget_space_f(dset_phi, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_phi, H5T_NATIVE_DOUBLE,&
     &                q2lr(:,:,kstartr:kendr), dims,&
     &                hdf_error, file_space_id = slabspace,&
     &                mem_space_id = memspace,&
     &                xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_phi, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)


      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam3, H5F_ACC_TRUNC_F, file_id, hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'wr', H5T_NATIVE_DOUBLE,filespace, dset_phi, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)

      call h5dget_space_f(dset_phi, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)
      call h5dwrite_f(dset_phi, H5T_NATIVE_DOUBLE,&
     &                q3lr(:,:,kstartr:kendr), dims,  &
     &                hdf_error, file_space_id = slabspace,&
     &                mem_space_id = memspace,&
     &                xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_phi, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)


      call h5sclose_f(filespace, hdf_error)
      call h5close_f(hdf_error)

      return
      end subroutine write_fine
