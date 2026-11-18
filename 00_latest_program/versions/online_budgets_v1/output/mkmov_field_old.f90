!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcWriteQ.F90                                 !
!    CONTAINS: subroutine CalcWriteQ                      !
!                                                         ! 
!    PURPOSE: Compute and write the 3D q-criteria field   !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      subroutine mkmov_field

      use local_arrays, only: q3,q2,q1,dens,dsal
      use mpi_param
      use mpih
      use hdf5
      use param
!      use stat_arrays
      implicit none

      real, allocatable, dimension(:,:,:) :: temp_mea
      real, allocatable, dimension(:,:,:) :: vz_mea
      real, allocatable, dimension(:,:,:) :: vy_mea
      real, allocatable, dimension(:,:,:) :: vx_mea
  
      real :: dvxx1,dvxx2,dvxx3
      real :: dvyx1,dvyx2,dvyx3
      real :: dvzx1,dvzx2,dvzx3
      real :: strn, omeg

      integer :: ic,jc,kc,ip,jp,kp,im,jm,km

      !CS   For h5 and xmf --------------------------
      integer hdf_error
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

!      integer(HID_T) :: file_qtens
      integer(HID_T) :: file_temp


!      integer(HID_T) :: dset_qtens
      integer(HID_T) :: dset_temp
      integer(HID_T) :: dset_VZ
      integer(HID_T) :: dset_VY
      integer(HID_T) :: dset_VX


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
      !----------------------------------------------

      call update_both_ghosts(n1,n2,q1,kstart,kend)
      call update_both_ghosts(n1,n2,q2,kstart,kend)
      call update_both_ghosts(n1,n2,q3,kstart,kend)
      call update_both_ghosts(n1,n2,dsal,kstart,kend)
     ! call update_both_ghosts(n1,n2,vx_me,kstart,kend)



      allocate(temp_mea(1:n1m,1:n2m,kstart:kend))
      allocate(vz_mea(1:n1m,1:n2m,kstart:kend))
      allocate(vy_mea(1:n1m,1:n2m,kstart:kend))
      allocate(vx_mea(1:n1m,1:n2m,kstart:kend))

      !nstatsamples = nstatsamples + 1.0
      
      do kc=kstart,kend
        do jc=1,n2m
          do ic=1,n1m
           vz_mea(ic,jc,kc)=q3(ic,jc,kc)
           vy_mea(ic,jc,kc)=q2(ic,jc,kc)
           vx_mea(ic,jc,kc)=q1(ic,jc,kc)
           temp_mea(ic,jc,kc)=dsal(ic,jc,kc)
          end do
        end do
      end do

      !CS   Begin h5 and xmf routine

      !RO   Sort out MPI definitions and file names

      tprfi = 1/tframe
      itime=nint(time*tprfi)
      write(ipfi,'(i5.5)')itime

      filnam1='movie/field'//ipfi//'.h5'
      filnamxdm = 'movie/field'//ipfi//'.xmf' 

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

      !RO   Set offsets and element counts

      ndims=3

      dims(1)=n1m
      dims(2)=n2m
      dims(3)=n3m

      data_count(1) = n1m
      data_count(2) = n2m
      data_count(3) = kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = kstart-1 


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

      !RO   Open first continua file for dsal
      call h5fcreate_f(filnam1, H5F_ACC_TRUNC_F, file_temp, hdf_error, access_prp=file_plist)

      !RO   Create dataset on file

      call h5dcreate_f(file_temp, 'T_me', H5T_NATIVE_DOUBLE,filespace, dset_temp, hdf_error)
      call h5dcreate_f(file_temp, 'VX_me', H5T_NATIVE_DOUBLE,filespace, dset_VX, hdf_error)
      call h5dcreate_f(file_temp, 'VY_me', H5T_NATIVE_DOUBLE,filespace, dset_VY, hdf_error)
      call h5dcreate_f(file_temp, 'VZ_me', H5T_NATIVE_DOUBLE,filespace, dset_VZ, hdf_error)


      !RO   Set hyperslab----temp
      call h5dget_space_f(dset_temp, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)

      call h5dwrite_f(dset_temp, H5T_NATIVE_DOUBLE,&
     &   temp_mea(1:n1m,1:n2m,kstart:kend), dims, &
     &   hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
     &   xfer_prp = slab_plist)

      !RO   Set hyperslab----vz
      call h5dget_space_f(dset_VZ, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)

      call h5dwrite_f(dset_VZ, H5T_NATIVE_DOUBLE,&
     &   vz_mea(1:n1m,1:n2m,kstart:kend), dims, &
     &   hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
     &   xfer_prp = slab_plist)


      !RO   Set hyperslab----vy
      call h5dget_space_f(dset_VY, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)

      call h5dwrite_f(dset_VY, H5T_NATIVE_DOUBLE,&
     &   vy_mea(1:n1m,1:n2m,kstart:kend), dims, &
     &   hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
     &   xfer_prp = slab_plist)


      !RO   Set hyperslab----vx
      call h5dget_space_f(dset_VX, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)

      call h5dwrite_f(dset_VX, H5T_NATIVE_DOUBLE,&
     &   vx_mea(1:n1m,1:n2m,kstart:kend), dims, &
     &   hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
     &   xfer_prp = slab_plist)


      !RO   Close dataset and file for dsal

      call h5dclose_f(dset_temp, hdf_error)
      call h5dclose_f(dset_VZ, hdf_error)
      call h5dclose_f(dset_VY, hdf_error)
      call h5dclose_f(dset_VX, hdf_error)

      call h5fclose_f(file_temp, hdf_error)

      !RO   Close all other stuff

      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5pclose_f(file_plist, hdf_error)
      call h5pclose_f(slab_plist, hdf_error)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      !EP   Write the xdm

      if (myid.eq.0) then

      open(45,file=filnamxdm,status='unknown')
      rewind(45)
      write(45,'("<?xml version=""1.0"" ?>")')
      write(45,'("<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>")')
      write(45,'("<Xdmf Version=""2.0"">")')
      write(45,'("<Domain>")')
      write(45,'("<Grid Name=""RB Cartesian"" GridType=""Uniform"">")')
      write(45,'("<Topology TopologyType=""3DRectMesh"" &
     &NumberOfElements=""",i4," ",i4," ",i4,"""/>")') n3m,n2m,n1m
      write(45,'("<Geometry GeometryType=""VXVYVZ"">")')
      write(45,'("<DataItem Dimensions=""",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1m
      write(45,'("cordin_info.h5:/x")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2m
      write(45,'("cordin_info.h5:/y")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m
      write(45,'("cordin_info.h5:/z")')
      write(45,'("</DataItem>")')
      write(45,'("</Geometry>")')

      write(45,'("<Attribute Name=""T_me""&
     & AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')&
     & n3m,n2m,n1m
      write(45,'("field",i5.5,".h5:/T_me")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')

      write(45,'("<Attribute Name=""VZ_me""&
     & AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')&
     & n3m,n2m,n1m
      write(45,'("field",i5.5,".h5:/VZ_me")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')

      write(45,'("<Attribute Name=""VY_me""&
     & AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')&
     & n3m,n2m,n1m
      write(45,'("field",i5.5,".h5:/VY_me")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')

      write(45,'("<Attribute Name=""VX_me""&
     & AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')&
     & n3m,n2m,n1m
      write(45,'("field",i5.5,".h5:/VX_me")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')


      write(45,'("<Time Value=""",e12.5,""" />")') time
      write(45,'("</Grid>")')
      write(45,'("</Domain>")')
      write(45,'("</Xdmf>")')
      close(45)

      endif

      call h5close_f(hdf_error)

      if(allocated(temp_mea)) deallocate(temp_mea)
      if(allocated(vz_mea)) deallocate(vz_mea)
      if(allocated(vy_mea)) deallocate(vy_mea)
      if(allocated(vx_mea)) deallocate(vx_mea)

      end subroutine mkmov_field

