!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcWriteQ.F90                                 !
!    CONTAINS: subroutine CalcWriteQ                      !
!                                                         ! 
!    PURPOSE: Compute and write the 3D q-criteria field   !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
      subroutine mkmov_qtens

      use local_arrays, only: q3,q2,q1
      use mpi_param
      use mpih
      use hdf5
      use param

      implicit none

      real, allocatable, dimension(:,:,:) :: qtens
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

      integer(HID_T) :: file_qtens

      integer(HID_T) :: dset_qtens

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

      allocate(qtens(1:n1m,1:n2m,kstart:kend))

      do kc=kstart,kend
        km=kc-1
        kp=kc+1
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
          do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)
      !do ic=xstart(3),xend(3)
      !im=ic-1
      !ip=ic+1
      !do jc=xstart(2),xend(2)
      !jm=jc-1
      !jp=jc+1
      !do kc=2,nxm
      !km=kc-1
      !kp=kc+1

      dvxx1=(q1(ip,jc,kc)-q1(ic,jc,kc))*dx1

      dvxx2=((q1(ic,jp,kc)+q1(ip,jp,kc))- &
            (q1(ic,jm,kc)+q1(ip,jm,kc)))*0.25d0*dx2

      dvxx3=((q1(ic,jc,kp)+q1(ip,jc,kp))- &
            (q1(ic,jc,km)+q1(ip,jc,km)))*0.25d0*udx3m(kc)


      dvyx1=((q2(ip,jc,kc)+q2(ip,jp,kc))- &
            (q2(im,jc,kc)+q2(im,jp,kc)))*0.25d0*dx1

      dvyx2=(q2(ic,jp,kc)-q2(ic,jc,kc))*dx2

      dvyx3=((q2(ic,jc,kp)+q2(ic,jp,kp))- &
            (q2(ic,jc,km)+q2(ic,jp,km)))*0.25d0*udx3m(kc)


      dvzx1=((q3(ip,jc,kc)+q3(ip,jc,kp))- &
            (q3(im,jc,kc)+q3(im,jc,kp)))*0.25d0*dx1

      dvzx2=((q3(ic,jp,kc)+q3(ic,jp,kp))- &
            (q3(ic,jm,kc)+q3(ic,jm,kp)))*0.25d0*dx2

      dvzx3=(q3(ic,jc,kp)-q3(ic,jc,kc))*udx3m(kc)

      strn=dvxx1**2.0+dvyx2**2.0+dvzx3**2.0+ &
           0.5*((dvyx1+dvxx2)**2.0+ &
                (dvzx1+dvxx3)**2.0+ &
                (dvyx3+dvzx2)**2.0)

      omeg=0.5*((dvyx1-dvxx2)**2.0+ &
                (dvzx1-dvxx3)**2.0+ &
                (dvzx2-dvyx3)**2.0)

      qtens(ic,jc,kc)=0.5*(omeg-strn)

      end do
      end do
      end do

      !CS   Begin h5 and xmf routine

      !RO   Sort out MPI definitions and file names

      tprfi = 1/tframe
      itime=nint(time*tprfi)
      write(ipfi,'(i5.5)')itime

      filnam1='movie/qcrit'//ipfi//'.h5'
      filnamxdm = 'movie/qcrit'//ipfi//'.xmf' 

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
      call h5fcreate_f(filnam1, H5F_ACC_TRUNC_F, file_qtens, hdf_error, access_prp=file_plist)

      !RO   Create dataset on file

      call h5dcreate_f(file_qtens, 'Q', H5T_NATIVE_DOUBLE,filespace, dset_qtens, hdf_error)

      !RO   Set hyperslab
      call h5dget_space_f(dset_qtens, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)

      call h5dwrite_f(dset_qtens, H5T_NATIVE_DOUBLE,&
     &   qtens(1:n1m,1:n2m,kstart:kend), dims, &
     &   hdf_error, file_space_id = slabspace, mem_space_id = memspace, &
     &   xfer_prp = slab_plist)

      !RO   Close dataset and file for dsal

      call h5dclose_f(dset_qtens, hdf_error)
      call h5fclose_f(file_qtens, hdf_error)

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
      write(45,'("<Attribute Name=""Q""&
     & AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,"""&
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')&
     & n3m,n2m,n1m
      write(45,'("qcrit",i5.5,".h5:/Q")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Time Value=""",e12.5,""" />")') time
      write(45,'("</Grid>")')
      write(45,'("</Domain>")')
      write(45,'("</Xdmf>")')
      close(45)

      endif

      call h5close_f(hdf_error)

      if(allocated(qtens)) deallocate(qtens)

      end subroutine mkmov_qtens

