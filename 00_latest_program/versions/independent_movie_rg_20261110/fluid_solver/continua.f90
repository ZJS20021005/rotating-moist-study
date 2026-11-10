!================================================
      subroutine mpi_write_continua
      use param
      use mpih
      use mpi_param, only: kstart,kend,kstartr,kendr
      use local_arrays, only: q2,q3,q1,dsal,qvap
      use hdf5
      implicit none

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_q1
      integer(HID_T) :: dset_q2
      integer(HID_T) :: dset_q3
      integer(HID_T) :: dset_dsal
      integer(HID_T) :: dset_qvap

      integer(HID_T) :: plist_id

      integer(HSIZE_T), dimension(3) :: dims
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info
      integer :: ndims

      character filnam1*30
      character filnam2*30
      character filnam3*30
      character filnam4*30
      character filnam5*30
      character filnam6*30
      character filnam9*30
      character filnam7*30
      character filnam8*30

      call h5open_f(hdf_error)

      ! Sort out MPI definitions
      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

      ! Form the name of the file

      filnam2 = 'continua_q1.h5'
      filnam3 = 'continua_q2.h5'
      filnam4 = 'continua_q3.h5'
      filnam5 = 'continua_dsal.h5'
      filnam6 = 'continua_pr.h5'
      filnam9 = 'continua_qvap.h5'

      ! Set offsets and element counts
   
      ndims = 3

      dims(1)=n1
      dims(2)=n2
      dims(3)=n3m

      data_count(1) = n1
      data_count(2) = n2
      data_count(3) = kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = kstart-1

      ! q1
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam2, H5F_ACC_TRUNC_F, file_id,hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'Vth', H5T_NATIVE_DOUBLE,filespace, dset_q1, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 
      call h5dget_space_f(dset_q1, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_q1, H5T_NATIVE_DOUBLE,q1(1:n1,1:n2,kstart:kend), dims, hdf_error, file_space_id = slabspace,mem_space_id = memspace, xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_q1, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      ! q2
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam3, H5F_ACC_TRUNC_F, file_id,hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'Vr', H5T_NATIVE_DOUBLE,filespace, dset_q2, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 
      call h5dget_space_f(dset_q2, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_q2, H5T_NATIVE_DOUBLE,q2(1:n1,1:n2,kstart:kend), dims,hdf_error, file_space_id = slabspace,mem_space_id = memspace, xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_q2, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      ! qvap (mirror of dsal write)
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam9, H5F_ACC_TRUNC_F, file_id, hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'qvap', H5T_NATIVE_DOUBLE,filespace, dset_qvap, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_qvap, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_qvap, H5T_NATIVE_DOUBLE,qvap(1:n1,1:n2,kstart:kend), dims,hdf_error, file_space_id = slabspace,mem_space_id = memspace,xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_qvap, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      ! q2
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam4, H5F_ACC_TRUNC_F, file_id,hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'Vz', H5T_NATIVE_DOUBLE,filespace, dset_q3, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 
      call h5dget_space_f(dset_q3, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_q3, H5T_NATIVE_DOUBLE,q3(1:n1,1:n2,kstart:kend), dims, hdf_error, file_space_id = slabspace,mem_space_id = memspace, xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_q3, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)


      ! Set offsets and element counts

      ndims = 3
      dims(1)=n1r
      dims(2)=n2r
      dims(3)=n3mr

      data_count(1) = n1r
      data_count(2) = n2r
      data_count(3) = kendr-kstartr+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = kstartr-1

      ! dsal

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fcreate_f(filnam5, H5F_ACC_TRUNC_F, file_id, hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)
      call h5dcreate_f(file_id, 'dsal', H5T_NATIVE_DOUBLE,filespace, dset_dsal, hdf_error)
      call h5sclose_f(filespace, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error)
      call h5dget_space_f(dset_dsal, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,hdf_error)
      call h5dwrite_f(dset_dsal, H5T_NATIVE_DOUBLE,dsal(1:n1r,1:n2r,kstartr:kendr), dims,hdf_error, file_space_id = slabspace,mem_space_id = memspace,xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)
      call h5dclose_f(dset_dsal, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)


      if (myid .eq. 0) then
        open(13,file='continua_grid.dat',status='unknown')
        rewind(13)                                                      
        write(13,'(3i8)') n1,n2,n3
        write(13,'(3f25.15)') rext1,rext2,time
        write(13,'(i8,f25.15)') istr3,str3
        write(13,'(3i8)') mref1, mref2, mref3
        close(13)
      endif

      return
      end subroutine mpi_write_continua

!================================================      
      subroutine mpi_read_continua(n1o,n2o,n3o,ks,ke,intvar,qua)
      use mpih
      use param
      use hdf5
      implicit none
      integer, intent(in) :: ks,ke,n2o,n1o,n3o
      real, dimension(1:n1o,1:n2o,ks-lvlhalo:ke+lvlhalo),intent(out)::qua
      integer k,j,i

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_qua

      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info
      integer :: ndims

      integer, intent(in) :: intvar
      character(70) :: filnam1
      character(10) :: dsetname

      call h5open_f(hdf_error)

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

      ! Select file and dataset based on intvar

      select case (intvar)
        case (1)
          dsetname = trim('Vth')
          filnam1 = trim('continua_q1.h5')
        case (2)
          dsetname = trim('Vr')
          filnam1 = trim('continua_q2.h5')
        case (3)
          dsetname = trim('Vz')
          filnam1 = trim('continua_q3.h5')
        case (5)
          dsetname = trim('dsal')
          filnam1 = trim('continua_dsal.h5')
        case (8)
          dsetname = trim('qvap')
          filnam1 = trim('continua_qvap.h5')
        case (6)
          dsetname = trim('dens_in')
          filnam1 = trim('continua_dens_in.h5')
        case (7)
          dsetname = trim('dsal_in')
          filnam1 = trim('continua_dsal_in.h5')
      end select

      do k=ks,ke
        do j=1,n2o
          do i=1,n1o
            qua(i,j,k)=0.d0
          enddo
        enddo
      enddo

      ! Set offsets and element counts
   
      ndims = 3

      dims(1)=n1o
      dims(2)=n2o
      dims(3)=n3o-1


      data_count(1) = n1o
      data_count(2) = n2o
      data_count(3) = ke-ks+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = ks-1

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf_error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, hdf_error)
      call h5fopen_f(filnam1, H5F_ACC_RDONLY_F, file_id, hdf_error, access_prp=plist_id)
      call h5pclose_f(plist_id,hdf_error)

      call h5dopen_f(file_id, dsetname, dset_qua, hdf_error)
      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 
      call h5dget_space_f(dset_qua, slabspace, hdf_error)
      call h5sselect_hyperslab_f(slabspace, H5S_SELECT_SET_F,data_offset, data_count, hdf_error)

      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf_error)

      call h5dread_f(dset_qua, H5T_NATIVE_DOUBLE,qua(1:n1o,1:n2o,ks:ke), dims, hdf_error, file_space_id = slabspace,mem_space_id = memspace, xfer_prp = plist_id)

      call h5dclose_f(dset_qua, hdf_error)
      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5pclose_f(plist_id, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      call h5close_f(hdf_error)

      if(myid.eq.0)write(*,'(5x,a)')'reading complete: '//filnam1

      return
      end subroutine mpi_read_continua
