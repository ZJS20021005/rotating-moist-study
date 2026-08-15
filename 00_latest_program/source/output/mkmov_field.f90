!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Refined three-dimensional movie output.
!
! Files are intentionally separated:
!   fieldNNNNN.h5               : DSAL_me, VZ_me, RH_me
!   horizontal_velocityNNNNN.h5 : VX_me, VY_me
!   mprimeNNNNN.h5              : MPRIME_me
!   condensationNNNNN.h5        : COND_me
!
! MPRIME_me is computed independently on every horizontal plane:
!   m' = b + gamma*qvap - <b + gamma*qvap>_xy(z,t).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mkmov_field

      use local_arrays, only: dsal,qvap,qs,cond_term
      use mgrd_arrays, only: q1lr,q2lr,q3lr
      use mpi_param
      use mpih
      use hdf5
      use param
      implicit none

      real, allocatable :: dsal_mea(:,:,:)
      real, allocatable :: vx_mea(:,:,:)
      real, allocatable :: vy_mea(:,:,:)
      real, allocatable :: vz_mea(:,:,:)
      real, allocatable :: rh_mea(:,:,:)
      real, allocatable :: mprime_mea(:,:,:)
      real, allocatable :: cond_mea(:,:,:)

      integer :: ic,jc,kc,ip,jp,kp
      integer :: hdf_error,comm,info,ndims,itime,ios
      integer(HID_T) :: filespace,slabspace,memspace
      integer(HID_T) :: file_id,file_plist,slab_plist
      integer(HID_T) :: dset_dsal,dset_vx,dset_vy,dset_vz,dset_rh
      integer(HID_T) :: dset_mprime,dset_cond
      integer(HSIZE_T) :: dims(3)
      integer(HSIZE_T) :: data_count(3)
      integer(HSSIZE_T) :: data_offset(3)
      real :: tprfi,qsval,mbar,mcell
      character(70) :: field_h5,hvel_h5,mprime_h5,cond_h5
      character(70) :: field_xmf,hvel_xmf,mprime_xmf,cond_xmf
      character(5) :: ipfi
      character(len=256) :: iomsg

      call update_both_ghosts(n1r,n2r,dsal,kstartr,kendr)
      call update_both_ghosts(n1r,n2r,qvap,kstartr,kendr)
      call update_both_ghosts(n1r,n2r,qs,kstartr,kendr)
      call update_both_ghosts(n1r,n2r,cond_term,kstartr,kendr)
      call update_both_ghosts(n1r,n2r,q1lr,kstartr,kendr)
      call update_both_ghosts(n1r,n2r,q2lr,kstartr,kendr)
      call update_both_ghosts(n1r,n2r,q3lr,kstartr,kendr)

      allocate(dsal_mea(1:n1mr,1:n2mr,kstartr:kendr))
      allocate(vx_mea(1:n1mr,1:n2mr,kstartr:kendr))
      allocate(vy_mea(1:n1mr,1:n2mr,kstartr:kendr))
      allocate(vz_mea(1:n1mr,1:n2mr,kstartr:kendr))
      allocate(rh_mea(1:n1mr,1:n2mr,kstartr:kendr))
      allocate(mprime_mea(1:n1mr,1:n2mr,kstartr:kendr))
      allocate(cond_mea(1:n1mr,1:n2mr,kstartr:kendr))

      do kc=kstartr,kendr
        mbar = 0.0d0
        do jc=1,n2mr
          do ic=1,n1mr
            mbar = mbar+dsal(ic,jc,kc)+gamma*qvap(ic,jc,kc)
          enddo
        enddo
        mbar = mbar/(dble(n1mr)*dble(n2mr))

        kp = kc+1
        do jc=1,n2mr
          jp = jpvr(jc)
          do ic=1,n1mr
            ip = ipvr(ic)
            dsal_mea(ic,jc,kc) = dsal(ic,jc,kc)
            vx_mea(ic,jc,kc) = 0.5d0*(q1lr(ic,jc,kc)+ &
     &                                q1lr(ip,jc,kc))
            vy_mea(ic,jc,kc) = 0.5d0*(q2lr(ic,jc,kc)+ &
     &                                q2lr(ic,jp,kc))
            vz_mea(ic,jc,kc) = 0.5d0*(q3lr(ic,jc,kc)+q3lr(ic,jc,kp))
            cond_mea(ic,jc,kc) = cond_term(ic,jc,kc)
            mcell = dsal(ic,jc,kc)+gamma*qvap(ic,jc,kc)
            mprime_mea(ic,jc,kc) = mcell-mbar

            qsval = qs(ic,jc,kc)
            if(qsval.gt.0.0d0 .and. qsval.eq.qsval)then
              rh_mea(ic,jc,kc) = qvap(ic,jc,kc)/qsval
            else
              rh_mea(ic,jc,kc) = 0.0d0
            endif
          enddo
        enddo
      enddo

      tprfi = 1.0d0/tframe
      itime = nint(time*tprfi)
      write(ipfi,'(i5.5)') itime

      field_h5 = 'movie/field'//ipfi//'.h5'
      field_xmf = 'movie/field'//ipfi//'.xmf'
      hvel_h5 = 'movie/horizontal_velocity'//ipfi//'.h5'
      hvel_xmf = 'movie/horizontal_velocity'//ipfi//'.xmf'
      mprime_h5 = 'movie/mprime'//ipfi//'.h5'
      mprime_xmf = 'movie/mprime'//ipfi//'.xmf'
      cond_h5 = 'movie/condensation'//ipfi//'.h5'
      cond_xmf = 'movie/condensation'//ipfi//'.xmf'

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL
      ndims = 3
      dims = (/n1mr,n2mr,n3mr/)
      data_count = (/n1mr,n2mr,kendr-kstartr+1/)
      data_offset = (/0,0,kstartr-1/)

      call h5open_f(hdf_error)
      call h5pcreate_f(H5P_FILE_ACCESS_F,file_plist,hdf_error)
      call h5pset_fapl_mpio_f(file_plist,comm,info,hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F,slab_plist,hdf_error)
      call h5pset_dxpl_mpio_f(slab_plist,H5FD_MPIO_COLLECTIVE_F, &
     &                        hdf_error)
      call h5screate_simple_f(ndims,dims,filespace,hdf_error)
      call h5screate_simple_f(ndims,data_count,memspace,hdf_error)

!     Main requested field file: b, vertical velocity, and RH.
      call h5fcreate_f(field_h5,H5F_ACC_TRUNC_F,file_id,hdf_error, &
     &                 access_prp=file_plist)
      call h5dcreate_f(file_id,'DSAL_me',H5T_NATIVE_DOUBLE,filespace, &
     &                 dset_dsal,hdf_error)
      call h5dcreate_f(file_id,'VZ_me',H5T_NATIVE_DOUBLE,filespace, &
     &                 dset_vz,hdf_error)
      call h5dcreate_f(file_id,'RH_me',H5T_NATIVE_DOUBLE,filespace, &
     &                 dset_rh,hdf_error)

      call h5dget_space_f(dset_dsal,slabspace,hdf_error)
      call h5sselect_hyperslab_f(slabspace,H5S_SELECT_SET_F, &
     &                           data_offset,data_count,hdf_error)
      call h5dwrite_f(dset_dsal,H5T_NATIVE_DOUBLE,dsal_mea,dims, &
     & hdf_error,file_space_id=slabspace,mem_space_id=memspace, &
     & xfer_prp=slab_plist)

      call h5dget_space_f(dset_vz,slabspace,hdf_error)
      call h5sselect_hyperslab_f(slabspace,H5S_SELECT_SET_F, &
     &                           data_offset,data_count,hdf_error)
      call h5dwrite_f(dset_vz,H5T_NATIVE_DOUBLE,vz_mea,dims, &
     & hdf_error,file_space_id=slabspace,mem_space_id=memspace, &
     & xfer_prp=slab_plist)

      call h5dget_space_f(dset_rh,slabspace,hdf_error)
      call h5sselect_hyperslab_f(slabspace,H5S_SELECT_SET_F, &
     &                           data_offset,data_count,hdf_error)
      call h5dwrite_f(dset_rh,H5T_NATIVE_DOUBLE,rh_mea,dims, &
     & hdf_error,file_space_id=slabspace,mem_space_id=memspace, &
     & xfer_prp=slab_plist)

      call h5dclose_f(dset_dsal,hdf_error)
      call h5dclose_f(dset_vz,hdf_error)
      call h5dclose_f(dset_rh,hdf_error)
      call h5fclose_f(file_id,hdf_error)

!     Horizontal velocity components in their own file.
      call h5fcreate_f(hvel_h5,H5F_ACC_TRUNC_F,file_id,hdf_error, &
     &                 access_prp=file_plist)
      call h5dcreate_f(file_id,'VX_me',H5T_NATIVE_DOUBLE,filespace, &
     &                 dset_vx,hdf_error)
      call h5dcreate_f(file_id,'VY_me',H5T_NATIVE_DOUBLE,filespace, &
     &                 dset_vy,hdf_error)
      call h5dget_space_f(dset_vx,slabspace,hdf_error)
      call h5sselect_hyperslab_f(slabspace,H5S_SELECT_SET_F, &
     &                           data_offset,data_count,hdf_error)
      call h5dwrite_f(dset_vx,H5T_NATIVE_DOUBLE,vx_mea,dims, &
     & hdf_error,file_space_id=slabspace,mem_space_id=memspace, &
     & xfer_prp=slab_plist)
      call h5dget_space_f(dset_vy,slabspace,hdf_error)
      call h5sselect_hyperslab_f(slabspace,H5S_SELECT_SET_F, &
     &                           data_offset,data_count,hdf_error)
      call h5dwrite_f(dset_vy,H5T_NATIVE_DOUBLE,vy_mea,dims, &
     & hdf_error,file_space_id=slabspace,mem_space_id=memspace, &
     & xfer_prp=slab_plist)
      call h5dclose_f(dset_vx,hdf_error)
      call h5dclose_f(dset_vy,hdf_error)
      call h5fclose_f(file_id,hdf_error)

!     Moist-static-energy anomaly in its own file.
      call h5fcreate_f(mprime_h5,H5F_ACC_TRUNC_F,file_id,hdf_error, &
     &                 access_prp=file_plist)
      call h5dcreate_f(file_id,'MPRIME_me',H5T_NATIVE_DOUBLE,filespace, &
     &                 dset_mprime,hdf_error)
      call h5dget_space_f(dset_mprime,slabspace,hdf_error)
      call h5sselect_hyperslab_f(slabspace,H5S_SELECT_SET_F, &
     &                           data_offset,data_count,hdf_error)
      call h5dwrite_f(dset_mprime,H5T_NATIVE_DOUBLE,mprime_mea,dims, &
     & hdf_error,file_space_id=slabspace,mem_space_id=memspace, &
     & xfer_prp=slab_plist)
      call h5dclose_f(dset_mprime,hdf_error)
      call h5fclose_f(file_id,hdf_error)

!     Condensation rate in its own file.
      call h5fcreate_f(cond_h5,H5F_ACC_TRUNC_F,file_id,hdf_error, &
     &                 access_prp=file_plist)
      call h5dcreate_f(file_id,'COND_me',H5T_NATIVE_DOUBLE,filespace, &
     &                 dset_cond,hdf_error)
      call h5dget_space_f(dset_cond,slabspace,hdf_error)
      call h5sselect_hyperslab_f(slabspace,H5S_SELECT_SET_F, &
     &                           data_offset,data_count,hdf_error)
      call h5dwrite_f(dset_cond,H5T_NATIVE_DOUBLE,cond_mea,dims, &
     & hdf_error,file_space_id=slabspace,mem_space_id=memspace, &
     & xfer_prp=slab_plist)
      call h5dclose_f(dset_cond,hdf_error)
      call h5fclose_f(file_id,hdf_error)

      call h5sclose_f(memspace,hdf_error)
      call h5sclose_f(slabspace,hdf_error)
      call h5sclose_f(filespace,hdf_error)
      call h5pclose_f(file_plist,hdf_error)
      call h5pclose_f(slab_plist,hdf_error)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if(myid.eq.0)then
        call write_main_field_xmf(field_xmf,itime,ios,iomsg)
        call write_horizontal_velocity_xmf(hvel_xmf,itime,ios,iomsg)
        call write_single_field_xmf(mprime_xmf,'mprime', &
     &       'MPRIME_me',itime,ios,iomsg)
        call write_single_field_xmf(cond_xmf,'condensation', &
     &       'COND_me',itime,ios,iomsg)
      endif

      call h5close_f(hdf_error)

      deallocate(dsal_mea,vx_mea,vy_mea,vz_mea,rh_mea)
      deallocate(mprime_mea,cond_mea)
      return
      end

!=============================================================================
      subroutine write_main_field_xmf(filename,itime,ios,iomsg)
      use param
      use mpih
      implicit none
      character(*), intent(in) :: filename
      integer, intent(in) :: itime
      integer, intent(out) :: ios
      character(*), intent(out) :: iomsg

      open(45,file=filename,status='replace',action='write', &
     &     iostat=ios,iomsg=iomsg)
      if(ios.ne.0)then
        write(*,'(a,1x,a)') 'Error opening XMF file:',trim(filename)
        call MPI_Abort(MPI_COMM_WORLD,1,ierr)
      endif

      call write_xmf_header(45)
      call write_xmf_attribute(45,'DSAL_me','field',itime)
      call write_xmf_attribute(45,'VZ_me','field',itime)
      call write_xmf_attribute(45,'RH_me','field',itime)
      call write_xmf_footer(45)
      close(45)
      return
      end

!=============================================================================
      subroutine write_horizontal_velocity_xmf(filename,itime,ios,iomsg)
      use param
      use mpih
      implicit none
      character(*), intent(in) :: filename
      integer, intent(in) :: itime
      integer, intent(out) :: ios
      character(*), intent(out) :: iomsg

      open(45,file=filename,status='replace',action='write', &
     &     iostat=ios,iomsg=iomsg)
      if(ios.ne.0)then
        write(*,'(a,1x,a)') 'Error opening XMF file:',trim(filename)
        call MPI_Abort(MPI_COMM_WORLD,1,ierr)
      endif

      call write_xmf_header(45)
      call write_xmf_attribute(45,'VX_me','horizontal_velocity',itime)
      call write_xmf_attribute(45,'VY_me','horizontal_velocity',itime)
      call write_xmf_footer(45)
      close(45)
      return
      end

!=============================================================================
      subroutine write_single_field_xmf(filename,h5prefix,attribute, &
     &                                  itime,ios,iomsg)
      use param
      use mpih
      implicit none
      character(*), intent(in) :: filename,h5prefix,attribute
      integer, intent(in) :: itime
      integer, intent(out) :: ios
      character(*), intent(out) :: iomsg

      open(45,file=filename,status='replace',action='write', &
     &     iostat=ios,iomsg=iomsg)
      if(ios.ne.0)then
        write(*,'(a,1x,a)') 'Error opening XMF file:',trim(filename)
        call MPI_Abort(MPI_COMM_WORLD,1,ierr)
      endif

      call write_xmf_header(45)
      call write_xmf_attribute(45,attribute,h5prefix,itime)
      call write_xmf_footer(45)
      close(45)
      return
      end

!=============================================================================
      subroutine write_xmf_header(iunit)
      use param
      implicit none
      integer, intent(in) :: iunit

      write(iunit,'("<?xml version=""1.0"" ?>")')
      write(iunit,'("<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>")')
      write(iunit,'("<Xdmf Version=""2.0""><Domain>")')
      write(iunit,'("<Grid Name=""RB Cartesian"" GridType=""Uniform"">")')
      write(iunit,'("<Topology TopologyType=""3DRectMesh"" ", &
     & "NumberOfElements=""",i0," ",i0," ",i0,"""/>")') &
     & n3mr,n2mr,n1mr
      write(iunit,'("<Geometry GeometryType=""VXVYVZ"">")')
      write(iunit,'("<DataItem Dimensions=""",i0, &
     & """ NumberType=""Float"" Precision=""8"" Format=""HDF"">")') n1mr
      write(iunit,'("cordin_info.h5:/x</DataItem>")')
      write(iunit,'("<DataItem Dimensions=""",i0, &
     & """ NumberType=""Float"" Precision=""8"" Format=""HDF"">")') n2mr
      write(iunit,'("cordin_info.h5:/y</DataItem>")')
      write(iunit,'("<DataItem Dimensions=""",i0, &
     & """ NumberType=""Float"" Precision=""8"" Format=""HDF"">")') n3mr
      write(iunit,'("cordin_info.h5:/z</DataItem>")')
      write(iunit,'("</Geometry>")')
      return
      end

!=============================================================================
      subroutine write_xmf_attribute(iunit,attribute,h5prefix,itime)
      use param
      implicit none
      integer, intent(in) :: iunit,itime
      character(*), intent(in) :: attribute,h5prefix

      write(iunit,'("<Attribute Name=""",a, &
     & """ AttributeType=""Scalar"" Center=""Node"">")') trim(attribute)
      write(iunit,'("<DataItem Dimensions=""",i0," ",i0," ",i0, &
     & """ NumberType=""Float"" Precision=""8"" Format=""HDF"">")') &
     & n3mr,n2mr,n1mr
      write(iunit,'(a,i5.5,a,a,a)') trim(h5prefix),itime, &
     & '.h5:/',trim(attribute),'</DataItem>'
      write(iunit,'("</Attribute>")')
      return
      end

!=============================================================================
      subroutine write_xmf_footer(iunit)
      use param
      implicit none
      integer, intent(in) :: iunit
      write(iunit,'("<Time Value=""",es20.10,""" />")') time
      write(iunit,'("</Grid></Domain></Xdmf>")')
      return
      end
