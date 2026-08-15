      subroutine mkmov_hdf_zcut

      use param
      use hdf5
      use mpih
      use mpi_param, only: kstart,kend
      use local_arrays, only: q1,q2,q3,pr,dens
      use mgrd_arrays, only: dsalc

      IMPLICIT NONE

      integer im,jm,km
      integer ic,jc,kc
      integer ip,jp,kp
      
      character(70) :: namfile,xdmnam
      character(5) :: ipfi
      integer hdf_error, itime
      integer zk

      integer(HID_T) :: file_id
      integer(HID_T) :: dspace
      integer(HSIZE_T) :: dims(2)

      integer(HID_T) :: dset_q1v
      integer(HID_T) :: dset_q2v
      integer(HID_T) :: dset_q3v
      integer(HID_T) :: dset_prv
      integer(HID_T) :: dset_sal
      integer(HID_T) :: dset_den
      
      real tprfi
      real, dimension(n1m,n2m) :: prx,v1,v2,v3,dsalx,densx

      zk = mov_zcut_k

      tprfi=1.d0/tframe
      itime=int(time/tframe)
      write(ipfi,'(i5.5)')itime

      namfile='flowmov/frame_'//ipfi//'_zcut.h5'
      xdmnam='flowmov/frame_'//ipfi//'_zcut.xmf'

      !-- find the proc owner 
      if(zk.ge.kstart.and.zk.le.kend) then
         kc=zk
         kp=kc+1
         do jc=1,n2m
            jp=jc+1
            do ic=1,n1m
               ip=ic+1
               v1(ic,jc) = (q1(ic,jc,kc)+q1(ip,jc,kc))*0.5d0
               v2(ic,jc) = (q2(ic,jc,kc)+q2(ic,jp,kc))*0.5d0
               v3(ic,jc) = (q3(ic,jc,kc)+q3(ic,jc,kp))*0.5d0
               prx(ic,jc) = pr(ic,jc,kc)
               dsalx(ic,jc) = dsalc(ic,jc,kc)
            end do
         end do
      
         call h5open_f(hdf_error)

         !-- set up mpi file properties
         call h5fcreate_f(namfile,H5F_ACC_TRUNC_F,file_id,hdf_error)

         dims(1)=n1m
         dims(2)=n2m
     
         call h5screate_simple_f(2,dims,dspace,hdf_error)
         call h5dcreate_f(file_id,'Vx',H5T_NATIVE_DOUBLE,dspace,dset_q1v,hdf_error)
         call h5dwrite_f(dset_q1v,H5T_NATIVE_DOUBLE,v1(1:n1m,1:n2m),dims,hdf_error)
         call h5dclose_f(dset_q1v,hdf_error)
         call h5sclose_f(dspace,hdf_error)

         call h5screate_simple_f(2,dims,dspace,hdf_error)
         call h5dcreate_f(file_id,'Vy',H5T_NATIVE_DOUBLE,dspace,dset_q2v,hdf_error)
         call h5dwrite_f(dset_q2v,H5T_NATIVE_DOUBLE,v2(1:n1m,1:n2m),dims,hdf_error)
         call h5dclose_f(dset_q2v,hdf_error)
         call h5sclose_f(dspace,hdf_error)

         call h5screate_simple_f(2,dims,dspace,hdf_error)
         call h5dcreate_f(file_id,'Vz',H5T_NATIVE_DOUBLE,dspace,dset_q3v,hdf_error)
         call h5dwrite_f(dset_q3v,H5T_NATIVE_DOUBLE,v3(1:n1m,1:n2m),dims,hdf_error)
         call h5dclose_f(dset_q3v,hdf_error)
         call h5sclose_f(dspace,hdf_error)

         call h5screate_simple_f(2,dims,dspace,hdf_error)
         call h5dcreate_f(file_id,'Pr',H5T_NATIVE_DOUBLE,dspace,dset_prv,hdf_error)
         call h5dwrite_f(dset_prv,H5T_NATIVE_DOUBLE,prx(1:n1m,1:n2m),dims,hdf_error)
         call h5dclose_f(dset_prv,hdf_error)
         call h5sclose_f(dspace,hdf_error)

         call h5screate_simple_f(2,dims,dspace,hdf_error)
         call h5dcreate_f(file_id,'Sal',H5T_NATIVE_DOUBLE,dspace,dset_sal,hdf_error)
         call h5dwrite_f(dset_sal,H5T_NATIVE_DOUBLE,dsalx(1:n1m,1:n2m),dims,hdf_error)
         call h5dclose_f(dset_sal,hdf_error)
         call h5sclose_f(dspace,hdf_error)

         !CS Sufficient for current proc to close file_id, otherwise error thrown
         call h5fclose_f(file_id,hdf_error)
         call h5close_f(hdf_error)

      end if

      if (myid.eq.0) then

         open(45,file=xdmnam,status='unknown')
         rewind(45)
         write(45,'("<?xml version=""1.0"" ?>")')
         write(45,'("<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>")')
         write(45,'("<Xdmf Version=""2.0"">")')
         write(45,'("<Domain>")')
         write(45,'("<Grid Name=""zcut"" GridType=""Uniform"">")')
         write(45,'("<Topology TopologyType=""2DSMesh"" NumberOfElements=""",i4," ",i4,"""/>")') n2m,n1m
         write(45,'("<Geometry GeometryType=""X_Y"">")')
         write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2m,n1m
         write(45,'("cordin_zcut_info.h5:/x")')
         write(45,'("</DataItem>")')
         write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2m,n1m
         write(45,'("cordin_zcut_info.h5:/y")')
         write(45,'("</DataItem>")')
         write(45,'("</Geometry>")')
         write(45,'("<Attribute Name=""X-Velocity"" AttributeType=""Scalar"" Center=""Node"">")')
         write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2m,n1m
         write(45,'("frame_",i5.5,"_zcut.h5:/Vx")') itime
         write(45,'("</DataItem>")')
         write(45,'("</Attribute>")')
         write(45,'("<Attribute Name=""Y-Velocity"" AttributeType=""Scalar"" Center=""Node"">")')
         write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2m,n1m
         write(45,'("frame_",i5.5,"_zcut.h5:/Vy")') itime
         write(45,'("</DataItem>")')
         write(45,'("</Attribute>")')
         write(45,'("<Attribute Name=""Z-Velocity"" AttributeType=""Scalar"" Center=""Node"">")')
         write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2m,n1m
         write(45,'("frame_",i5.5,"_zcut.h5:/Vz")') itime
         write(45,'("</DataItem>")')
         write(45,'("</Attribute>")')
         write(45,'("<Attribute Name=""Pressure"" AttributeType=""Scalar"" Center=""Node"">")')
         write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2m,n1m
         write(45,'("frame_",i5.5,"_zcut.h5:/Pr")') itime
         write(45,'("</DataItem>")')
         write(45,'("</Attribute>")')
         write(45,'("<Attribute Name=""Salinity"" AttributeType=""Scalar"" Center=""Node"">")')
         write(45,'("<DataItem Dimensions=""",i4," ",i4,""" NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2m,n1m
         write(45,'("frame_",i5.5,"_zcut.h5:/Sal")') itime
         write(45,'("</DataItem>")')
         write(45,'("</Attribute>")')
         write(45,'("<Time Value=""",e12.5""" />")')time
         write(45,'("</Grid>")')
         write(45,'("</Domain>")')
         write(45,'("</Xdmf>")')
         close(45)

      end if

      return
      end
