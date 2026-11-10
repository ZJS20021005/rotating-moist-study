      subroutine inimov
      use mpih
      use param
      use hdf5

      IMPLICIT none

      integer hdf_error

      integer(HID_T) :: file_grid
      integer(HID_T) :: dset_grid
      integer(HID_T) :: dspace_grid

      integer(HSIZE_T) :: dims_grid(1)

      character(70) namfile

      call h5open_f(hdf_error)

      if (myid.eq.0) then
      namfile='movie/cordin_info.h5'

      !RO   Write the grid information normally

      call h5fcreate_f(namfile,H5F_ACC_TRUNC_F,file_grid,hdf_error)

      dims_grid(1)=n1mr
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'x', H5T_NATIVE_DOUBLE,dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, xmr(1:n1mr),dims_grid, hdf_error)


      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=n2mr
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)

      call h5dcreate_f(file_grid, 'y', H5T_NATIVE_DOUBLE, dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, ymr(1:n2mr),dims_grid,hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=n3mr
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'z', H5T_NATIVE_DOUBLE,dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zmr(1:n3mr),dims_grid, hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5fclose_f(file_grid, hdf_error)
      endif

      call h5close_f(hdf_error)

      return
      end


