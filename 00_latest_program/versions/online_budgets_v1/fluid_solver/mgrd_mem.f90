      subroutine mgrd_mem_alloc
      use mpih
      use param
      use mpi_param
      use mgrd_arrays

      implicit none
      integer :: merr, merr_all

      integer :: i,j,k

      merr_all = 0

      !---------------------------------------------------
      ! Auxiliary arrays for inform trasfer between meshes
      !---------------------------------------------------

      allocate(q1lr(1:n1r,1:n2r,kstartr-lvlhalo:kendr+lvlhalo),stat=merr)
      merr_all = merr_all + merr
      allocate(q2lr(1:n1r,1:n2r,kstartr-lvlhalo:kendr+lvlhalo),stat=merr)
      merr_all = merr_all + merr
      allocate(q3lr(1:n1r,1:n2r,kstartr-lvlhalo:kendr+lvlhalo),stat=merr)
      merr_all = merr_all + merr

      allocate(dsalc(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo),stat=merr)
      merr_all = merr_all + merr
      allocate(qvapc(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo),stat=merr)
      merr_all = merr_all + merr

      if(merr_all.ne.0)then
        write(*,*)myid, 'multi grid memory alloc error'
        write(*,*)merr_all
      endif

!m===============================================
      do k=kstartr-1,kendr+1
        do j=1,n2r
          do i=1,n1r
            q1lr(i,j,k)=0.d0
            q2lr(i,j,k)=0.d0
            q3lr(i,j,k)=0.d0
          enddo
        enddo
      enddo

      do k=kstart-1,kend+1
        do j=1,n2
          do i=1,n1
            dsalc(i,j,k)=0.d0
            qvapc(i,j,k)=0.d0
          enddo
        enddo
      enddo

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end subroutine mgrd_mem_alloc

!==================================================      
      subroutine mgrd_mem_dealloc
      use mgrd_arrays

      implicit none

      if(allocated(q1lr)) deallocate(q1lr)
      if(allocated(q2lr)) deallocate(q2lr)
      if(allocated(q3lr)) deallocate(q3lr)

      if(allocated(dsalc)) deallocate(dsalc)
      if(allocated(qvapc)) deallocate(qvapc)

      return
      end subroutine mgrd_mem_dealloc
