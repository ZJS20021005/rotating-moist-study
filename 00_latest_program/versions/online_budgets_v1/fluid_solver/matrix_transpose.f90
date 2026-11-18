      subroutine PackZ_UnpackR(aa,bb)
      use mpih
      use mpi_param
      use param, only: n1m,n2m,n3m
      implicit none
      real, intent(in) :: aa(1:n1m,1:n2m,kstart:kend)
      real, intent(out) :: bb(1:n3m,1:n1m,jstart:jend)
      real, allocatable :: sbuf(:),rbuf(:)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:)
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr

      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))
      do i=0,numtasks-1
        aaj(i)= dk* countj(i)*n1m
        aak(i)= dj* countk(i)*n1m
      end do
      
      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1)) 
      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do
     
      if(.not. allocated(sbuf)) allocate(sbuf(0:n1m*n2m*dk-1),stat=merr)
      if(merr .ne. 0) then
        write(*,*)"process  ",myid," failed to allocate memory for sbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      if(.not. allocated(rbuf)) allocate(rbuf(0:n1m*n3m*dj-1),stat=merr)
      
      if(merr .ne. 0) then
        write(*,*)"process  ",myid," failed to allocate memory for sbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif
      do kk = 0, numtasks-1
        nc = dispj(kk)
        dr = countj(kk)
        offsetr = offsetj(kk)
        do k = kstart,kend
          do j=1,dr
            do i=1,n1m
              sbuf(nc) = aa(i,j+offsetr,k)
              nc=nc+1
            enddo
          enddo
        enddo
      enddo
      call MPI_ALLTOALLV(sbuf, aaj,dispj, MDP,rbuf, aak,dispk, MDP,MPI_COMM_WORLD, ierr)
     
      do kk = 0, numtasks-1
        nc = dispk(kk)
        dz = countk(kk)
        offsetz = offsetk(kk)
        do k = 1,dz
          do j=jstart,jend
            do i=1,n1m
              bb(k+offsetz,i,j) = rbuf(nc)
              nc=nc+1
            enddo
          enddo
        enddo
      enddo
      if(allocated(sbuf)) deallocate(sbuf)
      if(allocated(rbuf)) deallocate(rbuf)
     
      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)
     
      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
      end subroutine PackZ_UnpackR
 
!============================================
      subroutine PackR_UnpackZ(aa,bb)
      use mpih
      use mpi_param
      use param, only: n1m,n2m,n3m
      implicit none
      real, intent(in) :: aa(1:n3m,1:n1m,jstart:jend)
      real, intent(out) :: bb(1:n1m,1:n2m,kstart:kend)
      real, allocatable :: sbuf(:),rbuf(:)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:) 
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr

      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dk* countj(i)*n1m
        aak(i)= dj* countk(i)*n1m
      end do
       
      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1)) 
      
      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do
      
      if(.not. allocated(rbuf)) allocate(rbuf(0:n1m*n2m*dk-1),stat=merr)

      if(merr .ne. 0) then
        write(*,*)"process  ",myid," failed to allocate memory for sbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      
      if(.not. allocated(sbuf)) allocate(sbuf(0:n1m*n3m*dj-1),stat=merr)

      if(merr .ne. 0) then
        write(*,*)"process  ",myid," failed to allocate memory for rbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 

      do kk = 0, numtasks-1
        nc = dispk(kk)
        dz= countk(kk)
        offsetz = offsetk(kk)
        do k = 1,dz
          do j=jstart,jend
            do i=1,n1m
              sbuf(nc) = aa(k+offsetz,i,j) 
              nc=nc+1
            enddo
          enddo
        enddo
      enddo
 
      call MPI_ALLTOALLV(sbuf, aak,dispk, MDP,rbuf, aaj,dispj, MDP,MPI_COMM_WORLD, ierr)
     
      do kk = 0, numtasks-1
        nc = dispj(kk)
        dr= countj(kk)
        offsetr = offsetj(kk)
        do k = kstart,kend
          do j=1,dr
            do i=1,n1m
              bb(i,j+offsetr,k) = rbuf(nc)
              nc=nc+1
            enddo
          enddo
        enddo
      enddo

      if(allocated(sbuf)) deallocate(sbuf)
      if(allocated(rbuf)) deallocate(rbuf)
     
      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)
     
      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
      
      end subroutine PackR_UnpackZ

!===========================================
      subroutine PackZ_UnpackR_refi(aa,bb)
      use mpih
      use mpi_param
      use param, only: n1mr,n2mr,n3mr
      implicit none
      real, intent(in) :: aa(1:n1mr,1:n2mr,kstartr:kendr)
      real, intent(out) :: bb(1:n3mr,1:n1mr,jstartr:jendr)
      real, allocatable :: sbuf(:),rbuf(:)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:)
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr

      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dkr* countjr(i)*n1mr
        aak(i)= djr* countkr(i)*n1mr
      end do
      
      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1)) 

      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do
     
      if(.not. allocated(sbuf)) then
         allocate(sbuf(0:n1mr*n2mr*dkr-1),stat=merr)
      end if
 
      if(merr .ne. 0) then
        write(*,*)"process  ",myid," failed to allocate memory for sbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      
      if(.not. allocated(rbuf)) then
        allocate(rbuf(0:n1mr*n3mr*djr-1),stat=merr)
      end if
      
      if(merr .ne. 0) then
        write(*,*)"process  ",myid," failed to allocate memory for sbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 

      do kk = 0, numtasks-1
        nc = dispj(kk)
        dr = countjr(kk)
        offsetr = offsetjr(kk)
        do k = kstartr,kendr
          do j=1,dr
            do i=1,n1mr
              sbuf(nc) = aa(i,j+offsetr,k)
              nc=nc+1
            enddo
          enddo
        enddo
      enddo

      call MPI_ALLTOALLV(sbuf, aaj, dispj, MDP,rbuf, aak, dispk, MDP,MPI_COMM_WORLD, ierr)
     
      do kk = 0, numtasks-1
        nc = dispk(kk)
        dz = countkr(kk)
        offsetz = offsetkr(kk)
        do k = 1,dz
          do j=jstartr,jendr
            do i=1,n1mr
              bb(k+offsetz,i,j) = rbuf(nc)
              nc=nc+1
            enddo
          enddo
        enddo
      enddo

      if(allocated(sbuf)) deallocate(sbuf)
      if(allocated(rbuf)) deallocate(rbuf)
     
      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)
     
      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
      
      end subroutine PackZ_UnpackR_refi
 
!============================================
      subroutine PackR_UnpackZ_refi(aa,bb)
      use mpih
      use mpi_param
      use param, only: n1mr,n2mr,n3mr
      implicit none
      real, intent(in) ::  aa(1:n3mr,1:n1mr,jstartr:jendr)
      real, intent(out) :: bb(1:n1mr,1:n2mr,kstartr:kendr)
      real, allocatable :: sbuf(:),rbuf(:)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:) 
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr

      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dkr* countjr(i)*n1mr
        aak(i)= djr* countkr(i)*n1mr
      end do
       
      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1)) 
      
      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do
      
      if(.not. allocated(rbuf)) then 
       allocate(rbuf(0:n1mr*n2mr*dkr-1),stat=merr)
      end if

      if(merr .ne. 0) then
        write(*,*)"process  ",myid," failed to allocate memory for sbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      
      if(.not. allocated(sbuf)) then 
       allocate(sbuf(0:n1mr*n3mr*djr-1),stat=merr)
      end if

      if(merr .ne. 0) then
        write(*,*)"process  ",myid," failed to allocate memory for rbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 

      do kk = 0, numtasks-1
        nc = dispk(kk)
        dz = countkr(kk)
        offsetz = offsetkr(kk)
        do k = 1,dz
          do j=jstartr,jendr
            do i=1,n1mr
              sbuf(nc) = aa(k+offsetz,i,j)
              nc=nc+1
            enddo
          enddo
        enddo
      enddo
      
      call MPI_ALLTOALLV(sbuf, aak, dispk, MDP,rbuf, aaj, dispj, MDP,MPI_COMM_WORLD, ierr)
     
      do kk = 0, numtasks-1
        nc = dispj(kk)
        dr = countjr(kk)
        offsetr = offsetjr(kk)
        do k = kstartr,kendr
          do j=1,dr
            do i=1,n1mr
              bb(i,j+offsetr,k) = rbuf(nc)
              nc=nc+1
            enddo
          enddo
        enddo
      enddo

      if(allocated(sbuf)) deallocate(sbuf)
      if(allocated(rbuf)) deallocate(rbuf)
     
      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)
     
      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
      
      end subroutine PackR_UnpackZ_refi

!==========================================
     subroutine PackZ_UnpackRP(aa,bb)
      use mpih
      use mpi_param
      use param, only: n1m,n2m,n3m
      implicit none
      real,intent(in) :: aa(1:n1m,1:n2m+2,kstart:kend)
      real,intent(out) :: bb(1:n3m,1:n1m,jstartp:jendp)
      real, allocatable :: sbuf(:),rbuf(:)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:) 
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr

      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dk* countjp(i)*n1m
        aak(i)= djp* countk(i)*n1m
      end do
      
      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1)) 

      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do
 
     
      if(.not. allocated(sbuf)) allocate(sbuf(0:n1m*(n2m+2)*dk-1),stat=merr)

      if(.not. allocated(rbuf)) allocate(rbuf(0:n1m*n3m*djp-1),stat=merr)

      do kk = 0, numtasks-1
        nc = dispj(kk)
        dr = countjp(kk)
        offsetr = offsetjp(kk)
        do k = kstart,kend
          do j=1,dr
            do i=1,n1m
              sbuf(nc) = aa(i,j+offsetr,k)
              nc=nc+1
            enddo
          enddo
        enddo
      enddo
      
      call MPI_ALLTOALLV(sbuf, aaj,dispj, MDP,rbuf, aak,dispk, MDP,MPI_COMM_WORLD, ierr)
     
     
      do kk = 0, numtasks-1
        nc = dispk(kk)
        dz = countk(kk)
        offsetz = offsetk(kk)
        do k = 1,dz
          do j=jstartp,jendp
            do i=1,n1m
              bb(k+offsetz,i,j) = rbuf(nc)
              nc=nc+1
            enddo
          enddo
        enddo
      enddo

      if(allocated(sbuf)) deallocate(sbuf)
      if(allocated(rbuf)) deallocate(rbuf)
     
      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)
     
      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
      
      end subroutine PackZ_UnpackRP
 
!============================================
      subroutine PackR_UnpackZP(aa,bb)
      use mpih
      use mpi_param
      use param, only: n1m,n2m,n3m
      implicit none
      real, intent(in) :: aa(1:n3m,1:n1m,jstartp:jendp)
      real, intent(out) :: bb(1:n1m,1:n2m+2,kstart:kend)
      real, allocatable :: sbuf(:),rbuf(:)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:) 
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr

      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dk* countjp(i)*n1m
        aak(i)= djp* countk(i)*n1m
      end do
       
      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1)) 
      
      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do
      
      if(.not. allocated(rbuf)) allocate(rbuf(0:n1m*(n2m+2)*dk-1),stat=merr)

      if(merr .ne. 0) then
        write(*,*)"process  ",myid," failed to allocate memory for sbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      
      if(.not. allocated(sbuf)) allocate(sbuf(0:n1m*n3m*djp-1),stat=merr)

      if(merr .ne. 0) then
        write(*,*)"process  ",myid," failed to allocate memory for rbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      
      do kk = 0, numtasks-1
        nc = dispk(kk)
        dz = countk(kk)
        offsetz = offsetk(kk)
        do k = 1,dz
          do j=jstartp,jendp
            do i=1,n1m
              sbuf(nc) = aa(k+offsetz,i,j) 
              nc=nc+1
            enddo
          enddo
        enddo
      enddo
      
      call MPI_ALLTOALLV(sbuf, aak,dispk, MDP, rbuf, aaj,dispj, MDP,MPI_COMM_WORLD, ierr)
     
      do kk = 0, numtasks-1
        nc = dispj(kk)
        dr = countjp(kk)
        offsetr = offsetjp(kk)
        do k = kstart,kend
          do j=1,dr
            do i=1,n1m
              bb(i,j+offsetr,k) = rbuf(nc)
              nc=nc+1
            enddo
          enddo
        enddo
      enddo

      if(allocated(sbuf)) deallocate(sbuf)
      if(allocated(rbuf)) deallocate(rbuf)
     
      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)
     
      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
      
      end subroutine PackR_UnpackZP
 
!============================================
      subroutine PackZ_UnpackR_refi_withghost(aa,bb)
      use mpih
      use mpi_param
      use param, only: n1mr,n2mr,n3mr
      implicit none
      real, intent(in) :: aa(1:n1mr,1:n2mr,kstartr-lvlhalo:kendr+lvlhalo)
      real, intent(out) :: bb(1:n3mr,1:n1mr,jstartr-lvlhalo:jendr+lvlhalo)
      real, allocatable :: sbuf(:),rbuf(:)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:)
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr

      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dkr* countjr(i)*n1mr
        aak(i)= djr* countkr(i)*n1mr
      end do
      
      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1)) 

      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do
     
      if(.not. allocated(sbuf)) then
         allocate(sbuf(0:n1mr*n2mr*dkr-1),stat=merr)
      end if

      if(.not. allocated(rbuf)) then
        allocate(rbuf(0:n1mr*n3mr*(djr+2*lvlhalo)-1),stat=merr)
      end if

      do kk = 0, numtasks-1
        nc = dispj(kk)
        dr = countjr(kk)
        offsetr = offsetjr(kk)
        do k = kstartr,kendr
          do j=1,dr
            do i=1,n1mr
              sbuf(nc) = aa(i,j+offsetr,k)
              nc=nc+1
            enddo
          enddo
        enddo
      enddo

      call MPI_ALLTOALLV(sbuf, aaj, dispj, MDP,rbuf, aak, dispk, MDP,MPI_COMM_WORLD, ierr)
     
      do kk = 0, numtasks-1
        nc = dispk(kk)
        dz = countkr(kk)
        offsetz = offsetkr(kk)
        do k = 1,dz
          do j=jstartr-lvlhalo,jendr+lvlhalo
            do i=1,n1mr
              bb(k+offsetz,i,j) = rbuf(nc)
              nc=nc+1
            enddo
          enddo
        enddo
      enddo

      if(allocated(sbuf)) deallocate(sbuf)
      if(allocated(rbuf)) deallocate(rbuf)
     
      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)
     
      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
      
      end subroutine PackZ_UnpackR_refi_withghost
