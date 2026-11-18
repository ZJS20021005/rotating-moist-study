      subroutine solq3k
      use param
      use local_arrays, only : q3,rhs
      use mpi_param
      implicit none
      integer :: jc,kc,ic,info
      real :: betadx,fkl(m3,m1m)
      real :: amkT(m3-1),ackT(m3),apkT(m3-1)
      real, allocatable, dimension(:,:,:) :: rhst
      allocate(rhst(1:n3m,1:n1m,jstart:jend))

      call PackZ_UnpackR(rhs,rhst)

      betadx=beta*al

      ackT(1) = 1.0d0
      do kc=2,n3m
        ackT(kc) = (1.d0-ac3ck(kc)*betadx)
      enddo
      ackT(n3) = 1.0d0

      amkT(1) = 1.0d0
      amkT(2:n3m-1)=-betadx
      amkT(n3m) = 0.0d0

      apkT(1) = 0.0d0
      apkT(2:n3m-1)=-betadx
      apkT(n3m) = 1.0d0


      call ddttrfb(n3,amkT,ackT,apkT,info)

      do jc=jstart,jend
        do ic=1,n1m
          fkl(1,ic)= 0.d0
          fkl(2:n3m,ic) = rhst(2:n3m,ic,jc)
          fkl(n3,ic)= 0.d0
        enddo
          

        call ddttrsb('N',n3,n1m,amkT,ackT,apkT,fkl,n3,info)

        do ic=1,n1m
          rhst(1:n3m,ic,jc)= fkl(1:n3m,ic)
        enddo
      end do

      call PackR_UnpackZ(rhst,rhs) 

      do kc=kstart,kend
        do jc=1,n2m
          do ic=1,n1m
              q3(ic,jc,kc) =  q3(ic,jc,kc) + rhs(ic,jc,kc)
          enddo
        enddo
      enddo

      if(allocated(rhst)) deallocate(rhst)

      return
      end
