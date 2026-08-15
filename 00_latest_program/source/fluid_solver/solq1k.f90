      subroutine solq1k
      use param
      use local_arrays, only : q1, rhs
      use mpi_param
      use mpih
      implicit none
      integer :: ic,jc,kc,info,ipkv(m3m)
      real :: betadx, fkl(m3m,m1m)
      real :: amkT(m3m-1),ackT(m3m),apkT(m3m-1)
      real,allocatable :: rhst(:,:,:)

      allocate(rhst(1:n3m,1:n1m,jstart:jend))

      call PackZ_UnpackR(rhs,rhst)

      betadx=beta*al

      do kc=1,n3m
        ackT(kc)=(1.d0-ac3sk(kc)*betadx)
      enddo
    
      amkT(1:n3m-1) = -betadx*am3sk(2:n3m)
      apkT(1:n3m-1) = -betadx*ap3sk(1:n3m-1)

      call ddttrfb(n3m,amkT,ackT,apkT,info)

      do jc=jstart,jend
        do ic=1,n1m
          do kc=1,n3m
            fkl(kc,ic)=rhst(kc,ic,jc)
          end do
        enddo

        call ddttrsb('N',n3m,n1m,amkT,ackT,apkT,fkl,n3m,info)

        do ic=1,n1m
          do kc=1,n3m
            rhst(kc,ic,jc)=fkl(kc,ic)
          end do
        enddo
      end do

      call PackR_UnpackZ(rhst,rhs)

      do kc=kstart,kend
        do jc=1,n2m
          do ic=1,n1m
              q1(ic,jc,kc) =  q1(ic,jc,kc) + rhs(ic,jc,kc)
          enddo
        enddo
      enddo

      if(allocated(rhst)) deallocate(rhst)

      return
      end
