      subroutine solsak
      use param
      use local_arrays, only : dsal,rhsr
      use mpi_param
      use mpih
      implicit none
      integer :: jc,kc,ic,info,ipkv(m3mr)
      real :: betadx,fkl(m3mr,m1mr)
      real :: amkT(m3mr-1),ackT(m3mr),apkT(m3mr-1)
      real, allocatable, dimension(:,:,:) :: rhst

      allocate(rhst(1:n3mr,1:n1mr,jstartr:jendr))

      call PackZ_UnpackR_refi(rhsr,rhst)

      betadx=0.5d0*al*dts*kps
     
      do kc=1,n3mr
        ackT(kc)=(1.d0-ac3sskr(kc)*betadx)
      enddo

      amkT(1:n3mr-1)=-betadx*am3sskr(2:n3mr)
      apkT(1:n3mr-1)=-betadx*ap3sskr(1:n3mr-1)

      call ddttrfb(n3mr,amkT,ackT,apkT,info)

      do jc=jstartr,jendr
        do ic=1,n1mr
          do kc=1,n3mr
            fkl(kc,ic)=rhst(kc,ic,jc)
          end do
        enddo

       call ddttrsb('N',n3mr,n1mr,amkT,ackT,apkT,fkl,n3mr,info)
          
        do ic=1,n1mr
          do kc=1,n3mr
            rhst(kc,ic,jc)= fkl(kc,ic)
          end do
        enddo
      end do

      call PackR_UnpackZ_refi(rhst,rhsr)

      do kc=kstartr,kendr
        do jc=1,n2mr
          do ic=1,n1mr
            dsal(ic,jc,kc) = dsal(ic,jc,kc) + rhsr(ic,jc,kc)
          enddo
        enddo
      enddo

      if(allocated(rhst)) deallocate(rhst)

      return
      end
