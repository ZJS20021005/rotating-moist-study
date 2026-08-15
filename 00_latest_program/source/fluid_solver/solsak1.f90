      subroutine solsak
      use param
      use local_arrays, only : dsal,rhsr,kpa
      use mpi_param
      use mpih
      implicit none
      integer :: jc,kc,ic,info,ipkv(m3mr)
      real :: betadx,fkl(m3mr),ackl_b(m3m)
      real :: amkT(m3mr-1),ackT(m3mr),apkT(m3mr-1)
      real, allocatable, dimension(:,:,:) :: rhst,kpat
      real :: K_p,K_c,K_m

      allocate(rhst(1:n3mr,1:n1mr,jstartr:jendr))
      allocate(kpat(1:n3mr,1:n1mr,jstartr:jendr))

      call PackZ_UnpackR_refi(rhsr,rhst)
      call PackZ_UnpackR_refi(kpa,kpat)

      betadx=0.5d0*al*dts
     

      do jc=jstartr,jendr
        do ic=1,n1mr
          do kc=1,n3mr

            K_p = (kpat(kc+1,ic,jc) + kpat(kc,ic,jc))*0.5d0
            K_m = (kpat(kc,ic,jc) + kpat(kc-1,ic,jc))*0.5d0
            K_c = (kpat(kc+1,ic,jc) + 2.0d0*kpat(kc,ic,jc) + kpat(kc-1,ic,jc))*0.25d0
            ackl_b(kc)=1.d0/(1.d0-ac3sskr(kc)*betadx*K_c) ! question

          enddo

          do kc=1,n3mr-1
            amkT(kc)=-K_m*betadx*am3sskr(kc+1)*ackl_b(kc+1)
            apkT(kc)=-K_p*betadx*ap3sskr(kc)*ackl_b(kc)
          enddo

          ackT(1) = 1.d0+apkT(1)
          ackT(2:m3mr-1) = 1.d0
          ackT(m3mr) = 1.d0 + amkT(m3mr-1)

!          ackT(1:n3mr) = 1.d0

          call ddttrfb(n3mr,amkT,ackT,apkT,info)

          fkl(1:n3mr)=rhst(1:n3mr,ic,jc)*ackl_b(1:n3mr)

          call ddttrsb('N',n3mr,1,amkT,ackT,apkT,fkl,n3mr,info)

          rhst(1:n3mr,ic,jc)=fkl(1:n3mr)

        end do
      end do

      call PackR_UnpackZ_refi(rhst, rhsr)

      do kc=kstartr,kendr
        do jc=1,n2mr
          dsal(1:n1mr,jc,kc) = dsal(1:n1mr,jc,kc) + rhsr(1:n1mr,jc,kc)
        enddo
      enddo

      if(allocated(rhst)) deallocate(rhst)
      if(allocated(kpat)) deallocate(kpat)

      return
      end
