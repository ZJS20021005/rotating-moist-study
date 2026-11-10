      subroutine solxrj(betadx)
      use param
      use local_arrays, only : rhsr
      use mpi_param, only: kstartr,kendr
      implicit none
      integer :: jc,kc,ic,info
      real,intent(in) :: betadx
      real :: amjl(m2mr-1),apjl(m2mr-1),acjl(m2mr)
      real :: fjl(m2mr,m1mr)

      amjl(1:n2mr-1)=-betadx
      apjl(1:n2mr-1)=-betadx

      acjl(1) = (1.d0+1.d0*betadx)
      acjl(2:n2mr-1) = (1.d0+2.d0*betadx)
      acjl(n2mr) = (1.d0+1.d0*betadx)

      call ddttrfb(n2mr,amjl,acjl,apjl,info)

      do kc=kstartr,kendr
        do ic=1,n1mr
          do jc=1,n2mr
            fjl(jc,ic)=rhsr(ic,jc,kc)
          enddo
        end do

        call ddttrsb('N',n2mr,n1mr,amjl,acjl,apjl,fjl,n2mr,info)

        do ic=1,n1mr
          do jc=1,n2mr
            rhsr(ic,jc,kc) = fjl(jc,ic)
          enddo
        end do
      end do 

      return
      end
