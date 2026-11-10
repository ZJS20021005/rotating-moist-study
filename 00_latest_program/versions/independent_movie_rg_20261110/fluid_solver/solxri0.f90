      subroutine solxri(betadx)
      use param
      use local_arrays, only : rhsr
      use mpi_param, only: kstartr,kendr
      implicit none
      integer :: jc,kc,ic,info
      real,intent(in) :: betadx
      real :: amil(m1mr-1),apil(m1mr-1),acil(m1mr)
      real :: fil(m1mr,m2mr)

      amil(1:n1mr-1)=-betadx
      apil(1:n1mr-1)=-betadx

      acil(1) = (1.d0+1.d0*betadx)
      acil(2:n1mr-1) = (1.d0+2.d0*betadx)
      acil(n1mr) = (1.d0+1.d0*betadx)

      call ddttrfb(n1mr,amil,acil,apil,info)

      do kc=kstartr,kendr
        do jc=1,n2mr
          do ic=1,n1mr
            fil(ic,jc)=rhsr(ic,jc,kc)
          enddo
        end do

        call ddttrsb('N',n1mr,n2mr,amil,acil,apil,fil,n1mr,info)

        do jc=1,n2mr
          do ic=1,n1mr
            rhsr(ic,jc,kc) = fil(ic,jc)
          enddo
        end do
      end do 

      return
      end
