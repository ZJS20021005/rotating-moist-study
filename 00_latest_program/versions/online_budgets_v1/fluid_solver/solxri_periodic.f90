      subroutine solxri_periodic(betadx)
      use param
      use local_arrays, only : rhsr
      use mpi_param, only: kstartr,kendr
      implicit none
      integer :: jc,kc,ic,info
      real,intent(in) :: betadx
      real :: amil(m1mr-1),apil(m1mr-1),acil(m1mr)
      real :: ackl_b,aoffdiag,fil(m1mr)
      real :: qpd(m1mr), vtq, ytq

      ackl_b = 1.d0/(1.d0+2.d0*betadx)
      aoffdiag = -betadx*ackl_b

      amil(1:n1mr-1)=aoffdiag
      apil(1:n1mr-1)=aoffdiag

      acil(1) = 2.d0
      acil(2:n1mr-1) = 1.d0
      acil(n1mr) = 1.d0 + aoffdiag*aoffdiag

      call ddttrfb(n1mr,amil,acil,apil,info)

      qpd(1) = -1.d0
      qpd(2:n1mr-1) = 0.d0
      qpd(n1mr) = aoffdiag

      call ddttrsb('N',n1mr,1,amil,acil,apil,qpd,n1mr,info)

      vtq = 1.d0/(1.d0+qpd(1)-aoffdiag*qpd(n1mr))
      qpd(1:n1mr) = qpd(1:n1mr)*vtq

      do kc=kstartr,kendr
        do jc=1,n2mr
          fil(1:n1mr)=rhsr(1:n1mr,jc,kc)*ackl_b

          call ddttrsb('N',n1mr,1,amil,acil,apil,fil,n1mr,info)

          ytq = fil(1)-aoffdiag*fil(n1mr)
          rhsr(1:n1mr,jc,kc) = fil(1:n1mr) - ytq*qpd(1:n1mr)
        end do
      end do 

      return
      end
