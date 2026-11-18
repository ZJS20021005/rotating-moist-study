      subroutine solxrj_periodic(betadx)
      use param
      use local_arrays, only : rhsr
      use mpi_param, only: kstartr,kendr
      implicit none
      integer :: jc,kc,ic,info
      real,intent(in) :: betadx
      real :: amjl(m2mr-1),apjl(m2mr-1),acjl(m2mr)
      real :: ackl_b, aoffdiag, fjl(m2mr)
      real :: qpd(m2mr), vtq, ytq   ! for periodic modification

      ackl_b = 1.d0/(1.d0+2.d0*betadx)
      aoffdiag = -betadx*ackl_b

      amjl(1:n2mr-1)=aoffdiag
      apjl(1:n2mr-1)=aoffdiag

      acjl(1) = 2.d0
      acjl(2:n2mr-1) = 1.d0
      acjl(n2mr) = 1.d0 + aoffdiag*aoffdiag

      call ddttrfb(n2mr,amjl,acjl,apjl,info)

      qpd(1) = -1.d0
      qpd(2:n2mr-1) = 0.d0
      qpd(n2mr) = aoffdiag

      call ddttrsb('N',n2mr,1,amjl,acjl,apjl,qpd,n2mr,info)

      vtq = 1.d0/(1.d0+qpd(1)-aoffdiag*qpd(n2mr))
      qpd(1:n2mr) = qpd(1:n2mr)*vtq

      do kc=kstartr,kendr
        do ic=1,n1mr
          fjl(1:n2mr)=rhsr(ic,1:n2mr,kc)*ackl_b

          call ddttrsb('N',n2mr,1,amjl,acjl,apjl,fjl,n2mr,info)

          ytq = fjl(1) - aoffdiag*fjl(n2mr)
          rhsr(ic,1:n2mr,kc) = fjl(1:n2mr) - ytq*qpd(1:n2mr)
        end do
      end do 

      return
      end
