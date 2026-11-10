      subroutine solxi_periodic(betadx)
      use param
      use local_arrays, only : rhs
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc,ic,info
      real,intent(in) :: betadx
      real :: amil(m1m-1),apil(m1m-1),acil(m1m)
      real :: ackl_b, aoffdiag, fil(m1m)
      real :: qpd(m1m), vtq, ytq   ! for periodic modification

      ackl_b = 1.d0/(1.d0+2.d0*betadx)
      aoffdiag = -betadx*ackl_b

      amil(1:m1m-1)=aoffdiag
      apil(1:m1m-1)=aoffdiag

      acil(1) = 2.d0
      acil(2:m1m-1) = 1.d0
      acil(m1m) = 1.d0 + aoffdiag*aoffdiag

      call ddttrfb(m1m,amil,acil,apil,info)

      qpd(1) = -1.d0
      qpd(2:m1m-1) = 0.d0
      qpd(m1m) = aoffdiag

      call ddttrsb('N',m1m,1,amil,acil,apil,qpd,m1m,info)

      vtq = 1.d0/(1.d0+qpd(1)-aoffdiag*qpd(m1m))
      qpd(1:m1m) = qpd(1:m1m)*vtq

      do kc=kstart,kend
        do jc=1,m2m
          fil(1:m1m)=rhs(1:m1m,jc,kc)*ackl_b

          call ddttrsb('N',m1m,1,amil,acil,apil,fil,m1m,info)

          ytq = fil(1)-aoffdiag*fil(m1m)
          rhs(1:m1m,jc,kc) = fil(1:m1m) - ytq*qpd(1:m1m)
        end do
      end do 

      return
      end
