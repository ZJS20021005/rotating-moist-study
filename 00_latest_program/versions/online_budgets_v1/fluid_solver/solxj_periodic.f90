      subroutine solxj_periodic(betadx)
      use param
      use local_arrays, only : rhs
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc,ic,info
      real,intent(in) :: betadx
      real :: amjl(m2m-1),apjl(m2m-1),acjl(m2m)
      real :: ackl_b, aoffdiag, fjl(m2m)
      real :: qpd(m2m), vtq, ytq   ! for periodic modification

      ackl_b = 1.d0/(1.d0+2.d0*betadx)
      aoffdiag = -betadx*ackl_b

      amjl(1:n2m-1)=aoffdiag
      apjl(1:n2m-1)=aoffdiag

      acjl(1) = 2.d0
      acjl(2:n2m-1) = 1.d0
      acjl(n2m) = 1.d0 + aoffdiag*aoffdiag

      call ddttrfb(n2m,amjl,acjl,apjl,info)

      qpd(1) = -1.d0
      qpd(2:n2m-1) = 0.d0
      qpd(n2m) = aoffdiag

      call ddttrsb('N',n2m,1,amjl,acjl,apjl,qpd,n2m,info)

      vtq = 1.d0/(1.d0+qpd(1)-aoffdiag*qpd(n2m))
      qpd(1:n2m) = qpd(1:n2m)*vtq

      do kc=kstart,kend
        do ic=1,n1m
          fjl(1:n2m)=rhs(ic,1:n2m,kc)*ackl_b

          call ddttrsb('N',n2m,1,amjl,acjl,apjl,fjl,n2m,info)

          ytq = fjl(1) - aoffdiag*fjl(n2m)
          rhs(ic,1:n2m,kc) = fjl(1:n2m) - ytq*qpd(1:n2m)
        end do
      end do 

      return
      end
