      subroutine solxi_fs(betadx)
      use param
      use local_arrays, only : rhs
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc,ic,info
      real,intent(in) :: betadx
      real :: amil(m1m-1),apil(m1m-1),acil(m1m)
      real :: fil(m1m,m2m)

      amil(1:n1m-1)=-betadx
      apil(1:n1m-1)=-betadx

      acil(1) = (1.d0+1.d0*betadx)
      acil(2:n1m-1) = (1.d0+2.d0*betadx)
      acil(n1m) = (1.d0+1.d0*betadx)

      call ddttrfb(n1m,amil,acil,apil,info)

      do kc=kstart,kend
        do jc=1,n2m
          do ic=1,n1m
            fil(ic,jc)=rhs(ic,jc,kc)
          enddo
        end do

        call ddttrsb('N',n1m,n2m,amil,acil,apil,fil,n1m,info)

        do jc=1,n2m
          do ic=1,n1m
            rhs(ic,jc,kc) = fil(ic,jc)
          enddo
        end do
      end do 

      return
      end
