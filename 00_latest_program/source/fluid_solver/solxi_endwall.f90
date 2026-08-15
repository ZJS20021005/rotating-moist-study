      subroutine solxi_endwall(betadx)
      use param
      use local_arrays, only : rhs
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc,ic,info
      real,intent(in) :: betadx
      real :: amil(m1m),apil(m1m),acil(m1)
      real :: fil(m1,m2m)

      acil(1) = 1.0d0
      acil(2:n1m) = (1.d0+2.d0*betadx)
      acil(n1) = 1.0d0

      amil(1) = 1.0d0
      amil(2:n1m-1)=-betadx
      amil(n1m) = 0.0d0

      apil(1) = 0.0d0
      apil(2:n1m-1)=-betadx
      apil(n1m) = 1.0d0


      call ddttrfb(n1,amil,acil,apil,info)

      do kc=kstart,kend
        do jc=1,n2m
          do ic=1,n1
            if(ic.eq.1) then
               fil(ic,jc)=0.0d0
            elseif(ic.eq.n1) then
               fil(ic,jc)=0.0d0
            else
               fil(ic,jc)=rhs(ic,jc,kc)
            endif
          enddo
        end do

        call ddttrsb('N',n1,n2m,amil,acil,apil,fil,n1,info)

        do jc=1,n2m
          do ic=1,n1m
            rhs(ic,jc,kc) = fil(ic,jc)
          enddo
        end do
      end do 

      return
      end
