      subroutine solxj_endwall(betadx)
      use param
      use local_arrays, only : rhs
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc,ic,info
      real,intent(in) :: betadx
      real :: amjl(m2m-1),apjl(m2m-1),acjl(m2m)
      real :: fjl(m2m,m1m)

      amjl(1:n2m-1)=-betadx
      apjl(1:n2m-1)=-betadx
      apjl(1) = 0.0d0

      acjl(1) = 1.d0
      acjl(2:n2m-1) = (1.d0+2.d0*betadx)
      acjl(n2m) = (1.d0+2.d0*betadx)


      call ddttrfb(n2m,amjl,acjl,apjl,info)

      do kc=kstart,kend
        do ic=1,n1m
          do jc=1,n2m
            if(jc.eq.1) then
               fjl(jc,ic)=0.0d0
            else
               fjl(jc,ic)=rhs(ic,jc,kc)
            endif
          enddo
        end do

        call ddttrsb('N',n2m,n1m,amjl,acjl,apjl,fjl,n2m,info)

        do ic=1,n1m
          do jc=1,n2m
            rhs(ic,jc,kc) = fjl(jc,ic)
          enddo
        end do
      end do 

      return
      end
