      subroutine solxrj(betadx)
      use param
      use local_arrays, only : rhsr,kpa
      use mpi_param, only: kstartr,kendr
      implicit none
      integer :: jc,kc,ic,info
      real,intent(in) :: betadx
      real :: amil(m2mr-1),apil(m2mr-1),acil(m2mr)
      real :: ackl_b, amoffdiag(m2mr-1),apoffdiag(m2mr-1),fil(m2mr)
      real :: K_p,K_c,K_m
      integer :: jp,jm

      do kc=kstartr,kendr
         do ic=1,m1mr
           do jc=1,m2mr-1
            jm=jmvr(jc)
            jp=jpvr(jc)
            K_p = (kpa(ic,jp,kc) + kpa(ic,jc,kc))*0.5d0
            K_c = (kpa(ic,jm,kc) + 2.0d0*kpa(ic,jc,kc) + kpa(ic,jp,kc))*0.5
            K_m = (kpa(ic,jc,kc) + kpa(ic,jm,kc))*0.5d0

            ackl_b = 1.d0/(1.d0+betadx*K_c)
            amoffdiag(jc) = -K_m*betadx*ackl_b
            apoffdiag(jc) = -K_p*betadx*ackl_b

            amil(jc) = amoffdiag(jc)
            apil(jc) = apoffdiag(jc)

          end do

          acil(1) = 1.d0+apil(1)
          acil(2:m2mr-1) = 1.d0
          acil(m2mr) = 1.d0 + amil(m2mr-1)

          call ddttrfb(n2mr,amil,acil,apil,info)

          do jc=1,m2mr
            jm=jmvr(jc)
            jp=jpvr(jc)
            K_c = (kpa(ic,jm,kc) + 2.0d0*kpa(ic,jc,kc) + kpa(ic,jp,kc))*0.5
            ackl_b = 1.d0/(1.d0+betadx*K_c)
            fil(jc) = rhsr(ic,jc,kc)*ackl_b
          end do

          call ddttrsb('N',m2mr,1,amil,acil,apil,fil,m2mr,info)

          rhsr(ic,1:m2mr,kc) = fil(1:m2mr)

        end do
      end do

      return
      end
