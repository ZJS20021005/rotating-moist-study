      subroutine solxri(betadx)
      use param
      use local_arrays, only : rhsr,kpa
      use mpi_param, only: kstartr,kendr
      implicit none
      integer :: jc,kc,ic,info
      real,intent(in) :: betadx
      real :: amil(m1mr-1),apil(m1mr-1),acil(m1mr)
      real :: ackl_b, amoffdiag(m1mr-1),apoffdiag(m1mr-1),fil(m1mr)
      real :: K_p,K_c,K_m
      integer :: im,ip

      do kc=kstartr,kendr
        do jc=1,m2mr
          do ic=1,m1mr-1
            im=imvr(ic)
            ip=ipvr(ic)
            K_p = (kpa(ip,jc,kc) + kpa(ic,jc,kc))*0.5d0
            K_c = (kpa(im,jc,kc) + 2.0d0*kpa(ic,jc,kc) + kpa(ip,jc,kc))*0.5
            K_m = (kpa(ic,jc,kc) + kpa(im,jc,kc))*0.5d0

            ackl_b = 1.d0/(1.d0+betadx*K_c)
            amoffdiag(ic) = -K_m*betadx*ackl_b
            apoffdiag(ic) = -K_p*betadx*ackl_b

            amil(ic) = amoffdiag(ic)
            apil(ic) = apoffdiag(ic)

          end do

          acil(1) = 1.d0+apil(1)
          acil(2:m1mr-1) = 1.d0
          acil(m1mr) = 1.d0 + amil(m1mr-1)

          call ddttrfb(n1mr,amil,acil,apil,info)

          do ic=1,m1mr
            im=imvr(ic)
            ip=ipvr(ic)
            K_c = (kpa(im,jc,kc) + 2.0d0*kpa(ic,jc,kc) + kpa(ip,jc,kc))*0.5
            ackl_b = 1.d0/(1.d0+betadx*K_c)
            fil(ic) = rhsr(ic,jc,kc)*ackl_b
          end do

          call ddttrsb('N',m1mr,1,amil,acil,apil,fil,m1mr,info)

          rhsr(1:m1mr,jc,kc) = fil(1:m1mr)

        end do
      end do

      return
      end
