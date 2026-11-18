subroutine solqvap
      use param
      use local_arrays, only : qvap,rhsr,kpa
      use mpi_param
      use mpih
      implicit none
      integer :: jc,kc,ic,info,ipkv(m3mr)
      real :: betadx,fkl(m3mr),ackl_b(m3m)
      real :: amkT(m3mr-1),ackT(m3mr),apkT(m3mr-1)
      real, allocatable, dimension(:,:,:) :: rhst,kpat
      real :: K_p,K_c,K_m
      integer :: km,kp

      allocate(rhst(1:n3mr,1:n1mr,jstartr:jendr))
      allocate(kpat(1:n3mr,1:n1mr,jstartr:jendr))

      call PackZ_UnpackR_refi(rhsr,rhst)
      call PackZ_UnpackR_refi(kpa,kpat)

      betadx=0.5d0*al*dts
     

      do jc=jstartr,jendr
        do ic=1,n1mr
          do kc=1,n3mr
            km=kmvr(kc)
            kp=kpvr(kc)
            if(kc.eq.1) then        
               K_p = Sm*(kpat(kp,ic,jc)*(zcr(kp)-zmr(kc)) + kpat(kc,ic,jc)*(zmr(kp)-zcr(kp)))/(zmr(kp)-zmr(kc))
               K_m = Sm*(kpat(kc,ic,jc) + kpat(km,ic,jc))*0.5d0
            elseif( (kc.ge.2) .and. (kc.le.n3mr-1) ) then
               K_p = Sm*(kpat(kp,ic,jc)*(zcr(kp)-zmr(kc)) + kpat(kc,ic,jc)*(zmr(kp)-zcr(kp)))/(zmr(kp)-zmr(kc))
               K_m = Sm*(kpat(kc,ic,jc)*(zcr(kc)-zmr(km)) + kpat(km,ic,jc)*(zmr(kc)-zcr(kc)))/(zmr(kc)-zmr(km))
            elseif(kc.eq.n3mr) then
               K_p = Sm*(kpat(kp,ic,jc) + kpat(kc,ic,jc))*0.5d0
               K_m = Sm*(kpat(kc,ic,jc)*(zcr(kc)-zmr(km)) + kpat(km,ic,jc)*(zmr(kc)-zcr(kc)))/(zmr(kc)-zmr(km))
            endif
            K_c = (K_p+K_m)

            ackl_b(kc)=1.d0/(1.d0-ac3sskr(kc)*0.5*betadx*K_c)
          enddo

          do kc=1,n3mr-1
            km=kmvr(kc)
            kp=kpvr(kc)
            if(kc.eq.1) then        
               K_p = Sm*(kpat(kp,ic,jc)*(zcr(kp)-zmr(kc)) + kpat(kc,ic,jc)*(zmr(kp)-zcr(kp)))/(zmr(kp)-zmr(kc))
               K_m = Sm*(kpat(kc,ic,jc) + kpat(km,ic,jc))*0.5d0
            elseif( (kc.ge.2) .and. (kc.le.n3mr-1) ) then
               K_p = Sm*(kpat(kp,ic,jc)*(zcr(kp)-zmr(kc)) + kpat(kc,ic,jc)*(zmr(kp)-zcr(kp)))/(zmr(kp)-zmr(kc))
               K_m = Sm*(kpat(kc,ic,jc)*(zcr(kc)-zmr(km)) + kpat(km,ic,jc)*(zmr(kc)-zcr(kc)))/(zmr(kc)-zmr(km))
            elseif(kc.eq.n3mr) then
               K_p = Sm*(kpat(kp,ic,jc) + kpat(kc,ic,jc))*0.5d0
               K_m = Sm*(kpat(kc,ic,jc)*(zcr(kc)-zmr(km)) + kpat(km,ic,jc)*(zmr(kc)-zcr(kc)))/(zmr(kc)-zmr(km))
            endif
            amkT(kc)=-K_m*betadx*am3sskr(kp)*ackl_b(kp)
            apkT(kc)=-K_p*betadx*ap3sskr(kc)*ackl_b(kc)
          enddo

          ackT(1) = 1.d0+apkT(1)
          ackT(2:m3mr-1) = 1.d0
          ackT(m3mr) = 1.d0 + amkT(m3mr-1)

          call ddttrfb(n3mr,amkT,ackT,apkT,info)

          fkl(1:n3mr)=rhst(1:n3mr,ic,jc)*ackl_b(1:n3mr)

          call ddttrsb('N',n3mr,1,amkT,ackT,apkT,fkl,n3mr,info)

          rhst(1:n3mr,ic,jc)=fkl(1:n3mr)

        end do
      end do

      call PackR_UnpackZ_refi(rhst, rhsr)

      do kc=kstartr,kendr
        do jc=1,n2mr
          qvap(1:n1mr,jc,kc) = qvap(1:n1mr,jc,kc) + rhsr(1:n1mr,jc,kc)
        enddo
      enddo

      if(allocated(rhst)) deallocate(rhst)
      if(allocated(kpat)) deallocate(kpat)

      return
      end
