      subroutine invtrsa
      use param
      use local_arrays, only:q2,dsal,hsal,rusal,rhsr,kpa
      use mpi_param, only: kstartr,kendr
      use mpih
      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,ip,im
      real    :: dq32,dq33,dcq3,dq31
      real    :: app,acc,amm
      real    :: K_p, K_c, K_m
      real    :: alpecl,udx1qr,udx2qr
      real    :: del1, del2, fcder
      real    :: dsaltop_outflow

      alpecl=al !*kps
      udx1qr=dx1qr
      udx2qr=dx2qr
      del2 = zmr(2)-zmr(1)

      do kc=kstartr,kendr
      !-- inner
      if( (kc.ge.2) .and. (kc.le.n3mr-1) ) then
        km=kmvr(kc)
        kp=kpvr(kc)
        app=ap3sskr(kc)
        acc=ac3sskr(kc)
        amm=am3sskr(kc)
        do jc=1,n2mr
          jm=jmvr(jc)
          jp=jpvr(jc)
          do ic=1,n1mr
            im=imvr(ic)
            ip=ipvr(ic)
            K_p = 0.5d0*(kpa(ip,jc,kc)+kpa(ic,jc,kc))
            K_c = 0.5d0*(kpa(ip,jc,kc)+2.0d0*kpa(ic,jc,kc)+kpa(im,jc,kc))
            K_m = 0.5d0*(kpa(ic,jc,kc)+kpa(im,jc,kc))

            dq31=(K_p*dsal(ip,jc,kc)-K_c*dsal(ic,jc,kc) +K_m*dsal(im,jc,kc))*udx1qr

            K_p = 0.5d0*(kpa(ic,jp,kc)+kpa(ic,jc,kc))
            K_c = 0.5d0*(kpa(ic,jp,kc)+2.0d0*kpa(ic,jc,kc)+kpa(ic,jm,kc))
            K_m = 0.5d0*(kpa(ic,jc,kc)+kpa(ic,jm,kc))

            dq32=(K_p*dsal(ic,jp,kc)-K_c*dsal(ic,jc,kc)+K_m*dsal(ic,jm,kc))*udx2qr

!-----------only for even grid ----
            K_p = 0.5d0*(kpa(ic,jc,kp)+kpa(ic,jc,kc))
            K_c = 0.5d0*(kpa(ic,jc,kp)+2.0d0*kpa(ic,jc,kc)+kpa(ic,jc,km))
            K_m = 0.5d0*(kpa(ic,jc,kc)+kpa(ic,jc,km))

            dq33= (K_p*dsal(ic,jc,kp)-K_c*dsal(ic,jc,kc)+K_m*dsal(ic,jc,km))/del2/del2
!            dq33= dsal(ic,jc,kp)*app+dsal(ic,jc,kc)*acc+dsal(ic,jc,km)*amm

            dcq3=dq32+dq33+dq31

            rhsr(ic,jc,kc)=(ga*hsal(ic,jc,kc)+ro*rusal(ic,jc,kc)+alpecl*dcq3)*dts
            rusal(ic,jc,kc)=hsal(ic,jc,kc)
            enddo
        enddo
      endif

      !-- bottom
      if(kc.eq.1) then
        del1 = zmr(1)-zcr(1)
        del2 = zmr(2)-zmr(1)
        fcder = 2.d0/(del1*del2*(del1+del2))
        kp = kc + 1
        do jc=1,n2mr
          jm=jmvr(jc)
          jp=jpvr(jc)
          do ic=1,n1mr
            im=imvr(ic)
            ip=ipvr(ic)
            dq31=(dsal(ip,jc,kc)-2.d0*dsal(ic,jc,kc)+dsal(im,jc,kc))*udx1qr
            dq32=(dsal(ic,jp,kc)-2.d0*dsal(ic,jc,kc)+dsal(ic,jm,kc))*udx2qr
            dq33=(dsal(ic,jc,kp)*del1-dsal(ic,jc,kc)*(del1+del2*dble(sbcbot))+dsalbot*del2*dble(sbcbot))*fcder
            dcq3=dq32+dq33+dq31
            rhsr(ic,jc,kc)=(ga*hsal(ic,jc,kc)+ro*rusal(ic,jc,kc) +alpecl*kps*dcq3)*dts
            rusal(ic,jc,kc)=hsal(ic,jc,kc)
          enddo
        enddo
      endif

      !-- top
      if(kc.eq.n3mr) then
        del1 = zcr(n3r)-zmr(n3mr)
        del2 = zmr(n3mr)-zmr(n3mr-1)
        fcder = 2.d0/(del1*del2*(del1+del2))
        km = kc - 1
        do jc=1,n2mr
          jm=jmvr(jc)
          jp=jpvr(jc)
          do ic=1,n1mr
            im=imvr(ic)
            ip=ipvr(ic)
            dq31=( dsal(ip,jc,kc)-2.d0*dsal(ic,jc,kc)+dsal(im,jc,kc))*udx1qr
            dq32=(dsal(ic,jp,kc)-2.d0*dsal(ic,jc,kc)+dsal(ic,jm,kc))*udx2qr
            dq33=(dsal(ic,jc,km)*del1-dsal(ic,jc,kc)*(del1+del2*dble(sbctop))+dsaltop*del2*dble(sbctop))*fcder
            dcq3=dq32+dq33+dq31

            rhsr(ic,jc,kc)=(ga*hsal(ic,jc,kc)+ro*rusal(ic,jc,kc)+alpecl*kps*dcq3)*dts
            rusal(ic,jc,kc)=hsal(ic,jc,kc)
          enddo
        enddo
      endif
      enddo

!      do kc=kstartr,kendr
!        do jc=1,n2mr
!          do ic=1,n1mr
!            dsal(ic,jc,kc) = dsal(ic,jc,kc) + rhsr(ic,jc,kc)
!          enddo
!        enddo
!      enddo


!!      -- implicit matrix solver
      call solxri( alpecl*dts*0.5d0*dx1qr )
      call solxrj( alpecl*dts*0.5d0*dx2qr )
      call solsak

      return
      end
