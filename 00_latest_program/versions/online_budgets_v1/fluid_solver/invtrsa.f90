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
      real    :: del1, del2, del3
      real    :: cmup,cmum
      real    :: dsaltop_val,dsalbot_val

      alpecl=al
      udx1qr=dx1qr
      udx2qr=dx2qr
      del2 = zmr(2)-zmr(1)

      do kc=kstartr,kendr
      !-- inner
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

            !-- x
            K_p = 0.5d0*(kpa(ip,jc,kc)+kpa(ic,jc,kc))
            K_c = 0.5d0*(kpa(ip,jc,kc)+2.0d0*kpa(ic,jc,kc)+kpa(im,jc,kc))
            K_m = 0.5d0*(kpa(ic,jc,kc)+kpa(im,jc,kc))

            dq31=(K_p*dsal(ip,jc,kc)-K_c*dsal(ic,jc,kc) +K_m*dsal(im,jc,kc))*udx1qr

            !-- y
            K_p = 0.5d0*(kpa(ic,jp,kc)+kpa(ic,jc,kc))
            K_c = 0.5d0*(kpa(ic,jp,kc)+2.0d0*kpa(ic,jc,kc)+kpa(ic,jm,kc))
            K_m = 0.5d0*(kpa(ic,jc,kc)+kpa(ic,jm,kc))

            dq32=(K_p*dsal(ic,jp,kc)-K_c*dsal(ic,jc,kc)+K_m*dsal(ic,jm,kc))*udx2qr


            !-- z
            if(kc.eq.1) then        
               del1 = zmr(1)-zcr(1)
               del2 = zmr(2)-zmr(1)
               del3 = zcr(2)-zcr(1)
               cmup = (kpa(ic,jc,2)*(zcr(2)-zmr(1))+kpa(ic,jc,1)*(zmr(2)-zcr(2)))/(zmr(2)-zmr(1))
               cmum = (kpa(ic,jc,1)+kpa(ic,jc,1))*0.5 ! for adiabatic condition

               dsalbot_val = dsalbot +  A_sbotmod*sin(2.0*pi*k_sbotmod*ymr(jc))*sin(2.0*pi*k_sbotmod*xmr(ic))
               dq33=((dsal(ic,jc,kp) - dsal(ic,jc,kc))*cmup/del2 &
                  -  (dsal(ic,jc,kc) - dsalbot_val   )*cmum/del1)/del3
            elseif( (kc.ge.2) .and. (kc.le.n3mr-1) ) then
               cmup = (kpa(ic,jc,kp)*(zcr(kp)-zmr(kc))+kpa(ic,jc,kc)*(zmr(kp)-zcr(kp)))/(zmr(kp)-zmr(kc))
               cmum = (kpa(ic,jc,kc)*(zcr(kc)-zmr(km))+kpa(ic,jc,km)*(zmr(kc)-zcr(kc)))/(zmr(kc)-zmr(km))

               dq33=(dsal(ic,jc,kp)*app + dsal(ic,jc,kc)*acc*0.5d0)*cmup &
                 +  (dsal(ic,jc,km)*amm + dsal(ic,jc,kc)*acc*0.5d0)*cmum
            elseif(kc.eq.n3mr) then
               del1 = zcr(n3r)-zmr(n3mr)
               del2 = zmr(n3mr)-zmr(n3mr-1)
               del3 = zcr(n3r)-zcr(n3mr)
               cmup = (kpa(ic,jc,n3mr)+kpa(ic,jc,n3mr))*0.5  ! for adiabatic condition
               cmum = (kpa(ic,jc,n3mr)*(zcr(n3mr)-zmr(n3mr-1))+kpa(ic,jc,n3mr-1)*(zmr(n3mr)-zcr(n3mr)))/(zmr(n3mr)-zmr(n3mr-1))

               dsaltop_val = dsaltop +  A_stopmod*sin(2.0*pi*k_stopmod*ymr(jc))
               dq33=((dsaltop_val    - dsal(ic,jc,kc))*cmup/del1 &
                  -  (dsal(ic,jc,kc) - dsal(ic,jc,km))*cmum/del2)/del3
            endif

            !-- sum
            dcq3=dq32+dq33+dq31

            rhsr(ic,jc,kc)=(ga*hsal(ic,jc,kc)+ro*rusal(ic,jc,kc)+alpecl*dcq3)*dts
            rusal(ic,jc,kc)=hsal(ic,jc,kc)
            enddo
        enddo
      enddo

      !-- implicit matrix solver
      call solxri_periodic( alpecl*kps*dts*0.5d0*dx1qr ) ! no variable kpa
      call solxrj_periodic( alpecl*kps*dts*0.5d0*dx2qr ) ! no variable kpa
      call solsak

      return
      end
