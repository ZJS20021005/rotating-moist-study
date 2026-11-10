      subroutine invtrqvap
      use param
  use local_arrays, only:q2,qvap,ruqvap,rhsr,kpa,hqvap
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
      real    :: qvaptop_val,qvapbot_val

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
            K_p = 0.5d0*(Sm*kpa(ip,jc,kc)+Sm*kpa(ic,jc,kc))
            K_c = 0.5d0*(Sm*kpa(ip,jc,kc)+2.0d0*Sm*kpa(ic,jc,kc)+Sm*kpa(im,jc,kc))
            K_m = 0.5d0*(Sm*kpa(ic,jc,kc)+Sm*kpa(im,jc,kc))

            dq31=(K_p*qvap(ip,jc,kc)-K_c*qvap(ic,jc,kc) +K_m*qvap(im,jc,kc))*udx1qr

            !-- y
            K_p = 0.5d0*(Sm*kpa(ic,jp,kc)+Sm*kpa(ic,jc,kc))
            K_c = 0.5d0*(Sm*kpa(ic,jp,kc)+2.0d0*Sm*kpa(ic,jc,kc)+Sm*kpa(ic,jm,kc))
            K_m = 0.5d0*(Sm*kpa(ic,jc,kc)+Sm*kpa(ic,jm,kc))

            dq32=(K_p*qvap(ic,jp,kc)-K_c*qvap(ic,jc,kc)+K_m*qvap(ic,jm,kc))*udx2qr


            !-- z
            if(kc.eq.1) then        
               del1 = zmr(1)-zcr(1)
               del2 = zmr(2)-zmr(1)
               del3 = zcr(2)-zcr(1)
               cmup = (Sm*kpa(ic,jc,2)*(zcr(2)-zmr(1))+Sm*kpa(ic,jc,1)*(zmr(2)-zcr(2)))/(zmr(2)-zmr(1))
               cmum = (Sm*kpa(ic,jc,1)+Sm*kpa(ic,jc,1))*0.5 ! for adiabatic condition

               qvapbot_val = qvapbot
               dq33=((qvap(ic,jc,kp) - qvap(ic,jc,kc))*cmup/del2 &
                  -  (qvap(ic,jc,kc) - qvapbot_val   )*cmum/del1)/del3
            elseif( (kc.ge.2) .and. (kc.le.n3mr-1) ) then
               cmup = (Sm*kpa(ic,jc,kp)*(zcr(kp)-zmr(kc))+Sm*kpa(ic,jc,kc)*(zmr(kp)-zcr(kp)))/(zmr(kp)-zmr(kc))
               cmum = (Sm*kpa(ic,jc,kc)*(zcr(kc)-zmr(km))+Sm*kpa(ic,jc,km)*(zmr(kc)-zcr(kc)))/(zmr(kc)-zmr(km))

               dq33=(qvap(ic,jc,kp)*app + qvap(ic,jc,kc)*acc*0.5d0)*cmup &
                 +  (qvap(ic,jc,km)*amm + qvap(ic,jc,kc)*acc*0.5d0)*cmum
            elseif(kc.eq.n3mr) then
               del1 = zcr(n3r)-zmr(n3mr)
               del2 = zmr(n3mr)-zmr(n3mr-1)
               del3 = zcr(n3r)-zcr(n3mr)
               cmup = (Sm*kpa(ic,jc,n3mr)+Sm*kpa(ic,jc,n3mr))*0.5  ! for adiabatic condition
               cmum = (Sm*kpa(ic,jc,n3mr)*(zcr(n3mr)-zmr(n3mr-1))+Sm*kpa(ic,jc,n3mr-1)*(zmr(n3mr)-zcr(n3mr)))/(zmr(n3mr)-zmr(n3mr-1))

               qvaptop_val = qvaptop
               dq33=((qvaptop_val    - qvap(ic,jc,kc))*cmup/del1 &
                  -  (qvap(ic,jc,kc) - qvap(ic,jc,km))*cmum/del2)/del3
            endif

            !-- sum
            dcq3=dq32+dq33+dq31

            rhsr(ic,jc,kc)=(ga*hqvap(ic,jc,kc)+ro*ruqvap(ic,jc,kc)+alpecl*dcq3)*dts
            ruqvap(ic,jc,kc)=hqvap(ic,jc,kc)
            enddo
        enddo
      enddo

      !-- implicit matrix solver
  call solxri_periodic( alpecl*Sm*kps*dts*0.5d0*dx1qr ) ! no variable kpa
  call solxrj_periodic( alpecl*Sm*kps*dts*0.5d0*dx2qr ) ! no variable kpa
  call solqvap

      return
      end
