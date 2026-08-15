      subroutine avgvar
!EP   This routine calculates the avg of variables
      use param
      use local_arrays, only: q2,q3,q1,dsal,qvap,T,kpa
      use mgrd_arrays, only: q3lr
      use mpi_param, only: kstart,kend,kstartr,kendr
      use mpih
      implicit none
!
!     Online statistics written by this routine:
!
!     unit 96, data/avgvar.out
!       1 time
!       2 trms            : volume standard deviation of T
!       3 Re              : sqrt(Ra/Pr)*sqrt(<u^2+v^2+w^2>)
!       4 vrms            : volume mean <u^2+v^2+w^2>
!       5 vrm_horizontal  : volume mean <u^2+v^2>
!       6 moistenergy     : volume mean <m>, where m=b+gamma*qvap
!       7 wb              : volume mean <w*b>, kept as the old wb output
!       8 L0mid           : 2*pi weighted horizontal-velocity length at z~0.5
!       9 L0_075          : 2*pi weighted horizontal-velocity length at z~0.75
!      10 nuh             : volume mean T flux, convfac*w*T-dT/dz
!
!     unit 99, data/nu_profiles.out
!       1 time, 2 z_face, 3 NuT(z), 4 Num(z), 5 Nuq(z)
!       NuT = convfac*w*T - kappa*dT/dz
!       Num = convfac*w*m - kappa*db/dz - gamma*Sm*kappa*dq/dz
!       Nuq = convfac*w*q - Sm*kappa*dq/dz
!
!     unit 100, data/kh_energy.out
!       1 time, 2-13 shell energy for kh=1..12 at z~0.5
!       The spectrum uses horizontal velocity components q1,q2 only.
!
      integer, parameter :: nkh = 12
      integer :: jc,kc,kp,ic
      integer :: kface,kl,ku
      integer :: kxi,kyi,ikh
      integer :: kmid,k075,n1mh,n2mh
      real    :: cellvol,my_liquidcellvol,liquidcellvol
!     The speed average was kept in older versions but is not output now.
!      real    :: totalcellvol,my_spdavg,spdavg
!      real    :: spdval
      real    :: my_tavg,tavg
      real    :: wt
!     nuv was an unused volume convective contribution; nuh is the output flux.
!      real    :: my_nuv,nuv
      real    :: my_nuh,nuh
      real    :: my_trms,trms
      real    :: my_vrms,vrms,vrm,Re,my_vrm_horizontal,vrm_horizontal,my_moistenergy,moistenergy,wb,mywb
      real    :: my_l0mid,l0mid,my_l075,l075
      real    :: num_l0,den_l0,num_l075,den_l075
      real    :: coefnorm,kx,ky,khat
      real    :: dtdz,heatflux
      real    :: dsalbot_val,dsaltop_val
      real    :: qvapbot_val,qvaptop_val
      real    :: cellarea,convfac
      real    :: bcell,bup,bdown,bbot,btop
      real    :: qcell,qup,qdown
      real    :: tcell,tup,tdown,tbot,ttop
      real    :: mcell,bface,qface,tface,mface
      real    :: dbdz,dqdz,dtempdz
      real    :: dz,kfacecoef
      real    :: my_nut_prof(n3r),nut_prof(n3r)
      real    :: my_num_prof(n3r),num_prof(n3r)
      real    :: my_nuq_prof(n3r),nuq_prof(n3r)
      real    :: my_area_prof(n3r),area_prof(n3r)
      real    :: xr_u(m2m,m1m),xr_v(m2m,m1m)
      complex :: xa_u(m2mh,m1m),xa_v(m2mh,m1m)
      real    :: spec_e
      real    :: my_kh_energy(nkh),kh_energy(nkh)
!      totalcellvol = 0.0d0
!      spdavg       = 0.0d0
!      my_spdavg    = 0.0d0
      my_tavg      = 0.0d0
      my_vrms  = 0.0
      my_trms = 0.0
!      my_nuv= 0.0
      my_nuh= 0.0
      my_vrm_horizontal = 0.0d0
      my_moistenergy = 0.0d0
      my_liquidcellvol = 0.0d0
      mywb = 0.0d0
      my_l0mid = 0.0d0
      my_l075 = 0.0d0
      my_kh_energy = 0.0d0
      my_nut_prof = 0.0d0
      my_num_prof = 0.0d0
      my_nuq_prof = 0.0d0
      my_area_prof = 0.0d0
      convfac = dsqrt(Ra*Prs)
      kmid = 1
      k075 = 1
      do kc=2,n3m
        if(dabs(zm(kc)-0.5d0*alx3).lt.dabs(zm(kmid)-0.5d0*alx3)) kmid=kc
        if(dabs(zm(kc)-0.75d0*alx3).lt.dabs(zm(k075)-0.75d0*alx3)) k075=kc
      enddo
      do kc=kstart,kend
        kp = kc + 1
        do jc=1,n2m
          do ic=1,n1m
            cellarea = (xc(ic+1)-xc(ic))*(yc(jc+1)-yc(jc))
            cellvol  = (xc(ic+1)-xc(ic))*(yc(jc+1)-yc(jc))*(zc(kc+1)-zc(kc))
!             spdval = dsqrt(Ra/Prs) &
!     &              *dsqrt(q1(ic,jc,kc)**2.0+q2(ic,jc,kc)**2.0+q3(ic,jc,kc)**2.0)
             bcell = dsal(ic,jc,kc)
             qcell = qvap(ic,jc,kc)
             tcell = T(ic,jc,kc)
             mcell = bcell + gamma*qcell
             wt = convfac*tcell*q3(ic,jc,kc)
             qvapbot_val = qvapbot
             qvaptop_val = qvaptop
             if(kc.eq.1)then
               dsalbot_val = dsalbot + A_sbotmod*sin(2.0*pi*k_sbotmod*ym(jc))*sin(2.0*pi*k_sbotmod*xm(ic))
               bbot = dsalbot_val
               tbot = bbot - betaqs*zc(1)
               dbdz = (bcell-bbot)/(zm(kc)-zc(1))
               dqdz = (qcell-qvapbot_val)/(zm(kc)-zc(1))
               dtempdz = (tcell-tbot)/(zm(kc)-zc(1))
               dtdz = dtempdz
             elseif(kc.eq.n3m)then
               dsaltop_val = dsaltop + A_stopmod*sin(2.0*pi*k_stopmod*ym(jc))
               btop = dsaltop_val
               ttop = btop - betaqs*zc(n3)
               dbdz = (btop-bcell)/(zc(n3)-zm(kc))
               dqdz = (qvaptop_val-qcell)/(zc(n3)-zm(kc))
               dtempdz = (ttop-tcell)/(zc(n3)-zm(kc))
               dtdz = dtempdz
             else
               bup = dsal(ic,jc,kc+1)
               bdown = dsal(ic,jc,kc-1)
               qup = qvap(ic,jc,kc+1)
               qdown = qvap(ic,jc,kc-1)
               tup = T(ic,jc,kc+1)
               tdown = T(ic,jc,kc-1)
               dbdz = (bup-bdown)/(zm(kc+1)-zm(kc-1))
               dqdz = (qup-qdown)/(zm(kc+1)-zm(kc-1))
               dtempdz = (tup-tdown)/(zm(kc+1)-zm(kc-1))
               dtdz = dtempdz
             endif
             heatflux = wt - dtdz
             vrm     = q1(ic,jc,kc)**2.0+q2(ic,jc,kc)**2.0+q3(ic,jc,kc)**2.0
             vrm_horizontal =q1(ic,jc,kc)**2.0+q2(ic,jc,kc)**2.0
             moistenergy = mcell
             wb = q3(ic,jc,kc)*dsal(ic,jc,kc)
             my_vrm_horizontal = my_vrm_horizontal + vrm_horizontal*cellvol
             my_liquidcellvol = my_liquidcellvol + cellvol
!             my_spdavg       = my_spdavg       + spdval*cellvol
             my_vrms       = my_vrms       + vrm*cellvol
             my_tavg       = my_tavg       + tcell*cellvol
!              my_nuv  = my_nuv + wt*cellvol
              my_nuh  = my_nuh + heatflux*cellvol
              my_trms  = my_trms  + tcell**2*cellvol
              my_moistenergy = my_moistenergy + moistenergy*cellvol
             mywb = mywb + wb*cellvol
          enddo
        enddo
      enddo

      do kface=kstartr,kendr
        do jc=1,n2mr
          do ic=1,n1mr
            cellarea = (xcr(ic+1)-xcr(ic))*(ycr(jc+1)-ycr(jc))
            if(kface.eq.1)then
              dsalbot_val = dsalbot + A_sbotmod*sin(2.0*pi*k_sbotmod*ymr(jc))*sin(2.0*pi*k_sbotmod*xmr(ic))
              bbot = dsalbot_val
              qvapbot_val = qvapbot
              tbot = bbot - betaqs*zcr(1)
              bcell = dsal(ic,jc,1)
              qcell = qvap(ic,jc,1)
              tcell = T(ic,jc,1)
              dz = zmr(1)-zcr(1)
              kfacecoef = kpa(ic,jc,1)/kps
              bface = 0.5d0*(bbot+bcell)
              qface = 0.5d0*(qvapbot_val+qcell)
              tface = 0.5d0*(tbot+tcell)
              dbdz = (bcell-bbot)/dz
              dqdz = (qcell-qvapbot_val)/dz
              dtempdz = (tcell-tbot)/dz
            else
              kl = kface-1
              ku = kface
              dz = zmr(ku)-zmr(kl)
              kfacecoef = (kpa(ic,jc,ku)*(zcr(kface)-zmr(kl)) &
     &          + kpa(ic,jc,kl)*(zmr(ku)-zcr(kface)))/dz/kps
              bface = (dsal(ic,jc,ku)*g3rmr(kl)+dsal(ic,jc,kl)*g3rmr(ku)) &
     &          /(g3rmr(kl)+g3rmr(ku))
              qface = (qvap(ic,jc,ku)*g3rmr(kl)+qvap(ic,jc,kl)*g3rmr(ku)) &
     &          /(g3rmr(kl)+g3rmr(ku))
              tface = (T(ic,jc,ku)*g3rmr(kl)+T(ic,jc,kl)*g3rmr(ku)) &
     &          /(g3rmr(kl)+g3rmr(ku))
              dbdz = (dsal(ic,jc,ku)-dsal(ic,jc,kl))/dz
              dqdz = (qvap(ic,jc,ku)-qvap(ic,jc,kl))/dz
              dtempdz = (T(ic,jc,ku)-T(ic,jc,kl))/dz
            endif
            mface = bface + gamma*qface
            my_nut_prof(kface) = my_nut_prof(kface) &
     &        + (convfac*q3lr(ic,jc,kface)*tface - kfacecoef*dtempdz)*cellarea
            my_num_prof(kface) = my_num_prof(kface) &
     &        + (convfac*q3lr(ic,jc,kface)*mface &
     &        - kfacecoef*dbdz - gamma*Sm*kfacecoef*dqdz)*cellarea
            my_nuq_prof(kface) = my_nuq_prof(kface) &
     &        + (convfac*q3lr(ic,jc,kface)*qface - Sm*kfacecoef*dqdz)*cellarea
            my_area_prof(kface) = my_area_prof(kface) + cellarea
          enddo
        enddo
      enddo

      if(myid.eq.numtasks-1)then
        kface = n3r
        do jc=1,n2mr
          do ic=1,n1mr
            cellarea = (xcr(ic+1)-xcr(ic))*(ycr(jc+1)-ycr(jc))
            dsaltop_val = dsaltop + A_stopmod*sin(2.0*pi*k_stopmod*ymr(jc))
            btop = dsaltop_val
            qvaptop_val = qvaptop
            ttop = btop - betaqs*zcr(n3r)
            bcell = dsal(ic,jc,n3mr)
            qcell = qvap(ic,jc,n3mr)
            tcell = T(ic,jc,n3mr)
            dz = zcr(n3r)-zmr(n3mr)
            kfacecoef = kpa(ic,jc,n3mr)/kps
            bface = 0.5d0*(bcell+btop)
            qface = 0.5d0*(qcell+qvaptop_val)
            tface = 0.5d0*(tcell+ttop)
            dbdz = (btop-bcell)/dz
            dqdz = (qvaptop_val-qcell)/dz
            dtempdz = (ttop-tcell)/dz
            mface = bface + gamma*qface
            my_nut_prof(kface) = my_nut_prof(kface) &
     &        + (convfac*q3lr(ic,jc,kface)*tface - kfacecoef*dtempdz)*cellarea
            my_num_prof(kface) = my_num_prof(kface) &
     &        + (convfac*q3lr(ic,jc,kface)*mface &
     &        - kfacecoef*dbdz - gamma*Sm*kfacecoef*dqdz)*cellarea
            my_nuq_prof(kface) = my_nuq_prof(kface) &
     &        + (convfac*q3lr(ic,jc,kface)*qface - Sm*kfacecoef*dqdz)*cellarea
            my_area_prof(kface) = my_area_prof(kface) + cellarea
          enddo
        enddo
      endif

      if(kmid.ge.kstart .and. kmid.le.kend)then
        n1mh = n1m/2 + 1
        n2mh = n2m/2 + 1
        coefnorm = 1.d0/(dble(n1m)*dble(n2m))

        do jc=1,n2m
          do ic=1,n1m
            xr_u(jc,ic)=q1(ic,jc,kmid)
            xr_v(jc,ic)=q2(ic,jc,kmid)
          enddo
        enddo

        call dfftw_execute_dft_r2c(fwd_plan,xr_u,xa_u)
        call dfftw_execute_dft_r2c(fwd_plan,xr_v,xa_v)

        num_l0 = 0.d0
        den_l0 = 0.d0
        do jc=1,n2mh
          ky = 2.d0*pi*dble(jc-1)/rext2
          do ic=1,n1m
            if(ic.le.n1mh)then
              kx = 2.d0*pi*dble(ic-1)/rext1
              kxi = ic-1
            else
              kx = -2.d0*pi*dble(n1m-ic+1)/rext1
              kxi = -(n1m-ic+1)
            endif
            kyi = jc-1
            khat = dsqrt(kx**2 + ky**2)
            if(khat.gt.0.d0)then
              spec_e = coefnorm**2 * (cdabs(xa_u(jc,ic))**2 + cdabs(xa_v(jc,ic))**2)
              ikh = nint(dsqrt(dble(kxi*kxi + kyi*kyi)))
              if(ikh.ge.1 .and. ikh.le.nkh)then
                my_kh_energy(ikh) = my_kh_energy(ikh) + spec_e
              endif
              num_l0 = num_l0 + spec_e/khat
              den_l0 = den_l0 + spec_e
            endif
          enddo
        enddo

        if(den_l0.gt.0.d0) my_l0mid = 2.d0*pi*num_l0/den_l0
      endif

      if(k075.ge.kstart .and. k075.le.kend)then
        n1mh = n1m/2 + 1
        n2mh = n2m/2 + 1
        coefnorm = 1.d0/(dble(n1m)*dble(n2m))

        do jc=1,n2m
          do ic=1,n1m
            xr_u(jc,ic)=q1(ic,jc,k075)
            xr_v(jc,ic)=q2(ic,jc,k075)
          enddo
        enddo

        call dfftw_execute_dft_r2c(fwd_plan,xr_u,xa_u)
        call dfftw_execute_dft_r2c(fwd_plan,xr_v,xa_v)

        num_l075 = 0.d0
        den_l075 = 0.d0
        do jc=1,n2mh
          ky = 2.d0*pi*dble(jc-1)/rext2
          do ic=1,n1m
            if(ic.le.n1mh)then
              kx = 2.d0*pi*dble(ic-1)/rext1
            else
              kx = -2.d0*pi*dble(n1m-ic+1)/rext1
            endif
            khat = dsqrt(kx**2 + ky**2)
            if(khat.gt.0.d0)then
              spec_e = coefnorm**2 * (cdabs(xa_u(jc,ic))**2 + cdabs(xa_v(jc,ic))**2)
              num_l075 = num_l075 + spec_e/khat
              den_l075 = den_l075 + spec_e
            endif
          enddo
        enddo

        if(den_l075.gt.0.d0) my_l075 = 2.d0*pi*num_l075/den_l075
      endif

      call MPI_ALLREDUCE(my_liquidcellvol,liquidcellvol,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
!      call MPI_ALLREDUCE(my_spdavg,spdavg,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_tavg,tavg,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_vrms,vrms,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
!      call MPI_ALLREDUCE(my_nuv,nuv,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_nuh,nuh,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_trms,trms,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_vrm_horizontal,vrm_horizontal,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_moistenergy,moistenergy,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(mywb,wb,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_l0mid,l0mid,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_l075,l075,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_kh_energy,kh_energy,nkh,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_nut_prof,nut_prof,n3r,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_num_prof,num_prof,n3r,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_nuq_prof,nuq_prof,n3r,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_area_prof,area_prof,n3r,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
!      spdavg = spdavg/liquidcellvol
      tavg = tavg/liquidcellvol
!      nuv  = 1+nuv/liquidcellvol
      nuh  = nuh/liquidcellvol
      trms  = trms/liquidcellvol
      vrms  = vrms/liquidcellvol
      trms = dsqrt(dmax1(0.0d0, trms - tavg**2))
      Re = dsqrt(Ra/Prs)*dsqrt(vrms)
      moistenergy = moistenergy/liquidcellvol
      wb = wb/liquidcellvol
      if(myid.eq.0)then
        write(96,510) time, trms, Re, vrms, vrm_horizontal,moistenergy,wb,l0mid,l075,nuh
        write(100,512) time, kh_energy(1:nkh)
        do kc=1,n3r
          if(area_prof(kc).gt.0.0d0)then
            write(99,511) time, zcr(kc), nut_prof(kc)/area_prof(kc), &
     &        num_prof(kc)/area_prof(kc), nuq_prof(kc)/area_prof(kc)
          endif
        enddo
      endif
 510   format(1x,f10.4,9(1x,ES20.8))
 511   format(1x,f10.4,4(1x,ES20.8))
 512   format(1x,f10.4,12(1x,ES20.8))

      return   
      end
