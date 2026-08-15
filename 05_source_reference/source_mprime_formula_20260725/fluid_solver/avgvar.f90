      subroutine avgvar
!EP   This routine calculates the avg of variables
      use param
      use local_arrays, only: q2,q3,q1,dens,dsal,qvap,T,kpa
      use mgrd_arrays, only: q3lr
      use mpi_param, only: kstart,kend,kstartr,kendr
      use mpih
      implicit none
      integer, parameter :: nkh = 12
      integer :: jc,kc,kp,ic
      integer :: kface,kl,ku
      integer :: kxi,kyi,ikh
      integer :: kmid,k075,n1mh,n2mh
      integer :: kwlo,kwhi
      real    :: cellvol,totalcellvol,my_liquidcellvol,liquidcellvol
      real    :: my_spdavg,spdavg
      real    :: spdval
      real    :: my_tavg,tavg
      real    :: my_nuv,nuv,wt
      real    :: my_nuh,nuh
      real    :: my_trms,trms
      real    :: my_vrms,vrms,vrm,Re,my_vrm_horizontal,vrm_horizontal,my_moistenergy,moistenergy,wb,mywb
      real    :: my_l0mid,l0mid,my_l075,l075
      real    :: num_l0,den_l0,num_l075,den_l075
      real    :: my_wspec_num,wspec_num,my_wspec_den,wspec_den
      real    :: w_length_z075,w_wavelength_z075
      real    :: coefnorm,kx,ky,khat,spec_weight
      real    :: target_wz,wz_weight
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
      real    :: my_msum_prof(n3m),msum_prof(n3m)
      real    :: my_m2_prof(n3m),m2_prof(n3m)
      real    :: my_marea_prof(n3m),marea_prof(n3m)
      real    :: mbar_prof(n3m),mvar_prof(n3m)
      real    :: my_m2d(m2m,m1m),m2d(m2m,m1m)
      real    :: mse_var_horizontal,mse_rms_horizontal
      real    :: mprime2_global,mprime2_z075
      real    :: zweight,total_zweight
      real    :: xr_u(m2m,m1m),xr_v(m2m,m1m),xr_w(m2m,m1m)
      real    :: xr_m(m2m,m1m),xr_corr(m2m,m1m)
      complex :: xa_u(m2mh,m1m),xa_v(m2mh,m1m),xa_w(m2mh,m1m)
      complex :: xa_m(m2mh,m1m),xa_corr(m2mh,m1m)
      real    :: spec_e
      real    :: my_kh_energy(nkh),kh_energy(nkh)
      real    :: mse_shell(n1m+n2m)
      real    :: corr_sum(n1m+n2m),corr_count(n1m+n2m)
      real    :: corr_profile(n1m+n2m)
      real    :: mse_l_peak,mse_l_integral
      real    :: mse_l_moment,mse_lambda_moment
      real    :: mse_spec_num,mse_spec_den,peak_energy
      real    :: dk_shell,dr_corr,dx_m,dy_m
      real    :: rx,ry,radius,corr_center,corr_prev,corr_now
      real    :: r_prev,r_zero,m2d_mean,m2d_sigma
      real    :: cluster_mean0,cluster_max0
      real    :: cluster_mean_sigma,cluster_max_sigma
      integer :: shell_index,peak_shell,bin_index,max_corr_bin
      integer :: cluster_count0,cluster_count_sigma



      totalcellvol = 0.0d0
      spdavg       = 0.0d0
      my_spdavg       = 0.0d0
      my_tavg      = 0.0d0
      my_vrms  = 0.0
      my_trms = 0.0
      my_nuv= 0.0
      my_nuh= 0.0
      my_vrm_horizontal = 0.0d0
      my_moistenergy = 0.0d0
      my_liquidcellvol = 0.0d0
      mywb = 0.0d0
      my_l0mid = 0.0d0
      my_l075 = 0.0d0
      my_wspec_num = 0.0d0
      my_wspec_den = 0.0d0
      my_kh_energy = 0.0d0
      my_nut_prof = 0.0d0
      my_num_prof = 0.0d0
      my_nuq_prof = 0.0d0
      my_area_prof = 0.0d0
      my_msum_prof = 0.0d0
      my_m2_prof = 0.0d0
      my_marea_prof = 0.0d0
      my_m2d = 0.0d0
      m2d = 0.0d0
      mbar_prof = 0.0d0
      mvar_prof = 0.0d0
      convfac = dsqrt(Ra*Prs)
      kmid = 1
      k075 = 1
      do kc=2,n3m
        if(dabs(zm(kc)-0.5d0*alx3).lt.dabs(zm(kmid)-0.5d0*alx3)) kmid=kc
        if(dabs(zm(kc)-0.75d0*alx3).lt.dabs(zm(k075)-0.75d0*alx3)) k075=kc
      enddo
      target_wz = 0.75d0*alx3
      kwlo = 1
      do kc=1,n3m-1
        if(zm(kc).le.target_wz .and. zm(kc+1).ge.target_wz)then
          kwlo = kc
          exit
        endif
      enddo
      kwhi = min(kwlo+1,n3m)
      wz_weight = 0.0d0
      if(zm(kwhi).gt.zm(kwlo))then
        wz_weight = (target_wz-zm(kwlo))/(zm(kwhi)-zm(kwlo))
      endif
      do kc=kstart,kend
        kp = kc + 1
        do jc=1,n2m
          do ic=1,n1m
            cellarea = (xc(ic+1)-xc(ic))*(yc(jc+1)-yc(jc))
            cellvol  = (xc(ic+1)-xc(ic))*(yc(jc+1)-yc(jc))*(zc(kc+1)-zc(kc))
             spdval          = dsqrt(Ra/Prs)*dsqrt(q1(ic,jc,kc)**2.0+q2(ic,jc,kc)**2.0+q3(ic,jc,kc)**2.0)
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
            my_spdavg       = my_spdavg       + spdval*cellvol
            my_vrms       = my_vrms       + vrm*cellvol
             my_tavg       = my_tavg       + tcell*cellvol
             my_nuv  = my_nuv + wt*cellvol
             my_nuh  = my_nuh + heatflux*cellvol
             my_trms  = my_trms  + tcell**2*cellvol
              my_moistenergy = my_moistenergy + moistenergy*cellvol
              mywb = mywb + wb*cellvol
              my_msum_prof(kc) = my_msum_prof(kc) + mcell*cellarea
              my_m2_prof(kc) = my_m2_prof(kc) + mcell*mcell*cellarea
              my_marea_prof(kc) = my_marea_prof(kc) + cellarea
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

!     z/H=0.75 vertical-velocity spectral length requested for the
!     aggregation/merger study.  Interpolate w to z/H=0.75 first, then
!     identify the field spectrum before forming the two spectral moments:
!
!       ell_w = sum(|w_hat|^2) / sum(k_h |w_hat|^2)
!       lambda_w = 2*pi*ell_w
!
!     The Hermitian multiplicity restores the omitted negative-ky modes of
!     the real-to-complex transform.  The zero horizontal mode is excluded.
      if(kwlo.ge.kstart .and. kwlo.le.kend)then
        n1mh = n1m/2 + 1
        n2mh = n2m/2 + 1
        coefnorm = 1.d0/(dble(n1m)*dble(n2m))

        do jc=1,n2m
          do ic=1,n1m
            xr_w(jc,ic) = (1.d0-wz_weight)*q3(ic,jc,kwlo) &
     &                   + wz_weight*q3(ic,jc,kwhi)
          enddo
        enddo

        call dfftw_execute_dft_r2c(fwd_plan,xr_w,xa_w)

        my_wspec_num = 0.d0
        my_wspec_den = 0.d0
        do jc=1,n2mh
          ky = 2.d0*pi*dble(jc-1)/rext2
          spec_weight = 2.d0
          if(jc.eq.1) spec_weight = 1.d0
          if(mod(n2m,2).eq.0 .and. jc.eq.n2mh) spec_weight = 1.d0
          do ic=1,n1m
            if(ic.le.n1mh)then
              kx = 2.d0*pi*dble(ic-1)/rext1
            else
              kx = -2.d0*pi*dble(n1m-ic+1)/rext1
            endif
            khat = dsqrt(kx**2 + ky**2)
            if(khat.gt.0.d0)then
              spec_e = spec_weight*coefnorm**2*cdabs(xa_w(jc,ic))**2
              my_wspec_num = my_wspec_num + spec_e
              my_wspec_den = my_wspec_den + khat*spec_e
            endif
          enddo
        enddo
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
      call MPI_ALLREDUCE(my_spdavg,spdavg,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_tavg,tavg,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_vrms,vrms,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_nuv,nuv,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_nuh,nuh,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_trms,trms,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_vrm_horizontal,vrm_horizontal,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_moistenergy,moistenergy,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(mywb,wb,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_l0mid,l0mid,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_l075,l075,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_wspec_num,wspec_num,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_wspec_den,wspec_den,1,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_kh_energy,kh_energy,nkh,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_nut_prof,nut_prof,n3r,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_num_prof,num_prof,n3r,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_nuq_prof,nuq_prof,n3r,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_area_prof,area_prof,n3r,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_msum_prof,msum_prof,n3m,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_m2_prof,m2_prof,n3m,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_marea_prof,marea_prof,n3m,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      spdavg = spdavg/liquidcellvol
      tavg = tavg/liquidcellvol
      nuv  = 1+nuv/liquidcellvol
      nuh  = nuh/liquidcellvol
      trms  = trms/liquidcellvol
      vrms  = vrms/liquidcellvol
      trms = dsqrt(dmax1(0.0d0, trms - tavg**2))
      Re = dsqrt(Ra/Prs)*dsqrt(vrms)
      moistenergy = moistenergy/liquidcellvol
      wb = wb/liquidcellvol
      w_length_z075 = 0.0d0
      if(wspec_den.gt.0.0d0) w_length_z075 = wspec_num/wspec_den
      w_wavelength_z075 = 2.0d0*pi*w_length_z075

!     Moist-convective self-aggregation order parameter.
!     First remove the horizontal mean independently at every z level,
!     then average the horizontal variance in z.  This is deliberately not
!     the variance about one three-dimensional volume mean.
      mse_var_horizontal = 0.0d0
      total_zweight = 0.0d0
      do kc=1,n3m
        if(marea_prof(kc).gt.0.0d0)then
          mbar_prof(kc) = msum_prof(kc)/marea_prof(kc)
          mvar_prof(kc) = dmax1(0.0d0, &
     &      m2_prof(kc)/marea_prof(kc)-mbar_prof(kc)**2)
          zweight = zc(kc+1)-zc(kc)
          mse_var_horizontal = mse_var_horizontal+mvar_prof(kc)*zweight
          total_zweight = total_zweight+zweight
        endif
      enddo
      if(total_zweight.gt.0.0d0)then
        mse_var_horizontal = mse_var_horizontal/total_zweight
      endif
      mse_rms_horizontal = dsqrt(dmax1(0.0d0,mse_var_horizontal))
      mprime2_global = mse_var_horizontal
      mprime2_z075 = 0.0d0
      if(marea_prof(kwlo).gt.0.0d0 .and. &
     &   marea_prof(kwhi).gt.0.0d0)then
        mprime2_z075 = (1.0d0-wz_weight)*mvar_prof(kwlo) &
     &    +wz_weight*mvar_prof(kwhi)
      endif

!     Build the vertically weighted horizontal MSE-anomaly field used by
!     the existing offline self-aggregation analysis:
!
!       m2d(x,y) = integral [m(x,y,z)-<m>_xy(z)] dz / integral dz.
!
!     The operation is performed at every statistics time so the resulting
!     length diagnostics have the same continuous cadence as avgvar.out.
      my_m2d = 0.0d0
      do kc=kstart,kend
        zweight = zc(kc+1)-zc(kc)
        do jc=1,n2m
          do ic=1,n1m
            mcell = dsal(ic,jc,kc)+gamma*qvap(ic,jc,kc)
            my_m2d(jc,ic) = my_m2d(jc,ic) &
     &        +(mcell-mbar_prof(kc))*zweight
          enddo
        enddo
      enddo
      call MPI_ALLREDUCE(my_m2d,m2d,n1m*n2m,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)

      mse_l_peak = 0.0d0
      mse_l_integral = 0.0d0
      mse_l_moment = 0.0d0
      mse_lambda_moment = 0.0d0
      cluster_count0 = 0
      cluster_count_sigma = 0
      cluster_mean0 = 0.0d0
      cluster_max0 = 0.0d0
      cluster_mean_sigma = 0.0d0
      cluster_max_sigma = 0.0d0

      if(myid.eq.0 .and. total_zweight.gt.0.0d0)then
        m2d = m2d/total_zweight
        m2d_mean = sum(m2d)/(dble(n1m)*dble(n2m))
        m2d = m2d-m2d_mean
        m2d_sigma = dsqrt(sum(m2d*m2d)/(dble(n1m)*dble(n2m)))

        n1mh = n1m/2 + 1
        n2mh = n2m/2 + 1
        coefnorm = 1.d0/(dble(n1m)*dble(n2m))
        dx_m = rext1/dble(n1m)
        dy_m = rext2/dble(n2m)
        dk_shell = 2.d0*pi/dmax1(rext1,rext2)
        xr_m = m2d
        call dfftw_execute_dft_r2c(fwd_plan,xr_m,xa_m)

        mse_shell = 0.0d0
        mse_spec_num = 0.0d0
        mse_spec_den = 0.0d0
        do jc=1,n2mh
          ky = 2.d0*pi*dble(jc-1)/rext2
          spec_weight = 2.d0
          if(jc.eq.1) spec_weight = 1.d0
          if(mod(n2m,2).eq.0 .and. jc.eq.n2mh) spec_weight = 1.d0
          do ic=1,n1m
            if(ic.le.n1mh)then
              kx = 2.d0*pi*dble(ic-1)/rext1
            else
              kx = -2.d0*pi*dble(n1m-ic+1)/rext1
            endif
            khat = dsqrt(kx**2+ky**2)
            spec_e = spec_weight*coefnorm**2*cdabs(xa_m(jc,ic))**2
            shell_index = nint(khat/dk_shell)+1
            if(shell_index.ge.1 .and. shell_index.le.n1m+n2m)then
              mse_shell(shell_index) = mse_shell(shell_index)+spec_e
            endif
            if(khat.gt.0.0d0)then
              mse_spec_num = mse_spec_num+spec_e
              mse_spec_den = mse_spec_den+khat*spec_e
            endif
          enddo
        enddo

        peak_shell = 0
        peak_energy = -1.0d0
        do shell_index=2,n1m+n2m
          if(mse_shell(shell_index).gt.peak_energy)then
            peak_energy = mse_shell(shell_index)
            peak_shell = shell_index
          endif
        enddo
        if(peak_shell.gt.1 .and. peak_energy.gt.0.0d0)then
          mse_l_peak = 2.0d0*pi/(dble(peak_shell-1)*dk_shell)
        endif
        if(mse_spec_den.gt.0.0d0)then
          mse_l_moment = mse_spec_num/mse_spec_den
          mse_lambda_moment = 2.0d0*pi*mse_l_moment
        endif

!       Correlation integral scale: radially average the periodic 2-D
!       autocorrelation and integrate its positive lobe to the first zero.
        do jc=1,n2mh
          do ic=1,n1m
            xa_corr(jc,ic) = cmplx(cdabs(xa_m(jc,ic))**2,0.0d0)
          enddo
        enddo
        call dfftw_execute_dft_c2r(bck_plan,xa_corr,xr_corr)
        dr_corr = dmin1(dx_m,dy_m)
        max_corr_bin = int(0.5d0*dmin1(rext1,rext2)/dr_corr)
        corr_sum = 0.0d0
        corr_count = 0.0d0
        corr_profile = 0.0d0
        do jc=1,n2m
          ry = dble(min(jc-1,n2m-jc+1))*dy_m
          do ic=1,n1m
            rx = dble(min(ic-1,n1m-ic+1))*dx_m
            radius = dsqrt(rx**2+ry**2)
            bin_index = nint(radius/dr_corr)+1
            if(bin_index.le.max_corr_bin+1)then
              corr_sum(bin_index) = corr_sum(bin_index)+xr_corr(jc,ic)
              corr_count(bin_index) = corr_count(bin_index)+1.0d0
            endif
          enddo
        enddo
        if(corr_count(1).gt.0.0d0)then
          corr_center = corr_sum(1)/corr_count(1)
        else
          corr_center = 0.0d0
        endif
        if(corr_center.gt.0.0d0)then
          do bin_index=1,max_corr_bin+1
            if(corr_count(bin_index).gt.0.0d0)then
              corr_profile(bin_index) = corr_sum(bin_index) &
     &          /corr_count(bin_index)/corr_center
            endif
          enddo
          mse_l_integral = 0.0d0
          do bin_index=2,max_corr_bin+1
            corr_prev = corr_profile(bin_index-1)
            corr_now = corr_profile(bin_index)
            r_prev = dble(bin_index-2)*dr_corr
            if(corr_now.le.0.0d0)then
              if(corr_prev.gt.0.0d0)then
                r_zero = r_prev+corr_prev/(corr_prev-corr_now)*dr_corr
                mse_l_integral = mse_l_integral &
     &            +0.5d0*corr_prev*(r_zero-r_prev)
              endif
              exit
            else
              mse_l_integral = mse_l_integral &
     &          +0.5d0*(corr_prev+corr_now)*dr_corr
            endif
          enddo
        endif

        call periodic_cluster_stats(m2d,n1m,n2m,dx_m,dy_m,0.0d0, &
     &       cluster_count0,cluster_mean0,cluster_max0)
        call periodic_cluster_stats(m2d,n1m,n2m,dx_m,dy_m,m2d_sigma, &
     &       cluster_count_sigma,cluster_mean_sigma,cluster_max_sigma)
      endif

      if(myid.eq.0)then
        write(96,510) time, trms, Re, vrms, vrm_horizontal,moistenergy,wb,l0mid,l075,nuh
        write(100,512) time, kh_energy(1:nkh)
        write(101,513) time,mse_var_horizontal,mse_rms_horizontal
        write(103,515) time,w_length_z075,w_wavelength_z075
        write(105,513) time,mprime2_global,mprime2_z075
        write(104,516) time,mse_l_peak,mse_l_integral,mse_l_moment, &
     &    mse_lambda_moment,dble(cluster_count0),cluster_mean0, &
     &    cluster_max0,dble(cluster_count_sigma),cluster_mean_sigma, &
     &    cluster_max_sigma
        do kc=1,n3m
          write(102,514) time,zm(kc),mbar_prof(kc),mvar_prof(kc)
        enddo
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
 513   format(1x,f10.4,2(1x,ES20.8))
 514   format(1x,f10.4,3(1x,ES20.8))
 515   format(1x,f10.4,2(1x,ES20.8))
 516   format(1x,f10.4,10(1x,ES20.8))

      return   
      end


!=======================================================================
! Periodic 8-neighbour connected components for thresholded m'_2D.
! Returns the number of clusters and their mean/maximum equivalent radius.
      subroutine periodic_cluster_stats(field,nx,ny,dx,dy,threshold, &
     &                                   ncluster,mean_radius,max_radius)
      implicit none
      integer, intent(in) :: nx,ny
      real, intent(in) :: field(ny,nx),dx,dy,threshold
      integer, intent(out) :: ncluster
      real, intent(out) :: mean_radius,max_radius
      logical :: mask(ny,nx),visited(ny,nx)
      integer :: queue_x(nx*ny),queue_y(nx*ny)
      integer :: ix,iy,jx,jy,nx_neighbor,ny_neighbor
      integer :: di,dj,head,tail,count
      real :: radius,radius_sum,pi_local

      pi_local = 4.0d0*datan(1.0d0)
      mask = field.gt.threshold
      visited = .false.
      ncluster = 0
      radius_sum = 0.0d0
      max_radius = 0.0d0

      do iy=1,ny
        do ix=1,nx
          if(mask(iy,ix) .and. .not.visited(iy,ix))then
            ncluster = ncluster+1
            head = 1
            tail = 1
            queue_x(1) = ix
            queue_y(1) = iy
            visited(iy,ix) = .true.
            count = 0
            do while(head.le.tail)
              jx = queue_x(head)
              jy = queue_y(head)
              head = head+1
              count = count+1
              do dj=-1,1
                do di=-1,1
                  if(di.ne.0 .or. dj.ne.0)then
                    nx_neighbor = modulo(jx-1+di,nx)+1
                    ny_neighbor = modulo(jy-1+dj,ny)+1
                    if(mask(ny_neighbor,nx_neighbor) .and. &
     &                 .not.visited(ny_neighbor,nx_neighbor))then
                      tail = tail+1
                      queue_x(tail) = nx_neighbor
                      queue_y(tail) = ny_neighbor
                      visited(ny_neighbor,nx_neighbor) = .true.
                    endif
                  endif
                enddo
              enddo
            enddo
            radius = dsqrt(dble(count)*dx*dy/pi_local)
            radius_sum = radius_sum+radius
            max_radius = dmax1(max_radius,radius)
          endif
        enddo
      enddo

      mean_radius = 0.0d0
      if(ncluster.gt.0) mean_radius = radius_sum/dble(ncluster)
      return
      end
