      subroutine aggregation_scale_diagnostics
!     Complete time-dependent scale and thermodynamic diagnostics.
!
!     convective_scale.dat:
!       time, lc_w(z~0.5), vertically averaged lc_w(z)
!
!     moist_integral_scale.dat:
!       time, Lm, Lq, Lzeta
!       These use stretched-grid-weighted vertical means of the
!       per-height horizontal anomalies.
!
!     peak_scale.dat:
!       time, Lpeak_w, Lpeak_m, Lpeak_q
!       The w spectrum is the stretched-grid-weighted sum of the
!       per-height spectra. The m and q spectra use their 2-D fields.
!
!     condensation.dat:
!       time, number of cells with cond_term>0, cell fraction
!
!     mprime_square.dat:
!       time, <m'^2>xy at z~0.25, z~0.50, z~0.75,
!             and stretched-grid volume mean <m'^2>
      use param
      use local_arrays, only: q1,q2,q3,dsal,qvap,cond_term
      use mpi_param, only: kstart,kend,kstartr,kendr
      use mpih
      implicit none

      integer :: ic,jc,kc,ip,im,jp,jm,kp
      integer :: kmid,k025,k050,k075,nshell
      integer :: my_ncond,ncond,total_cells
      real :: zweight,my_zweight,total_zweight
      real :: mbar,qbar,wbar,mcell,qcell,wcell
      real :: vip,vim,ujp,ujm,zeta,zbar
      real :: pnum,pden,plane_length
      real :: my_lcw_sum,lcw_sum,lcw_mid_num,lcw_mid_den
      real :: my_lcw_mid_num,my_lcw_mid_den,lcw_mid,lcw_volume
      real :: lm,lq,lzeta,lpeak_w,lpeak_m,lpeak_q
      real :: m2_plane,my_m2_volume,m2_volume
      real :: my_m2_target(3),m2_target(3)
      real :: cond_fraction
      real :: field(n2m,n1m)
      real :: my_m2d(n2m,n1m),m2d(n2m,n1m)
      real :: my_q2d(n2m,n1m),q2d(n2m,n1m)
      real :: my_zeta2d(n2m,n1m),zeta2d(n2m,n1m)
      real :: my_w_shell(n1m+n2m),w_shell(n1m+n2m)
      real :: dummy_shell(n1m+n2m)
      real :: peak_wavelength
      integer :: nearest_base_level,nearest_refined_level

      call update_both_ghosts(n1,n2,q1,kstart,kend)
      call update_both_ghosts(n1,n2,q2,kstart,kend)
      call update_both_ghosts(n1,n2,q3,kstart,kend)
      call update_both_ghosts(n1r,n2r,dsal,kstartr,kendr)
      call update_both_ghosts(n1r,n2r,qvap,kstartr,kendr)
      call update_both_ghosts(n1r,n2r,cond_term,kstartr,kendr)

      nshell = n1m+n2m
      my_zweight = 0.0d0
      my_lcw_sum = 0.0d0
      my_lcw_mid_num = 0.0d0
      my_lcw_mid_den = 0.0d0
      my_m2_volume = 0.0d0
      my_m2_target = 0.0d0
      my_m2d = 0.0d0
      my_q2d = 0.0d0
      my_zeta2d = 0.0d0
      my_w_shell = 0.0d0

      kmid = nearest_base_level(0.5d0*alx3)
      k025 = nearest_refined_level(0.25d0*alx3)
      k050 = nearest_refined_level(0.50d0*alx3)
      k075 = nearest_refined_level(0.75d0*alx3)

!     Base-grid spectra and the vertically averaged anomaly fields.
      do kc=kstart,kend
        kp = kc+1
        zweight = zc(kc+1)-zc(kc)
        mbar = 0.0d0
        qbar = 0.0d0
        wbar = 0.0d0
        zbar = 0.0d0

        do jc=1,n2m
          jp = jpv(jc)
          jm = jmv(jc)
          do ic=1,n1m
            ip = ipv(ic)
            im = imv(ic)
            mbar = mbar+dsal(ic,jc,kc)+gamma*qvap(ic,jc,kc)
            qbar = qbar+qvap(ic,jc,kc)
            wbar = wbar+0.5d0*(q3(ic,jc,kc)+q3(ic,jc,kp))
            vip = 0.5d0*(q2(ip,jc,kc)+q2(ip,jp,kc))
            vim = 0.5d0*(q2(im,jc,kc)+q2(im,jp,kc))
            ujp = 0.5d0*(q1(ic,jp,kc)+q1(ip,jp,kc))
            ujm = 0.5d0*(q1(ic,jm,kc)+q1(ip,jm,kc))
            zbar = zbar+0.5d0*dx1*(vip-vim) &
     &                  -0.5d0*dx2*(ujp-ujm)
          enddo
        enddo
        mbar = mbar/(dble(n1m)*dble(n2m))
        qbar = qbar/(dble(n1m)*dble(n2m))
        wbar = wbar/(dble(n1m)*dble(n2m))
        zbar = zbar/(dble(n1m)*dble(n2m))

        do jc=1,n2m
          jp = jpv(jc)
          jm = jmv(jc)
          do ic=1,n1m
            ip = ipv(ic)
            im = imv(ic)
            wcell = 0.5d0*(q3(ic,jc,kc)+q3(ic,jc,kp))
            mcell = dsal(ic,jc,kc)+gamma*qvap(ic,jc,kc)
            qcell = qvap(ic,jc,kc)
            vip = 0.5d0*(q2(ip,jc,kc)+q2(ip,jp,kc))
            vim = 0.5d0*(q2(im,jc,kc)+q2(im,jp,kc))
            ujp = 0.5d0*(q1(ic,jp,kc)+q1(ip,jp,kc))
            ujm = 0.5d0*(q1(ic,jm,kc)+q1(ip,jm,kc))
            zeta = 0.5d0*dx1*(vip-vim) &
     &             -0.5d0*dx2*(ujp-ujm)
            field(jc,ic) = wcell-wbar
            my_m2d(jc,ic) = my_m2d(jc,ic)+(mcell-mbar)*zweight
            my_q2d(jc,ic) = my_q2d(jc,ic)+(qcell-qbar)*zweight
            my_zeta2d(jc,ic) = my_zeta2d(jc,ic)+(zeta-zbar)*zweight
          enddo
        enddo

        dummy_shell = 0.0d0
        call plane_spectral_moments(field,pnum,pden,dummy_shell, &
     &                              nshell,zweight)
        if(pden.gt.0.0d0)then
          plane_length = pnum/pden
          my_lcw_sum = my_lcw_sum+plane_length*zweight
        endif
        my_w_shell = my_w_shell+dummy_shell

        if(kc.eq.kmid)then
          my_lcw_mid_num = pnum
          my_lcw_mid_den = pden
        endif
        my_zweight = my_zweight+zweight
      enddo

      call MPI_ALLREDUCE(my_zweight,total_zweight,1,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_lcw_sum,lcw_sum,1,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_lcw_mid_num,lcw_mid_num,1,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_lcw_mid_den,lcw_mid_den,1,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_m2d,m2d,n1m*n2m,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_q2d,q2d,n1m*n2m,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_zeta2d,zeta2d,n1m*n2m,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_w_shell,w_shell,nshell,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)

      lcw_mid = 0.0d0
      lcw_volume = 0.0d0
      if(lcw_mid_den.gt.0.0d0) lcw_mid = lcw_mid_num/lcw_mid_den
      if(total_zweight.gt.0.0d0)then
        lcw_volume = lcw_sum/total_zweight
        m2d = m2d/total_zweight
        q2d = q2d/total_zweight
        zeta2d = zeta2d/total_zweight
      endif
      m2d = m2d-sum(m2d)/(dble(n1m)*dble(n2m))
      q2d = q2d-sum(q2d)/(dble(n1m)*dble(n2m))
      zeta2d = zeta2d-sum(zeta2d)/(dble(n1m)*dble(n2m))

      dummy_shell = 0.0d0
      call plane_spectral_moments(m2d,pnum,pden,dummy_shell, &
     &                            nshell,1.0d0)
      lm = 0.0d0
      if(pden.gt.0.0d0) lm = pnum/pden
      lpeak_m = peak_wavelength(dummy_shell,nshell)

      dummy_shell = 0.0d0
      call plane_spectral_moments(q2d,pnum,pden,dummy_shell, &
     &                            nshell,1.0d0)
      lq = 0.0d0
      if(pden.gt.0.0d0) lq = pnum/pden
      lpeak_q = peak_wavelength(dummy_shell,nshell)

      dummy_shell = 0.0d0
      call plane_spectral_moments(zeta2d,pnum,pden,dummy_shell, &
     &                            nshell,1.0d0)
      lzeta = 0.0d0
      if(pden.gt.0.0d0) lzeta = pnum/pden
      lpeak_w = peak_wavelength(w_shell,nshell)

!     Refined-grid condensation count and mprime-squared statistics.
      my_ncond = 0
      do kc=kstartr,kendr
        zweight = zcr(kc+1)-zcr(kc)
        mbar = 0.0d0
        do jc=1,n2mr
          do ic=1,n1mr
            mbar = mbar+dsal(ic,jc,kc)+gamma*qvap(ic,jc,kc)
            if(cond_term(ic,jc,kc).gt.0.0d0) my_ncond=my_ncond+1
          enddo
        enddo
        mbar = mbar/(dble(n1mr)*dble(n2mr))

        m2_plane = 0.0d0
        do jc=1,n2mr
          do ic=1,n1mr
            mcell = dsal(ic,jc,kc)+gamma*qvap(ic,jc,kc)
            m2_plane = m2_plane+(mcell-mbar)**2
          enddo
        enddo
        m2_plane = m2_plane/(dble(n1mr)*dble(n2mr))
        my_m2_volume = my_m2_volume+m2_plane*zweight
        if(kc.eq.k025) my_m2_target(1) = m2_plane
        if(kc.eq.k050) my_m2_target(2) = m2_plane
        if(kc.eq.k075) my_m2_target(3) = m2_plane
      enddo

      call MPI_ALLREDUCE(my_ncond,ncond,1,MPI_INTEGER,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_m2_volume,m2_volume,1,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_m2_target,m2_target,3,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)

      total_cells = n1mr*n2mr*n3mr
      cond_fraction = dble(ncond)/dble(total_cells)
      if(total_zweight.gt.0.0d0) m2_volume=m2_volume/total_zweight

      if(myid.eq.0)then
        write(109,510) time,lcw_mid,lcw_volume
        write(110,511) time,lm,lq,lzeta
        write(111,511) time,lpeak_w,lpeak_m,lpeak_q
        write(112,512) time,ncond,cond_fraction
        write(113,513) time,m2_target(1),m2_target(2), &
     &                 m2_target(3),m2_volume
      endif

510   format(1x,f10.4,2(1x,ES20.8))
511   format(1x,f10.4,3(1x,ES20.8))
512   format(1x,f10.4,1x,i14,1x,ES20.8)
513   format(1x,f10.4,4(1x,ES20.8))
      return
      end

!=============================================================================
      subroutine plane_spectral_moments(field,numerator,denominator, &
     &                                  shell_power,nshell,weight)
!     numerator=sum(E), denominator=sum(k_h E), excluding k_h=0.
      use param
      implicit none
      integer, intent(in) :: nshell
      real, intent(in) :: field(n2m,n1m),weight
      real, intent(out) :: numerator,denominator
      real, intent(out) :: shell_power(nshell)
      integer :: ic,jc,kxi,ish,n1mh_local,n2mh_local
      real :: xr(n2m,n1m)
      complex :: xa(m2mh,m1m)
      real :: kx,ky,kh,energy,multiplicity,coefnorm,dk

      xr = field
      call dfftw_execute_dft_r2c(fwd_plan,xr,xa)
      numerator = 0.0d0
      denominator = 0.0d0
      shell_power = 0.0d0
      n1mh_local = n1m/2+1
      n2mh_local = n2m/2+1
      coefnorm = 1.0d0/(dble(n1m)*dble(n2m))
      dk = dmin1(2.0d0*pi/rext1,2.0d0*pi/rext2)

      do jc=1,n2mh_local
        ky = 2.0d0*pi*dble(jc-1)/rext2
        multiplicity = 2.0d0
        if(jc.eq.1) multiplicity = 1.0d0
        if(mod(n2m,2).eq.0 .and. jc.eq.n2mh_local) &
     &     multiplicity = 1.0d0
        do ic=1,n1m
          if(ic.le.n1mh_local)then
            kxi = ic-1
          else
            kxi = -(n1m-ic+1)
          endif
          kx = 2.0d0*pi*dble(kxi)/rext1
          kh = dsqrt(kx*kx+ky*ky)
          if(kh.gt.0.0d0)then
            energy = multiplicity*coefnorm**2*cdabs(xa(jc,ic))**2
            numerator = numerator+energy
            denominator = denominator+kh*energy
            ish = nint(kh/dk)
            if(ish.ge.1 .and. ish.le.nshell) &
     &         shell_power(ish)=shell_power(ish)+weight*energy
          endif
        enddo
      enddo
      return
      end

!=============================================================================
      real function peak_wavelength(shell_power,nshell)
      use param
      implicit none
      integer, intent(in) :: nshell
      real, intent(in) :: shell_power(nshell)
      integer :: ish,ipeak
      real :: dk

      peak_wavelength = 0.0d0
      ipeak = 0
      do ish=1,nshell
        if(shell_power(ish).gt.0.0d0)then
          if(ipeak.eq.0)then
            ipeak = ish
          elseif(shell_power(ish).gt.shell_power(ipeak))then
            ipeak = ish
          endif
        endif
      enddo
      if(ipeak.gt.0)then
        dk = dmin1(2.0d0*pi/rext1,2.0d0*pi/rext2)
        peak_wavelength = 2.0d0*pi/(dble(ipeak)*dk)
      endif
      return
      end

!=============================================================================
      integer function nearest_base_level(target_z)
      use param
      implicit none
      real, intent(in) :: target_z
      integer :: kc
      nearest_base_level = 1
      do kc=2,n3m
        if(dabs(zm(kc)-target_z).lt. &
     &     dabs(zm(nearest_base_level)-target_z)) nearest_base_level=kc
      enddo
      return
      end

!=============================================================================
      integer function nearest_refined_level(target_z)
      use param
      implicit none
      real, intent(in) :: target_z
      integer :: kc
      nearest_refined_level = 1
      do kc=2,n3mr
        if(dabs(zmr(kc)-target_z).lt. &
     &     dabs(zmr(nearest_refined_level)-target_z)) &
     &     nearest_refined_level=kc
      enddo
      return
      end
