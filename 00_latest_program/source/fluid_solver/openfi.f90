      subroutine openfi
      use param
      use mpih
      implicit none

      IF(myid.eq.0)then

!       Create the diagnostic tree before opening time-series files.
        call execute_command_line('mkdir -p diagnostics/scale diagnostics/vortex diagnostics/thermo')

!       Unit 95: time, max|q1|, max|q2|, max|q3|.
        open(95, file='./data/vmax.out',status='unknown',access='sequential',position='append')

!       Unit 96: time, trms, Re, <u^2+v^2+w^2>, <u^2+v^2>,
!                <m>, <w*b>, L0(z~0.5), L0(z~0.75), volume NuT.
        open(96, file='./data/avgvar.out',status='unknown',access='sequential',position='append')

!       Unit 97: wall Nusselt diagnostics from nusse_walls.
        open(97, file='./data/nusse_walls.out',status='unknown',access='sequential',position='append')

!       Unit 98: layer diagnostics from CalcStats.
        open(98, file='./data/layer.out',status='unknown',access='sequential',position='append')

!       Unit 99: time, z_face, NuT(z), Num(z), Nuq(z).
        open(99, file='./data/nu_profiles.out',status='unknown',access='sequential',position='append')

!       Unit 100: time, horizontal kinetic-energy shells kh=1..12 at z~0.5.
        open(100, file='./data/kh_energy.out',status='unknown',access='sequential',position='append')

!       Unit 105: MSE positive-structure tracking, m2d>0.
!                 time, id, radius_positive, area_positive, status.
        open(105, file='./diagnostics/vortex/vortex_positive_history.dat', &
     &       status='unknown',access='sequential',position='append')

!       Unit 106: MSE core-structure tracking, m2d>sigma_m2d.
!                 time, id, radius_core, area_core, status.
        open(106, file='./diagnostics/vortex/vortex_core_history.dat', &
     &       status='unknown',access='sequential',position='append')

!       Unit 107: time, Nv_positive, Rmax_positive, Nv_core, Rmax_core.
        open(107, file='./diagnostics/vortex/maximum_radius.dat', &
     &       status='unknown',access='sequential',position='append')

!       Unit 108: time, Rmax_positive/Rmax_core, total_area_positive/total_area_core.
        open(108, file='./diagnostics/vortex/structure_ratio.dat', &
     &       status='unknown',access='sequential',position='append')

!       Unit 109: time, lc_w(z~0.5), vertically averaged lc_w(z).
        open(109, file='./diagnostics/scale/convective_scale.dat', &
     &       status='unknown',access='sequential',position='append')

!       Unit 110: time, Lm, Lq, Lzeta.
        open(110, file='./diagnostics/scale/moist_integral_scale.dat', &
     &       status='unknown',access='sequential',position='append')

!       Unit 111: time, Lpeak_w, Lpeak_m, Lpeak_q.
        open(111, file='./diagnostics/scale/peak_scale.dat', &
     &       status='unknown',access='sequential',position='append')

!       Unit 112: time, N_condensation_cells, condensation_fraction.
        open(112, file='./diagnostics/thermo/condensation.dat', &
     &       status='unknown',access='sequential',position='append')

!       Unit 113: time, m2_z025, m2_z050, m2_z075, m2_volume.
        open(113, file='./diagnostics/thermo/mprime_square.dat', &
     &       status='unknown',access='sequential',position='append')

!       Unit 114: time, threshold, id, mean(u,v,w), mean/rms(zeta),
!                 mean(mprime) inside each current structure.
        open(114, file='./diagnostics/vortex/vortex_kinematics.dat', &
     &       status='unknown',access='sequential',position='append')

!       Unit 115: time, horizontal pressure-gradient RMS Fp.
        open(115, file='./data/pressure_force.out',status='unknown', &
     &       access='sequential',position='append')

!       Unit 116: strict online force-balance bulk time series.
        open(116, file='./data/force_balance.out',status='unknown', &
     &       access='sequential',position='append')

!       Unit 117: strict force-balance profiles at every z and time.
        open(117, file='./data/force_balance_z.out',status='unknown', &
     &       access='sequential',position='append')

      ! reset the time history
      if(ireset.eq.1 .or. nread.eq.0)then
          rewind(95)
          rewind(96)
          rewind(97)
          rewind(98)
          rewind(99)
          rewind(100)
          rewind(105)
          rewind(106)
          rewind(107)
          rewind(108)
          rewind(109)
          rewind(110)
          rewind(111)
          rewind(112)
          rewind(113)
          rewind(114)
          rewind(115)
          rewind(116)
          rewind(117)
      endif

      ENDIF

      return
      end   
      
!==============================================

      subroutine closefi
      use mpih
      implicit none
      
      if(myid.eq.0)then

      close(95)
      close(96)
      close(97)
      close(98)
      close(99)
      close(100)
      close(105)
      close(106)
      close(107)
      close(108)
      close(109)
      close(110)
      close(111)
      close(112)
      close(113)
      close(114)
      close(115)
      close(116)
      close(117)
      endif
     
      return      
      end
