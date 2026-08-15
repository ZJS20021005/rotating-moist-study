subroutine tschem
      use param
      use local_arrays
      use mgrd_arrays
      use mpih
      use mpi_param, only: kstart,kend,kstartr,kendr
      implicit none
      integer :: ns,inp,ntr
      integer :: j,k,i,mstep

      dts = dt/dble(lmax)
!====== Begin RK scheme ======!
      do ns=1,nsst
        al=alm(ns)
        ga=gam(ns)
        ro=rom(ns)

        !-- Set var to zero
        dph(:,:,:)=0.0d0
        dq(:,:,:)=0.0d0
        qcap(:,:,:)=0.0d0

        !-- cal properties
        call caldiffusivity
        call update_both_ghosts(n1mr,n2mr,kpa,kstartr,kendr)

        !-- Explicit terms
        call compute_qsat_from_b  ! 计算qs和cond_term
        call hdnlq1
        call hdnlq2
        call hdnlq3
        call hdnlsa
        call hdnlqvap

        !-- Implicit terms
        call invtrq1
        call invtrq2
        call invtrq3
        call update_upper_ghost(n1,n2,q3)
        call invtrsa
        call update_both_ghosts(n1r,n2r,dsal,kstartr,kendr)
        call invtrqvap
        call update_both_ghosts(n1r,n2r,qvap,kstartr,kendr)
        ! Keep diagnostic temperature/moisture fields synchronized with
        ! the newly updated scalars before any statistics or output.
        call compute_qsat_from_b
        call update_both_ghosts(n1r,n2r,T,kstartr,kendr)
        call update_both_ghosts(n1r,n2r,qs,kstartr,kendr)
        call update_both_ghosts(n1r,n2r,cond_term,kstartr,kendr)

        !-- velbc
        call velbc

        !-- Swap before pressure solver !-- Important
        call update_both_ghosts(n1,n2,q1,kstart,kend)
        call update_both_ghosts(n1,n2,q2,kstart,kend)
        call update_both_ghosts(n1,n2,q3,kstart,kend)
        call update_both_ghosts(n1r,n2r,dsal,kstartr,kendr)
        call update_both_ghosts(n1r,n2r,qvap,kstartr,kendr)
        call update_both_ghosts(n1r,n2r,dens,kstartr,kendr)

        !-- Velocity pressure correction
        dph(:,:,:)=0.0d0
        call update_upper_ghost(n1,n2,q3)
        call divg
        call phcalc
        call update_both_ghosts(n1m,n2m+2,dph,kstart,kend)
        call updvp                 ! SOLENOIDAL VEL FIELD
        call update_both_ghosts(n1,n2,q1,kstart,kend)
        call update_both_ghosts(n1,n2,q2,kstart,kend)
        call update_both_ghosts(n1,n2,q3,kstart,kend)
        call prcalc                         !! PRESSURE FIELD
        call update_both_ghosts(n1,n2,pr,kstart,kend)

        !==== Multiresol info exchange ====!
        call mgrd_velitp
        call mgrd_dsalc
        call mgrd_qvapc
        call update_both_ghosts(n1,n2,dsalc,kstart,kend)
        call update_both_ghosts(n1,n2,qvapc,kstart,kend)
      enddo
!====== End RK scheme ======!

      return                                                            
      end
