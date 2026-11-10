!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: StatRoutines.F90                               !
!    CONTAINS: subroutine CalcStats,WriteStats            !
!                                                         !
!    PURPOSE: Calculates and writes out statistics for    !
!     the flow field. All quantities are averaged in the  !
!     two horizontal (homogeneous) directions             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcStats
      use param
      use local_arrays, only: kpa,q3,q2,q1,dens,dsal
      use mpi_param
      use stat_arrays



      implicit none
      integer :: jc,kc,ic

      nstatsamples = nstatsamples + 1.0

      do kc=kstart,kend
        do jc=1,n2
          do ic=1,n1
            vz_me(ic,jc,kc)=vz_me(ic,jc,kc)+q3(ic,jc,kc)
            vy_me(ic,jc,kc)=vy_me(ic,jc,kc)+q2(ic,jc,kc)
            vx_me(ic,jc,kc)=vx_me(ic,jc,kc)+q1(ic,jc,kc)
          enddo
        enddo
      enddo


      do kc=kstartr,kendr
        do jc=1,n2r
          do ic=1,n1r
            temp_me(ic,jc,kc)=temp_me(ic,jc,kc)+dsal(ic,jc,kc)
          enddo
        enddo
      enddo

      return
      end

