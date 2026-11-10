!
!     this subroutine calculates divg(q).
!     q are the "fluxes"
!
      subroutine divg
      use param
      use local_arrays, only: q1,q2,q3,dph
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,jp,kc,kp,ic,ip
      real    :: usdtal,dqcap
      real    :: q1up,q1low,q2up,q2low,q3up,q3low

      usdtal = 1.d0/(dt*al)

      do kc=kstart,kend

        kp=kc+1
        do jc=1,n2m
          jp=jpv(jc)
          do ic=1,n1m
            ip=ipv(ic)

               dqcap= (q1(ip,jc,kc)-q1(ic,jc,kc))*dx1&
        &            +(q2(ic,jp,kc)-q2(ic,jc,kc))*dx2&
        &            +(q3(ic,jc,kp)-q3(ic,jc,kc))*udx3m(kc)
            dph(ic,jc,kc)=dqcap*usdtal
          enddo
        enddo
      enddo

      return
      end
