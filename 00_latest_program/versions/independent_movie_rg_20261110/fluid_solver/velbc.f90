      subroutine velbc
      use param
      use mpih
      use local_arrays, only: q2,q3,q1
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc,ic

      real dl1q, dl2q, dlf1, dlf2

      do kc=kstart,kend
        do jc=1,n2
          q1(n1,jc,kc) = q1(1,jc,kc)
        enddo

        do ic=1,n1
          q2(ic,n2,kc) = q2(ic,1,kc)
        enddo
      enddo

      if(kstart.eq.1) then
        do jc=1,n2
            do ic=1,n1
               q3(ic,jc,1) = 0.0d0
            enddo
        enddo
      endif

      if(kend.eq.n3m) then
        do jc=1,n2
            do ic=1,n1
               q3(ic,jc,n3) = 0.0d0
            enddo
        enddo
      endif

      return
      end

