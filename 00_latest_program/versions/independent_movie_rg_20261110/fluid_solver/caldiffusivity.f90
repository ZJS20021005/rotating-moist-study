      subroutine caldiffusivity
      use param
      use local_arrays, only: kpa
      use mpi_param, only: kstartr,kendr
      implicit none
      integer :: jc,kc,ic

      do kc=kstartr,kendr
        do jc=1,n2mr
          do ic=1,n1mr
            kpa(ic,jc,kc) = kps
          enddo
        enddo
      enddo

      return
      end

