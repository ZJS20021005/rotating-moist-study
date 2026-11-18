subroutine compute_qsat_from_b()
      use param
      use local_arrays, only: T, qs, cond_term, dsal, qvap
     use mpi_param, only: kstartr,kendr
      implicit none
      integer :: ic,jc,kc
      real    :: x,H 
      do kc=kstartr,kendr
        do jc=1,n2mr
          do ic=1,n1mr
            T(ic,jc,kc) = dsal(ic,jc,kc) - betaqs * zmr(kc)
            qs(ic,jc,kc) = exp(alphaqs * T(ic,jc,kc))
            x = qvap(ic,jc,kc)-qs(ic,jc,kc)
            H = 0.5d0*(1.d0 + tanh(1d8*x))
            cond_term(ic,jc,kc) = H*(qvap(ic,jc,kc)-qs(ic,jc,kc))/tau_cond
!            cond_term(ic,jc,kc) = 0.0
          end do
        end do
      end do
      end subroutine
