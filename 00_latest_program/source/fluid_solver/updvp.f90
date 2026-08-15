      subroutine updvp
      use param
      use local_arrays, only: q2,q3,dph,q1
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,jm,kc,km,ic,im
      real    :: usukm,udx2,udx1,locdph
      real    :: q1nend,q2nend,q3nend

      udx1 = al*dt*dx1
      udx2 = al*dt*dx2
      do kc=kstart,kend
        usukm = al*dt*udx3c(kc)
        do jc=1,n2m
          do ic=1,n1m
             im=imv(ic)
             jm=jmv(jc)
             km=kc-1
             q1(ic,jc,kc)=q1(ic,jc,kc)-(dph(ic,jc,kc)-dph(im,jc,kc))*udx1
             q2(ic,jc,kc)=q2(ic,jc,kc)-(dph(ic,jc,kc)-dph(ic,jm,kc))*udx2
             if(kc.eq.1 .or. kc.eq.n3) then
               q3(ic,jc,kc)=q3(ic,jc,kc)
             else
               q3(ic,jc,kc)=q3(ic,jc,kc)-(dph(ic,jc,kc)-dph(ic,jc,km))*usukm
             endif
        enddo 
       enddo
      enddo
      return
      end

