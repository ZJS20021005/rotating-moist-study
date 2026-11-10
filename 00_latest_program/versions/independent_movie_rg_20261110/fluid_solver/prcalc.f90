      subroutine prcalc
      use param
      use local_arrays, only: pr,dph
      use mpi_param, only: kstart,kend
      implicit none
      integer :: kp,km,jm,jp,jc,kc,ic,ip,im
      real    :: be,amm,acc,app

      !    the pressure is evaluated at the center of the box.
      !     p^{n+1} = p^{n} + phi^{n+1} - b * Nabla^2 phi^{n+1}
      be=al*beta
      do kc=kstart,kend
        kp=kpv(kc)
        if(kc.eq.n3m) kp=kc
        km=kmv(kc)
        if(kc.eq.1) km=kc
        amm=amphk(kc)
        acc=acphk(kc)
        app=apphk(kc)
        do jc=1,n2m
          jm=jmv(jc)
          if(jc.eq.1) jm=jc
          jp=jpv(jc)
          if(jc.eq.n2m) jp=jc
          do ic=1,n1m
            im=imv(ic)
            if(ic.eq.1) im=ic
            ip=ipv(ic)
            if(ic.eq.n1m) ip=ic
            pr(ic,jc,kc)=pr(ic,jc,kc)+dph(ic,jc,kc)-be*(&
     &        (dph(ip,jc,kc)-2.d0*dph(ic,jc,kc)+dph(im,jc,kc))*dx1q&
     &       +(dph(ic,jp,kc)-2.d0*dph(ic,jc,kc)+dph(ic,jm,kc))*dx2q&
     &       +(dph(ic,jc,kp)*app+dph(ic,jc,kc)*acc+dph(ic,jc,km)*amm))
          enddo
        enddo
      enddo
      return
      end
