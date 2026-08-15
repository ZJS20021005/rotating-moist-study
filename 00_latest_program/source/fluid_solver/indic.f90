      subroutine indic
      use param
      use mpih
      implicit none
      integer :: jc,kc,ic,k

      !  k-dir
      do kc=1,n3m
        kmv(kc)=kc-1
        kpv(kc)=kc+1
        if(kc.eq.1) kmv(kc)=kc
        if(kc.eq.n3m) kpv(kc)=kc
      enddo
      do kc=1,n3m
        kpc(kc)=kpv(kc)-kc
        kmc(kc)=kc-kmv(kc)
        kup(kc)=1-kpc(kc)
        kum(kc)=1-kmc(kc)
      enddo

      !  i-dir
      do ic=1,n1m
        imv(ic)=ic-1
        ipv(ic)=ic+1
        if(ic.eq.1) imv(ic)=n1m
        if(ic.eq.n1m) ipv(ic)=1
      enddo

      !  j-dir
      do jc=1,n2m
        jmv(jc)=jc-1
        jpv(jc)=jc+1
        if(jc.eq.1) jmv(jc)=n2m
        if(jc.eq.n2m) jpv(jc)=1
      enddo

      do jc = 1,n2+1
       jmhv(jc) = mod(jc,n2m/2+1)
       if(jmhv(jc).eq.0) jmhv(jc) = n2m/2 + 1
      enddo

      return
      end
