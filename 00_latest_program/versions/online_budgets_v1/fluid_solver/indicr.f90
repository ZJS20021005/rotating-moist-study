      subroutine indicr
      use param
      use mpih
      implicit none
      integer :: jc,kc,ic
   
      !  i-dir
      do ic=1,n1mr
        imvr(ic)=ic-1
        ipvr(ic)=ic+1
        if(ic.eq.1) imvr(ic)=n1mr
        if(ic.eq.n1mr) ipvr(ic)=1
      enddo

      !  j-dir
      do jc=1,n2mr
        jmvr(jc)=jc-1
        jpvr(jc)=jc+1
        if(jc.eq.1) jmvr(jc)=n2mr
        if(jc.eq.n2mr) jpvr(jc)=1
      enddo

      !  k-dir
      do kc=1,n3mr
        kmvr(kc)=kc-1
        kpvr(kc)=kc+1
        if(kc.eq.1) kmvr(kc)=kc
        if(kc.eq.n3mr) kpvr(kc)=kc
      end do

      return
      end
