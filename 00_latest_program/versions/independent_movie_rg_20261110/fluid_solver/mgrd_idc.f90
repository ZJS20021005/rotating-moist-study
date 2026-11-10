      subroutine mgrd_idc
      use param
      use mgrd_arrays
      use mpih
      implicit none
      integer :: jc,kc,ic,i,j,k
      integer :: icr, jcr, kcr   

      real xxl(0:m1r), yyl(0:m2r), zzl(0:m3r)
      real xxs(-1:m1+1), yys(-1:m2+1), zzs(-1:m3+1)
      real h00, h01, h10, h11
      real lxm,lxp, lym,lyp, lzm,lzp, lxa,lya,lza
      real dlc, dlm, dlp

      !m============================================================
      ! for second order interpolate velocity
      IF(mref1.eq.1)then
      irangs(0) = 1
      do ic=1,n1m
        irangs(ic) = ic
      enddo
      irangs(n1) = n1

      ELSE

      irangs(0) = 1
      do ic=1,n1m
        do i=irangs(ic-1),n1mr
          if(xmr(i).lt.xm(ic) .and. xmr(i+1).ge.xm(ic))then
            irangs(ic) = i+1
          endif
        enddo
      enddo
      irangs(n1) = n1r

      ENDIF

      IF(mref2.eq.1)then

      jrangs(0) = 1
      do jc=1,n2m
        jrangs(jc) = jc
      enddo
      jrangs(n2) = n2

      ELSE

      jrangs(0) = 1
      do jc=1,n2m
        do j=jrangs(jc-1),n2mr
          if(ymr(j).lt.ym(jc) .and. ymr(j+1).ge.ym(jc))then
            jrangs(jc) = j+1
          endif
        enddo
      enddo
      jrangs(n2) = n2r

      ENDIF

      IF(mref3.eq.1)then

      krangs(0) = 1
      do kc=1,n3m
        krangs(kc) = kc
      enddo
      krangs(n3) = n3

      ELSE

      krangs(0) = 1
      do kc=1,n3m
        do k=krangs(kc-1),n3mr
          if(zmr(k).lt.zm(kc) .and. zmr(k+1).ge.zm(kc))then
            krangs(kc) = k+1
          endif         
        enddo
      enddo
      krangs(n3) = n3r

      ENDIF


      ! interpolate coefficients
      !  for q1

      xxl(1:n1r) = xcr(1:n1r)
      xxl(0) = 2.d0*xxl(1) - xxl(2)
      xxs(1:n1) = xc(1:n1)
      xxs(0) = 2.d0*xxs(1) - xxs(2)
      xxs(-1) = 2.d0*xxs(0) - xxs(1)
      xxs(n1+1) = 2.d0*xxs(n1) - xxs(n1m)

      yyl(1:n2mr) = ymr(1:n2mr)
      yyl(0) = 2.d0*yyl(1) - yyl(2)
      yyl(n2r) = 2.d0*yyl(n2mr) - yyl(n2mr-1)
      yys(1:n2m) = ym(1:n2m)
      yys(0) = 2.d0*yys(1) - yys(2)
      yys(-1) = 2.d0*yys(0) - yys(1)
      yys(n2) = 2.d0*yys(n2m) - yys(n2m-1)
      yys(n2+1) = 2.d0*yys(n2) - yys(n2m)

      zzl(0) = 0.d0
      zzl(1:n3mr) = zmr(1:n3mr)
      zzl(n3r) = alx3
      zzs(0) = 0.d0
      zzs(1:n3m) = zm(1:n3m)
      zzs(n3) = alx3

      do kc=0,n3m
       if(kc.eq.0)then
        dlc = zzs(kc+1)-zzs(kc)
        dlp = zzs(kc+2)-zzs(kc+1)
        lza = 1.d0/dlc
        do kcr=max(krangs(kc),1),min(krangs(kc+1)-1,n3mr)
         lzm = (zzl(kcr)-zzs(kc))*lza
         lzp = 1.d0 - lzm
         h00=(1.d0+2.d0*lzm)*lzp*lzp
         h10=lzm*lzp*lzp
         h01=(1.d0+2.d0*lzp)*lzm*lzm
         h11=-lzp*lzm*lzm
         czq1(1,kcr)=0.d0                                  
         czq1(2,kcr)=h00-h11*dlp/(dlp+dlc)-dble(ubcbot)*h10
         czq1(3,kcr)=h01+h11*(dlp-dlc)/dlp+dble(ubcbot)*h10
         czq1(4,kcr)=h11*dlc*dlc/dlp/(dlp+dlc)
        enddo
       elseif(kc.eq.n3m)then
        dlc = zzs(kc+1)-zzs(kc)
        dlm = zzs(kc)-zzs(kc-1)
        lza = 1.d0/dlc
        do kcr=max(krangs(kc),1),min(krangs(kc+1)-1,n3mr)
         lzm = (zzl(kcr)-zzs(kc))*lza
         lzp = 1.d0 - lzm
         h00=(1.d0+2.d0*lzm)*lzp*lzp
         h10=lzm*lzp*lzp
         h01=(1.d0+2.d0*lzp)*lzm*lzm
         h11=-lzp*lzm*lzm
         czq1(1,kcr)=-h10*dlc*dlc/dlm/(dlc+dlm)
         czq1(2,kcr)=h00+h10*(dlc-dlm)/dlm-dble(ubctop)*h11
         czq1(3,kcr)=h01+h10*dlm/(dlc+dlm)+dble(ubctop)*h11
         czq1(4,kcr)=0.d0
        enddo
       else
        dlc = zzs(kc+1)-zzs(kc)
        dlm = zzs(kc)-zzs(kc-1)
        dlp = zzs(kc+2)-zzs(kc+1)
        lza = 1.d0/dlc
        do kcr=max(krangs(kc),1),min(krangs(kc+1)-1,n3mr)
         lzm = (zzl(kcr) - zzs(kc))*lza
         lzp = 1.d0 - lzm
         h00=(1.d0+2.d0*lzm)*lzp*lzp
         h10=lzm*lzp*lzp
         h01=(1.d0+2.d0*lzp)*lzm*lzm
         h11=-lzp*lzm*lzm
         czq1(1,kcr)=-h10*dlc*dlc/dlm/(dlc+dlm)
         czq1(2,kcr)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
         czq1(3,kcr)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
         czq1(4,kcr)=h11*dlc*dlc/dlp/(dlp+dlc)
        enddo
       endif
      enddo
      do jc=0,n2m
        dlc = yys(jc+1)-yys(jc)
        dlm = yys(jc)-yys(jc-1)
        dlp = yys(jc+2)-yys(jc+1)
        lya=1.d0/dlc
        do jcr=max(jrangs(jc),1),min(jrangs(jc+1)-1,n2mr)
         lym = (yyl(jcr) - yys(jc))*lya
         lyp = 1.d0 - lym
         h00=(1.d0+2.d0*lym)*lyp*lyp
         h10=lym*lyp*lyp
         h01=(1.d0+2.d0*lyp)*lym*lym
         h11=-lyp*lym*lym
         cyq1(1,jcr)=-h10*dlc*dlc/dlm/(dlc+dlm)
         cyq1(2,jcr)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
         cyq1(3,jcr)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
         cyq1(4,jcr)=h11*dlc*dlc/dlp/(dlp+dlc)
        enddo
      enddo
      do ic=1,n1m
        dlc = xxs(ic+1)-xxs(ic)
        dlm = xxs(ic)-xxs(ic-1)
        dlp = xxs(ic+2)-xxs(ic+1)
        lxa=1.d0/dlc
        do icr=max((ic-1)*mref1+1,1),min(ic*mref1,n1mr)
         lxm = (xxl(icr) - xxs(ic))*lxa
         lxp = 1.d0 - lxm
         h00=(1.d0+2.d0*lxm)*lxp*lxp
         h10=lxm*lxp*lxp
         h01=(1.d0+2.d0*lxp)*lxm*lxm
         h11=-lxp*lxm*lxm
         cxq1(1,icr)=-h10*dlc*dlc/dlm/(dlc+dlm)
         cxq1(2,icr)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
         cxq1(3,icr)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
         cxq1(4,icr)=h11*dlc*dlc/dlp/(dlp+dlc)
        enddo
      enddo

      !  for q2

      xxl(1:n1mr) = xmr(1:n1mr)
      xxl(0) = 2.d0*xxl(1) - xxl(2)
      xxl(n1r) = 2.d0*xxl(n1mr) - xxl(n1mr-1)
      xxs(1:n1m) = xm(1:n1m)
      xxs(0) = 2.d0*xxs(1) - xxs(2)
      xxs(-1) = 2.d0*xxs(0) - xxs(1)
      xxs(n1) = 2.d0*xxs(n1m) - xxs(n1m-1)
      xxs(n1+1) = 2.d0*xxs(n1) - xxs(n1m)

      yyl(1:n2r) = ycr(1:n2r)
      yyl(0) = 2.d0*yyl(1) - yyl(2)
      yys(1:n2) = yc(1:n2)
      yys(0) = 2.d0*yys(1) - yys(2)
      yys(-1) = 2.d0*yys(0) - yys(1)
      yys(n2+1) = 2.d0*yys(n2) - yys(n2m)


      zzl(0) = zcr(1)
      zzl(1:n3mr) = zmr(1:n3mr)
      zzl(n3r) = zcr(n3r)
      zzs(0) = zc(1)
      zzs(1:n3m) = zm(1:n3m)
      zzs(n3) = zc(n3)

      do kc=0,n3m
       if(kc.eq.0)then
        dlc = zzs(kc+1)-zzs(kc)
        dlp = zzs(kc+2)-zzs(kc+1)
        lza = 1.d0/dlc
        do kcr=max(krangs(kc),1),min(krangs(kc+1)-1,n3mr)
         lzm = (zzl(kcr) - zzs(kc))*lza
         lzp = 1.d0 - lzm
         h00=(1.d0+2.d0*lzm)*lzp*lzp
         h10=lzm*lzp*lzp
         h01=(1.d0+2.d0*lzp)*lzm*lzm
         h11=-lzp*lzm*lzm
         czq2(1,kcr)=0.d0
         czq2(2,kcr)=h00-h11*dlp/(dlp+dlc)-dble(ubcbot)*h10
         czq2(3,kcr)=h01+h11*(dlp-dlc)/dlp+dble(ubcbot)*h10
         czq2(4,kcr)=h11*dlc*dlc/dlp/(dlp+dlc)
        enddo
       elseif(kc.eq.n3m)then
        dlc = zzs(kc+1)-zzs(kc)
        dlm = zzs(kc)-zzs(kc-1)
        lza = 1.d0/dlc
        do kcr=max(krangs(kc),1),min(krangs(kc+1)-1,n3mr)
         lzm = (zzl(kcr) - zzs(kc))*lza
         lzp = 1.d0 - lzm
         h00=(1.d0+2.d0*lzm)*lzp*lzp
         h10=lzm*lzp*lzp
         h01=(1.d0+2.d0*lzp)*lzm*lzm
         h11=-lzp*lzm*lzm
         czq2(1,kcr)=-h10*dlc*dlc/dlm/(dlc+dlm)
         czq2(2,kcr)=h00+h10*(dlc-dlm)/dlm-dble(ubctop)*h11
         czq2(3,kcr)=h01+h10*dlm/(dlm+dlc)+dble(ubctop)*h11
         czq2(4,kcr)=0.d0
        enddo
       else
        dlc = zzs(kc+1)-zzs(kc)
        dlm = zzs(kc)-zzs(kc-1)
        dlp = zzs(kc+2)-zzs(kc+1)
        lza = 1.d0/dlc
        do kcr=max(krangs(kc),1),min(krangs(kc+1)-1,n3mr)
         lzm = (zzl(kcr) - zzs(kc))*lza
         lzp = 1.d0 - lzm
         h00=(1.d0+2.d0*lzm)*lzp*lzp
         h10=lzm*lzp*lzp
         h01=(1.d0+2.d0*lzp)*lzm*lzm
         h11=-lzp*lzm*lzm
         czq2(1,kcr)=-h10*dlc*dlc/dlm/(dlc+dlm)
         czq2(2,kcr)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
         czq2(3,kcr)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
         czq2(4,kcr)=h11*dlc*dlc/dlp/(dlp+dlc)
        enddo
       endif
      enddo
      do jc=1,n2m
        dlc = yys(jc+1)-yys(jc)
        dlm = yys(jc)-yys(jc-1)
        dlp = yys(jc+2)-yys(jc+1)
        lya=1.d0/dlc
        do jcr=max((jc-1)*mref2+1,1),min(jc*mref2,n2mr)
         lym = (yyl(jcr) - yys(jc))*lya
         lyp = 1.d0 - lym
         h00=(1.d0+2.d0*lym)*lyp*lyp
         h10=lym*lyp*lyp
         h01=(1.d0+2.d0*lyp)*lym*lym
         h11=-lyp*lym*lym
         cyq2(1,jcr)=-h10*dlc*dlc/dlm/(dlc+dlm)
         cyq2(2,jcr)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
         cyq2(3,jcr)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
         cyq2(4,jcr)=h11*dlc*dlc/dlp/(dlp+dlc)
        enddo
      enddo
      do ic=0,n1m
        dlc = xxs(ic+1)-xxs(ic)
        dlm = xxs(ic)-xxs(ic-1)
        dlp = xxs(ic+2)-xxs(ic+1)
        lxa=1.d0/dlc
        do icr=max(irangs(ic),1),min(irangs(ic+1)-1,n1mr)
         lxm = (xxl(icr) - xxs(ic))*lxa
         lxp = 1.d0 - lxm
         h00=(1.d0+2.d0*lxm)*lxp*lxp
         h10=lxm*lxp*lxp
         h01=(1.d0+2.d0*lxp)*lxm*lxm
         h11=-lxp*lxm*lxm
         cxq2(1,icr)=-h10*dlc*dlc/dlm/(dlc+dlm)
         cxq2(2,icr)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
         cxq2(3,icr)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
         cxq2(4,icr)=h11*dlc*dlc/dlp/(dlp+dlc)
        enddo
      enddo


      ! for q3

      xxl(1:n1mr) = xmr(1:n1mr)
      xxl(0) = 2.d0*xxl(1) - xxl(2)
      xxl(n1r) = 2.d0*xxl(n1mr) - xxl(n1mr-1)
      xxs(1:n1m) = xm(1:n1m)
      xxs(0) = 2.d0*xxs(1) - xxs(2)
      xxs(-1) = 2.d0*xxs(0) - xxs(1)
      xxs(n1) = 2.d0*xxs(n1m) - xxs(n1m-1)
      xxs(n1+1) = 2.d0*xxs(n1) - xxs(n1m)

      yyl(1:n2mr) = ymr(1:n2mr)
      yyl(0) = 2.d0*yyl(1) - yyl(2)
      yyl(n2r) = 2.d0*yyl(n2mr) - yyl(n2mr-1)
      yys(1:n2m) = ym(1:n2m)
      yys(0) = 2.d0*yys(1) - yys(2)
      yys(-1) = 2.d0*yys(0) - yys(1)
      yys(n2) = 2.d0*yys(n2m) - yys(n2m-1)
      yys(n2+1) = 2.d0*yys(n2) - yys(n2m)

      zzl(1:n3r) = zcr(1:n3r)
      zzs(1:n3) = zc(1:n3)

      do kc=1,n3m
       if(kc.eq.1)then
        dlc = zzs(kc+1)-zzs(kc)
        dlp = zzs(kc+2)-zzs(kc+1)
        lza = 1.d0/dlc
        do kcr=max((kc-1)*mref3+1,1),min(kc*mref3,n3mr)
         lzm = (zzl(kcr) - zzs(kc))*lza
         lzp = 1.d0 - lzm
         h00=(1.d0+2.d0*lzm)*lzp*lzp
         h10=lzm*lzp*lzp
         h01=(1.d0+2.d0*lzp)*lzm*lzm
         h11=-lzp*lzm*lzm
         czq3(1,kcr)=0.d0                                               
         czq3(2,kcr)=h00-h10-h11*dlp/(dlp+dlc)
         czq3(3,kcr)=h10+h01+h11*(dlp-dlc)/dlp
         czq3(4,kcr)=h11*dlc*dlc/dlp/(dlp+dlc)
        enddo
       elseif(kc.eq.n3m)then
        dlc = zzs(kc+1)-zzs(kc)
        dlm = zzs(kc)-zzs(kc-1)
        lza = 1.d0/dlc
        do kcr=max((kc-1)*mref3+1,1),min(kc*mref3,n3mr)
         lzm = (zzl(kcr) - zzs(kc))*lza
         lzp = 1.d0 - lzm
         h00=(1.d0+2.d0*lzm)*lzp*lzp
         h10=lzm*lzp*lzp
         h01=(1.d0+2.d0*lzp)*lzm*lzm
         h11=-lzp*lzm*lzm
         czq3(1,kcr)=-h10*dlc*dlc/dlm/(dlc+dlm)
         czq3(2,kcr)=h00-h11+h10*(dlc-dlm)/dlm
         czq3(3,kcr)=h01+h11+h10*dlm/(dlm+dlc)
         czq3(4,kcr)=0.d0
        enddo
       else
        dlc = zzs(kc+1)-zzs(kc)
        dlm = zzs(kc)-zzs(kc-1)
        dlp = zzs(kc+2)-zzs(kc+1)
        lza = 1.d0/dlc
        do kcr=max((kc-1)*mref3+1,1),min(kc*mref3,n3mr)
         lzm = (zzl(kcr) - zzs(kc))*lza
         lzp = 1.d0 - lzm
         h00=(1.d0+2.d0*lzm)*lzp*lzp
         h10=lzm*lzp*lzp
         h01=(1.d0+2.d0*lzp)*lzm*lzm
         h11=-lzp*lzm*lzm
         czq3(1,kcr)=-h10*dlc*dlc/dlm/(dlc+dlm)
         czq3(2,kcr)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
         czq3(3,kcr)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
         czq3(4,kcr)=h11*dlc*dlc/dlp/(dlp+dlc)
        enddo
       endif
      enddo
      do jc=0,n2m
        dlc = yys(jc+1)-yys(jc)
        dlm = yys(jc)-yys(jc-1)
        dlp = yys(jc+2)-yys(jc+1)
        lya=1.d0/dlc
        do jcr=max(jrangs(jc),1),min(jrangs(jc+1)-1,n2mr)
         lym = (yyl(jcr) - yys(jc))*lya
         lyp = 1.d0 - lym
         h00=(1.d0+2.d0*lym)*lyp*lyp
         h10=lym*lyp*lyp
         h01=(1.d0+2.d0*lyp)*lym*lym
         h11=-lyp*lym*lym
         cyq3(1,jcr)=-h10*dlc*dlc/dlm/(dlc+dlm)
         cyq3(2,jcr)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
         cyq3(3,jcr)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
         cyq3(4,jcr)=h11*dlc*dlc/dlp/(dlp+dlc)
        enddo
      enddo
      do ic=0,n1m
        dlc = xxs(ic+1)-xxs(ic)
        dlm = xxs(ic)-xxs(ic-1)
        dlp = xxs(ic+2)-xxs(ic+1)
        lxa=1.d0/dlc
        do icr=max(irangs(ic),1),min(irangs(ic+1)-1,n1mr)
         lxm = (xxl(icr) - xxs(ic))*lxa
         lxp = 1.d0 - lxm
         h00=(1.d0+2.d0*lxm)*lxp*lxp
         h10=lxm*lxp*lxp
         h01=(1.d0+2.d0*lxp)*lxm*lxm
         h11=-lxp*lxm*lxm
         cxq3(1,icr)=-h10*dlc*dlc/dlm/(dlc+dlm)
         cxq3(2,icr)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
         cxq3(3,icr)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
         cxq3(4,icr)=h11*dlc*dlc/dlp/(dlp+dlc)
        enddo
      enddo

      !m=========================================================
      !  for second order interpolation du_i/dx_i 

      xxl(1:n1mr) = xmr(1:n1mr)
      xxl(0) = 2.d0*xxl(1) - xxl(2)
      xxl(n1r) = 2.d0*xxl(n1mr) - xxl(n1mr-1)
      xxs(1:n1m) = xm(1:n1m)
      xxs(0) = 2.d0*xxs(1) - xxs(2)
      xxs(-1) = 2.d0*xxs(0) - xxs(1)
      xxs(n1) = 2.d0*xxs(n1m) - xxs(n1m-1)
      xxs(n1+1) = 2.d0*xxs(n1) - xxs(n1m)

      yyl(1:n2mr) = ymr(1:n2mr)
      yyl(0) = 2.d0*yyl(1) - yyl(2)
      yyl(n2r) = 2.d0*yyl(n2mr) - yyl(n2mr-1)
      yys(1:n2m) = ym(1:n2m)
      yys(0) = 2.d0*yys(1) - yys(2)
      yys(-1) = 2.d0*yys(0) - yys(1)
      yys(n2) = 2.d0*yys(n2m) - yys(n2m-1)
      yys(n2+1) = 2.d0*yys(n2) - yys(n2m)

      zzl(0) = zcr(1)
      zzl(1:n3mr) = zmr(1:n3mr)
      zzl(n3r) = zcr(n3r)
      zzs(0) = zc(1)
      zzs(1:n3m) = zm(1:n3m)
      zzs(n3) = zc(n3)

      do kc=0,n3m
       if(kc.eq.0)then
        dlc = zzs(kc+1)-zzs(kc)
        dlp = zzs(kc+2)-zzs(kc+1)
        lza = 1.d0/dlc
        do kcr=max(krangs(kc),1),min(krangs(kc+1)-1,n3mr)
         lzm = (zzl(kcr) - zzs(kc))*lza
         lzp = 1.d0 - lzm
         h00=(1.d0+2.d0*lzm)*lzp*lzp
         h10=lzm*lzp*lzp
         h01=(1.d0+2.d0*lzp)*lzm*lzm
         h11=-lzp*lzm*lzm
         czrs(1,kcr)=0.d0
         czrs(2,kcr)=lzp  !h00-h11*dlp/(dlp+dlc)-dble(ubcbot)*h10
         czrs(3,kcr)=lzm  !h01+h11*(dlp-dlc)/dlp+dble(ubcbot)*h10
         czrs(4,kcr)=0.d0 !h11*dlc*dlc/dlp/(dlp+dlc)
        enddo
       elseif(kc.eq.n3m)then
        dlc = zzs(kc+1)-zzs(kc)
        dlm = zzs(kc)-zzs(kc-1)
        lza = 1.d0/dlc
        do kcr=max(krangs(kc),1),min(krangs(kc+1)-1,n3mr)
         lzm = (zzl(kcr) - zzs(kc))*lza
         lzp = 1.d0 - lzm
         h00=(1.d0+2.d0*lzm)*lzp*lzp
         h10=lzm*lzp*lzp
         h01=(1.d0+2.d0*lzp)*lzm*lzm
         h11=-lzp*lzm*lzm
         czrs(1,kcr)=0.d0 !-h10*dlc*dlc/dlm/(dlc+dlm)
         czrs(2,kcr)=lzp  !h00+h10*(dlc-dlm)/dlm-dble(ubctop)*h11
         czrs(3,kcr)=lzm  !h01+h10*dlm/(dlm+dlc)+dble(ubctop)*h11
         czrs(4,kcr)=0.d0
        enddo
       else
        dlc = zzs(kc+1)-zzs(kc)
        dlm = zzs(kc)-zzs(kc-1)
        dlp = zzs(kc+2)-zzs(kc+1)
        lza = 1.d0/dlc
        do kcr=max(krangs(kc),1),min(krangs(kc+1)-1,n3mr)
         lzm = (zzl(kcr) - zzs(kc))*lza
         lzp = 1.d0 - lzm
         h00=(1.d0+2.d0*lzm)*lzp*lzp
         h10=lzm*lzp*lzp
         h01=(1.d0+2.d0*lzp)*lzm*lzm
         h11=-lzp*lzm*lzm
         czrs(1,kcr)=-h10*dlc*dlc/dlm/(dlc+dlm)
         czrs(2,kcr)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
         czrs(3,kcr)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
         czrs(4,kcr)=h11*dlc*dlc/dlp/(dlp+dlc)
        enddo
       endif
      enddo
      do jc=0,n2m
        dlc = yys(jc+1)-yys(jc)
        dlm = yys(jc)-yys(jc-1)
        dlp = yys(jc+2)-yys(jc+1)
        lya=1.d0/dlc
        do jcr=max(jrangs(jc),1),min(jrangs(jc+1)-1,n2mr)
         lym = (yyl(jcr) - yys(jc))*lya
         lyp = 1.d0 - lym
         h00=(1.d0+2.d0*lym)*lyp*lyp
         h10=lym*lyp*lyp
         h01=(1.d0+2.d0*lyp)*lym*lym
         h11=-lyp*lym*lym
         cyrs(1,jcr)=-h10*dlc*dlc/dlm/(dlc+dlm)
         cyrs(2,jcr)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
         cyrs(3,jcr)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
         cyrs(4,jcr)=h11*dlc*dlc/dlp/(dlp+dlc)
        enddo
      enddo
      do ic=0,n1m
        dlc = xxs(ic+1)-xxs(ic)
        dlm = xxs(ic)-xxs(ic-1)
        dlp = xxs(ic+2)-xxs(ic+1)
        lxa=1.d0/dlc
        do icr=max(irangs(ic),1),min(irangs(ic+1)-1,n1mr)
         lxm = (xxl(icr) - xxs(ic))*lxa
         lxp = 1.d0 - lxm
         h00=(1.d0+2.d0*lxm)*lxp*lxp
         h10=lxm*lxp*lxp
         h01=(1.d0+2.d0*lxp)*lxm*lxm
         h11=-lxp*lxm*lxm
         cxrs(1,icr)=-h10*dlc*dlc/dlm/(dlc+dlm)
         cxrs(2,icr)=h00-h11*dlp/(dlp+dlc)+h10*(dlc-dlm)/dlm
         cxrs(3,icr)=h01+h10*dlm/(dlm+dlc)+h11*(dlp-dlc)/dlp
         cxrs(4,icr)=h11*dlc*dlc/dlp/(dlp+dlc)
        enddo
      enddo

      return
      end subroutine mgrd_idc
!
