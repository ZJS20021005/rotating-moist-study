      subroutine hdnlq2
      use param
      use local_arrays, only: q1,q2,q3,dph,dens
      use mgrd_arrays, only: dsalc
      use mpi_param, only: kstart,kend
      implicit none
      integer :: kc,kp,jp,jm,jc,ic,im,ip
      integer :: kmm,kpp
      real    :: h22,h23,udx1,udx2,h21
      real    :: dsalit,coriol
      real    :: uxc,uzc
      real    :: porous_term,friction_term
      udx1=dx1*0.25d0
      udx2=dx2*0.25d0

    do kc=kstart,kend
        kmm=kc-1
        kpp=kc+1
        kp=kc+1
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
          do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)

               h21=( (q2(ip,jc,kc)+q2(ic,jc,kc))*(q1(ip,jc,kc)+q1(ip,jm,kc))&
        &           -(q2(ic,jc,kc)+q2(im,jc,kc))*(q1(ic,jc,kc)+q1(ic,jm,kc)))*udx1

               h22=( (q2(ic,jp,kc)+q2(ic,jc,kc))*(q2(ic,jp,kc)+q2(ic,jc,kc))&
        &           -(q2(ic,jm,kc)+q2(ic,jc,kc))*(q2(ic,jm,kc)+q2(ic,jc,kc)))*udx2

            if(kc.eq.1) then
               h23=( (q3(ic,jc,kp)+q3(ic,jm,kp))*(q2(ic,jc,kpp)+q2(ic,jc,kc)))*udx3m(kc)*0.25d0
            elseif(kc.eq.n3m) then
               h23=(-(q3(ic,jc,kc)+q3(ic,jm,kc))*(q2(ic,jc,kc)+q2(ic,jc,kmm)))*udx3m(kc)*0.25d0
            else
               h23=( (q3(ic,jc,kp)+q3(ic,jm,kp))*(q2(ic,jc,kpp)+q2(ic,jc,kc))&
        &           -(q3(ic,jc,kc)+q3(ic,jm,kc))*(q2(ic,jc,kc)+q2(ic,jc,kmm)))*udx3m(kc)*0.25d0
            endif
            uxc = ( q1(ic,jc,kc)+q1(ip,jc,kc)  &
        &          +q1(ic,jm,kc)+q1(ip,jm,kc))*0.25d0
            uzc = ( q3(ic,jc,kc)+q3(ic,jc,kp)  &
        &          +q3(ic,jm,kc)+q3(ic,jm,kp))*0.25d0
            friction_term = -alpha*q2(ic,jc,kc) 
            coriol = - ( invRo*uxc )
!            porous_term = -invDa*dsqrt(Prs/Ra)*1.0*q2(ic,jc,kc)
 
            dph(ic,jc,kc)=-(h21+h22+h23)+friction_term + coriol

          enddo
        enddo
      enddo
      
      return
      end

