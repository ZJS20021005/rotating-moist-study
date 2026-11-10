      subroutine hdnlq1
      use param
      use local_arrays, only: q1,q2,q3,dq
      use mpi_param, only: kstart,kend
      implicit none
      integer :: kc,kp,jp,jm,jc,ic,im,ip,js
      integer :: kmm,kpp
      real    :: h11,h12,h13,udx1,udx2
      real    :: porous_term,friction_term      
      real    :: coriol
      real    :: uyc,uzc

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

               h11=( (q1(ip,jc,kc)+q1(ic,jc,kc))*(q1(ip,jc,kc)+q1(ic,jc,kc))&
        &           -(q1(im,jc,kc)+q1(ic,jc,kc))*(q1(im,jc,kc)+q1(ic,jc,kc)))*udx1

               h12=( (q2(ic,jp,kc)+q2(im,jp,kc))*(q1(ic,jp,kc)+q1(ic,jc,kc))&
        &           -(q2(ic,jc,kc)+q2(im,jc,kc))*(q1(ic,jc,kc)+q1(ic,jm,kc)))*udx2

            if(kc.eq.1) then
               h13=( (q3(ic,jc,kp)+q3(im,jc,kp))*(q1(ic,jc,kpp)+q1(ic,jc,kc)))*udx3m(kc)*0.25d0
            elseif(kc.eq.n3m) then
               h13=( -(q3(ic,jc,kc)+q3(im,jc,kc))*(q1(ic,jc,kc)+q1(ic,jc,kmm)))*udx3m(kc)*0.25d0
            else
               h13=( (q3(ic,jc,kp)+q3(im,jc,kp))*(q1(ic,jc,kpp)+q1(ic,jc,kc))&
        &           -(q3(ic,jc,kc)+q3(im,jc,kc))*(q1(ic,jc,kc)+q1(ic,jc,kmm)))*udx3m(kc)*0.25d0
            endif

            uyc = ( q2(ic,jc,kc)+q2(im,jc,kc) +  &
     &              q2(ic,jp,kc)+q2(im,jp,kc) )*0.25d0
            uzc = ( q3(ic,jc,kc)+q3(im,jc,kc) +  &
     &              q3(ic,jc,kp)+q3(im,jc,kp) )*0.25d0
            coriol = - ( - invRo*uyc )
             friction_term = -alpha*q1(ic,jc,kc)
 !           porous_term = -invDa*dsqrt(Prs/Ra)*1.0*q1(ic,jc,kc)

            dq(ic,jc,kc)=-(h11+h12+h13)+coriol+friction_term

          enddo
        enddo
      enddo

      return
      end

