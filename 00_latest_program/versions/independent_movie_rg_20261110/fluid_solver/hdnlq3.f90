      subroutine hdnlq3
      use param
      use local_arrays, only: q1,q2,q3,qcap,dens
      use local_arrays, only:dsal
      use mgrd_arrays, only:dsalc
      use mpi_param, only: kstart,kend
      implicit none
      integer :: ic,jc,kc,js
      integer :: km,kp,jm,jp,im,ip
      real    :: h32,h33,h31
      real    :: densit,dsalit,udx1,udx2
      integer :: kstartp
      real    :: porous_term,friction_term

      if(kstart.eq.1) then
        kstartp=2
      else
        kstartp=kstart
      endif

      udx1=dx1*0.25d0
      udx2=dx2*0.25d0

      do kc=kstartp,kend
        km=kc-1
        kp=kc+1
        do jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
          do ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)

               h31=( (q1(ip,jc,kc)+q1(ip,jc,km))*(q3(ip,jc,kc)+q3(ic,jc,kc))&
        &           -(q1(ic,jc,kc)+q1(ic,jc,km))*(q3(ic,jc,kc)+q3(im,jc,kc)))*udx1

               h32=( (q2(ic,jp,kc)+q2(ic,jp,km))*(q3(ic,jp,kc)+q3(ic,jc,kc))&
        &           -(q2(ic,jc,kc)+q2(ic,jc,km))*(q3(ic,jc,kc)+q3(ic,jm,kc)))*udx2

            if(kc.eq.2) then
               h33=( (q3(ic,jc,kp)+q3(ic,jc,kc))*(q3(ic,jc,kp)+q3(ic,jc,kc))&
        &           -(q3(ic,jc,kc)+       0.0d0)*(q3(ic,jc,kc)+      0.0d0))*udx3c(kc)*0.25d0
            elseif(kc.eq.n3m) then
               h33=( (0.0d0       +q3(ic,jc,kc))*(0.0d0       +q3(ic,jc,kc))&
        &           -(q3(ic,jc,kc)+q3(ic,jc,km))*(q3(ic,jc,kc)+q3(ic,jc,km)))*udx3c(kc)*0.25d0
            else
               h33=( (q3(ic,jc,kp)+q3(ic,jc,kc))*(q3(ic,jc,kp)+q3(ic,jc,kc))&
        &           -(q3(ic,jc,kc)+q3(ic,jc,km))*(q3(ic,jc,kc)+q3(ic,jc,km)))*udx3c(kc)*0.25d0
            endif

            dsalit = (dsalc(ic,jc,kc)+dsalc(ic,jc,kc-1))*0.5d0
            friction_term = -alpha*q3(ic,jc,kc)
!            porous_term = -invDa*dsqrt(Prs/Ra)*1.0*q3(ic,jc,kc)

            qcap(ic,jc,kc) = -(h31+h32+h33) + dsalit + friction_term
!              qcap(ic,jc,kc) = -(h31+h32+h33) +  friction_term
          enddo
        enddo
      enddo
      return
      end
