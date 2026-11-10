      subroutine invtrq3
      use param
      use local_arrays, only: q3,qcap,pr,ru3,rhs
      use mpi_param, only: kstart,kend

      implicit none
      integer :: jc,kc,km,kp,jp,jm,ic,ip,im
      real    :: udx3
      real    :: dq32,dq33,dcq3,dpx33,dq31
      real    :: app,acc,amm
      real    :: alre,udx1q,udx2q
      integer :: kstartp
      
      if(kstart.eq.1) then
      kstartp=2
      else
      kstartp=kstart
      endif
      alre=al*nu

      udx1q=dx1q
      udx2q=dx2q

      do kc=kstartp,kend
         km=kmv(kc)
         kp=kc+1
         udx3=al*udx3c(kc)
         amm=am3ck(kc)
         acc=ac3ck(kc)
         app=ap3ck(kc)
         do jc=1,n2m
            jm=jmv(jc)
            jp=jpv(jc)
            do ic=1,n1m
               im=imv(ic)
               ip=ipv(ic)

               dq31=(q3(im,jc,kc)-2.d0*q3(ic,jc,kc)+q3(ip,jc,kc))*udx1q
               dq32=(q3(ic,jm,kc)-2.d0*q3(ic,jc,kc)+q3(ic,jp,kc))*udx2q

               if(kc.eq.2) then
                  dq33=q3(ic,jc,kp)*app+q3(ic,jc,kc)*acc+0.0d0*amm
               elseif(kc.eq.n3m) then
                  dq33=0.0d0*app+q3(ic,jc,kc)*acc+q3(ic,jc,km)*amm
               else
                  dq33=q3(ic,jc,kp)*app+q3(ic,jc,kc)*acc+q3(ic,jc,km)*amm
               endif

               dcq3=dq32+dq33+dq31
               dpx33=(pr(ic,jc,kc)-pr(ic,jc,km))*udx3
 
               rhs(ic,jc,kc)=(ga*qcap(ic,jc,kc)+ro*ru3(ic,jc,kc)+alre*dcq3-dpx33)*dt 
               ru3(ic,jc,kc)=qcap(ic,jc,kc)
            enddo
         enddo
      enddo

      !-- implicit matrix solver
      call solxi_periodic(beta*al*dx1q)
      call solxj_periodic(beta*al*dx2q)
      call solq3k

      return
      end
