      subroutine hdnlqvap
      use param
      use local_arrays, only: hqvap,qvap,T,qs,cond_term
      use mgrd_arrays, only: q1lr,q2lr,q3lr
      use mpi_param, only: kstartr,kendr
      implicit none
      integer :: jc,kc,ic
      integer :: kpp,kmm,kp,jp,jmm,jpp,ip,imm,ipp
      real    :: h31,h32,h33,udx2r,udx1r
      real    :: stepfunction_term

       udx1r=dx1r*0.5d0
      udx2r=dx2r*0.5d0

      do kc=kstartr,kendr
        kmm=kc-1
        kpp=kc+1
        kp=kc+1
        do jc=1,n2mr
          jp=jc+1
          jmm=jmvr(jc)
          jpp=jpvr(jc)
          do ic=1,n1mr
            ip=ic+1
            ipp=ipvr(ic)
            imm=imvr(ic)

               h31=( q1lr(ip,jc,kc)*(qvap(ipp,jc,kc)+qvap(ic,jc,kc)) &
        &           -q1lr(ic,jc,kc)*(qvap(ic,jc,kc)+qvap(imm,jc,kc)))*udx1r

               h32=( q2lr(ic,jp,kc)*(qvap(ic,jpp,kc)+qvap(ic,jc,kc)) &
        &           -q2lr(ic,jc,kc)*(qvap(ic,jc,kc)+qvap(ic,jmm,kc)))*udx2r

            if(kc.eq.1) then
               h33=( (qvap(ic,jc,kpp)*g3rmr(kc)+qvap(ic,jc,kc)*g3rmr(kpp))&
        &              /(g3rmr(kc)+g3rmr(kpp))*q3lr(ic,jc,kp))*udx3mr(kc)
            elseif(kc.eq.n3m) then
               h33=(-(qvap(ic,jc,kc)*g3rmr(kmm)+qvap(ic,jc,kmm)*g3rmr(kc))&
        &              /(g3rmr(kc)+g3rmr(kmm))*q3lr(ic,jc,kc)&
        &          )*udx3mr(kc)
            else
               h33=( (qvap(ic,jc,kpp)*g3rmr(kc)+qvap(ic,jc,kc)*g3rmr(kpp))&
        &              /(g3rmr(kc)+g3rmr(kpp))*q3lr(ic,jc,kp)&
        &           -(qvap(ic,jc,kc)*g3rmr(kmm)+qvap(ic,jc,kmm)*g3rmr(kc))&
        &              /(g3rmr(kc)+g3rmr(kmm))*q3lr(ic,jc,kc)&
        &          )*udx3mr(kc)
            endif
            stepfunction_term = -cond_term(ic,jc,kc)
            hqvap(ic,jc,kc)= -(h31+h32+h33)+stepfunction_term
     
          enddo
        enddo
      enddo

      return
      end

