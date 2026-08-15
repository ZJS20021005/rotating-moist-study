!***********************************************************************
      subroutine divgck(qmaxc,qmaxr)
      use param
      use local_arrays, only: q2,q3,q1
      use mgrd_arrays, only: q2lr,q3lr,q1lr
      use mpih
      use mpi_param, only: kstart,kend,kstartr,kendr
      implicit none
      real,intent(out) :: qmaxc, qmaxr
      integer :: jc,kc,kp,jp,ic,ip
      real    :: dqcap,my_qmaxc,my_qmaxr

      !-- Base grid
      my_qmaxc =-100.d0
      do kc=kstart,kend
        kp=kc+1
        if(kc.ge.2 .and. kc.le.n3m-1)then
        do jc=1,n2m
          jp=jpv(jc)
          do ic=1,n1m
            ip=ipv(ic)
            dqcap= (q1(ip,jc,kc)-q1(ic,jc,kc))*dx1&
     &            +(q2(ic,jp,kc)-q2(ic,jc,kc))*dx2&
     &            +(q3(ic,jc,kp)-q3(ic,jc,kc))*udx3m(kc)
            my_qmaxc = dmax1(dabs(dqcap),my_qmaxc)          
          enddo
        enddo
        endif
      enddo
     
      qmaxc=0.d0
      call MPI_ALLREDUCE(my_qmaxc,qmaxc,1,MDP,MPI_MAX,MPI_COMM_WORLD,ierr)


      !-- Refined grid
      my_qmaxr =-100.d0
      do kc=kstartr,kendr
        kp=kc+1
        if(kc.ge.mref3+1 .and. kc.le.(n3m)*mref3-1)then
        do jc=mref2,(n2m)*mref2
          jp=jpvr(jc)
          do ic=mref1,(n1m)*mref1
            ip=ipvr(ic)
            dqcap= (q1lr(ip,jc,kc)-q1lr(ic,jc,kc))*dx1r&
     &            +(q2lr(ic,jp,kc)-q2lr(ic,jc,kc))*dx2r&
     &            +(q3lr(ic,jc,kp)-q3lr(ic,jc,kc))*udx3mr(kc)
            my_qmaxr = dmax1(dabs(dqcap),my_qmaxr)
          enddo
        enddo
        endif
      enddo

      qmaxr=0.d0
      call MPI_ALLREDUCE(my_qmaxr,qmaxr,1,MDP,MPI_MAX,MPI_COMM_WORLD,ierr)
      
      return     
      end         
