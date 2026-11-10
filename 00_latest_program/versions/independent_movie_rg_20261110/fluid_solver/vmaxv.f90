      subroutine vmaxv
!EP   This routine calculates the maximum velocities
      use param
      use local_arrays, only: q2,q3,q1,dens,dsal
      use mgrd_arrays, only: q3lr
      use mpi_param, only: kstart,kend,kstartr,kendr
      use mpih
      implicit none
      real    :: my_vmax2,my_vmax3,my_vmax1
      integer :: jc,kc,kp,ic

      my_vmax1=-100.d0
      my_vmax2=-100.d0
      my_vmax3=-100.d0
      do kc=kstart,kend
        kp = kc + 1
        do jc=1,n2m
          do ic=1,n1m
            my_vmax1 = dmax1(my_vmax1,dabs(q1(ic,jc,kc)))
            my_vmax2 = dmax1(my_vmax2,dabs(q2(ic,jc,kc)))
            my_vmax3 = dmax1(my_vmax3,dabs(q3(ic,jc,kc)))
          enddo
        enddo
      enddo

      call MPI_ALLREDUCE(my_vmax1,vmax(1),1,MDP,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_vmax2,vmax(2),1,MDP,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_vmax3,vmax(3),1,MDP,MPI_MAX,MPI_COMM_WORLD,ierr)

      if(myid.eq.0)then
        write(95,510) time, dabs(vmax(1)), dabs(vmax(2)), dabs(vmax(3))
      endif
510   format(1x,f10.4,3(1x,ES20.8))

      return   
      end
