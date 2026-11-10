      subroutine cflr(cflmr)
      use param
      use mgrd_arrays, only: q2lr,q3lr,q1lr
      use mpih
      use mpi_param, only: kstartr,kendr
      implicit none
      real,intent(inout)    :: cflmr
      real    :: my_cflmr
      integer :: j,k,jp,kp,i,ip
      real :: qcf,udx3
      
      my_cflmr=1.d-8
                                                                       
      do k=kstartr,kendr
        udx3=udx3mr(k)
        kp=k+1
        do j=1,n2mr
          jp=j+1
          do i=1,n1mr
            ip=i+1
            qcf=( dabs((q1lr(i,j,k)+q1lr(ip,j,k))*0.5d0*dx1r)&
     &           +dabs((q2lr(i,j,k)+q2lr(i,jp,k))*0.5d0*dx2r)&
     &           +dabs((q3lr(i,j,k)+q3lr(i,j,kp))*0.5d0*udx3))

            my_cflmr = dmax1(my_cflmr,qcf)
          enddo
        enddo
      enddo
            
      call MPI_ALLREDUCE(my_cflmr,cflmr,1,MDP,MPI_MAX,MPI_COMM_WORLD,ierr)

      return  
      end                                                               
