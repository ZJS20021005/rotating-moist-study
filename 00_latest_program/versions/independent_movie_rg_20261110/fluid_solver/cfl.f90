      subroutine cfl(cflm)
      use param
      use local_arrays, only: q2,q3,q1
      use mpih
      use mpi_param, only: kstart,kend
      implicit none
      real,intent(inout)    :: cflm
      real    :: my_cflm
      integer :: j,k,jp,kp,i,ip
      real :: qcf,udx3
      
      my_cflm=1.d-8
                                                                      
      do k=kstart,kend
        udx3=udx3m(k)
        kp=k+1
        do j=1,n2m
          jp=j+1
          do i=1,n1m
            ip=i+1
            qcf=( dabs((q1(i,j,k)+q1(ip,j,k))*0.5d0*dx1)&
     &           +dabs((q2(i,j,k)+q2(i,jp,k))*0.5d0*dx2)&
     &           +dabs((q3(i,j,k)+q3(i,j,kp))*0.5d0*udx3))

            my_cflm = dmax1(my_cflm,qcf)
          enddo
        enddo
      enddo
            
      call MPI_ALLREDUCE(my_cflm,cflm,1,MDP,MPI_MAX,MPI_COMM_WORLD,ierr)

      return  
      end                                                               
