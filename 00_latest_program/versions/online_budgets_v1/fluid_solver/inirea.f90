      subroutine inirea
      use mpih
      use mpi_param, only: kstart,kend,kstartr,kendr
      use local_arrays, only: dens,q2,q3,q1,dsal,qvap
      use param
      IMPLICIT NONE
      real :: dummyr 
      character(70) :: filcnw2
      integer :: i,j
      integer, parameter :: ghosts = 4
      
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      !  One to one HDF read
      call mpi_read_continua(n1,n2,n3,kstart,kend,1,q1)
      call mpi_read_continua(n1,n2,n3,kstart,kend,2,q2)
      call mpi_read_continua(n1,n2,n3,kstart,kend,3,q3)
      call mpi_read_continua(n1r,n2r,n3r,kstartr,kendr,5,dsal)
      call mpi_read_continua(n1r,n2r,n3r,kstartr,kendr,8,qvap)

      if (ireset.eq.1) then                                             
        time=0.d0
      endif                                                             

!     Reading old grid information by rank0
      if (myid .eq. 0) then
        filcnw2 = 'continua_grid.dat'
        open(13,file=filcnw2,status='old')
        rewind(13)    
        read(13,*) dummyr,dummyr,dummyr                                                
        read(13,*) dummyr,dummyr,time
        close(13)
        write(*,'(5x,a)') 'old grid info read'
      endif
      call MPI_BCAST(time,1,MDP,0,MPI_COMM_WORLD,ierr)

      return                                                            
      end                                                                                                                                  
