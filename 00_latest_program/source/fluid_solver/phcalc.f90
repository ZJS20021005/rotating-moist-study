!  this subroutine perform the calculation of dph , periodic direction
!  along x3 and x1to use the real fourier transform
      subroutine phcalc
      use param
      use local_arrays, only: dph
      use mpi_param
      use mpih
      implicit none
      !include "mkl_service.fi"
      integer :: i,j,k,info
      real :: coefnorm,acphT_b(m3m)
      real :: xr(m2m,m1m)
      complex :: xa(m2mh,m1m)
      real :: amphT(m3m-1),acphT(m3m),apphT(m3m-1),drhs(m3m)
      integer :: n2mh,jmh
      real,allocatable,dimension(:,:,:) :: dpht,dpho

      allocate(dpht(1:n3m,1:n1m,jstartp:jendp))
      allocate(dpho(1:n1m,1:n2m+2,kstart:kend))

      n2mh=n2m/2+1

      coefnorm = 1.d0/(dble(n1m)*dble(n2m))

      do k=kstart,kend
        do j=1,n2m
          do i=1,n1m
           xr(j,i)=dph(i,j,k)
          enddo
        enddo
        
        call dfftw_execute_dft_r2c(fwd_plan,xr,xa)

        do j=1,n2mh
         do i=1,n1m
         dpho(i,j,k)=dreal(xa(j,i))*coefnorm
         dpho(i,j+n2mh,k)=dimag(xa(j,i))*coefnorm
        enddo
        enddo
      end do

!=================================================
      call PackZ_UnpackRP(dpho,dpht)
!=================================================
      call mkl_set_num_threads(1)

      do j=jstartp,jendp
        jmh=jmhv(j)
        do i=1,n1m
          do k = 1,n3m
            acphT_b(k)=1.d0/(acphk(k)-ak2(jmh)-ak1(i))
            drhs(k)=dpht(k,i,j)*acphT_b(k)
          enddo

          amphT(1:n3m-1)=amphk(2:n3m)*acphT_b(2:n3m)
          apphT(1:n3m-1)=apphk(1:n3m-1)*acphT_b(1:n3m-1)
          acphT(1:n3m)=1.d0

          call ddttrfb(n3m,amphT,acphT,apphT,info)
          call ddttrsb('N',n3m,1,amphT,acphT,apphT,drhs,n3m,info)

          do k=1,n3m
            dpht(k,i,j) = drhs(k)
          enddo
        enddo
      enddo

      call mkl_set_num_threads(numthreads)
!=============================================
      call PackR_UnpackZP(dpht,dpho)
!=============================================
      do k=kstart,kend
       do j=1,n2mh
        do i=1,n1m
          xa(j,i)=dcmplx(dpho(i,j,k),dpho(i,j+n2mh,k))
        enddo
       end do

      call dfftw_execute_dft_c2r(bck_plan,xa,xr)

       do j=1,n2m
         do i=1,n1m
           dph(i,j,k)=xr(j,i)
         enddo
       end do
      end do

      if(allocated(dpht)) deallocate(dpht)
      if(allocated(dpho)) deallocate(dpho)

      return
      end
