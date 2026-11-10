      subroutine fftqua
      use param
      integer :: n2mh,n2mp,j,i,n1mh,n1mp
      real    :: ao

      n1mh=n1m/2+1
      n2mh=n2m/2+1
      !  Eigenvalue
      do i=1,n1m
        if(i.le.n1mh) then
         ao=(i-1)*2.d0*pi
        else
         ao=-(n1m-i+1)*2.d0*pi
        endif
        ak1(i)=2.d0*(1.d0-dcos(ao/n1m))*(dble(n1m)/rext1)**2
      enddo

      do j=1,n2m
        if(j.le.n2mh) then
         ao=(j-1)*2.d0*pi
        else
         ao=-(n2m-j+1)*2.d0*pi
        endif
        ak2(j)=2.d0*(1.d0-dcos(ao/n2m))*(dble(n2m)/rext2)**2
      enddo

      return
      end

!=====================================================
      subroutine phini
      use param
      use mpi_param
      use mpih
      implicit none
      integer,parameter ::  fftw_es =64

      real, dimension(m2m,m1m) :: xr
      real, dimension(m2m,m1m) :: xa
    
      !RO   Initialize tridiag matrices
      call tridiag_matrices   

      !m    Initialize FFTW
      call fftqua

      call dfftw_plan_dft_r2c_2d(fwd_plan,m2m,m1m,xr,xa,fftw_es)
      call dfftw_plan_dft_c2r_2d(bck_plan,m2m,m1m,xa,xr,fftw_es)

      return
      end
      
!=======================================================================
      subroutine tridiag_matrices
      use param
      implicit none
      integer  :: kc,km,kp
      real :: ugmmm,a33icc,a33icp

      !   tridiagonal matrix coefficients at each k and i
      !   x1 and x3 cartesian coordinates
      do kc=1,n3m
        km=kmv(kc)
        kp=kpv(kc)
        a33icc=kmc(kc)*dx3q/g3rc(kc)
        a33icp=kpc(kc)*dx3q/g3rc(kp)
        ugmmm=1.0d0/g3rm(kc)
        amphk(kc)=a33icc*ugmmm
        apphk(kc)=a33icp*ugmmm
        if(kc.eq.1)then
        acphk(kc)=-(apphk(kc))
        elseif(kc.eq.n3m)then
        acphk(kc)=-(amphk(kc))
        else
        acphk(kc)=-(amphk(kc)+apphk(kc))
        endif
      enddo

      end subroutine tridiag_matrices

!==================================================
      subroutine phend
      use param
      implicit none

      call dfftw_destroy_plan(fwd_plan)
      call dfftw_destroy_plan(bck_plan)

      call dfftw_cleanup_threads()

      return
      end subroutine phend
