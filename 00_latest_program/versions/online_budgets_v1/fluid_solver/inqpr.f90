!***********************************************************************
!                       INITIAL CONDITION                              *
!***********************************************************************
      subroutine inqpr
      use param
      use local_arrays, only: q2,q3,dens,q1,dsal,qvap
      use mpi_param, only: kstart,kend,kstartr,kendr
      use mpih
      implicit none

      integer :: j,k,i,p,ios,nprof,seed_size
      integer, allocatable :: random_seed_values(:)
      real :: perturb_amp,varptb,r,u1,u2,phi,env,zeta
      real :: bbase,qbase,frac,tolerance
      real :: profile_params(13),expected_params(12)
      real, allocatable :: zprof(:),bprof(:),qprof(:)
      logical :: profile_exists

!     The profile is generated before submission by the user-supplied
!     linear-stability drizzle solver. File format:
!       row 1: nprofile, buoyancy perturbation amplitude
!       row 2: Ra,Pr,invRo,gamma,Sm,alphaqs,betaqs,tau_cond,
!              bbot,btop,qbot,qtop,saturation_width
!       rows 3+: z/H, drizzle b(z), drizzle q(z)
!
!     Only buoyancy is perturbed:
!       b(x,y,z)=b_drizzle(z)+amp*sin(pi*z/H)*N(0,1)
!       q(x,y,z)=q_drizzle(z)
!     Velocity starts from zero.

      tolerance = 2.0d-10
      profile_exists = .false.
      nprof = 0
      perturb_amp = 0.0d0
      profile_params = 0.0d0

      if(myid.eq.0)then
        inquire(file='drizzle_init.dat',exist=profile_exists)
        if(.not.profile_exists)then
          write(*,*) 'ERROR: nread=0 requires drizzle_init.dat'
          call MPI_ABORT(MPI_COMM_WORLD,71,ierr)
        endif

        open(unit=97,file='drizzle_init.dat',status='old',iostat=ios)
        if(ios.ne.0)then
          write(*,*) 'ERROR: cannot open drizzle_init.dat, iostat=',ios
          call MPI_ABORT(MPI_COMM_WORLD,72,ierr)
        endif
        read(97,*,iostat=ios) nprof,perturb_amp
        if(ios.ne.0 .or. nprof.lt.33)then
          write(*,*) 'ERROR: invalid drizzle profile header'
          call MPI_ABORT(MPI_COMM_WORLD,73,ierr)
        endif
        read(97,*,iostat=ios) profile_params
        if(ios.ne.0)then
          write(*,*) 'ERROR: invalid drizzle parameter row'
          call MPI_ABORT(MPI_COMM_WORLD,74,ierr)
        endif

        allocate(zprof(nprof),bprof(nprof),qprof(nprof))
        do p=1,nprof
          read(97,*,iostat=ios) zprof(p),bprof(p),qprof(p)
          if(ios.ne.0)then
            write(*,*) 'ERROR: invalid drizzle profile row ',p
            call MPI_ABORT(MPI_COMM_WORLD,75,ierr)
          endif
        enddo
        close(97)

        if(perturb_amp.le.0.0d0 .or. perturb_amp.gt.1.0d-3)then
          write(*,*) 'ERROR: invalid drizzle perturbation amplitude ', &
     &               perturb_amp
          call MPI_ABORT(MPI_COMM_WORLD,76,ierr)
        endif
        if(dabs(A_stopmod).gt.1.0d-14 .or. &
     &     dabs(A_sbotmod).gt.1.0d-14)then
          write(*,*) 'ERROR: 1-D drizzle initialization requires ', &
     &               'A_stopmod=A_sbotmod=0'
          call MPI_ABORT(MPI_COMM_WORLD,77,ierr)
        endif
        if(dabs(zprof(1)).gt.tolerance .or. &
     &     dabs(zprof(nprof)-1.0d0).gt.tolerance)then
          write(*,*) 'ERROR: drizzle profile must cover 0<=z/H<=1'
          call MPI_ABORT(MPI_COMM_WORLD,78,ierr)
        endif
        do p=1,nprof-1
          if(zprof(p+1).le.zprof(p))then
            write(*,*) 'ERROR: drizzle z grid is not increasing at ',p
            call MPI_ABORT(MPI_COMM_WORLD,79,ierr)
          endif
        enddo

        expected_params(1) = Ra
        expected_params(2) = Prs
        expected_params(3) = invRo
        expected_params(4) = gamma
        expected_params(5) = Sm
        expected_params(6) = alphaqs
        expected_params(7) = betaqs
        expected_params(8) = tau_cond
        expected_params(9) = dsalbot
        expected_params(10) = dsaltop
        expected_params(11) = qvapbot
        expected_params(12) = qvaptop
        do p=1,12
          if(dabs(profile_params(p)-expected_params(p)).gt. &
     &       tolerance*dmax1(1.0d0,dabs(expected_params(p))))then
            write(*,*) 'ERROR: drizzle parameter mismatch, index=',p, &
     &                 ' profile=',profile_params(p), &
     &                 ' bou.in=',expected_params(p)
            call MPI_ABORT(MPI_COMM_WORLD,80,ierr)
          endif
        enddo
        if(dabs(profile_params(13)-1.0d-8).gt.1.0d-12)then
          write(*,*) 'ERROR: drizzle saturation width must match DNS ', &
     &               'width=1e-8, profile=',profile_params(13)
          call MPI_ABORT(MPI_COMM_WORLD,81,ierr)
        endif

        if(dabs(bprof(1)-dsalbot).gt.tolerance .or. &
     &     dabs(bprof(nprof)-dsaltop).gt.tolerance .or. &
     &     dabs(qprof(1)-qvapbot).gt.tolerance .or. &
     &     dabs(qprof(nprof)-qvaptop).gt.tolerance)then
          write(*,*) 'ERROR: drizzle endpoints do not match bou.in'
          call MPI_ABORT(MPI_COMM_WORLD,82,ierr)
        endif
      endif

      call MPI_BCAST(nprof,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(perturb_amp,1,MDP,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(profile_params,13,MDP,0,MPI_COMM_WORLD,ierr)
      if(myid.ne.0) allocate(zprof(nprof),bprof(nprof),qprof(nprof))
      call MPI_BCAST(zprof,nprof,MDP,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(bprof,nprof,MDP,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(qprof,nprof,MDP,0,MPI_COMM_WORLD,ierr)

!     Deterministic rank-dependent random seeds avoid repeating the same
!     horizontal pattern on every MPI z slab.
      call random_seed(size=seed_size)
      allocate(random_seed_values(seed_size))
      do p=1,seed_size
        random_seed_values(p) = mod(13579+104729*(myid+1)+7919*p, &
     &                              2147483000)
      enddo
      call random_seed(put=random_seed_values)
      deallocate(random_seed_values)

      do k=kstart,kend
        do j=1,n2m
          do i=1,n1m
            dens(i,j,k) = 0.0d0
          enddo
        enddo
      enddo

      do k=kstartr,kendr
        zeta = zmr(k)/alx3
        if(zeta.le.zprof(1))then
          bbase = bprof(1)
          qbase = qprof(1)
        elseif(zeta.ge.zprof(nprof))then
          bbase = bprof(nprof)
          qbase = qprof(nprof)
        else
          p = 1
          do while(p.lt.nprof-1 .and. zprof(p+1).lt.zeta)
            p = p+1
          enddo
          frac = (zeta-zprof(p))/(zprof(p+1)-zprof(p))
          bbase = bprof(p)+frac*(bprof(p+1)-bprof(p))
          qbase = qprof(p)+frac*(qprof(p+1)-qprof(p))
        endif

        env = dsin(pi*zeta)
        do j=1,n2mr
          do i=1,n1mr
            call random_number(u1)
            call random_number(u2)
            if(u1.le.1.0d-12) u1 = 1.0d-12
            r = dsqrt(-2.0d0*dlog(u1))
            phi = 2.0d0*pi*u2
            varptb = r*dcos(phi)
            dsal(i,j,k) = bbase+perturb_amp*env*varptb
            qvap(i,j,k) = qbase
          enddo
        enddo
      enddo

!     Initialize saturation and condensation from the final perturbed b and
!     unperturbed drizzle q fields.
      call compute_qsat_from_b

      do k=kstart-1,kend+1
        do j=1,n2
          do i=1,n1
            q1(i,j,k) = 0.0d0
            q2(i,j,k) = 0.0d0
            q3(i,j,k) = 0.0d0
          enddo
        enddo
      enddo

      if(myid.eq.0)then
        write(*,'(3x,a,i6)') 'Drizzle profile points: ',nprof
        write(*,'(3x,a,es12.4)') 'Buoyancy perturbation amplitude: ', &
     &                            perturb_amp
        write(*,'(3x,a,es12.4)') 'Drizzle saturation width: ', &
     &                            profile_params(13)
      endif

      deallocate(zprof,bprof,qprof)
      return
      end
