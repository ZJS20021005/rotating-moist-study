!===================================================
!     this code is made for cartesian coordinates     
!     walls in all directions  
!===================================================                                                                                                     
      program papero
      use mpih
      use mpi_param
      use param
      implicit none
      character(len=4) :: dummy
      character(len=256) :: input_line
      integer :: tfini,tfin,n,ns
      integer :: read_ios
      integer :: stat_me_dummy
      real :: ts,te
      real :: tframe_me_dummy
      real :: movie_values(8), diagnostic_values(3)
      integer :: omp_get_num_threads

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)

      numthreads = 1

      if (myid .eq. 0) then
        write(*,'(3x,a,i8)') 'MPI tasks=', numtasks
        write(*,'(3x,a,i8)') 'OMP threads per task=', numthreads
        write(*,'(3x,a,i8)') 'No. of processors=', numthreads*numtasks
      endif
      
      !*******************************************************
      !******* Read input file bou.in by all processes********
      !*******************************************************
      open(unit=15,file='bou.in',status='old')
        read(15,301) dummy;read(15,301) dummy 
        read(15,*) n1,n2,n3,nsst
        read(15,301) dummy;read(15,301) dummy 
        read(15,*) ntst,trestart,tpin,tframe,tmax,idtv
        read(15,301) dummy;read(15,301) dummy 
        read(15,*) nread,ireset,pread
        read(15,301) dummy;read(15,301) dummy 
        read(15,*) alx3,rext1,rext2,istr3,str3,lmax
        read(15,301) dummy;read(15,301) dummy 
        read(15,*) ubcbot,ubctop
        read(15,301) dummy;read(15,301) dummy 
        read(15,*) Ra,Prs,invRo,alpha,gamma,Sm,alphaqs,betaqs,tau_cond
        read(15,301) dummy;read(15,301) dummy 
        read(15,*) A_stopmod,k_stopmod,A_sbotmod,k_sbotmod,dsaltop,dsalbot,qvaptop,qvapbot
        read(15,301) dummy;read(15,301) dummy 
        read(15,*) dtmax,resid,cflmax,dt,cfllim
        read(15,301) dummy;read(15,301) dummy 
       ! The legacy movie/control row may contain one, four, five, or eight
       ! values.  Columns 6--8 remain a compatibility fallback only.  New
       ! bou.in files use a separate Diagnostic output intervals section.
       ! Defaults are chosen for the multi-day t=1000 production workflow.
       tslice = tframe
       vortex_lc = 0.0d0
       vortex_q_z = 0.7d0
       tstat_diag = 0.5d0
       tscale_diag = 1.0d0
       tforce_diag = 5.0d0
       movie_values = -1.0d0
       read(15,'(A)') input_line
       read(input_line,*,iostat=read_ios) movie_values
       if(movie_values(1).ge.0.d0) mov_zcut_k=nint(movie_values(1))
       if(movie_values(2).gt.0.d0) tslice=movie_values(2)
       if(movie_values(4).ge.0.d0) vortex_lc=movie_values(4)
       if(movie_values(5).gt.0.d0) vortex_q_z=movie_values(5)
       if(movie_values(6).gt.0.d0) tstat_diag=movie_values(6)
       if(movie_values(7).gt.0.d0) tscale_diag=movie_values(7)
       if(movie_values(8).gt.0.d0) tforce_diag=movie_values(8)
       ! Optional explicit section following the movie row:
       ! !---------- Diagnostic output intervals --------------!
       ! TSTAT_DIAG  TSCALE_DIAG  TFORCE_DIAG
       ! 0.5d0       1.0d0        5.0d0
       read(15,'(A)',iostat=read_ios) input_line
       if(read_ios.eq.0) then
         if(index(input_line,'Diagnostic output intervals').gt.0) then
           read(15,'(A)',iostat=read_ios) input_line
           if(read_ios.ne.0) then
             if(myid.eq.0) write(*,*) &
     &         'Missing diagnostic interval column header in bou.in'
             call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
           endif
           diagnostic_values = -1.0d0
           read(15,*,iostat=read_ios) diagnostic_values
           if(read_ios.ne.0 .or. any(diagnostic_values.le.0.d0)) then
             if(myid.eq.0) write(*,*) &
     &         'Invalid diagnostic intervals in bou.in: require 3 positive values'
             call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
           endif
           tstat_diag=diagnostic_values(1)
           tscale_diag=diagnostic_values(2)
           tforce_diag=diagnostic_values(3)
         endif
       endif
       if(tstat_diag.le.0.d0) tstat_diag=0.5d0
       if(tscale_diag.le.0.d0) tscale_diag=1.0d0
       if(tforce_diag.le.0.d0) tforce_diag=5.0d0
       if(myid.eq.0) then
         write(*,'(3x,a,f8.3)') 'ordinary statistics interval = ',tstat_diag
         write(*,'(3x,a,f8.3)') 'scale diagnostics interval   = ',tscale_diag
         write(*,'(3x,a,f8.3)') 'force-balance interval       = ',tforce_diag
       endif
301     format(a4)                
      close(15)
      if(ubcbot.ne.0 .and. ubcbot.ne.1)then
        if(myid.eq.0) write(*,*) 'Invalid UBCBOT in bou.in: use 0=free-slip or 1=no-slip'
        call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
      endif
      if(ubctop.ne.0 .and. ubctop.ne.1)then
        if(myid.eq.0) write(*,*) 'Invalid UBCTOP in bou.in: use 0=free-slip or 1=no-slip'
        call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      !============================================
      !    basic parameters
      pi=2.d0*dasin(1.d0)                          
      tfini=dt*ntst                                                     
      tfin=tfini                                                        
      n1m=n1-1                                                          
      n2m=n2-1                                                          
      n3m=n3-1
      
      n1mr=n1m*mref1
      n2mr=n2m*mref2
      n3mr=n3m*mref3
      n1r=n1mr+1
      n2r=n2mr+1
      n3r=n3mr+1      

      usref1=1.d0/dble(mref1)
      usref2=1.d0/dble(mref2)
      usref3=1.d0/dble(mref3)

      !============================================
      !    allocate basic variables
      call initParams
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      !============================================
      !    coefficients in equations
      nu = dsqrt(Prs/Ra)
      kps = 1.d0 / dsqrt(Ra*Prs)
      byct = 1.0d0
      bycs = 1.0d0

      !============================================      
      !    boundary values for scalars
      sbctop = 1
      sbcbot = 1                    

      !============================================      
      !    open file headers
      call openfi

      !====================================================                                                                          
      !     assign coefficients for time marching schemes                     
      if(nsst.gt.1) then   
        gam(1)=8.d0/15.d0
        gam(2)=5.d0/12.d0
        gam(3)=3.d0/4.d0
        rom(1)=0.d0
        rom(2)=-17.d0/60.d0
        rom(3)=-5.d0/12.d0

        if(myid.eq.0) then
          write(*,100) (gam(n),n=1,nsst),(rom(n),n=1,nsst)
        endif
100     format(3x,'The time scheme is a III order Runge-Kutta',/,3x,'gam= ',3f8.3,4x,'ro= ',3f8.3)
      else                                                              
        gam(1)=1.5d0
        gam(2)=0.d0
        gam(3)=0.d0
        rom(1)=-0.5d0
        rom(2)=0.d0
        rom(3)=0.d0
        if(myid.eq.0) then
          write(*,110) gam(1),rom(1)
        endif
110     format(3x,'The time scheme is the Adams-Bashfort',/,3x,'gam= ',f8.3,4x,'ro= ',f8.3)                    
      endif    
                                                               
      do ns=1,nsst
        alm(ns)=(gam(ns)+rom(ns))
      enddo

      !====================================================                                                                          
      !     MPI start
      call mpi_workdistribution
      call mem_alloc
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      !m======================================================
      !     the solution of the problem starts
      !m======================================================

      ts=MPI_WTIME()
      call gcurv

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)      
      te=MPI_WTIME()
      
      if(myid.eq.0) then
        open(27,file="fact/Total_time.out")
        write(27,*)"Total simulation time in sec.: ", te-ts
        close(27)
      endif

      !============================================      
      !    close file headers
      call closefi

      !====================================================                                                                          
      !     MPI finalize
      call mem_dealloc
      call mgrd_mem_dealloc
      call MPI_FINALIZE(ierr)

      stop                                                              
      end
