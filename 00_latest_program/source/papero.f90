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
        vortex_lc = 0.0d0
        read(15,'(A)') input_line
        read(input_line,*,iostat=read_ios) mov_zcut_k,tframe_me_dummy, &
     &       stat_me_dummy,vortex_lc
        if(read_ios.ne.0)then
          vortex_lc = 0.0d0
          read(input_line,*) mov_zcut_k
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
