!===========================================================
! Declaration of global variables
!***********************************************************      
      module param
        implicit none
        !===========================================================
        !      grid size
        integer,parameter :: mref1=1
        integer,parameter :: mref2=1
        integer,parameter :: mref3=1
        integer :: m1,m2,m3

        integer :: m2m,m3m,m2mh,m1m
        integer :: m1r,m2r,m3r
        integer :: m1mr,m2mr,m3mr,m2mhr

       !==============================================================
       !      inital condition
        integer,parameter :: tag_ini = 1
        integer,parameter :: Nitf = 1

       !==============================================================
       !      bou.in
        integer   :: n1, n2, n3, nsst, nwrit, nread
        integer   :: ntst, ireset
        real      :: tpin, trestart, tmax
        real      :: alx3, str3, rext1, rext2
        integer   :: istr3, lmax
        real      :: Ra,Prs,invRo,alpha,gamma,Sm,alphaqs,betaqs,tau_cond
        real      :: A_stopmod,k_stopmod,A_sbotmod,k_sbotmod,dsaltop,dsalbot,qvaptop,qvapbot
        integer   :: idtv
        real      :: dtmax, resid, cflmax, dt,dt_o, cfllim, dts
        real      :: tframe, tslice
        integer   :: mov_zcut_k
        real      :: vortex_lc

        integer   :: ubctop, ubcbot
        real      :: sbcbot,sbctop
       !==============================================================
       !      part.in
        integer   :: imlsfor,imlssca,pread
        real      :: sclf, objinitheight
        character*50  gtsfx

       !==============================================================
       !     simulation parameters
        real      :: time
        real :: dx2,dx3,dx1
        real :: dx2q,dx3q,dx1q

        real :: dx2r,dx3r,dx1r
        real :: dx2qr,dx3qr,dx1qr
    
        real, allocatable, dimension(:) :: xc,xm
        real, allocatable, dimension(:) :: yc,ym
        real, allocatable, dimension(:) :: zc,zm,g3rc,g3rm
    
        real, allocatable, dimension(:) :: xcr,xmr
        real, allocatable, dimension(:) :: ycr,ymr
        real, allocatable, dimension(:) :: zcr,zmr,g3rcr,g3rmr

       !==============================================================
       !     quantities for derivatives
        real, allocatable, dimension(:) :: udx3c,udx3m
        real, allocatable, dimension(:) :: udx3cr,udx3mr

       !==============================================================
       !     grid indices
        integer, allocatable, dimension(:) :: jmv,jpv
        integer, allocatable, dimension(:) :: imv,ipv
        integer, allocatable, dimension(:) :: jmhv
        integer, allocatable, dimension(:) :: kmc,kpc,kmv,kpv,kup,kum

        integer, allocatable, dimension(:) :: imvr,ipvr
        integer, allocatable, dimension(:) :: jmvr,jpvr
        integer, allocatable, dimension(:) :: kmvr,kpvr
  
       !==============================================================
       !     metric coefficients
        real, allocatable, dimension(:) :: ap3j,ac3j,am3j
        real, allocatable, dimension(:) :: ap3ck,ac3ck,am3ck
        real, allocatable, dimension(:) :: ap3sk,ac3sk,am3sk
        real, allocatable, dimension(:) :: ap3ssk,ac3ssk,am3ssk   
        real, allocatable, dimension(:) :: ap3sskr,ac3sskr,am3sskr

       !==============================================================
       !     variables for FFTW and Poisson solver
        real, dimension(13) :: ifx1
        real, allocatable, dimension(:) :: trigx1
        real, allocatable, dimension(:) :: ak2
        real, allocatable, dimension(:) :: ak1
        real, allocatable, dimension(:) :: amphk,acphk,apphk
        integer*8 :: fwd_plan,bck_plan
        
       !==============================================================
       !     other variables
        integer  :: n2m, n3m, n1m
        integer :: n1r,n2r,n3r,n1mr,n2mr,n3mr
        real :: byct, bycs, nu,kps
        real :: pi
        real :: al,ga,ro
        real :: beta
        integer :: ntime
        integer, parameter:: ndv=3
        real, dimension(1:ndv) :: vmax
        real, dimension(1:3) :: gam,rom,alm
        real :: usref1, usref2, usref3

      end module param
      

!===========================================================
! 3D arrays, dynamically allocated by each process
!***********************************************************
      module local_arrays
        use param
        implicit none
        real,allocatable,dimension(:,:,:) :: q1,q2,q3,dens,dsal,qvap,kpa,T,qs,cond_term
        real,allocatable,dimension(:,:,:) :: hro,hsal,hqvap,rhs,rhsr
        real,allocatable,dimension(:,:,:) :: ru1,ru2,ru3,ruro,rusal,ruqvap
        real,allocatable,dimension(:,:,:) :: pr,qcap,dph,dq
      end module local_arrays
    
      module mpih
        implicit none
        include 'mpif.h'
        integer :: myid, numtasks, numthreads, ierr
        integer, parameter :: master=0
        integer, parameter :: lvlhalo=3
        integer :: MDP = MPI_DOUBLE_PRECISION
        integer :: STATUS(MPI_STATUS_SIZE,4)
        integer :: req(1:4)
        integer(kind=MPI_OFFSET_KIND) :: disp, offset
      end module mpih
      
!===========================================================
! mpi param
!***********************************************************
      module mpi_param
        implicit none
        integer :: istart,iend, jstart,jend, kstart,kend
        integer :: jstartr,jendr, kstartr,kendr
        integer :: jstartp,jendp
        integer :: dj,dk,mydata,mydatam
        integer :: djp,djr,dkr
        integer, allocatable, dimension(:) :: offsetj,offsetk
        integer, allocatable, dimension(:) :: offsetjr,offsetkr
        integer, allocatable, dimension(:) :: offsetjp
        integer, allocatable, dimension(:) :: countj,countk
        integer, allocatable, dimension(:) :: countjr,countkr
        integer, allocatable, dimension(:) :: countjp
        integer, allocatable, dimension(:) :: countf
        integer(8), allocatable, dimension(:) :: offsetf 
      end module mpi_param

!===========================================================
! multiresolution
!***********************************************************
      module mgrd_arrays
        use param
        implicit none
        integer, parameter :: mrefa=mref1*mref2*mref3
        integer, allocatable, dimension(:) :: irangs,jrangs,krangs
        real, allocatable, dimension(:,:) :: cxq1, cxq2, cxq3, cxrs
        real, allocatable, dimension(:,:) :: cyq1, cyq2, cyq3, cyrs
        real, allocatable, dimension(:,:) :: czq1, czq2, czq3, czrs
        real, allocatable,dimension(:,:,:) :: q1lr,q2lr,q3lr,dsalc,qvapc
      end module mgrd_arrays

!***********************************************************
! statistics cbzhao
!***********************************************************
      module stat_arrays
       implicit none
       real,allocatable, dimension(:,:,:) :: vz_me,vz_rms
       real,allocatable, dimension(:,:,:) :: vy_me,vx_me,vy_rms,vx_rms
       real,allocatable, dimension(:,:,:) :: vz_fi_me,vz_fi_rms
       real,allocatable, dimension(:,:,:) :: vy_fi_me,vx_fi_me,vy_fi_rms,vx_fi_rms
       real,allocatable, dimension(:,:,:) :: temp_me,temp_rms,qtens
       real, allocatable,dimension(:,:,:) :: disste,dissth,tempvx_me
       integer :: nstatsamples
      end module stat_arrays

