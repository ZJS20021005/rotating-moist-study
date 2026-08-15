!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    FILE: initParams.f90                              !
!    CONTAINS: subroutine InitVariables                   !
!    PURPOSE: Initialization routine. Sets to zero all    !
!     variables used in the code                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initParams
use param
use mgrd_arrays

!===========================================================
! global variables
!***********************************************************   
  m1 = n1
  m2 = n2
  m3 = n3

  !-- var sizes
  m1m   = m1 - 1
  m2m   = m2 - 1
  m3m   = m3 - 1
  m1mr  = m1m * mref1
  m2mr  = m2m * mref2
  m3mr  = m3m * mref3
  m1r   = m1mr + 1
  m2r   = m2mr + 1
  m3r   = m3mr + 1
  m2mh  = m2m / 2 + 1
  m2mhr = m2mr / 2 + 1

  !==============================================================
  !     simulation parameters: geometry
  allocate( xc(1:m1), xm(1:m1) )
  allocate( yc(1:m2), ym(1:m2) )
  allocate( zc(1:m3), zm(1:m3), g3rc(1:m3), g3rm(1:m3) )

  allocate( xcr(1:m1r), xmr(1:m1r) )
  allocate( ycr(1:m2r), ymr(1:m2r) )
  allocate( zcr(1:m3r), zmr(1:m3r), g3rcr(1:m3r), g3rmr(1:m3r) )

  !==============================================================
  !     quantities for derivatives
  allocate( udx3c(1:m3),   udx3m(1:m3) )
  allocate( udx3cr(1:m3r), udx3mr(1:m3r) )

  !==============================================================
  !     grid indices
  allocate( imv(1:m1), ipv(1:m1) )
  allocate( jmv(1:m2), jpv(1:m2) )
  allocate( jmhv(1:m2+1) )
  allocate( kmc(1:m3), kpc(1:m3), kmv(1:m3), kpv(1:m3), kup(1:m3), kum(1:m3) )

  allocate( imvr(1:m1r), ipvr(1:m1r) ) 
  allocate( jmvr(1:m2r), jpvr(1:m2r) )
  allocate( kmvr(1:m3r), kpvr(1:m3r) )


  !==============================================================
  !     metric coefficients
  allocate( ap3j(1:m2),     ac3j(1:m2),     am3j(1:m2)     )
  allocate( ap3ck(1:m3),    ac3ck(1:m3),    am3ck(1:m3)    )
  allocate( ap3sk(1:m3),    ac3sk(1:m3),    am3sk(1:m3)    )
  allocate( ap3ssk(1:m3),   ac3ssk(1:m3),   am3ssk(1:m3)   )
  allocate( ap3sskr(1:m3r), ac3sskr(1:m3r), am3sskr(1:m3r) )

  !==============================================================
  !     Variables for FFTW and Poisson solver
  allocate( trigx1(3*m2/2+1) )
  allocate( ak2(1:m2))
  allocate( ak1(1:m1))
  allocate( amphk(1:m3), acphk(1:m3), apphk(1:m3) )

!===========================================================
! multiresolution
!***********************************************************
  allocate( irangs(0:m1),  jrangs(0:m2),  krangs(0:m3) )
  allocate( cxq1(4,0:m1r), cxq2(4,0:m1r), cxq3(4,0:m1r), cxrs(4,0:m1r) )
  allocate( cyq1(4,0:m2r), cyq2(4,0:m2r), cyq3(4,0:m2r), cyrs(4,0:m2r) )
  allocate( czq1(4,0:m3r), czq2(4,0:m3r), czq3(4,0:m3r), czrs(4,0:m3r) )

end
