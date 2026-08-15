!************************************************************************
!                                                                       *
! ****************************** subrout coetar  ********************** *
!                                                                       *
!    this subroutine calculates the coefficients for the                *
!    integration in the radial direction with non-uniform coor. trasf.  *
!                                                                       *
!************************************************************************
      subroutine coetarr
      use param
      use mpih
      implicit none
      integer :: km,kc,kp,k
      real:: a33,a33m,a33p

!  **coefficients for q1, q2 ap3sk,am3sk,ac3sk
!   ap3ssk,am3ssk,ac3ssk, psc   for x3 differentation
!  s means staggered that is at k+1/2 location
      do kc=2,n3mr-1
        kp=kc+1
        km=kc-1
        a33 = dx3qr/g3rmr(kc)
        a33p= +a33/g3rcr(kp)
        a33m= +a33/g3rcr(kc)
        ap3sskr(kc)=a33p
        am3sskr(kc)=a33m
        ac3sskr(kc)=-(a33p+a33m)
      enddo

!    lower wall  bound. conditions 
!    differemtiation of sec. derivative at 1+1/2
      kc=1
      kp=kc+1
      a33 = dx3qr/g3rmr(kc)
      a33p= +a33/g3rcr(kp)
      a33m= +a33/g3rcr(kc)
      ap3sskr(kc)=a33p
      am3sskr(kc)=2.0*a33m
      ac3sskr(kc)=-(a33p+2.0*a33m)
 
!    upper wall  bound. conditions
!    differentiation of sec. derivative at n3-1/2
      kc=n3mr
      kp=kc+1
      a33=dx3qr/g3rmr(kc)
      a33p= +a33/g3rcr(kp)
      a33m= +a33/g3rcr(kc)
      ap3sskr(kc)=2.0*a33p
      am3sskr(kc)=a33m
      ac3sskr(kc)=-(a33m+2.0*a33p)

      return
      end

