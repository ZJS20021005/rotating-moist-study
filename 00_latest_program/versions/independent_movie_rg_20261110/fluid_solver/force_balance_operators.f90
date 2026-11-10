      module force_balance_operators
!     Shared DNS momentum-equation force operator.
!     The transient state is sampled at the ordinary-statistics cadence.
!     Full-field and vortex-conditioned force diagnostics can then be written
!     less frequently without turning F_T into a long-interval derivative.
      use param
      use local_arrays, only: q1,q2,q3,pr
      use mgrd_arrays, only: dsalc
      implicit none
      real, allocatable, save :: q1_prev(:,:,:),q2_prev(:,:,:),q3_prev(:,:,:)
      real, allocatable, save :: ft1(:,:,:),ft2(:,:,:),ft3(:,:,:)
      real, save :: force_time_prev=0.d0
      real, save :: force_dt=0.d0
      logical, save :: force_history_ready=.false.
      logical, save :: force_step_ready=.false.
      contains

      subroutine force_balance_sample
!     Save a backward-difference estimate of du/dt at the ordinary
!     statistics cadence.  The arrays are updated only after the estimate
!     has been formed, so force diagnostics at this same time use the
!     correct interval.
      if(.not.force_history_ready)then
        allocate(q1_prev(lbound(q1,1):ubound(q1,1), &
     &                   lbound(q1,2):ubound(q1,2), &
     &                   lbound(q1,3):ubound(q1,3)))
        allocate(q2_prev(lbound(q2,1):ubound(q2,1), &
     &                   lbound(q2,2):ubound(q2,2), &
     &                   lbound(q2,3):ubound(q2,3)))
        allocate(q3_prev(lbound(q3,1):ubound(q3,1), &
     &                   lbound(q3,2):ubound(q3,2), &
     &                   lbound(q3,3):ubound(q3,3)))
        allocate(ft1(lbound(q1,1):ubound(q1,1), &
     &               lbound(q1,2):ubound(q1,2), &
     &               lbound(q1,3):ubound(q1,3)))
        allocate(ft2(lbound(q2,1):ubound(q2,1), &
     &               lbound(q2,2):ubound(q2,2), &
     &               lbound(q2,3):ubound(q2,3)))
        allocate(ft3(lbound(q3,1):ubound(q3,1), &
     &               lbound(q3,2):ubound(q3,2), &
     &               lbound(q3,3):ubound(q3,3)))
        q1_prev=q1
        q2_prev=q2
        q3_prev=q3
        ft1=0.d0
        ft2=0.d0
        ft3=0.d0
        force_time_prev=time
        force_history_ready=.true.
        return
      endif
      force_dt=time-force_time_prev
      if(force_dt.le.0.d0)then
        q1_prev=q1
        q2_prev=q2
        q3_prev=q3
        force_time_prev=time
        force_step_ready=.false.
        return
      endif
      ft1=(q1-q1_prev)/force_dt
      ft2=(q2-q2_prev)/force_dt
      ft3=(q3-q3_prev)/force_dt
      q1_prev=q1
      q2_prev=q2
      q3_prev=q3
      force_time_prev=time
      force_step_ready=.true.
      return
      end subroutine force_balance_sample

      subroutine force_balance_commit
!     Kept as a no-op compatibility hook for older callers.
      return
      end subroutine force_balance_commit

      subroutine compute_force_at_cell(i,j,k,force)
      integer, intent(in) :: i,j,k
      real, intent(out) :: force(3,6)
      integer :: im,ip,jm,jp,km,kp
      real :: h11,h12,h13,h21,h22,h23,h31,h32,h33
      real :: uxc,uyc,amm,acc,app,d11,d22,d33

      im=imv(i)
      ip=ipv(i)
      jm=jmv(j)
      jp=jpv(j)
      km=kmv(k)
      kp=kpv(k)
      force=0.d0

!     F_T = du/dt from the most recent ordinary-statistics interval.
!     There is no additional damping term in this model.
      force(1,6)=ft1(i,j,k)
      force(2,6)=ft2(i,j,k)
      force(3,6)=ft3(i,j,k)

!     Conservative staggered-grid advection stencil used by the solver.
      h11=((q1(ip,j,k)+q1(i,j,k))**2 &
     &    -(q1(im,j,k)+q1(i,j,k))**2)*dx1*0.25d0
      h12=((q2(i,jp,k)+q2(im,jp,k))*(q1(i,jp,k)+q1(i,j,k)) &
     &    -(q2(i,j,k)+q2(im,j,k))*(q1(i,j,k)+q1(i,jm,k))) &
     &    *dx2*0.25d0
      if(k.eq.1)then
        h13=(q3(i,j,k+1)+q3(im,j,k+1)) &
     &      *(q1(i,j,k+1)+q1(i,j,k))*udx3m(k)*0.25d0
      elseif(k.eq.n3m)then
        h13=-(q3(i,j,k)+q3(im,j,k)) &
     &      *(q1(i,j,k)+q1(i,j,k-1))*udx3m(k)*0.25d0
      else
        h13=((q3(i,j,k+1)+q3(im,j,k+1))*(q1(i,j,k+1)+q1(i,j,k)) &
     &      -(q3(i,j,k)+q3(im,j,k))*(q1(i,j,k)+q1(i,j,k-1))) &
     &      *udx3m(k)*0.25d0
      endif
      force(1,1)=-(h11+h12+h13)
      uyc=(q2(i,j,k)+q2(im,j,k)+q2(i,jp,k)+q2(im,jp,k))*0.25d0
      force(1,2)=invRo*uyc

      h21=((q2(ip,j,k)+q2(i,j,k))*(q1(ip,j,k)+q1(ip,jm,k)) &
     &    -(q2(i,j,k)+q2(im,j,k))*(q1(i,j,k)+q1(i,jm,k))) &
     &    *dx1*0.25d0
      h22=((q2(i,jp,k)+q2(i,j,k))**2 &
     &    -(q2(i,jm,k)+q2(i,j,k))**2)*dx2*0.25d0
      if(k.eq.1)then
        h23=(q3(i,j,k+1)+q3(i,jm,k+1)) &
     &      *(q2(i,j,k+1)+q2(i,j,k))*udx3m(k)*0.25d0
      elseif(k.eq.n3m)then
        h23=-(q3(i,j,k)+q3(i,jm,k)) &
     &      *(q2(i,j,k)+q2(i,j,k-1))*udx3m(k)*0.25d0
      else
        h23=((q3(i,j,k+1)+q3(i,jm,k+1))*(q2(i,j,k+1)+q2(i,j,k)) &
     &      -(q3(i,j,k)+q3(i,jm,k))*(q2(i,j,k)+q2(i,j,k-1))) &
     &      *udx3m(k)*0.25d0
      endif
      force(2,1)=-(h21+h22+h23)
      uxc=(q1(i,j,k)+q1(ip,j,k)+q1(i,jm,k)+q1(ip,jm,k))*0.25d0
      force(2,2)=-invRo*uxc

      if(k.gt.1)then
        h31=((q1(ip,j,k)+q1(ip,j,k-1))*(q3(ip,j,k)+q3(i,j,k)) &
     &      -(q1(i,j,k)+q1(i,j,k-1))*(q3(i,j,k)+q3(im,j,k))) &
     &      *dx1*0.25d0
        h32=((q2(i,jp,k)+q2(i,jp,k-1))*(q3(i,jp,k)+q3(i,j,k)) &
     &      -(q2(i,j,k)+q2(i,j,k-1))*(q3(i,j,k)+q3(i,jm,k))) &
     &      *dx2*0.25d0
        if(k.eq.2)then
          h33=((q3(i,j,k+1)+q3(i,j,k))**2-q3(i,j,k)**2) &
     &        *udx3c(k)*0.25d0
        elseif(k.eq.n3m)then
          h33=(q3(i,j,k)**2-(q3(i,j,k)+q3(i,j,k-1))**2) &
     &        *udx3c(k)*0.25d0
        else
          h33=((q3(i,j,k+1)+q3(i,j,k))**2 &
     &        -(q3(i,j,k)+q3(i,j,k-1))**2)*udx3c(k)*0.25d0
        endif
        force(3,1)=-(h31+h32+h33)
      endif

!     Pressure includes the vertical derivative, which is essential for the
!     vertical residual after the horizontally averaged hydrostatic balance
!     has been removed.
      force(1,3)=-(pr(i,j,k)-pr(im,j,k))*dx1
      force(2,3)=-(pr(i,j,k)-pr(i,jm,k))*dx2
      if(k.gt.1) force(3,3)=-(pr(i,j,k)-pr(i,j,k-1))*udx3c(k)

!     Viscous terms use the solver's stretched-grid coefficients.
      amm=am3sk(k)
      acc=ac3sk(k)
      app=ap3sk(k)
      d11=(q1(ip,j,k)-2.d0*q1(i,j,k)+q1(im,j,k))*dx1q
      d22=(q1(i,jp,k)-2.d0*q1(i,j,k)+q1(i,jm,k))*dx2q
      if(k.eq.1)then
        d33=q1(i,j,kp)*app+q1(i,j,k)*acc
      elseif(k.eq.n3m)then
        d33=q1(i,j,k)*acc+q1(i,j,km)*amm
      else
        d33=q1(i,j,kp)*app+q1(i,j,k)*acc+q1(i,j,km)*amm
      endif
      force(1,4)=nu*(d11+d22+d33)

      d11=(q2(ip,j,k)-2.d0*q2(i,j,k)+q2(im,j,k))*dx1q
      d22=(q2(i,jp,k)-2.d0*q2(i,j,k)+q2(i,jm,k))*dx2q
      if(k.eq.1)then
        d33=q2(i,j,kp)*app+q2(i,j,k)*acc
      elseif(k.eq.n3m)then
        d33=q2(i,j,k)*acc+q2(i,j,km)*amm
      else
        d33=q2(i,j,kp)*app+q2(i,j,k)*acc+q2(i,j,km)*amm
      endif
      force(2,4)=nu*(d11+d22+d33)

      if(k.gt.1)then
        amm=am3ck(k)
        acc=ac3ck(k)
        app=ap3ck(k)
        d11=(q3(im,j,k)-2.d0*q3(i,j,k)+q3(ip,j,k))*dx1q
        d22=(q3(i,jm,k)-2.d0*q3(i,j,k)+q3(i,jp,k))*dx2q
        if(k.eq.2)then
          d33=q3(i,j,k+1)*app+q3(i,j,k)*acc
        elseif(k.eq.n3m)then
          d33=q3(i,j,k)*acc+q3(i,j,k-1)*amm
        else
          d33=q3(i,j,k+1)*app+q3(i,j,k)*acc+q3(i,j,k-1)*amm
        endif
        force(3,4)=nu*(d11+d22+d33)
        force(3,5)=0.5d0*(dsalc(i,j,k)+dsalc(i,j,k-1))
      endif
      return
      end subroutine compute_force_at_cell
      end module force_balance_operators
