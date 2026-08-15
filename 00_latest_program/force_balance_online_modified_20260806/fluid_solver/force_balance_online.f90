      module force_balance_state
      implicit none
      real, allocatable :: q1_prev(:,:,:),q2_prev(:,:,:),q3_prev(:,:,:)
      real :: force_time_prev=0.d0
      logical :: force_history_ready=.false.
      end module force_balance_state

      subroutine force_balance_online
!     Online force-balance diagnostics using the DNS momentum stencils.
!
!     Quantity order: inertia, Coriolis, pressure, viscous, buoyancy,
!     and the finite-difference transient term.
!     Each force component is horizontally demeaned at every z before RMS.
!     Unit 116 writes bulk values over 0.1H <= z <= 0.9H:
!       time,
!       F_I F_C F_P F_V F_B F_T,
!       F_Ih F_Ch F_Ph F_Vh F_Bh F_Th,
!       F_Iz F_Cz F_Pz F_Vz F_Bz F_Tz,
!       R_G R_z R_momentum,
!       F_Ix F_Cx F_Px F_Vx F_Bx F_Tx,
!       F_Iy F_Cy F_Py F_Vy F_By F_Ty.
!     Unit 117 writes the same quantities at every z and time.
      use param
      use local_arrays, only: q1,q2,q3,pr
      use mgrd_arrays, only: dsalc
      use mpi_param, only: kstart,kend
      use mpih
      use force_balance_state
      implicit none

      integer, parameter :: nf=6
      integer :: ic,jc,kc,n
      real :: cellarea,dz,bulkdepth
      real :: f(3,nf),fp(3,nf),res(3)
      real :: my_area(n3m),area(n3m)
      real :: my_mean(3,nf,n3m),meanf(3,nf,n3m)
      real :: my_total_z(nf,n3m),total_z(nf,n3m),total(nf)
      real :: my_h_z(nf,n3m),h_z(nf,n3m),fh(nf)
      real :: my_x_z(nf,n3m),x_z(nf,n3m),fx(nf)
      real :: my_y_z(nf,n3m),y_z(nf,n3m),fy(nf)
      real :: my_z_z(nf,n3m),z_z(nf,n3m),fz(nf)
      real :: my_rg_z(n3m),rg_z(n3m),rg
      real :: my_rz_z(n3m),rz_z(n3m),rz
      real :: my_rmomentum_z(n3m),rmomentum_z(n3m),rmomentum
      real :: total_prof(nf),h_prof(nf),x_prof(nf),y_prof(nf)
      real :: z_prof(nf),rg_prof,rz_prof,rmomentum_prof
      real :: force_dt

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
        q1_prev=q1
        q2_prev=q2
        q3_prev=q3
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
        return
      endif

      my_area=0.d0
      my_mean=0.d0
      do kc=kstart,kend
        do jc=1,n2m
          do ic=1,n1m
            call force_at_cell(ic,jc,kc,f)
            cellarea=(xc(ic+1)-xc(ic))*(yc(jc+1)-yc(jc))
            my_area(kc)=my_area(kc)+cellarea
            my_mean(:,:,kc)=my_mean(:,:,kc)+f(:,:)*cellarea
          enddo
        enddo
      enddo

      call MPI_ALLREDUCE(my_area,area,n3m,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_mean,meanf,3*nf*n3m,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      do kc=1,n3m
        if(area(kc).gt.0.d0) meanf(:,:,kc)=meanf(:,:,kc)/area(kc)
      enddo

      my_total_z=0.d0
      my_h_z=0.d0
      my_x_z=0.d0
      my_y_z=0.d0
      my_z_z=0.d0
      my_rg_z=0.d0
      my_rz_z=0.d0
      my_rmomentum_z=0.d0
      do kc=kstart,kend
        if(zm(kc).lt.0.1d0*alx3 .or. zm(kc).gt.0.9d0*alx3) cycle
        do jc=1,n2m
          do ic=1,n1m
            call force_at_cell(ic,jc,kc,f)
            fp(:,:)=f(:,:)-meanf(:,:,kc)
            cellarea=(xc(ic+1)-xc(ic))*(yc(jc+1)-yc(jc))
            do n=1,nf
              my_total_z(n,kc)=my_total_z(n,kc)+sum(fp(:,n)**2)*cellarea
              my_h_z(n,kc)=my_h_z(n,kc)+(fp(1,n)**2+fp(2,n)**2)*cellarea
              my_x_z(n,kc)=my_x_z(n,kc)+fp(1,n)**2*cellarea
              my_y_z(n,kc)=my_y_z(n,kc)+fp(2,n)**2*cellarea
              my_z_z(n,kc)=my_z_z(n,kc)+fp(3,n)**2*cellarea
            enddo
            my_rg_z(kc)=my_rg_z(kc)+((fp(1,2)+fp(1,3))**2 &
     &                  +(fp(2,2)+fp(2,3))**2)*cellarea
            my_rz_z(kc)=my_rz_z(kc)+(fp(3,1)+fp(3,3)+fp(3,4) &
     &                  +fp(3,5))**2*cellarea
            res(:)=fp(:,6)-sum(fp(:,1:5),dim=2)
            my_rmomentum_z(kc)=my_rmomentum_z(kc)+sum(res(:)**2)*cellarea
          enddo
        enddo
      enddo

      call MPI_ALLREDUCE(my_total_z,total_z,nf*n3m,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_h_z,h_z,nf*n3m,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_x_z,x_z,nf*n3m,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_y_z,y_z,nf*n3m,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_z_z,z_z,nf*n3m,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_rg_z,rg_z,n3m,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_rz_z,rz_z,n3m,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_rmomentum_z,rmomentum_z,n3m,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)

      total=0.d0
      fh=0.d0
      fx=0.d0
      fy=0.d0
      fz=0.d0
      rg=0.d0
      rz=0.d0
      rmomentum=0.d0
      bulkdepth=0.d0
      do kc=1,n3m
        if(zm(kc).lt.0.1d0*alx3 .or. zm(kc).gt.0.9d0*alx3) cycle
        if(area(kc).le.0.d0) cycle
        dz=zc(kc+1)-zc(kc)
        bulkdepth=bulkdepth+dz
        do n=1,nf
          total(n)=total(n)+dsqrt(dmax1(0.d0,total_z(n,kc)/area(kc)))*dz
          fh(n)=fh(n)+dsqrt(dmax1(0.d0,h_z(n,kc)/area(kc)))*dz
          fx(n)=fx(n)+dsqrt(dmax1(0.d0,x_z(n,kc)/area(kc)))*dz
          fy(n)=fy(n)+dsqrt(dmax1(0.d0,y_z(n,kc)/area(kc)))*dz
          fz(n)=fz(n)+dsqrt(dmax1(0.d0,z_z(n,kc)/area(kc)))*dz
        enddo
        rg=rg+dsqrt(dmax1(0.d0,rg_z(kc)/area(kc)))*dz
        rz=rz+dsqrt(dmax1(0.d0,rz_z(kc)/area(kc)))*dz
        rmomentum=rmomentum+dsqrt(dmax1(0.d0,rmomentum_z(kc)/area(kc)))*dz
      enddo
      if(bulkdepth.gt.0.d0)then
        total=total/bulkdepth
        fh=fh/bulkdepth
        fx=fx/bulkdepth
        fy=fy/bulkdepth
        fz=fz/bulkdepth
        rg=rg/bulkdepth
        rz=rz/bulkdepth
        rmomentum=rmomentum/bulkdepth
      endif
      if(myid.eq.0)then
        write(116,510) time,total(:),fh(:),fz(:),rg,rz,rmomentum,fx(:),fy(:)
        do kc=1,n3m
          if(area(kc).le.0.d0) cycle
          total_prof(:)=dsqrt(dmax1(0.d0,total_z(:,kc)/area(kc)))
          h_prof(:)=dsqrt(dmax1(0.d0,h_z(:,kc)/area(kc)))
          x_prof(:)=dsqrt(dmax1(0.d0,x_z(:,kc)/area(kc)))
          y_prof(:)=dsqrt(dmax1(0.d0,y_z(:,kc)/area(kc)))
          z_prof(:)=dsqrt(dmax1(0.d0,z_z(:,kc)/area(kc)))
          rg_prof=dsqrt(dmax1(0.d0,rg_z(kc)/area(kc)))
          rz_prof=dsqrt(dmax1(0.d0,rz_z(kc)/area(kc)))
          rmomentum_prof=dsqrt(dmax1(0.d0,rmomentum_z(kc)/area(kc)))
          write(117,511) time,zm(kc),total_prof(:),h_prof(:), &
     &        z_prof(:),rg_prof,rz_prof,rmomentum_prof,x_prof(:),y_prof(:)
        enddo
      endif
      q1_prev=q1
      q2_prev=q2
      q3_prev=q3
      force_time_prev=time
 510  format(1x,f12.4,33(1x,ES20.8))
 511  format(1x,f12.4,34(1x,ES20.8))
      return

      contains

      subroutine force_at_cell(i,j,k,force)
      integer, intent(in) :: i,j,k
      real, intent(out) :: force(3,nf)
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
      force(1,6)=(q1(i,j,k)-q1_prev(i,j,k))/force_dt
      force(2,6)=(q2(i,j,k)-q2_prev(i,j,k))/force_dt
      force(3,6)=(q3(i,j,k)-q3_prev(i,j,k))/force_dt

!     Inertial and Coriolis terms at the staggered velocity locations.
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

!     Pressure gradients, including the vertical component.
      force(1,3)=-(pr(i,j,k)-pr(im,j,k))*dx1
      force(2,3)=-(pr(i,j,k)-pr(i,jm,k))*dx2
      if(k.gt.1) force(3,3)=-(pr(i,j,k)-pr(i,j,k-1))*udx3c(k)

!     Viscous terms use the same metric coefficients as invtrq1/2/3.
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
      end subroutine force_at_cell

      end subroutine force_balance_online
