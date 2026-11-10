      subroutine force_balance_online(write_outputs)
!     Full-field online force-balance diagnostic.
!
!     The force operator is shared with the conditional vortex diagnostic in
!     force_balance_operators.f90.  This is important: both diagnostics use
!     the same staggered-grid stencil, including the vertical pressure
!     gradient and the DNS viscosity coefficients.
!
!     Force order: inertia, Coriolis, pressure, viscous, buoyancy, transient.
!     Unit 116 writes the bulk values over 0.1H <= z <= 0.9H.
!     Unit 117 writes the same quantities at every z.
      use param
      use local_arrays, only: q1,q2,q3,pr
      use mgrd_arrays, only: dsalc
      use mpi_param, only: kstart,kend
      use mpih
      use force_balance_operators, only: compute_force_at_cell, force_step_ready
      implicit none

      logical, intent(in) :: write_outputs
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

!     Synchronise the fields needed by the spatial force operators.  F_T has
!     already been sampled by force_balance_sample at the ordinary-statistics
!     cadence.
      call update_both_ghosts(n1,n2,q1,kstart,kend)
      call update_both_ghosts(n1,n2,q2,kstart,kend)
      call update_both_ghosts(n1,n2,q3,kstart,kend)
      call update_both_ghosts(n1,n2,pr,kstart,kend)
      call update_both_ghosts(n1,n2,dsalc,kstart,kend)
      geostrophic_ratio=-1.d0
      if(.not.force_step_ready) return

      my_area=0.d0
      my_mean=0.d0
      do kc=kstart,kend
        do jc=1,n2m
          do ic=1,n1m
            call compute_force_at_cell(ic,jc,kc,f)
            cellarea=(xc(ic+1)-xc(ic))*(yc(jc+1)-yc(jc))
            my_area(kc)=my_area(kc)+cellarea
            my_mean(:,:,kc)=my_mean(:,:,kc)+f(:,:)*cellarea
          enddo
        enddo
      enddo

      call MPI_ALLREDUCE(my_area,area,n3m,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_mean,meanf,3*nf*n3m,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
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
            call compute_force_at_cell(ic,jc,kc,f)
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

      call MPI_ALLREDUCE(my_total_z,total_z,nf*n3m,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_h_z,h_z,nf*n3m,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_x_z,x_z,nf*n3m,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_y_z,y_z,nf*n3m,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_z_z,z_z,nf*n3m,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_rg_z,rg_z,n3m,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_rz_z,rz_z,n3m,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_rmomentum_z,rmomentum_z,n3m,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)

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
      if(fh(3).gt.tiny(1.d0)) geostrophic_ratio=rg/fh(3)
      if(myid.eq.0.and.write_outputs)then
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
510   format(1x,f12.4,33(1x,ES20.8))
511   format(1x,f12.4,34(1x,ES20.8))
      return
      end
