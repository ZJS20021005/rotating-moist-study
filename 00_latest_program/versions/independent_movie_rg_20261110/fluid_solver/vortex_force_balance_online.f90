      subroutine vortex_force_balance_online
!     Force balance conditioned on the Q_2D>Q_rms vortex regions.
!
!     The regions are 2-D connected components on the configured z-cut.  The
!     reported forces are therefore conditional plane RMS values, not a claim
!     about a three-dimensional vortex volume.  A horizontal-plane mean is
!     removed from every force component before the conditional RMS, matching
!     the full-field diagnostic's removal of the height-wise mean.
      use param
      use local_arrays, only: q1,q2,q3,pr
      use mpi_param, only: kstart,kend
      use mpih
      use force_balance_operators, only: compute_force_at_cell,force_step_ready
      use dynamical_vortex_state
      implicit none
      integer, parameter :: nf=6
      integer :: ic,jc,kc,kq,c,n
      real :: f(3,nf),fp(3,nf),res(3)
      real :: dx,dy,cellarea,plane_area
      real :: my_plane_area,all_plane_area
      real :: my_mean(3,nf),meanf(3,nf)
      real, allocatable :: my_area(:),all_area(:)
      real, allocatable :: my_total(:,:,:),all_total(:,:,:)
      real, allocatable :: my_h(:,:,:),all_h(:,:,:)
      real, allocatable :: my_z(:,:,:),all_z(:,:,:)
      real, allocatable :: my_rg(:),all_rg(:),my_rz(:),all_rz(:)
      real, allocatable :: my_rm(:),all_rm(:)
      real :: total_rms(nf),h_rms(nf),z_rms(nf),rg,rz,rm
      real :: all_total_sum(nf),all_h_sum(nf),all_z_sum(nf)

      if(.not.q_ready) return
      if(.not.force_step_ready) return
      kq=q_k
      dx=rext1/dble(n1m)
      dy=rext2/dble(n2m)
      cellarea=dx*dy
      allocate(my_area(max(1,q_ncomp)),all_area(max(1,q_ncomp)))
      allocate(my_total(3,nf,max(1,q_ncomp)),all_total(3,nf,max(1,q_ncomp)))
      allocate(my_h(nf,max(1,q_ncomp),1),all_h(nf,max(1,q_ncomp),1))
      allocate(my_z(nf,max(1,q_ncomp),1),all_z(nf,max(1,q_ncomp),1))
      allocate(my_rg(max(1,q_ncomp)),all_rg(max(1,q_ncomp)))
      allocate(my_rz(max(1,q_ncomp)),all_rz(max(1,q_ncomp)))
      allocate(my_rm(max(1,q_ncomp)),all_rm(max(1,q_ncomp)))
      my_area=0.d0
      my_total=0.d0
      my_h=0.d0
      my_z=0.d0
      my_rg=0.d0
      my_rz=0.d0
      my_rm=0.d0
      my_mean=0.d0
      my_plane_area=0.d0

      if(kq.ge.kstart.and.kq.le.kend)then
        do jc=1,n2m
          do ic=1,n1m
            call compute_force_at_cell(ic,jc,kq,f)
            my_mean=my_mean+f
            my_plane_area=my_plane_area+cellarea
          enddo
        enddo
      endif
      call MPI_ALLREDUCE(my_mean,meanf,3*nf,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_plane_area,all_plane_area,1,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      if(all_plane_area.gt.0.d0)meanf=meanf/all_plane_area

      if(q_ncomp.gt.0.and.kq.ge.kstart.and.kq.le.kend)then
        do jc=1,n2m
          do ic=1,n1m
            c=q_label_core(jc,ic)
            if(c.le.0.or.c.gt.q_ncomp)cycle
            call compute_force_at_cell(ic,jc,kq,f)
            fp=f-meanf
            my_area(c)=my_area(c)+cellarea
            do n=1,nf
              my_total(1,n,c)=my_total(1,n,c)+fp(1,n)**2*cellarea
              my_total(2,n,c)=my_total(2,n,c)+fp(2,n)**2*cellarea
              my_total(3,n,c)=my_total(3,n,c)+fp(3,n)**2*cellarea
              my_h(n,c,1)=my_h(n,c,1)+(fp(1,n)**2+fp(2,n)**2)*cellarea
              my_z(n,c,1)=my_z(n,c,1)+fp(3,n)**2*cellarea
            enddo
            my_rg(c)=my_rg(c)+((fp(1,2)+fp(1,3))**2 &
     &                   +(fp(2,2)+fp(2,3))**2)*cellarea
            my_rz(c)=my_rz(c)+(fp(3,1)+fp(3,3)+fp(3,4) &
     &                   +fp(3,5))**2*cellarea
            res(:)=fp(:,6)-sum(fp(:,1:5),dim=2)
            my_rm(c)=my_rm(c)+sum(res(:)**2)*cellarea
          enddo
        enddo
      endif
      call MPI_ALLREDUCE(my_area,all_area,max(1,q_ncomp),MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_total,all_total,3*nf*max(1,q_ncomp),MDP, &
     &                   MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_h,all_h,nf*max(1,q_ncomp),MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_z,all_z,nf*max(1,q_ncomp),MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_rg,all_rg,max(1,q_ncomp),MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_rz,all_rz,max(1,q_ncomp),MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_rm,all_rm,max(1,q_ncomp),MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)

      if(myid.eq.0)then
        do c=1,q_ncomp
          if(all_area(c).le.0.d0)cycle
          do n=1,nf
            total_rms(n)=dsqrt(dmax1(0.d0, &
     &          (all_total(1,n,c)+all_total(2,n,c)+all_total(3,n,c))/all_area(c)))
            h_rms(n)=dsqrt(dmax1(0.d0,all_h(n,c,1)/all_area(c)))
            z_rms(n)=dsqrt(dmax1(0.d0,all_z(n,c,1)/all_area(c)))
          enddo
          rg=dsqrt(dmax1(0.d0,all_rg(c)/all_area(c)))
          rz=dsqrt(dmax1(0.d0,all_rz(c)/all_area(c)))
          rm=dsqrt(dmax1(0.d0,all_rm(c)/all_area(c)))
          write(122,610)time,q_id_core(c),q_status_core(c),all_area(c), &
     &       q_radius_core(c),total_rms(:),h_rms(:),z_rms(:),rg,rz,rm
        enddo

        all_total_sum=0.d0
        all_h_sum=0.d0
        all_z_sum=0.d0
        plane_area=0.d0
        rg=0.d0
        rz=0.d0
        rm=0.d0
        do c=1,q_ncomp
          if(all_area(c).le.0.d0)cycle
          plane_area=plane_area+all_area(c)
          do n=1,nf
            all_total_sum(n)=all_total_sum(n)+all_total(1,n,c) &
     &                         +all_total(2,n,c)+all_total(3,n,c)
            all_h_sum(n)=all_h_sum(n)+all_h(n,c,1)
            all_z_sum(n)=all_z_sum(n)+all_z(n,c,1)
          enddo
          rg=rg+all_rg(c)
          rz=rz+all_rz(c)
          rm=rm+all_rm(c)
        enddo
        total_rms=0.d0
        h_rms=0.d0
        z_rms=0.d0
        if(plane_area.gt.0.d0)then
          do n=1,nf
            total_rms(n)=dsqrt(dmax1(0.d0,all_total_sum(n)/plane_area))
            h_rms(n)=dsqrt(dmax1(0.d0,all_h_sum(n)/plane_area))
            z_rms(n)=dsqrt(dmax1(0.d0,all_z_sum(n)/plane_area))
          enddo
          rg=dsqrt(dmax1(0.d0,rg/plane_area))
          rz=dsqrt(dmax1(0.d0,rz/plane_area))
          rm=dsqrt(dmax1(0.d0,rm/plane_area))
        endif
        write(123,611)time,q_z,q_ncomp,plane_area,total_rms(:), &
     &       h_rms(:),z_rms(:),rg,rz,rm
      endif
      deallocate(my_area,all_area,my_total,all_total,my_h,all_h,my_z,all_z)
      deallocate(my_rg,all_rg,my_rz,all_rz,my_rm,all_rm)
610   format(1x,f12.4,2(1x,i12),2(1x,ES20.8),23(1x,ES20.8))
611   format(1x,f12.4,1x,ES20.8,1x,i8,1x,ES20.8,22(1x,ES20.8))
      return
      end
