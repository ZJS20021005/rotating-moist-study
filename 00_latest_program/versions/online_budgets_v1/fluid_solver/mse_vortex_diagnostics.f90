      subroutine mse_vortex_diagnostics
!     MSE self-aggregation structure diagnostics.
!
!     This routine is intentionally separate from avgvar:
!       avgvar                 -> bulk averages, Nu profiles, velocity spectra
!       mse_vortex_diagnostics -> MSE anomaly masks and structure tracking
!
!     The two-dimensional identification field is the stretched-grid-weighted
!     vertical mean of the per-height horizontal MSE anomaly,
!
!       m' (x,y,z) = m(x,y,z)-<m>_xy(z),  m=b+gamma*qvap
!       m2d(x,y)   = integral m'(x,y,z) dz / integral dz.
!
!     Two masks are handled independently:
!       positive structures: m2d > 0
!       core structures    : m2d > sigma_m2d
!
!     Output files:
!       unit 105, diagnostics/vortex/vortex_positive_history.dat:
!         time, id, radius_positive, area_positive, status
!       unit 106, diagnostics/vortex/vortex_core_history.dat:
!         time, id, radius_core, area_core, status
!       unit 107, diagnostics/vortex/maximum_radius.dat:
!         time, Nv_positive, Rmax_positive, Nv_core, Rmax_core
!       unit 108, diagnostics/vortex/structure_ratio.dat:
!         time, Rmax_positive/Rmax_core,
!               total_area_positive/total_area_core
!       unit 114, diagnostics/vortex/vortex_kinematics.dat:
!         time, threshold kind, id, mean(u,v,w), mean/rms(zeta),
!               mean(mprime)
!
!     status: 0 alive, 1 birth, 2 merge, 3 breakup, 4 death
      use param
      use local_arrays, only: q1,q2,q3,dsal,qvap
      use mpi_param, only: kstart,kend,kstartr,kendr
      use mpih
      implicit none
      integer :: ic,jc,kc,ip,im,jp,jm,kp
      integer :: npos,ncore
      real :: zweight,my_total_zweight,total_zweight
      real :: mcell,mbar,sigma_m2d
      real :: ucell,vcell,wcell,zeta
      real :: dx,dy
      real :: rmax_pos,rmax_core,total_area_pos,total_area_core
      real :: ratio_radius,ratio_area
      real :: my_m2d(n2m,n1m),m2d(n2m,n1m)
      real :: my_u2d(n2m,n1m),u2d(n2m,n1m)
      real :: my_v2d(n2m,n1m),v2d(n2m,n1m)
      real :: my_w2d(n2m,n1m),w2d(n2m,n1m)
      real :: my_zeta2d(n2m,n1m),zeta2d(n2m,n1m)
      integer :: label_pos(n2m,n1m),label_core(n2m,n1m)
      integer :: id_pos(n1m*n2m),id_core(n1m*n2m)
      real :: area_pos(n1m*n2m),radius_pos(n1m*n2m)
      real :: area_core(n1m*n2m),radius_core(n1m*n2m)

      call update_both_ghosts(n1,n2,q1,kstart,kend)
      call update_both_ghosts(n1,n2,q2,kstart,kend)
      call update_both_ghosts(n1,n2,q3,kstart,kend)
      call update_both_ghosts(n1r,n2r,dsal,kstartr,kendr)
      call update_both_ghosts(n1r,n2r,qvap,kstartr,kendr)

      my_m2d = 0.0d0
      my_u2d = 0.0d0
      my_v2d = 0.0d0
      my_w2d = 0.0d0
      my_zeta2d = 0.0d0
      m2d = 0.0d0
      my_total_zweight = 0.0d0

      do kc=kstart,kend
        kp = kc+1
        zweight = zc(kc+1)-zc(kc)
        mbar = 0.0d0
        do jc=1,n2m
          do ic=1,n1m
            mbar = mbar+dsal(ic,jc,kc)+gamma*qvap(ic,jc,kc)
          enddo
        enddo
        mbar = mbar/(dble(n1m)*dble(n2m))

        do jc=1,n2m
          jp = jpv(jc)
          jm = jmv(jc)
          do ic=1,n1m
            ip = ipv(ic)
            im = imv(ic)
            mcell = dsal(ic,jc,kc)+gamma*qvap(ic,jc,kc)
            my_m2d(jc,ic) = my_m2d(jc,ic)+(mcell-mbar)*zweight
            ucell = 0.5d0*(q1(ic,jc,kc)+q1(ip,jc,kc))
            vcell = 0.5d0*(q2(ic,jc,kc)+q2(ic,jp,kc))
            wcell = 0.5d0*(q3(ic,jc,kc)+q3(ic,jc,kp))
            zeta = 0.25d0*dx1*(q2(ip,jc,kc)+q2(ip,jp,kc) &
     &                         -q2(im,jc,kc)-q2(im,jp,kc)) &
     &            -0.25d0*dx2*(q1(ic,jp,kc)+q1(ip,jp,kc) &
     &                         -q1(ic,jm,kc)-q1(ip,jm,kc))
            my_u2d(jc,ic) = my_u2d(jc,ic)+ucell*zweight
            my_v2d(jc,ic) = my_v2d(jc,ic)+vcell*zweight
            my_w2d(jc,ic) = my_w2d(jc,ic)+wcell*zweight
            my_zeta2d(jc,ic) = my_zeta2d(jc,ic)+zeta*zweight
          enddo
        enddo
        my_total_zweight = my_total_zweight+zweight
      enddo

      call MPI_ALLREDUCE(my_m2d,m2d,n1m*n2m,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_u2d,u2d,n1m*n2m,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_v2d,v2d,n1m*n2m,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_w2d,w2d,n1m*n2m,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_zeta2d,zeta2d,n1m*n2m,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_total_zweight,total_zweight,1,MDP, &
     &                   MPI_SUM,MPI_COMM_WORLD,ierr)

      if(myid.eq.0 .and. total_zweight.gt.0.0d0)then
        m2d = m2d/total_zweight
        u2d = u2d/total_zweight
        v2d = v2d/total_zweight
        w2d = w2d/total_zweight
        zeta2d = zeta2d/total_zweight
        m2d = m2d-sum(m2d)/(dble(n1m)*dble(n2m))
        sigma_m2d = dsqrt(sum(m2d*m2d)/(dble(n1m)*dble(n2m)))

        dx = rext1/dble(n1m)
        dy = rext2/dble(n2m)

        call identify_periodic_components(m2d,n1m,n2m,dx,dy,0.0d0, &
     &       npos,label_pos,area_pos,radius_pos,rmax_pos,total_area_pos)
        call identify_periodic_components(m2d,n1m,n2m,dx,dy,sigma_m2d, &
     &       ncore,label_core,area_core,radius_core,rmax_core, &
     &       total_area_core)

        if(vortex_lc.gt.0.0d0)then
          call filter_components_by_radius(n1m,n2m,2.0d0*vortex_lc, &
     &         npos,label_pos,area_pos,radius_pos,rmax_pos, &
     &         total_area_pos)
          call filter_components_by_radius(n1m,n2m,2.0d0*vortex_lc, &
     &         ncore,label_core,area_core,radius_core,rmax_core, &
     &         total_area_core)
        endif

        call track_mse_components(1,time,npos,label_pos,area_pos, &
     &       radius_pos,105,id_pos)
        call track_mse_components(2,time,ncore,label_core,area_core, &
     &       radius_core,106,id_core)

        call write_component_kinematics(time,1,npos,label_pos,id_pos, &
     &       u2d,v2d,w2d,zeta2d,m2d)
        call write_component_kinematics(time,2,ncore,label_core,id_core, &
     &       u2d,v2d,w2d,zeta2d,m2d)

        write(107,510) time,npos,rmax_pos,ncore,rmax_core

        ratio_radius = 0.0d0
        ratio_area = 0.0d0
        if(rmax_core.gt.0.0d0) ratio_radius = rmax_pos/rmax_core
        if(total_area_core.gt.0.0d0) ratio_area = &
     &     total_area_pos/total_area_core
        write(108,511) time,ratio_radius,ratio_area
      endif

510   format(1x,f10.4,1x,i8,1x,ES20.8,1x,i8,1x,ES20.8)
511   format(1x,f10.4,2(1x,ES20.8))
      return
      end

!=============================================================================
      subroutine filter_components_by_radius(nx,ny,min_radius,ncomp,label, &
     &     area,radius,rmax,total_area)
!     Drop components smaller than the linear-instability scale criterion.
      implicit none
      integer, intent(in) :: nx,ny
      real, intent(in) :: min_radius
      integer, intent(inout) :: ncomp,label(ny,nx)
      real, intent(inout) :: area(nx*ny),radius(nx*ny),rmax,total_area
      integer :: old_to_new(nx*ny)
      integer :: ix,iy,c,new_n

      if(min_radius.le.0.0d0) return

      old_to_new = 0
      new_n = 0
      rmax = 0.0d0
      total_area = 0.0d0

      do c=1,ncomp
        if(radius(c).gt.min_radius)then
          new_n = new_n+1
          old_to_new(c) = new_n
          area(new_n) = area(c)
          radius(new_n) = radius(c)
          rmax = dmax1(rmax,radius(new_n))
          total_area = total_area+area(new_n)
        endif
      enddo

      do c=new_n+1,ncomp
        area(c) = 0.0d0
        radius(c) = 0.0d0
      enddo

      do iy=1,ny
        do ix=1,nx
          c = label(iy,ix)
          if(c.gt.0)then
            label(iy,ix) = old_to_new(c)
          endif
        enddo
      enddo

      ncomp = new_n
      return
      end

!=============================================================================
      subroutine identify_periodic_components(field,nx,ny,dx,dy,threshold, &
     &     ncomp,label,area,radius,rmax,total_area)
!     Periodic 8-neighbour connected components of field > threshold.
      implicit none
      integer, intent(in) :: nx,ny
      real, intent(in) :: field(ny,nx),dx,dy,threshold
      integer, intent(out) :: ncomp,label(ny,nx)
      real, intent(out) :: area(nx*ny),radius(nx*ny),rmax,total_area
      logical :: mask(ny,nx),visited(ny,nx)
      integer :: queue_x(nx*ny),queue_y(nx*ny)
      integer :: ix,iy,jx,jy,nx_neighbor,ny_neighbor
      integer :: di,dj,head,tail,count
      real :: cell_area,pi_local

      pi_local = 4.0d0*datan(1.0d0)
      cell_area = dx*dy
      mask = field.gt.threshold
      visited = .false.
      label = 0
      area = 0.0d0
      radius = 0.0d0
      ncomp = 0
      rmax = 0.0d0
      total_area = 0.0d0

      do iy=1,ny
        do ix=1,nx
          if(mask(iy,ix) .and. .not.visited(iy,ix))then
            ncomp = ncomp+1
            head = 1
            tail = 1
            count = 0
            queue_x(1) = ix
            queue_y(1) = iy
            visited(iy,ix) = .true.
            label(iy,ix) = ncomp

            do while(head.le.tail)
              jx = queue_x(head)
              jy = queue_y(head)
              head = head+1
              count = count+1

              do dj=-1,1
                do di=-1,1
                  if(di.ne.0 .or. dj.ne.0)then
                    nx_neighbor = modulo(jx-1+di,nx)+1
                    ny_neighbor = modulo(jy-1+dj,ny)+1
                    if(mask(ny_neighbor,nx_neighbor) .and. &
     &                 .not.visited(ny_neighbor,nx_neighbor))then
                      tail = tail+1
                      queue_x(tail) = nx_neighbor
                      queue_y(tail) = ny_neighbor
                      visited(ny_neighbor,nx_neighbor) = .true.
                      label(ny_neighbor,nx_neighbor) = ncomp
                    endif
                  endif
                enddo
              enddo
            enddo

            area(ncomp) = dble(count)*cell_area
            radius(ncomp) = dsqrt(area(ncomp)/pi_local)
            rmax = dmax1(rmax,radius(ncomp))
            total_area = total_area+area(ncomp)
          endif
        enddo
      enddo

      return
      end

!=============================================================================
      subroutine track_mse_components(track_kind,time_now,ncur,cur_label, &
     &     cur_area,cur_radius,out_unit,cur_id)
!     Independent ID/status tracking for the two MSE thresholds.
!       track_kind=1: positive mask, m2d>0
!       track_kind=2: core mask,     m2d>sigma_m2d
      use param
      implicit none
      integer, intent(in) :: track_kind,ncur,out_unit
      integer, intent(in) :: cur_label(n2m,n1m)
      integer, intent(out) :: cur_id(n1m*n2m)
      real, intent(in) :: time_now,cur_area(n1m*n2m),cur_radius(n1m*n2m)

      logical, save :: initialized_pos=.false.,initialized_core=.false.
      logical, save :: load_attempted_pos=.false.
      logical, save :: load_attempted_core=.false.
      integer, save :: next_id_pos=1,next_id_core=1
      integer, save :: nprev_pos=0,nprev_core=0
      integer, allocatable, save :: prev_label_pos(:,:),prev_label_core(:,:)
      integer, allocatable, save :: prev_id_pos(:),prev_id_core(:)
      integer, allocatable :: prev_label(:,:),prev_id(:)
      integer, allocatable :: overlap(:,:),cur_prev_count(:)
      integer, allocatable :: prev_cur_count(:),cur_best_prev(:)
      integer, allocatable :: prev_best_cur(:)
      logical, allocatable :: prev_used(:)
      integer :: initialized_int,next_id,nprev
      integer :: ix,iy,c,p,best_count,status
      integer :: best_child,best_child_count

      if(.not.allocated(prev_label_pos))then
        allocate(prev_label_pos(n2m,n1m),prev_label_core(n2m,n1m))
        allocate(prev_id_pos(n1m*n2m),prev_id_core(n1m*n2m))
        prev_label_pos = 0
        prev_label_core = 0
        prev_id_pos = 0
        prev_id_core = 0
      endif

!     Recover IDs and masks across a continuation run. The state files are
!     overwritten after every diagnostic output, independently by threshold.
      if(track_kind.eq.1 .and. .not.load_attempted_pos)then
        load_attempted_pos = .true.
        if(nread.ne.0 .and. ireset.ne.1)then
          call load_tracker_state( &
     &      'diagnostics/vortex/tracker_positive.restart', &
     &      initialized_pos,next_id_pos,nprev_pos, &
     &      prev_label_pos,prev_id_pos)
        endif
      elseif(track_kind.eq.2 .and. .not.load_attempted_core)then
        load_attempted_core = .true.
        if(nread.ne.0 .and. ireset.ne.1)then
          call load_tracker_state( &
     &      'diagnostics/vortex/tracker_core.restart', &
     &      initialized_core,next_id_core,nprev_core, &
     &      prev_label_core,prev_id_core)
        endif
      endif

      allocate(prev_label(n2m,n1m),prev_id(n1m*n2m))
      cur_id = 0

      if(track_kind.eq.1)then
        initialized_int = merge(1,0,initialized_pos)
        next_id = next_id_pos
        nprev = nprev_pos
        prev_label = prev_label_pos
        prev_id = prev_id_pos
      else
        initialized_int = merge(1,0,initialized_core)
        next_id = next_id_core
        nprev = nprev_core
        prev_label = prev_label_core
        prev_id = prev_id_core
      endif

      if(initialized_int.eq.0)then
        do c=1,ncur
          cur_id(c) = next_id
          next_id = next_id+1
          write(out_unit,510) time_now,cur_id(c),cur_radius(c), &
     &       cur_area(c),1
        enddo
      else
        allocate(overlap(max(1,nprev),max(1,ncur)))
        allocate(cur_prev_count(max(1,ncur)),prev_cur_count(max(1,nprev)))
        allocate(cur_best_prev(max(1,ncur)),prev_best_cur(max(1,nprev)))
        allocate(prev_used(max(1,nprev)))
        overlap = 0
        cur_prev_count = 0
        prev_cur_count = 0
        cur_best_prev = 0
        prev_best_cur = 0
        prev_used = .false.
        cur_id = 0

        do iy=1,n2m
          do ix=1,n1m
            c = cur_label(iy,ix)
            p = prev_label(iy,ix)
            if(c.gt.0 .and. p.gt.0) overlap(p,c) = overlap(p,c)+1
          enddo
        enddo

        do c=1,ncur
          best_count = 0
          do p=1,nprev
            if(overlap(p,c).gt.0) cur_prev_count(c)=cur_prev_count(c)+1
            if(overlap(p,c).gt.best_count)then
              best_count = overlap(p,c)
              cur_best_prev(c) = p
            endif
          enddo
        enddo

        do p=1,nprev
          best_child = 0
          best_child_count = 0
          do c=1,ncur
            if(overlap(p,c).gt.0) prev_cur_count(p)=prev_cur_count(p)+1
            if(overlap(p,c).gt.best_child_count)then
              best_child_count = overlap(p,c)
              best_child = c
            endif
          enddo
          prev_best_cur(p) = best_child
        enddo

        do c=1,ncur
          p = cur_best_prev(c)
          if(cur_prev_count(c).eq.0)then
            status = 1
            cur_id(c) = next_id
            next_id = next_id+1
          else
            if(cur_prev_count(c).gt.1)then
              status = 2
            elseif(prev_cur_count(p).gt.1)then
              status = 3
            else
              status = 0
            endif

            if(p.gt.0 .and. .not.prev_used(p) .and. &
     &         (status.ne.3 .or. prev_best_cur(p).eq.c))then
              cur_id(c) = prev_id(p)
              prev_used(p) = .true.
            else
              cur_id(c) = next_id
              next_id = next_id+1
            endif
          endif

          write(out_unit,510) time_now,cur_id(c),cur_radius(c), &
     &       cur_area(c),status
        enddo

        do p=1,nprev
          if(prev_cur_count(p).eq.0)then
            write(out_unit,510) time_now,prev_id(p),0.0d0,0.0d0,4
          endif
        enddo

        deallocate(overlap,cur_prev_count,prev_cur_count)
        deallocate(cur_best_prev,prev_best_cur,prev_used)
      endif

      prev_label = cur_label
      prev_id = 0
      do c=1,ncur
        prev_id(c) = cur_id(c)
      enddo

      if(track_kind.eq.1)then
        initialized_pos = .true.
        next_id_pos = next_id
        nprev_pos = ncur
        prev_label_pos = prev_label
        prev_id_pos = prev_id
      else
        initialized_core = .true.
        next_id_core = next_id
        nprev_core = ncur
        prev_label_core = prev_label
        prev_id_core = prev_id
      endif

      if(track_kind.eq.1)then
        call save_tracker_state( &
     &    'diagnostics/vortex/tracker_positive.restart', &
     &    initialized_pos,next_id_pos,nprev_pos, &
     &    prev_label_pos,prev_id_pos)
      else
        call save_tracker_state( &
     &    'diagnostics/vortex/tracker_core.restart', &
     &    initialized_core,next_id_core,nprev_core, &
     &    prev_label_core,prev_id_core)
      endif

      deallocate(prev_label,prev_id)

510   format(1x,f10.4,1x,i12,2(1x,ES20.8),1x,i2)
      return
      end

!=============================================================================
      subroutine write_component_kinematics(time_now,threshold_kind,ncomp, &
     &     label,component_id,u2d,v2d,w2d,zeta2d,mprime2d)
!     Area averages needed for later structure-based Rossby diagnostics.
!     Output unit 114:
!       time, threshold_kind, id, mean_u, mean_v, mean_w,
!       mean_zeta, rms_zeta, mean_mprime
      use param
      implicit none
      integer, intent(in) :: threshold_kind,ncomp
      integer, intent(in) :: label(n2m,n1m),component_id(n1m*n2m)
      real, intent(in) :: time_now
      real, intent(in) :: u2d(n2m,n1m),v2d(n2m,n1m)
      real, intent(in) :: w2d(n2m,n1m),zeta2d(n2m,n1m)
      real, intent(in) :: mprime2d(n2m,n1m)
      integer :: ic,jc,c,count
      real :: su,sv,sw,sz,sz2,sum_mprime,rinv

      do c=1,ncomp
        count = 0
        su = 0.0d0
        sv = 0.0d0
        sw = 0.0d0
        sz = 0.0d0
        sz2 = 0.0d0
        sum_mprime = 0.0d0
        do jc=1,n2m
          do ic=1,n1m
            if(label(jc,ic).eq.c)then
              count = count+1
              su = su+u2d(jc,ic)
              sv = sv+v2d(jc,ic)
              sw = sw+w2d(jc,ic)
              sz = sz+zeta2d(jc,ic)
              sz2 = sz2+zeta2d(jc,ic)**2
              sum_mprime = sum_mprime+mprime2d(jc,ic)
            endif
          enddo
        enddo
        if(count.gt.0)then
          rinv = 1.0d0/dble(count)
          write(114,510) time_now,threshold_kind,component_id(c), &
     &      su*rinv,sv*rinv,sw*rinv,sz*rinv,dsqrt(sz2*rinv), &
     &      sum_mprime*rinv
        endif
      enddo

510   format(1x,f10.4,2(1x,i12),6(1x,ES20.8))
      return
      end

!=============================================================================
      subroutine save_tracker_state(filename,initialized,next_id,nprev, &
     &                              prev_label,prev_id)
      use param
      implicit none
      character(*), intent(in) :: filename
      logical, intent(in) :: initialized
      integer, intent(in) :: next_id,nprev
      integer, intent(in) :: prev_label(n2m,n1m)
      integer, intent(in) :: prev_id(n1m*n2m)
      integer :: iunit,ios

      open(newunit=iunit,file=filename,status='replace', &
     &     form='unformatted',access='stream',action='write',iostat=ios)
      if(ios.ne.0) return
      write(iunit) n1m,n2m,initialized,next_id,nprev
      write(iunit) prev_label
      write(iunit) prev_id
      close(iunit)
      return
      end

!=============================================================================
      subroutine load_tracker_state(filename,initialized,next_id,nprev, &
     &                              prev_label,prev_id)
      use param
      implicit none
      character(*), intent(in) :: filename
      logical, intent(inout) :: initialized
      integer, intent(inout) :: next_id,nprev
      integer, intent(inout) :: prev_label(n2m,n1m)
      integer, intent(inout) :: prev_id(n1m*n2m)
      logical :: exists,stored_initialized
      integer :: iunit,ios,stored_n1m,stored_n2m
      integer :: stored_next_id,stored_nprev

      inquire(file=filename,exist=exists)
      if(.not.exists) return
      open(newunit=iunit,file=filename,status='old', &
     &     form='unformatted',access='stream',action='read',iostat=ios)
      if(ios.ne.0) return
      read(iunit,iostat=ios) stored_n1m,stored_n2m,stored_initialized, &
     &                     stored_next_id,stored_nprev
      if(ios.ne.0 .or. stored_n1m.ne.n1m .or. stored_n2m.ne.n2m &
     &   .or. stored_nprev.lt.0 .or. stored_nprev.gt.n1m*n2m)then
        close(iunit)
        return
      endif
      read(iunit,iostat=ios) prev_label
      if(ios.eq.0) read(iunit,iostat=ios) prev_id
      close(iunit)
      if(ios.ne.0)then
        prev_label = 0
        prev_id = 0
        return
      endif
      initialized = stored_initialized
      next_id = stored_next_id
      nprev = stored_nprev
      return
      end
