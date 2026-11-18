      module dynamical_vortex_state
!     State shared by Q/omega_z identification and conditional force balance.
      use param
      implicit none
      integer, allocatable, save :: q_label_core(:,:),q_id_core(:)
      integer, allocatable, save :: q_status_core(:)
      real, allocatable, save :: q_area_core(:),q_radius_core(:)
      integer, save :: q_ncomp=0
      integer, save :: q_k=1
      real, save :: q_z=0.7d0
      logical, save :: q_ready=.false.

      logical, save :: q_initialized_pos=.false.,q_initialized_core=.false.
      logical, save :: q_load_attempted_pos=.false.,q_load_attempted_core=.false.
      integer, save :: q_next_id_pos=1,q_next_id_core=1
      integer, save :: q_nprev_pos=0,q_nprev_core=0
      integer, allocatable, save :: q_prev_label_pos(:,:),q_prev_label_core(:,:)
      integer, allocatable, save :: q_prev_id_pos(:),q_prev_id_core(:)
      end module dynamical_vortex_state

!=============================================================================
      subroutine dynamical_vortex_diagnostics
!     Identify dynamical vortices on the z=z(mov_zcut_k) horizontal plane.
!
!     The definition follows the supplied notebook:
!       Q_2D = 1/2 (||Omega_h||^2 - ||S_h||^2)
!       omega_z = dv/dx - du/dy
!     The notebook's 2-D Q cannot be called a full 3-D Q criterion because a
!     single z cut contains no vertical derivatives.  We therefore label all
!     outputs explicitly as Q_2D.  Two thresholds are retained:
!       Q_2D > 0 and Q_2D > Q_rms, both with periodic 8-connectivity.
!     The Q_rms set is the default mask passed to the vortex force diagnostic.
      use param
      use local_arrays, only: q1,q2
      use mpi_param, only: kstart,kend
      use mpih
      use dynamical_vortex_state
      implicit none
      integer :: ic,jc,ip,im,jp,jm,kq
      integer :: npos,ncore
      real :: ucell,vcell,dudx,dudy,dvdx,dvdy
      real :: s11,s22,s12,o12,snorm2,onorm2
      real :: q_rms,dx,dy,pi_local
      real :: my_q(n2m,n1m),q2d(n2m,n1m)
      real :: my_zeta(n2m,n1m),zeta2d(n2m,n1m)
      integer :: label_pos(n2m,n1m),label_core(n2m,n1m)
      integer :: id_pos(n1m*n2m),id_core(n1m*n2m)
      integer :: status_pos(n1m*n2m),status_core(n1m*n2m)
      real :: area_pos(n1m*n2m),radius_pos(n1m*n2m)
      real :: area_core(n1m*n2m),radius_core(n1m*n2m)
      real :: rmax_pos,rmax_core,total_area_pos,total_area_core

      if(.not.allocated(q_label_core))then
        allocate(q_label_core(n2m,n1m),q_id_core(n1m*n2m))
        allocate(q_status_core(n1m*n2m))
        allocate(q_area_core(n1m*n2m),q_radius_core(n1m*n2m))
        q_label_core=0
        q_id_core=0
        q_status_core=0
        q_area_core=0.d0
        q_radius_core=0.d0
      endif

      call update_both_ghosts(n1,n2,q1,kstart,kend)
      call update_both_ghosts(n1,n2,q2,kstart,kend)
      my_q=0.d0
      my_zeta=0.d0
      q_z=vortex_q_z
      q_k=1
      do ic=2,n3m
        if(dabs(zm(ic)-q_z).lt.dabs(zm(q_k)-q_z)) q_k=ic
      enddo
      kq=q_k
      q_z=zm(q_k)

      if(kq.ge.kstart .and. kq.le.kend)then
        dx=rext1/dble(n1m)
        dy=rext2/dble(n2m)
        do jc=1,n2m
          jp=jpv(jc)
          jm=jmv(jc)
          do ic=1,n1m
            ip=ipv(ic)
            im=imv(ic)
            ucell=0.5d0*(q1(ic,jc,kq)+q1(ip,jc,kq))
            vcell=0.5d0*(q2(ic,jc,kq)+q2(ic,jp,kq))
            dudx=0.5d0/dx*( &
     &        0.5d0*(q1(ip,jc,kq)+q1(ipv(ip),jc,kq)) &
     &       -0.5d0*(q1(im,jc,kq)+q1(ic,jc,kq)))
            dudy=0.5d0/dy*( &
     &        0.5d0*(q1(ic,jp,kq)+q1(ip,jp,kq)) &
     &       -0.5d0*(q1(ic,jm,kq)+q1(ip,jm,kq)))
            dvdx=0.5d0/dx*( &
     &        0.5d0*(q2(ip,jc,kq)+q2(ip,jp,kq)) &
     &       -0.5d0*(q2(im,jc,kq)+q2(im,jp,kq)))
            dvdy=0.5d0/dy*( &
     &        0.5d0*(q2(ic,jp,kq)+q2(ic,jpv(jp),kq)) &
     &       -0.5d0*(q2(ic,jm,kq)+q2(ic,jc,kq)))
            s11=dudx
            s22=dvdy
            s12=0.5d0*(dudy+dvdx)
            o12=0.5d0*(dudy-dvdx)
            snorm2=s11*s11+s22*s22+2.d0*s12*s12
            onorm2=2.d0*o12*o12
            my_q(jc,ic)=0.5d0*(onorm2-snorm2)
            my_zeta(jc,ic)=dvdx-dudy
          enddo
        enddo
      endif

      call MPI_ALLREDUCE(my_q,q2d,n1m*n2m,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(my_zeta,zeta2d,n1m*n2m,MDP,MPI_SUM, &
     &                   MPI_COMM_WORLD,ierr)
      dx=rext1/dble(n1m)
      dy=rext2/dble(n2m)
      pi_local=4.d0*datan(1.d0)
      q_rms=dsqrt(sum(q2d*q2d)/(dble(n1m)*dble(n2m)))

      if(myid.eq.0)then
        call identify_periodic_q_components(q2d,n1m,n2m,dx,dy,0.d0, &
     &      npos,label_pos,area_pos,radius_pos,rmax_pos,total_area_pos)
        call identify_periodic_q_components(q2d,n1m,n2m,dx,dy,q_rms, &
     &      ncore,label_core,area_core,radius_core,rmax_core, &
     &      total_area_core)
        if(vortex_lc.gt.0.d0)then
          call filter_q_components(n1m,n2m,2.d0*vortex_lc,npos,label_pos, &
     &       area_pos,radius_pos,rmax_pos,total_area_pos)
          call filter_q_components(n1m,n2m,2.d0*vortex_lc,ncore,label_core, &
     &       area_core,radius_core,rmax_core,total_area_core)
        endif
        call track_q_components(1,time,npos,label_pos,area_pos,radius_pos, &
     &       118,id_pos,status_pos)
        call track_q_components(2,time,ncore,label_core,area_core,radius_core, &
     &       119,id_core,status_core)
        write(120,510) time,zm(kq),q_rms,npos,rmax_pos,ncore,rmax_core
        q_ncomp=ncore
        q_label_core=label_core
        q_id_core=0
        q_status_core=0
        q_area_core=0.d0
        q_radius_core=0.d0
        do ic=1,ncore
          q_id_core(ic)=id_core(ic)
          q_status_core(ic)=status_core(ic)
          q_area_core(ic)=area_core(ic)
          q_radius_core(ic)=radius_core(ic)
        enddo
      endif

      call MPI_BCAST(q_ncomp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(q_label_core,n1m*n2m,MPI_INTEGER,0, &
     &               MPI_COMM_WORLD,ierr)
      call MPI_BCAST(q_id_core,n1m*n2m,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(q_status_core,n1m*n2m,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(q_area_core,n1m*n2m,MDP,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(q_radius_core,n1m*n2m,MDP,0,MPI_COMM_WORLD,ierr)
      q_ready=.true.

510   format(1x,f12.4,1x,ES20.8,1x,ES20.8,1x,i8,1x,ES20.8,1x,i8,1x,ES20.8)
      return
      end

!=============================================================================
      subroutine identify_periodic_q_components(field,nx,ny,dx,dy,threshold, &
     &     ncomp,label,area,radius,rmax,total_area)
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

      pi_local=4.d0*datan(1.d0)
      cell_area=dx*dy
      mask=field.gt.threshold
      visited=.false.
      label=0
      area=0.d0
      radius=0.d0
      ncomp=0
      rmax=0.d0
      total_area=0.d0
      do iy=1,ny
        do ix=1,nx
          if(mask(iy,ix).and..not.visited(iy,ix))then
            ncomp=ncomp+1
            head=1
            tail=1
            count=0
            queue_x(1)=ix
            queue_y(1)=iy
            visited(iy,ix)=.true.
            label(iy,ix)=ncomp
            do while(head.le.tail)
              jx=queue_x(head)
              jy=queue_y(head)
              head=head+1
              count=count+1
              do dj=-1,1
                do di=-1,1
                  if(di.ne.0.or.dj.ne.0)then
                    nx_neighbor=modulo(jx-1+di,nx)+1
                    ny_neighbor=modulo(jy-1+dj,ny)+1
                    if(mask(ny_neighbor,nx_neighbor).and. &
     &                 .not.visited(ny_neighbor,nx_neighbor))then
                      tail=tail+1
                      queue_x(tail)=nx_neighbor
                      queue_y(tail)=ny_neighbor
                      visited(ny_neighbor,nx_neighbor)=.true.
                      label(ny_neighbor,nx_neighbor)=ncomp
                    endif
                  endif
                enddo
              enddo
            enddo
            area(ncomp)=dble(count)*cell_area
            radius(ncomp)=dsqrt(area(ncomp)/pi_local)
            rmax=dmax1(rmax,radius(ncomp))
            total_area=total_area+area(ncomp)
          endif
        enddo
      enddo
      return
      end

!=============================================================================
      subroutine filter_q_components(nx,ny,min_radius,ncomp,label,area,radius, &
     &     rmax,total_area)
      implicit none
      integer, intent(in) :: nx,ny
      real, intent(in) :: min_radius
      integer, intent(inout) :: ncomp,label(ny,nx)
      real, intent(inout) :: area(nx*ny),radius(nx*ny),rmax,total_area
      integer :: old_to_new(nx*ny)
      integer :: ix,iy,c,new_n
      if(min_radius.le.0.d0)return
      old_to_new=0
      new_n=0
      rmax=0.d0
      total_area=0.d0
      do c=1,ncomp
        if(radius(c).gt.min_radius)then
          new_n=new_n+1
          old_to_new(c)=new_n
          area(new_n)=area(c)
          radius(new_n)=radius(c)
          rmax=dmax1(rmax,radius(new_n))
          total_area=total_area+area(new_n)
        endif
      enddo
      do c=new_n+1,ncomp
        area(c)=0.d0
        radius(c)=0.d0
      enddo
      do iy=1,ny
        do ix=1,nx
          c=label(iy,ix)
          if(c.gt.0)label(iy,ix)=old_to_new(c)
        enddo
      enddo
      ncomp=new_n
      return
      end

!=============================================================================
      subroutine track_q_components(kind,time_now,ncur,cur_label,cur_area, &
     &     cur_radius,out_unit,cur_id,cur_status)
      use param
      use dynamical_vortex_state
      implicit none
      integer, intent(in) :: kind,ncur,out_unit
      integer, intent(in) :: cur_label(n2m,n1m)
      integer, intent(out) :: cur_id(n1m*n2m),cur_status(n1m*n2m)
      real, intent(in) :: time_now,cur_area(n1m*n2m),cur_radius(n1m*n2m)
      logical :: initialized,load_attempted
      integer :: next_id,nprev,ix,iy,c,p,best_count,status
      integer :: best_child,best_child_count
      integer, allocatable :: prev_label(:,:),prev_id(:)
      integer, allocatable :: overlap(:,:),cur_prev_count(:),prev_cur_count(:)
      integer, allocatable :: cur_best_prev(:),prev_best_cur(:)
      logical, allocatable :: prev_used(:)

      if(.not.allocated(q_prev_label_pos))then
        allocate(q_prev_label_pos(n2m,n1m),q_prev_label_core(n2m,n1m))
        allocate(q_prev_id_pos(n1m*n2m),q_prev_id_core(n1m*n2m))
        q_prev_label_pos=0
        q_prev_label_core=0
        q_prev_id_pos=0
        q_prev_id_core=0
      endif
      if(kind.eq.1)then
        initialized=q_initialized_pos
        load_attempted=q_load_attempted_pos
        next_id=q_next_id_pos
        nprev=q_nprev_pos
        prev_label=q_prev_label_pos
        prev_id=q_prev_id_pos
        if(.not.load_attempted)then
          q_load_attempted_pos=.true.
          if(nread.ne.0.and.ireset.ne.1)call load_q_tracker(1, &
     &       initialized,next_id,nprev,prev_label,prev_id)
        endif
      else
        initialized=q_initialized_core
        load_attempted=q_load_attempted_core
        next_id=q_next_id_core
        nprev=q_nprev_core
        prev_label=q_prev_label_core
        prev_id=q_prev_id_core
        if(.not.load_attempted)then
          q_load_attempted_core=.true.
          if(nread.ne.0.and.ireset.ne.1)call load_q_tracker(2, &
     &       initialized,next_id,nprev,prev_label,prev_id)
        endif
      endif
      cur_id=0
      cur_status=4
      if(.not.initialized)then
        do c=1,ncur
          cur_id(c)=next_id
          cur_status(c)=1
          next_id=next_id+1
          write(out_unit,510)time_now,cur_id(c),cur_radius(c),cur_area(c),1
        enddo
      else
        allocate(overlap(max(1,nprev),max(1,ncur)))
        allocate(cur_prev_count(max(1,ncur)),prev_cur_count(max(1,nprev)))
        allocate(cur_best_prev(max(1,ncur)),prev_best_cur(max(1,nprev)))
        allocate(prev_used(max(1,nprev)))
        overlap=0
        cur_prev_count=0
        prev_cur_count=0
        cur_best_prev=0
        prev_best_cur=0
        prev_used=.false.
        do iy=1,n2m
          do ix=1,n1m
            c=cur_label(iy,ix)
            p=prev_label(iy,ix)
            if(c.gt.0.and.p.gt.0)overlap(p,c)=overlap(p,c)+1
          enddo
        enddo
        do c=1,ncur
          best_count=0
          do p=1,nprev
            if(overlap(p,c).gt.0)cur_prev_count(c)=cur_prev_count(c)+1
            if(overlap(p,c).gt.best_count)then
              best_count=overlap(p,c)
              cur_best_prev(c)=p
            endif
          enddo
        enddo
        do p=1,nprev
          best_child=0
          best_child_count=0
          do c=1,ncur
            if(overlap(p,c).gt.0)prev_cur_count(p)=prev_cur_count(p)+1
            if(overlap(p,c).gt.best_child_count)then
              best_child_count=overlap(p,c)
              best_child=c
            endif
          enddo
          prev_best_cur(p)=best_child
        enddo
        do c=1,ncur
          p=cur_best_prev(c)
          if(cur_prev_count(c).eq.0)then
            status=1
            cur_id(c)=next_id
            next_id=next_id+1
          else
            if(cur_prev_count(c).gt.1)then
              status=2
            elseif(prev_cur_count(p).gt.1)then
              status=3
            else
              status=0
            endif
            if(p.gt.0.and..not.prev_used(p).and. &
     &         (status.ne.3.or.prev_best_cur(p).eq.c))then
              cur_id(c)=prev_id(p)
              prev_used(p)=.true.
            else
              cur_id(c)=next_id
              next_id=next_id+1
            endif
          endif
          cur_status(c)=status
          write(out_unit,510)time_now,cur_id(c),cur_radius(c),cur_area(c),status
        enddo
        do p=1,nprev
          if(prev_cur_count(p).eq.0) &
     &      write(out_unit,510)time_now,prev_id(p),0.d0,0.d0,4
        enddo
        deallocate(overlap,cur_prev_count,prev_cur_count)
        deallocate(cur_best_prev,prev_best_cur,prev_used)
      endif
      prev_label=cur_label
      prev_id=0
      do c=1,ncur
        prev_id(c)=cur_id(c)
      enddo
      initialized=.true.
      if(kind.eq.1)then
        q_initialized_pos=initialized
        q_next_id_pos=next_id
        q_nprev_pos=ncur
        q_prev_label_pos=prev_label
        q_prev_id_pos=prev_id
        q_load_attempted_pos=.true.
      else
        q_initialized_core=initialized
        q_next_id_core=next_id
        q_nprev_core=ncur
        q_prev_label_core=prev_label
        q_prev_id_core=prev_id
        q_load_attempted_core=.true.
      endif
      call save_q_tracker(kind,initialized,next_id,ncur,prev_label,prev_id)
      deallocate(prev_label,prev_id)
510   format(1x,f12.4,1x,i12,2(1x,ES20.8),1x,i2)
      return
      end

!=============================================================================
      subroutine save_q_tracker(kind,initialized,next_id,nprev,prev_label,prev_id)
      use param
      implicit none
      integer, intent(in) :: kind,next_id,nprev
      logical, intent(in) :: initialized
      integer, intent(in) :: prev_label(n2m,n1m),prev_id(n1m*n2m)
      character(len=120) :: filename
      integer :: iunit,ios
      if(kind.eq.1)then
        filename='diagnostics/vortex_q/tracker_q_positive.restart'
      else
        filename='diagnostics/vortex_q/tracker_q_core.restart'
      endif
      open(newunit=iunit,file=filename,status='replace',form='unformatted', &
     &     access='stream',action='write',iostat=ios)
      if(ios.ne.0)return
      write(iunit)n1m,n2m,initialized,next_id,nprev
      write(iunit)prev_label
      write(iunit)prev_id
      close(iunit)
      return
      end

!=============================================================================
      subroutine load_q_tracker(kind,initialized,next_id,nprev,prev_label,prev_id)
      use param
      implicit none
      integer, intent(in) :: kind
      logical, intent(inout) :: initialized
      integer, intent(inout) :: next_id,nprev
      integer, intent(inout) :: prev_label(n2m,n1m),prev_id(n1m*n2m)
      character(len=120) :: filename
      logical :: exists,stored_initialized
      integer :: iunit,ios,stored_n1m,stored_n2m,stored_next_id,stored_nprev
      if(kind.eq.1)then
        filename='diagnostics/vortex_q/tracker_q_positive.restart'
      else
        filename='diagnostics/vortex_q/tracker_q_core.restart'
      endif
      inquire(file=filename,exist=exists)
      if(.not.exists)return
      open(newunit=iunit,file=filename,status='old',form='unformatted', &
     &     access='stream',action='read',iostat=ios)
      if(ios.ne.0)return
      read(iunit,iostat=ios)stored_n1m,stored_n2m,stored_initialized, &
     &     stored_next_id,stored_nprev
      if(ios.ne.0.or.stored_n1m.ne.n1m.or.stored_n2m.ne.n2m.or. &
     &   stored_nprev.lt.0.or.stored_nprev.gt.n1m*n2m)then
        close(iunit)
        return
      endif
      read(iunit,iostat=ios)prev_label
      if(ios.eq.0)read(iunit,iostat=ios)prev_id
      close(iunit)
      if(ios.ne.0)then
        prev_label=0
        prev_id=0
        return
      endif
      initialized=stored_initialized
      next_id=stored_next_id
      nprev=stored_nprev
      return
      end
