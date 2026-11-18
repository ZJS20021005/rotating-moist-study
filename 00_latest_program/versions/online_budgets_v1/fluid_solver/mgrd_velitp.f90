      subroutine mgrd_velitp
      use param
      use local_arrays, only: q1,q2,q3
      use mgrd_arrays
      use mpi_param
      use mpih
      implicit none
       
      integer ic,jc,kc, ip,jp, icr,jcr,kcr
      integer ic0,icr0, jc0,jcr0, ipr,jpr

      real,dimension(4,4,4) :: qv3 
      real,dimension(4,4) :: qv2
      real,dimension(4) :: qv1

      real,dimension(1:m1mr,1:m2mr) :: dwloc,dwbot,dwall
      real,allocatable,dimension(:,:,:) :: tpdv, tpdvr, tpvr
      real,allocatable,dimension(:,:) :: q1yzc,q2xzc

      real udx1r, udx2r
      real dl1q, dl2q, dlf1, dlf2, lza

      !m=========================================================
      !  single mesh

      IF(mrefa.eq.1)then
        do kc=kstart,kend
         do jc=1,n2
         do ic=1,n1
           q1lr(ic,jc,kc) = q1(ic,jc,kc)
           q2lr(ic,jc,kc) = q2(ic,jc,kc)
           q3lr(ic,jc,kc) = q3(ic,jc,kc)
         enddo
         enddo
        enddo
      call update_both_ghosts(n1r,n2r,q1lr,kstartr,kendr)
      call update_both_ghosts(n1r,n2r,q2lr,kstartr,kendr)
      call update_both_ghosts(n1r,n2r,q3lr,kstartr,kendr)

      !m=========================================================
      !  Multiple mesh
      ELSE

      allocate(tpdv(-1:n1+1,-1:n2+1,kstart-2:kend+2))
      allocate(tpdvr(1:n1mr,1:n2mr,kstartr:kendr))
      allocate(tpvr(1:n1mr,1:n2mr,kstartr:kendr))
      allocate(q1yzc(-1:n2+1,kstart-2:kend+2))
      allocate(q2xzc(-1:n1+1,kstart-2:kend+2))

      tpdv = 0.d0
      tpdvr = 0.d0
      tpvr = 0.d0
      q1yzc = 0.d0
      q2xzc = 0.d0

      !m  interpolation of q1
      !  first interpolate dudx at cell center
      do kc=kstart-1,kend+1
       do jc=1,n2m
        do ic=1,n1m
         ip=ipv(ic)
         tpdv(ic,jc,kc)=(q1(ip,jc,kc)-q1(ic,jc,kc))*dx1
        enddo
       enddo
      enddo

      do kc=kstart-1,kend+1
!$OMP  PARALLEL DO
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(jc,ic,ip)
       do jc=1,n2m
        tpdv(-1,jc,kc)=tpdv(n1m-1,jc,kc)
        tpdv(0,jc,kc)=tpdv(n1m,jc,kc)
        tpdv(n1,jc,kc)=tpdv(1,jc,kc)
        tpdv(n1+1,jc,kc)=tpdv(2,jc,kc)
       enddo

       do ic=-1,n1+1
        tpdv(ic,-1,kc)=tpdv(ic,n2m-1,kc)
        tpdv(ic,0,kc)=tpdv(ic,n2m,kc)
        tpdv(ic,n2,kc)=tpdv(ic,1,kc)
        tpdv(ic,n2+1,kc)=tpdv(ic,2,kc)
       enddo
      enddo
      call update_itp_ghosts(n1,n2,tpdv,kstart,kend)

      do kc=kstart-1,kend
      do jc=0,n2m
      do ic=0,n1m
       qv3=tpdv(ic-1:ic+2,jc-1:jc+2,kc-1:kc+2)
       do kcr=max(krangs(kc),kstartr),min(krangs(kc+1)-1,kendr)
        qv2(:,:) = qv3(:,:,1)*czrs(1,kcr)+qv3(:,:,2)*czrs(2,kcr)&
     &             +qv3(:,:,3)*czrs(3,kcr)+qv3(:,:,4)*czrs(4,kcr)
        do jcr=max(jrangs(jc),1),min(jrangs(jc+1)-1,n2mr)
         qv1(:) = qv2(:,1)*cyrs(1,jcr)+qv2(:,2)*cyrs(2,jcr)&
     &           +qv2(:,3)*cyrs(3,jcr)+qv2(:,4)*cyrs(4,jcr)
         do icr=max(irangs(ic),1),min(irangs(ic+1)-1,n1mr)
          tpdvr(icr,jcr,kcr) = sum(qv1(:)*cxrs(:,icr))
         enddo
        enddo
       enddo
      enddo
      enddo
      enddo

      ! reconstruct q1lr from dudx
      ! interpolate q1 at an arbitrary y-z plane
      ic0 = 1
      icr0 = (ic0-1)*mref1+1
      do kc=kstart-1,kend+1
       do jc=1,n2m
        q1yzc(jc,kc)=q1(ic0,jc,kc)
       enddo
       q1yzc(0,kc)=q1yzc(n2m,kc)
       q1yzc(-1,kc)=q1yzc(n2m-1,kc)
       q1yzc(n2,kc)=q1yzc(1,kc)
       q1yzc(n2+1,kc)=q1yzc(2,kc)
      enddo
      call update_itp_ghosts_2d(n2,q1yzc,kstart,kend)
      do kc=kstart-1,kend
      do jc=0,n2m
       qv2=q1yzc(jc-1:jc+2,kc-1:kc+2)
       do kcr=max(krangs(kc),kstartr),min(krangs(kc+1)-1,kendr)
        qv1(:) = qv2(:,1)*czq1(1,kcr)+qv2(:,2)*czq1(2,kcr)&
     &          +qv2(:,3)*czq1(3,kcr)+qv2(:,4)*czq1(4,kcr)
        do jcr=max(jrangs(jc),1),min(jrangs(jc+1)-1,n2mr)
         tpvr(icr0,jcr,kcr) = sum(qv1(:)*cyq1(:,jcr)) 
        enddo
       enddo
      enddo
      enddo

      ! integrate along each x-line
      udx1r = 1.d0/dx1r
      do kcr=kstartr,kendr
       do jcr=1,n2mr
        icr = icr0
        ipr = ipvr(icr)
        do ic=1,n1mr-1
         tpvr(ipr,jcr,kcr) = tpvr(icr,jcr,kcr) &
     &           + tpdvr(icr,jcr,kcr)*udx1r
         icr = ipr
         ipr = ipvr(icr)
        enddo
       enddo
      enddo

      ! interpolate q1 at an arbitrary y-z plane
      ic0 = n1m/2 
      icr0 = (ic0-1)*mref1+1
      do kc=kstart-1,kend+1
       do jc=1,n2m
        q1yzc(jc,kc)=q1(ic0,jc,kc)
       enddo
       q1yzc(0,kc)=q1yzc(n2m,kc)
       q1yzc(-1,kc)=q1yzc(n2m-1,kc)
       q1yzc(n2,kc)=q1yzc(1,kc)
       q1yzc(n2+1,kc)=q1yzc(2,kc)
      enddo
      call update_itp_ghosts_2d(n2,q1yzc,kstart,kend)
      do kc=kstart-1,kend
      do jc=0,n2m
       qv2=q1yzc(jc-1:jc+2,kc-1:kc+2)
       do kcr=max(krangs(kc),kstartr),min(krangs(kc+1)-1,kendr)
        qv1(:) = qv2(:,1)*czq1(1,kcr)+qv2(:,2)*czq1(2,kcr)&
     &          +qv2(:,3)*czq1(3,kcr)+qv2(:,4)*czq1(4,kcr)
        do jcr=max(jrangs(jc),1),min(jrangs(jc+1)-1,n2mr)
         q1lr(icr0,jcr,kcr) = sum(qv1(:)*cyq1(:,jcr))
        enddo
       enddo
      enddo
      enddo

      udx1r = 1.d0/dx1r
      do kcr=kstartr,kendr
       do jcr=1,n2mr
        icr = icr0
        ipr = ipvr(icr)
        do ic=1,n1mr-1
         q1lr(ipr,jcr,kcr) = q1lr(icr,jcr,kcr)&
     &           + tpdvr(icr,jcr,kcr)*udx1r
         icr = ipr
         ipr = ipvr(icr)
        enddo
       enddo
      enddo

      ! average of two integrations
      do kcr=kstartr,kendr
       do jcr=1,n2mr
        do icr=1,n1mr
         q1lr(icr,jcr,kcr) = (q1lr(icr,jcr,kcr)+tpvr(icr,jcr,kcr))*0.5d0
        enddo
       enddo
      enddo

      !m interpolation of q2
      ! interpolate dvdy at cell center

      do kc=kstart-1,kend+1
       do jc=1,n2m
        jp=jpv(jc)
        do ic=1,n1m
         tpdv(ic,jc,kc)=(q2(ic,jp,kc)-q2(ic,jc,kc))*dx2
        enddo
       enddo
      enddo

      do kc=kstart-1,kend+1
       do jc=1,n2m
        tpdv(-1,jc,kc)=tpdv(n1m-1,jc,kc)
        tpdv(0,jc,kc)=tpdv(n1m,jc,kc)
        tpdv(n1,jc,kc)=tpdv(1,jc,kc)
        tpdv(n1+1,jc,kc)=tpdv(2,jc,kc)
       enddo

       do ic=-1,n1+1
        tpdv(ic,-1,kc)=tpdv(ic,n2m-1,kc)
        tpdv(ic,0,kc)=tpdv(ic,n2m,kc)
        tpdv(ic,n2,kc)=tpdv(ic,1,kc)
        tpdv(ic,n2+1,kc)=tpdv(ic,2,kc)
       enddo
      enddo
      call update_itp_ghosts(n1,n2,tpdv,kstart,kend)

      do kc=kstart-1,kend
      do jc=0,n2m
      do ic=0,n1m
       qv3=tpdv(ic-1:ic+2,jc-1:jc+2,kc-1:kc+2)
       do kcr=max(krangs(kc),kstartr),min(krangs(kc+1)-1,kendr)
        qv2(:,:) = qv3(:,:,1)*czrs(1,kcr)+qv3(:,:,2)*czrs(2,kcr)&
     &             +qv3(:,:,3)*czrs(3,kcr)+qv3(:,:,4)*czrs(4,kcr)
        do jcr=max(jrangs(jc),1),min(jrangs(jc+1)-1,n2mr)
         qv1(:) = qv2(:,1)*cyrs(1,jcr)+qv2(:,2)*cyrs(2,jcr)&
     &           +qv2(:,3)*cyrs(3,jcr)+qv2(:,4)*cyrs(4,jcr)
         do icr=max(irangs(ic),1),min(irangs(ic+1)-1,n1mr)
          tpdvr(icr,jcr,kcr) = sum(qv1(:)*cxrs(:,icr))
         enddo
        enddo
       enddo
      enddo
      enddo
      enddo

      ! construct q2lr from dvdy
      ! interpolate v and an arbitrary x-z plane
      jc0 = 1
      jcr0 = (jc0-1)*mref2+1
      do kc=kstart-1,kend+1
       do ic=1,n1m
        q2xzc(ic,kc)=q2(ic,jc0,kc)
       enddo
       q2xzc(0,kc)=q2xzc(n1m,kc)
       q2xzc(-1,kc)=q2xzc(n1m-1,kc)
       q2xzc(n1,kc)=q2xzc(1,kc)
       q2xzc(n1+1,kc)=q2xzc(2,kc)
      enddo
      call update_itp_ghosts_2d(n1,q2xzc,kstart,kend)
      do kc=kstart-1,kend
      do ic=0,n1m
       qv2=q2xzc(ic-1:ic+2,kc-1:kc+2)
       do kcr=max(krangs(kc),kstartr),min(krangs(kc+1)-1,kendr)
        qv1(:) = qv2(:,1)*czq2(1,kcr)+qv2(:,2)*czq2(2,kcr)&
     &          +qv2(:,3)*czq2(3,kcr)+qv2(:,4)*czq2(4,kcr)
        do icr=max(irangs(ic),1),min(irangs(ic+1)-1,n1mr)
         tpvr(icr,jcr0,kcr) = sum(qv1(:)*cxq2(:,icr))
        enddo
       enddo
      enddo
      enddo

      ! integrate along each y-line
      udx2r=1.0/dx2r
      do kcr=kstartr,kendr
       do icr=1,n1mr
        jcr = jcr0
        jpr = jpvr(jcr)
        do jc=1,n2mr-1
         tpvr(icr,jpr,kcr) = tpvr(icr,jcr,kcr) &
     &       + tpdvr(icr,jcr,kcr)*udx2r
         jcr = jpr
         jpr = jpvr(jcr)
        enddo
       enddo
      enddo

      ! interpolate v and an arbitrary x-z plane
      jc0 = n2m/2 
      jcr0 = (jc0-1)*mref2+1
      do kc=kstart-1,kend+1
       do ic=1,n1m
        q2xzc(ic,kc)=q2(ic,jc0,kc)
       enddo
       q2xzc(0,kc)=q2xzc(n1m,kc)
       q2xzc(-1,kc)=q2xzc(n1m-1,kc)
       q2xzc(n1,kc)=q2xzc(1,kc)
       q2xzc(n1+1,kc)=q2xzc(2,kc)
      enddo
      call update_itp_ghosts_2d(n1,q2xzc,kstart,kend)
      do kc=kstart-1,kend
      do ic=0,n1m
       qv2=q2xzc(ic-1:ic+2,kc-1:kc+2)
       do kcr=max(krangs(kc),kstartr),min(krangs(kc+1)-1,kendr)
        qv1(:) = qv2(:,1)*czq2(1,kcr)+qv2(:,2)*czq2(2,kcr)&
     &          +qv2(:,3)*czq2(3,kcr)+qv2(:,4)*czq2(4,kcr)
        do icr=max(irangs(ic),1),min(irangs(ic+1)-1,n1mr)
         q2lr(icr,jcr0,kcr) = sum(qv1(:)*cxq2(:,icr))
        enddo
       enddo
      enddo
      enddo
    
      ! integrate along each y-line
      udx2r=1.0/dx2r
      do kcr=kstartr,kendr
       do icr=1,n1mr
        jcr = jcr0
        jpr = jpvr(jcr)
        do jc=1,n2mr-1
         q2lr(icr,jpr,kcr) = q2lr(icr,jcr,kcr)&
     &       + tpdvr(icr,jcr,kcr)*udx2r
         jcr = jpr
         jpr = jpvr(jcr)
        enddo
       enddo
      enddo

      ! average of two integrations

      do kcr=kstartr,kendr
       do jcr=1,n2mr
        do icr=1,n1mr
         q2lr(icr,jcr,kcr) = (q2lr(icr,jcr,kcr)+tpvr(icr,jcr,kcr))*0.5d0
        enddo
       enddo
      enddo
 
      !   interpolation of q3
      ! interpolate dwdz at cell center after transpose

      do kc=kstart,kend
       do jc=1,n2m
        do ic=1,n1m
         tpdv(ic,jc,kc)=(q3(ic,jc,kc+1)-q3(ic,jc,kc))*udx3m(kc)
        enddo
       enddo
      enddo
      if(kstart.eq.1)then
       do jc=1,n2m
        do ic=1,n1m
         tpdv(ic,jc,0)=-(q1(ic+1,jc,0)-q1(ic,jc,0))*dx1&
     &                 -(q2(ic,jc+1,0)-q2(ic,jc,0))*dx2
        enddo
       enddo
      endif
      if(kend.eq.n3m)then
       do jc=1,n2m
        do ic=1,n1m
         tpdv(ic,jc,n3)=-(q1(ic+1,jc,n3)-q1(ic,jc,n3))*dx1&
     &                  -(q2(ic,jc+1,n3)-q2(ic,jc,n3))*dx2
        enddo
       enddo
      endif

      call update_both_ghosts_itp(-1,n1+1,-1,n2+1,kstart,kend,tpdv)
     
      do kc=kstart-1,kend+1
       do jc=1,n2m
        tpdv(0,jc,kc)=tpdv(n1m,jc,kc)
        tpdv(n1,jc,kc)=tpdv(1,jc,kc)
        tpdv(-1,jc,kc)=tpdv(n1m-1,jc,kc)
        tpdv(n1+1,jc,kc)=tpdv(2,jc,kc)
       enddo
       do ic=-1,n1+1
        tpdv(ic,0,kc)=tpdv(ic,n2m,kc)
        tpdv(ic,n2,kc)=tpdv(ic,1,kc)
        tpdv(ic,-1,kc)=tpdv(ic,n2m-1,kc)
        tpdv(ic,n2+1,kc)=tpdv(ic,2,kc)
       enddo
      enddo
      call update_itp_ghosts(n1,n2,tpdv,kstart,kend)

      do kc=kstart-1,kend
      do jc=0,n2m
      do ic=0,n1m
       qv3=tpdv(ic-1:ic+2,jc-1:jc+2,kc-1:kc+2)
       do kcr=max(krangs(kc),kstartr),min(krangs(kc+1)-1,kendr)
        qv2(:,:) = qv3(:,:,1)*czrs(1,kcr)+qv3(:,:,2)*czrs(2,kcr)&
     &             +qv3(:,:,3)*czrs(3,kcr)+qv3(:,:,4)*czrs(4,kcr)
        do jcr=max(jrangs(jc),1),min(jrangs(jc+1)-1,n2mr)
         qv1(:) = qv2(:,1)*cyrs(1,jcr)+qv2(:,2)*cyrs(2,jcr)&
     &           +qv2(:,3)*cyrs(3,jcr)+qv2(:,4)*cyrs(4,jcr)
         do icr=max(irangs(ic),1),min(irangs(ic+1)-1,n1mr)
          tpdvr(icr,jcr,kcr) = sum(qv1*cxrs(:,icr))
         enddo
        enddo
       enddo
      enddo
      enddo
      enddo

      q3lr(:,:,kstartr)=0.d0
      do kcr=kstartr,kendr
       lza = zcr(kcr+1)-zcr(kcr)
       do jcr=1,n2mr
       do icr=1,n1mr
        q3lr(icr,jcr,kcr+1) = q3lr(icr,jcr,kcr)+ tpdvr(icr,jcr,kcr)*lza
       enddo
       enddo
      enddo
      dwloc=q3lr(1:n1mr,1:n2mr,kendr+1)

      call MPI_SCAN(dwloc,dwbot,n1mr*n2mr,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myid.eq.numtasks-1)then
        dwall = dwbot
      endif
      call MPI_BCAST(dwall,n1mr*n2mr,MDP, numtasks-1,MPI_COMM_WORLD,ierr)

      do kcr=kstartr,kendr+1
       q3lr(1:n1mr,1:n2mr,kcr)=q3lr(1:n1mr,1:n2mr,kcr)&
     &     + dwbot(1:n1mr,1:n2mr)-dwloc(1:n1mr,1:n2mr)&
     &     - dwall(1:n1mr,1:n2mr)*0.5d0
      enddo

      if(allocated(tpdv)) deallocate(tpdv)
      if(allocated(tpdvr)) deallocate(tpdvr)
      if(allocated(tpvr)) deallocate(tpvr)
      if(allocated(q1yzc)) deallocate(q1yzc)
      if(allocated(q2xzc)) deallocate(q2xzc)

      !   periodic B.C.
      do kc=kstartr,kendr
        do jc=1,n2mr
          q1lr(n1r,jc,kc) = q1lr(1,jc,kc)
          q2lr(n1r,jc,kc) = q2lr(1,jc,kc)
          q3lr(n1r,jc,kc) = q3lr(1,jc,kc)
        enddo

        do ic=1,n1r
          q1lr(ic,n2r,kc) = q1lr(ic,1,kc)
          q2lr(ic,n2r,kc) = q2lr(ic,1,kc)
          q3lr(ic,n2r,kc) = q3lr(ic,1,kc)
        enddo
      enddo

      call update_both_ghosts(n1r,n2r,q1lr,kstartr,kendr)
      call update_both_ghosts(n1r,n2r,q2lr,kstartr,kendr)
      call update_both_ghosts(n1r,n2r,q3lr,kstartr,kendr)

      !   free slip condition for q1 and q2

      if(kstartr.eq.1)then
       if(ubcbot.eq.0)then
        dl1q = (zmr(1)-zcr(1))**2
        dl2q = (zmr(2)-zcr(1))**2
        dlf1 = dl2q/(dl2q-dl1q)
        dlf2 = dl1q/(dl2q-dl1q)
        do jc=1,n2r
         do ic=1,n1r
          q1lr(ic,jc,0)&
     &      = dlf1*q1lr(ic,jc,1)-dlf2*q1lr(ic,jc,2)
          q2lr(ic,jc,0)&
     &      = dlf1*q2lr(ic,jc,1)-dlf2*q2lr(ic,jc,2)
         enddo
        enddo
       else
        do jc=1,n2r
         do ic=1,n1r
          q1lr(ic,jc,0) = 0.d0
          q2lr(ic,jc,0) = 0.d0
         enddo
        enddo
       endif
      endif

      if(kendr.eq.n3mr)then
       if(ubctop.eq.0)then
        dl1q = (zmr(n3mr)-zcr(n3r))**2
        dl2q = (zmr(n3mr-1)-zcr(n3r))**2
        dlf1 = dl2q/(dl2q-dl1q)
        dlf2 = dl1q/(dl2q-dl1q)
        do jc=1,n2r
         do ic=1,n1r
          q1lr(ic,jc,n3r)&
     &      = dlf1*q1lr(ic,jc,n3mr) - dlf2*q1lr(ic,jc,n3mr-1)
          q2lr(ic,jc,n3r)&
     &      = dlf1*q2lr(ic,jc,n3mr) - dlf2*q2lr(ic,jc,n3mr-1)
         enddo
        enddo
       else
        do jc=1,n2r
         do ic=1,n1r
          q1lr(ic,jc,n3r) = 0.d0
          q2lr(ic,jc,n3r) = 0.d0
         enddo
        enddo
       endif
      endif

      ENDIF

      return
      end  subroutine mgrd_velitp

!====================================================
      subroutine update_itp_ghosts(n1,n2,q1,ks,ke)
      use mpih
      implicit none
      integer, intent(in) :: ks,ke
      real,intent(inout) :: q1(-1:n1+1,-1:n2+1,ks-2:ke+2)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag

      mydata= (n1+3)*(n2+3)

      my_down=myid-1
      my_up=myid+1

      if(myid .eq. 0) my_down=MPI_PROC_NULL
      if(myid .eq. numtasks-1) my_up=MPI_PROC_NULL

      tag=1
      call MPI_ISEND(q1(-1,-1,ke-1), mydata, MDP,my_up,tag,MPI_COMM_WORLD,req(1),ierr)
      call MPI_ISEND(q1(-1,-1,ks+1), mydata,  MDP,my_down,tag,MPI_COMM_WORLD,req(2), ierr)
      call MPI_IRECV(q1(-1,-1,ks-2), mydata,  MDP,my_down,tag,MPI_COMM_WORLD,req(3),ierr)
      call MPI_IRECV(q1(-1,-1,ke+2), mydata,  MDP,my_up, tag,MPI_COMM_WORLD,req(4),ierr)
      call MPI_Waitall(4,req,status,ierr)

      end subroutine update_itp_ghosts

!===============================================
      subroutine update_both_ghosts_itp(is,ie,js,je,ks,ke,var)
      use mpih
      implicit none
      integer, intent(in) :: ks,ke,is,ie,js,je
      real,intent(inout) :: var(is:ie,js:je,ks-2:ke+2)
      integer :: mydata
      integer :: my_down, my_up,tag

      mydata= (ie-is+1)*(je-js+1)

      my_down=myid-1
      my_up=myid+1

      if(myid .eq. 0) my_down=MPI_PROC_NULL
      if(myid .eq. numtasks-1) my_up=MPI_PROC_NULL

      tag=1
      call MPI_ISEND(var(is,js,ke), mydata, MDP,my_up,tag,MPI_COMM_WORLD,req(1),ierr)
      call MPI_ISEND(var(is,js,ks), mydata,  MDP,my_down,tag,MPI_COMM_WORLD,req(2), ierr)
      call MPI_IRECV(var(is,js,ks-1), mydata,  MDP,my_down,tag,MPI_COMM_WORLD,req(3),ierr)
      call MPI_IRECV(var(is,js,ke+1), mydata,  MDP, my_up, tag,MPI_COMM_WORLD,req(4),ierr)
      call MPI_Waitall(4,req,status,ierr)

      end subroutine update_both_ghosts_itp

!===============================================
      subroutine update_itp_ghosts_2d(n1,q1,ks,ke)
      use mpih
      implicit none
      integer,intent(in) :: n1,ks,ke
      real,intent(inout) :: q1(-1:n1+1,ks-2:ke+2)
      integer :: mydata
      integer :: my_down, my_up,tag

      mydata= n1+3

      my_down=myid-1
      my_up=myid+1

      if(myid .eq. 0) my_down=MPI_PROC_NULL
      if(myid .eq. numtasks-1) my_up=MPI_PROC_NULL

      tag=1
      call MPI_ISEND(q1(-1,ke-1), mydata, MDP,my_up,tag,MPI_COMM_WORLD,req(1),ierr)
      call MPI_ISEND(q1(-1,ks+1), mydata,  MDP,my_down,tag,MPI_COMM_WORLD,req(2), ierr)
      call MPI_IRECV(q1(-1,ks-2), mydata,  MDP,my_down,tag,MPI_COMM_WORLD,req(3),ierr)
      call MPI_IRECV(q1(-1,ke+2), mydata,  MDP,my_up, tag,MPI_COMM_WORLD,req(4),ierr)
      call MPI_Waitall(4,req,status,ierr)

      end subroutine update_itp_ghosts_2d
