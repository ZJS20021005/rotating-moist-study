
!**************************************************************
!
      subroutine cordin
      use param                                                 
      use mpih
      implicit none
      integer  :: j,k,i,l
      real, dimension(1:m3r) :: etaz
      real, dimension(0:m3r*2) :: etazm
      integer nclipr, n3mor
      real :: tstr3, z2dp
      real :: x1,x2,x3,delet

      !   x direction
      do i=1,n1r
        x1=dble(i-1)/dble(n1mr)
        xcr(i) = rext1*x1
      enddo
      do i=1,n1
        xc(i) = xcr((i-1)*mref1+1)
      enddo

      do i=1,n1m
        xm(i)=(xc(i)+xc(i+1))*0.5d0
      end do
      xm(n1) = 2.d0*xm(n1m)-xm(n1m-1)

      do i=1,n1mr
       xmr(i) = (xcr(i+1)+xcr(i))*0.5d0
      end do
      xmr(n1r) = 2.d0*xmr(n1mr)-xmr(n1mr-1)

      !  y direction
      do j=1,n2r
        x2=dble(j-1)/dble(n2mr)
        ycr(j) = rext2*x2
      enddo
      do j=1,n2
        yc(j) = ycr((j-1)*mref2+1)
      enddo

      do j=1,n2m
        ym(j)=(yc(j)+yc(j+1))*0.5d0
      end do
      ym(n2) = 2.d0*ym(n2m)-ym(n2m-1)

      do j=1,n2mr
       ymr(j) = (ycr(j+1)+ycr(j))*0.5d0
      end do
      ymr(n2r) = 2.d0*ymr(n2mr)-ymr(n2mr-1)

      !  z direction
         if (istr3.eq.0) then
           do k=1,n3
             x3=dble(k-1)/dble(n3m)
             etaz(k)=alx3*x3
             zc(k)=etaz(k)
           enddo
           do k=1,n3m
             do l=1,mref3
               zcr(mref3*(k-1)+l) = zc(k)+(zc(k+1)-zc(k))*dble(l-1)*usref3
             end do
           end do
           zcr(n3r)=zc(n3)
         endif
         !   chebyshev
         if(istr3.eq.1) then
           nclipr = str3*mref3
           n3mor = n3r+nclipr+nclipr
           do k=1,n3mor
             etazm(k) = dcos(pi*(dble(k)-0.5d0)/dble(n3mor))
           end do
           do k=1,n3r
             etaz(k)=etazm(k+nclipr)
           end do
           delet = etaz(1)-etaz(n3r)
           do k=1,n3r
             etaz(k)=etaz(k)/(0.50d0*delet)
           end do
           zcr(1) = 0.d0
           do k=2,n3mr
             zcr(k) = alx3*(1.d0-etaz(k))*0.5d0
           end do
           zcr(n3r) = alx3
           
           do k=1,n3
             zc(k) = zcr((k-1)*mref3+1)
           end do
         endif

      do k=1,n3m
        zm(k)=(zc(k)+zc(k+1))*0.5d0
      enddo
      zm(n3) = alx3
      do k=1,n3mr
       zmr(k) = (zcr(k+1)+zcr(k))*0.5d0
      end do
      zmr(n3r) = alx3

      !m===================================================
      !     METRIC QUANTITIES
      do k=1,n3m
        g3rm(k)=(zc(k+1)-zc(k))*dx3
      enddo
      do k=2,n3m
        g3rc(k)=(zc(k+1)-zc(k-1))*dx3*0.5d0
      enddo
      g3rc(1)=(zc(2)-zc(1))*dx3
      g3rc(n3)= (zc(n3)-zc(n3m))*dx3

      do k=1,n3m
        udx3m(k) = dx3/g3rm(k)
        udx3c(k) = dx3/g3rc(k)
      end do
      udx3c(n3) = dx3/g3rc(n3)


      do k=1,n3mr
        g3rmr(k)=(zcr(k+1)-zcr(k))*dx3r
      enddo
      do k=2,n3mr
        g3rcr(k)=(zcr(k+1)-zcr(k-1))*dx3r*0.5d0
      enddo
      g3rcr(1)=(zcr(2)-zcr(1))*dx3r
      g3rcr(n3r)= (zcr(n3r)-zcr(n3mr))*dx3r

      do k=1,n3mr
        udx3mr(k) = dx3r/g3rmr(k)
        udx3cr(k) = dx3r/g3rcr(k)
      end do
      udx3cr(n3r) = dx3r/g3rcr(n3r)

      !m===================================================
      !     METRIC QUANTITIES
      if(myid.eq.0) then
        open(unit=78,file='fact/xcorrefi.out',status='unknown')
        do k=1,n1r
          write(78,'(i5,2(2x,es16.8))') k,xcr(k),xmr(k)
        end do
        close(78)
        open(unit=78,file='fact/ycorrefi.out',status='unknown')
        do k=1,n2r
          write(78,'(i5,2(2x,f16.8))') k,ycr(k),ymr(k)
        end do
        close(78)
        open(unit=78,file='fact/zcorrefi.out',status='unknown')
        do k=1,n3r
          write(78,345) k,zcr(k),zmr(k),g3rcr(k),g3rmr(k)
        end do
        close(78)
345     format(i5,4(2x,es16.8))

        !     QUANTITIES FOR DERIVATIVES
        open(unit=78,file='fact/udx3r.out',status='unknown')
        do k=1,n3mr
          write(78,'(i5,2(2x,es16.8))') k,udx3mr(k),udx3cr(k)
        end do
        write(78,'(i5,2(2x,es16.8))') n3r,udx3mr(n3r),udx3cr(n3r)
        close(78)

        open(unit=78,file='fact/xcorbase.out',status='unknown')
        do k=1,n1
          write(78,'(i5,2(2x,es16.8))') k,xc(k),xm(k)
        end do
        close(78)
        open(unit=78,file='fact/ycorbase.out',status='unknown')
        do k=1,n2
          write(78,'(i5,2(2x,es16.8))') k,yc(k),ym(k)
        end do
        close(78)
        open(unit=78,file='fact/zcorbase.out',status='unknown')
        do k=1,n3
          write(78,345) k,zc(k),zm(k),g3rc(k),g3rm(k)
        end do
        close(78)

        !     QUANTITIES FOR DERIVATIVES
        open(unit=78,file='fact/udx3.out',status='unknown')
        do k=1,n3m
          write(78,'(i5,2(2x,es16.8))') k,udx3m(k),udx3c(k)
        end do
        write(78,'(i5,2(2x,es16.8))') n3,udx3m(n3m),udx3c(n3)
        close(78)

      endif

      return                                                            
      end                                                               
