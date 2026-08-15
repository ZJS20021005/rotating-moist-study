      subroutine nusse_walls
      use param
      use local_arrays, only: dens, dsal
      use mpih
      implicit none
      integer :: j,i
      real :: fcder,fcdern
      real :: surf,fac2
      real :: anutlow, anutupp,anuslow,anusupp
      real :: my_anusupp
      real :: adsal0,adsal1
      real :: del1,del1n
      real :: dsaltop_val,dsalbot_val

      anuslow = 0.d0
      anusupp = 0.d0

      my_anusupp = 0.d0

      if(myid.eq.0) then
        !  salinity flux at lower plate
        del1 = zmr(1)-zcr(1)
        surf = rext1*rext2
        fac2 = 1.d0/(dx1r*dx2r)
        adsal0 = 0.d0
        adsal1 = 0.d0
        do j=1,n2mr
          do i=1,n1mr
            dsalbot_val = dsalbot +  A_sbotmod*sin(2.0*pi*k_sbotmod*ymr(j))*sin(2.0*pi*k_sbotmod*xmr(i))
            adsal0 = adsal0 + dsalbot_val*fac2
            adsal1 = adsal1 + dsal(i,j,1)*fac2
          enddo
        end do
        adsal0 = adsal0 / surf
        adsal1 = adsal1 / surf
        anuslow = (adsal0-adsal1)/del1
      endif

      if(myid.eq.numtasks-1) then
        ! salinity flux at upper plate
        del1n = zcr(n3r) - zmr(n3mr)
        surf = rext1*rext2
        fac2 = 1.d0/(dx1r*dx2r)
        adsal0 = 0.d0
        adsal1 = 0.d0
        do j=1,n2mr
          do i=1,n1mr
            dsaltop_val = dsaltop +  A_stopmod*sin(2.0*pi*k_stopmod*ymr(j))
            adsal0 = adsal0 + dsaltop_val*fac2
            adsal1 = adsal1 + dsal(i,j,n3mr)*fac2
          enddo
        end do
        adsal0 = adsal0 / surf
        adsal1 = adsal1 / surf
        my_anusupp = (adsal1-adsal0)/del1n
      endif

      call MPI_REDUCE(my_anusupp,anusupp,1,MDP,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      if(myid.eq.0)then
        write(97,546) time,anuslow, anusupp
      endif
546   format(1x,f10.4,2(1x,ES20.8))

      return         
      end               
