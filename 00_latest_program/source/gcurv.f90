      subroutine gcurv
      use mpih  
      use mpi_param,only: kstart,kend,kstartr,kendr
      use param
      use local_arrays
      use mgrd_arrays

      implicit none
      integer :: ntstf,l
      real    :: cflm,cflmr
      real    :: dmaxc,dmaxr
      real    :: aengc,aengr
      real    :: ti(2), tin(3)
      character*5 ipfi
      real    :: tprfi
      integer :: inp,i

      tin(1) = MPI_WTIME()
       
      !========================================================
      call initia
      dt_o=dt
      call meshes
      call cordin
      call coetar
      call coetarr
      call indic
      call indicr 

      !========================================================
      if(myid.eq.0) then
        write(*,754)n1,n2,n3
        write(*,756)mref1,mref2,mref3,lmax
      endif
754   format(3x,'grid resolution: ','   n1=',i5,' n2=',i5,'   n3=',i5)
756   format(3x,'refine scales:  ','    mref1=',i2,'    mref2=',i2,' mref3=',i2,'    lmax=',i2)
      call write_grid_info
      call phini
      do l=1,ndv                                                   
        vmax(l)=0.d0
      enddo

      !==============================================================
      !     create the initial conditions
      if(nread.eq.0) then
        if(myid.eq.0) then
          write(*,'(6x,a,/)')'nread=0 ---> new initialization'
        endif

        ntime=0                                                         
        time=0.d0
        cflm=0.d0
        cflmr=0.d0
        call inqpr
      else
        if(myid.eq.0) then
          write(*,'(6x,a,/)')'nread=1 ---> reading files'
        endif
        call inirea
        tmax = tmax + time
      endif

      call caldiffusivity
      call update_both_ghosts(n1,n2,q1,kstart,kend)
      call update_both_ghosts(n1,n2,q2,kstart,kend)
      call update_both_ghosts(n1,n2,q3,kstart,kend)
      call update_both_ghosts(n1,n2,dens,kstart,kend)
      call update_both_ghosts(n1r,n2r,dsal,kstartr,kendr)
      call update_both_ghosts(n1r,n2r,qvap,kstartr,kendr)
      call update_both_ghosts(n1,n2,pr,kstart,kend)
      call update_both_ghosts(n1mr,n2mr,kpa,kstartr,kendr)
      !============================================================
      !     the boundary conditions
      call velbc

      !=============================================================
      !    for multi grid 
      call mgrd_mem_alloc
      if(myid.eq.0)write(*,'(/,3x,a,/)')'multigrid mem allocated'
      call mgrd_idc
      call mgrd_velitp
      call mgrd_dsalc
      call update_both_ghosts(n1,n2,dsalc,kstart,kend)
      call mgrd_qvapc
      call update_both_ghosts(n1,n2,qvapc,kstart,kend)

      !=============================================================
      !    for movies initialization
#ifdef MOVIE
      call inimov
      call mkmov_field
#endif
#ifdef PLANEMOVIE
      call inimov
      !call inimov_hdf_ycut
      call inimov_hdf_xcut
      call inimov_hdf_zcut
#endif

      !=============================================================
      !     check divergence
      call divgck(dmaxc,dmaxr)
      if(myid.eq.0) then
       write(*,'(3x,a,es10.3,a,es10.3)')'Max div base mesh:',dmaxc,'    Max div refined mesh:',dmaxr
       write(*,'(3x,a,es10.3,a,es10.3,/)')'tot KTE base mesh:',aengc,'    tot KTE refined mesh:',aengr
      endif

      !===========================================================
      !     main loop for time integration 
      ntstf=ntst
      tin(2) = MPI_WTIME()
      if(myid.eq.0) write(*,'(3x,a,f10.3,a,/)') 'Initialization Time = ', tin(2)-tin(1), ' sec.'

      DO ntime=1,ntstf                                           
        call cfl(cflm)        
        call cflr(cflmr)
        if(idtv.eq.1) then
          dt=cflmax/cflmr
          if(dt.gt.dtmax) dt=dtmax
          if(dt.lt.1.d-8) go to 166
        else
          cflm=cflm*dt
          if(cflm.gt.cfllim) go to 165
        endif

        beta=dt*nu*0.5d0
        ti(1) = MPI_WTIME()
        call tschem            ! tscheme
        time=time+dt           ! simulation time
        ti(2) = MPI_WTIME()    ! CPU time

            !=============================================================
            !       computing statistics
            if( dabs(time-dble(nint(time/tpin))*tpin) .lt. dt*0.5d0 ) then
                call vmaxv
                call avgvar
                call force_balance_online
                call aggregation_scale_diagnostics
                call mse_vortex_diagnostics
                call nusse_walls
                call CalcStats

                call divgck(dmaxc,dmaxr)
                if(myid.eq.0) then
                  write(*,*) '---------------------------------------- '
                  write(*,'(a,f12.4,a,i8,a,es10.2)')'  T = ',time,'   NTIME = ',ntime,'   DT = ',dt
                  write(*,'(a,f6.2,a)') '  Iteration Time = ', ti(2)-ti(1), ' sec.'
                  write(*,'(a,es10.3,a,es10.3)')'  Max div base mesh:',dmaxc,'  Max div refined mesh:',dmaxr
                  write(*,'(a,es10.3,a,es10.3)')'  tot KTE base mesh:',aengc,'  tot KTE refined mesh:',aengr
                endif
            endif

            !===============================================================
            !      for movie
#ifdef MOVIE
              if(dabs(time-dble(nint(time/tframe))*tframe).lt.0.5d0*dt)then
                call mkmov_field
              endif
#endif

#ifdef PLANEMOVIE
              if(dabs(time-dble(nint(time/tframe))*tframe).lt.0.5d0*dt)then
                ! Q-criterion movie output is intentionally disabled.
                !call mkmov_hdf_ycut
                call mkmov_hdf_xcut
                call mkmov_hdf_zcut
                call mkmov_field
                !call mkmov_mean
              endif
#endif

            !==============================================================
            !   check if reaching time limitations
            if(time.ge.tmax) go to 333

            !   save data every 24h
            if( mod(ti(2)-tin(1), 86400.d0) .lt. (ti(2)-ti(1)) )then
                if(myid.eq.0)then
                  write(*,*)' '
                  write(*,*)'******************************'
                  write(*,*)'  write data every 24 hours'
                  write(*,*)' '
                endif
                call mpi_write_continua
                call MPI_BARRIER(MPI_COMM_WORLD,ierr)
                if(myid.eq.0)write(*,*)'  continua files written'
            endif

             !   save data at fixed iterations
             if(dabs(time-dble(nint(time/trestart))*trestart).lt.0.5d0*dt)then
                 call mpi_write_continua
             endif
      !     end loop for time integration 
      !===========================================================
      ENDDO


!===============================================================
!      Exit
!===============================================================
333   continue
  
      tin(3) = MPI_WTIME()
      if(myid.eq.0) then
        write(*,'(a,f9.2,a)')'  Total Iteration Time = ',(tin(3)-tin(2))/3600.0,'h'
      endif

      call mpi_write_continua
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(myid.eq.0)write(*,*)'  continua files written'

      go to 167  

!m==============================================================    
165   continue 
      if(myid.eq.0) then
        write(*,164) 
      endif 
164   format(10x,'cfl too large!  ') 
      go to 167 

!m===============================================================       
166   continue
      if(myid.eq.0) then
        write(*,168) dt 
      endif
168   format(10x,'dt too small, DT= ',e14.7)
      go to 167 
                       
!m================================================================          
266   continue
      if(myid.eq.0) then
        write(*,268)  
      endif 
268   format(10x,'velocities diverged!  ') 
      go to 167 
      
!m================================================================          
169   continue
      if(myid.eq.0) then
        write(*,178) dmaxc
      endif  
178   format(10x,'too large local residue for mass conservation : ',e12.5,' at ')     
!      call divgloc
 
!m==============================================================                               
167   continue

      call phend

      return  
      end
