      subroutine block(n, p, irank, istart, iend, blcsz)
      implicit none
      integer,intent(in) :: n,p,irank
      integer,intent(out) :: istart,iend
      integer :: i
      integer,dimension(0:p-1),intent(out) :: blcsz
      
      do i=0,p-1
        blcsz(i) = floor(real((n+p-i-1)/p))
      enddo
      istart = sum(blcsz(0:irank))-blcsz(irank)+1
      iend = istart+blcsz(irank)-1

      end subroutine block

!=================================================           
      subroutine mpi_workdistribution
      use param
      use mpih 
      use mpi_param
      implicit none
      integer :: i
      
      if(.not. allocated(countj)) allocate(countj(0:numtasks-1))
      if(.not. allocated(countjp)) allocate(countjp(0:numtasks-1))
      if(.not. allocated(countk)) allocate(countk(0:numtasks-1))

      if(.not. allocated(countjr)) allocate(countjr(0:numtasks-1))
      if(.not. allocated(countkr)) allocate(countkr(0:numtasks-1))

      !EP   For PERIODIC pressure solver
      call block(n2m+2, numtasks, myid, jstartp, jendp, countjp)
      djp=jendp-jstartp+1

      call block(n2m, numtasks, myid, jstart, jend, countj)
      dj=jend-jstart+1
      
      call block(n3m, numtasks, myid, kstart, kend, countk)
      dk=kend-kstart+1

      countjr = countj*mref2
      countkr = countk*mref3
      jstartr = sum(countjr(0:myid))-countjr(myid)+1
      jendr = jstartr+countjr(myid)-1
      kstartr = sum(countkr(0:myid))-countkr(myid)+1
      kendr = kstartr+countkr(myid)-1

      djr=jendr-jstartr+1
      dkr=kendr-kstartr+1


      if( dj .lt. 1 ) then            
       write(*,*)'process ',myid,' has work load <1 cell in j direction'
       write(*,*)"Check grid dimensions and number of processes"
       
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif

      if( dk .lt. 1 ) then            
       write(*,*)'process ',myid,' has work load <1 cell in k direction'
       write(*,*)"Check grid dimensions and number of processes"
       
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif
  
      if(.not. allocated(offsetjp)) allocate(offsetjp(0:numtasks-1))
      if(.not. allocated(offsetj)) allocate(offsetj(0:numtasks-1))
      if(.not. allocated(offsetk)) allocate(offsetk(0:numtasks-1))
      
      offsetjp(:)=0
      offsetj(:)=0
      offsetk(:)=0
      do i=1,numtasks-1
        offsetjp(i)= offsetjp(i-1) + countjp(i-1)
        offsetj(i)= offsetj(i-1) + countj(i-1)
        offsetk(i)= offsetk(i-1) + countk(i-1)
      end do
      
      if(.not. allocated(offsetjr)) allocate(offsetjr(0:numtasks-1))
      if(.not. allocated(offsetkr)) allocate(offsetkr(0:numtasks-1))
      
      offsetjr(:)=0
      offsetkr(:)=0
      do i=1,numtasks-1
        offsetjr(i)= offsetjr(i-1) + countjr(i-1)
        offsetkr(i)= offsetkr(i-1) + countkr(i-1)
      end do

      !-------For MPI-IO--------------------------------
      mydata= n2*dk*n1
      mydatam = n2m*dk*n1m

      if(myid .eq. numtasks-1) mydata = n2*(dk+1)*n1
      
      if(.not. allocated(countf)) allocate(countf(0:numtasks-1))
      if(.not. allocated(offsetf)) allocate(offsetf(0:numtasks-1))
       
      call MPI_ALLGATHER(mydata, 1, MPI_INTEGER, countf, 1, MPI_INTEGER,&
     & MPI_COMM_WORLD,ierr)
    
      offsetf(:)=0
      do i=1,numtasks-1
        offsetf(i)= offsetf(i-1) + countf(i-1)
      end do
      
      end subroutine mpi_workdistribution 

!===============================================
      subroutine update_both_ghosts(n1,n2,q1,ks,ke)
      use mpih
      implicit none
      integer, intent(in) :: ks,ke
      real,intent(inout) :: q1(n1,n2,ks-lvlhalo:ke+lvlhalo)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag
      
      mydata= n1*n2*lvlhalo
      
      my_down=myid-1
      
      my_up=myid+1

      if(myid .eq. 0) my_down=MPI_PROC_NULL
      if(myid .eq. numtasks-1) my_up=MPI_PROC_NULL

      tag=1
      call MPI_ISEND(q1(1,1,ke-lvlhalo+1), mydata, MDP,my_up,tag,MPI_COMM_WORLD,req(1),ierr)
      call MPI_ISEND(q1(1,1,ks), mydata,  MDP,my_down,tag,MPI_COMM_WORLD,req(2), ierr) 
      call MPI_IRECV(q1(1,1,ks-lvlhalo), mydata,  MDP, my_down,tag,MPI_COMM_WORLD,req(3),ierr)
      call MPI_IRECV(q1(1,1,ke+1), mydata,  MDP,my_up, tag,MPI_COMM_WORLD,req(4),ierr)
      call MPI_Waitall(4,req,status,ierr)

      end subroutine update_both_ghosts

!=========================================
      subroutine update_upper_ghost(n1,n2,q1)
      use mpih
      use mpi_param, only: kstart,kend,dk
      implicit none
      real,intent(inout) :: q1(n1,n2,kstart-lvlhalo:kend+lvlhalo)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag
       
      mydata= n1*n2*lvlhalo
      
      my_down= myid-1
      
      my_up= myid+1

      if(myid .eq. 0) my_down=MPI_PROC_NULL
      if(myid .eq. numtasks-1) my_up=MPI_PROC_NULL
     
      tag=1
      call MPI_ISEND(q1(1,1,kstart), mydata, MDP,my_down, tag, MPI_COMM_WORLD, req(1), ierr) 
      call MPI_IRECV(q1(1,1,kend+1), mydata, MDP,my_up,tag, MPI_COMM_WORLD, req(2), ierr)
      call MPI_Waitall(2,req,status,ierr)
    
      end subroutine update_upper_ghost

!=========================================
      subroutine update_lower_ghost(n1,n2,q1)
      use mpih
      use mpi_param, only: kstart,kend,dk
      implicit none
      real,intent(inout) :: q1(n1,n2,kstart-lvlhalo:kend+lvlhalo)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag
       
      mydata= n1*n2*lvlhalo
      
      my_down= myid-1
      my_up= myid+1

      if(myid .eq. 0) my_down=MPI_PROC_NULL
      if(myid .eq. numtasks-1) my_up=MPI_PROC_NULL
      
      tag=1
      
      call MPI_ISEND(q1(1,1,kend-lvlhalo+1), mydata,  MDP,my_up, tag, MPI_COMM_WORLD, req(1), ierr)
      call MPI_IRECV(q1(1,1,kstart-lvlhalo), mydata,  MDP,my_down,tag, MPI_COMM_WORLD, req(2), ierr)
      call MPI_Waitall(2,req,status,ierr)
    
      end subroutine update_lower_ghost

!=========================================
      subroutine update_add_upper_ghost(n1,n2,q1)
      use mpih
      use mpi_param, only: kstart,kend,dk
      implicit none
      real,intent(inout) :: q1(n1,n2,kstart-lvlhalo:kend+lvlhalo-1)
      real :: buf(n1,n2,lvlhalo)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag

      integer :: ic,jc
      real :: cksum

       
      mydata= n1*n2*lvlhalo
      
      my_down= myid-1
      my_up= myid+1

      buf=0.0d0

      if(myid .eq. 0) my_down= MPI_PROC_NULL
      if(myid .eq. numtasks-1) my_up= MPI_PROC_NULL
     
      tag=1

      call MPI_ISEND(q1(1,1,kstart-lvlhalo),mydata,MDP,my_down, tag, MPI_COMM_WORLD, req(1), ierr)
      call MPI_IRECV(buf(1,1,1), mydata, MDP,my_up,tag, MPI_COMM_WORLD, req(2), ierr)
      call MPI_Waitall(2,req,status,ierr)

      do ic=1,lvlhalo 
         jc=kend-lvlhalo+ic
         if(jc.eq.kend)then
            q1(:,:,jc) = buf(:,:,ic)
         else
            q1(:,:,jc) = q1(:,:,jc) + buf(:,:,ic)
         end if
      end do

      end subroutine update_add_upper_ghost

!=========================================
      subroutine update_add_lower_ghost(n1,n2,q1)
      use mpih
      use mpi_param, only: kstart,kend,dk
      implicit none
      real,intent(inout) :: q1(n1,n2,kstart-lvlhalo:kend+lvlhalo-1)
      real :: buf(n1,n2,lvlhalo)
      integer,intent(in) :: n1,n2
      integer :: mydata
      integer :: my_down, my_up,tag
      integer :: ic,jc
       
      mydata= n1*n2*(lvlhalo)
      
      my_down= myid-1
      my_up= myid+1

      buf=0.0d0
      if(myid .eq. 0) my_down= MPI_PROC_NULL
      if(myid .eq. numtasks-1) my_up= MPI_PROC_NULL
      
      tag=1
      call MPI_ISEND(q1(1,1,kend),mydata,MDP,my_up, tag, MPI_COMM_WORLD, req(1), ierr)
      call MPI_IRECV(buf(1,1,1), mydata, MDP,my_down,tag, MPI_COMM_WORLD, req(2), ierr)
      call MPI_Waitall(2,req,status,ierr)

      do ic=1,lvlhalo
         jc=kstart+ic-2
         q1(:,:,jc) = q1(:,:,jc) + buf(:,:,ic)
      end do

      end subroutine update_add_lower_ghost

!==============================================
        subroutine mpi_globalsum_double_arr(var,nvar)
        use mpih
          implicit none
          real,intent(inout),dimension(nvar) :: var
          real,dimension(nvar) :: var2
          integer,intent(in) :: nvar

          call MPI_ALLREDUCE(var,var2,nvar,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)          
          var = var2

        end subroutine mpi_globalsum_double_arr

!==============================================
      subroutine mem_alloc
      use mpih
      use param
      use mpi_param
      use local_arrays
      use stat_arrays

      implicit none
      integer :: merr, merr_all
     
      merr_all = 0 
      !-------------------------------------------------
      ! Arrays with ghost cells
      !-------------------------------------------------
      allocate(q1(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr)
      merr_all = merr_all + merr
      allocate(q2(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr)
      merr_all = merr_all + merr
      allocate(q3(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr)
      merr_all = merr_all + merr
      allocate(pr(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr)
      merr_all = merr_all + merr
      allocate(dens(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr)
      merr_all = merr_all + merr
      allocate(dph(1:n1m,1:n2m+2,kstart-lvlhalo:kend+lvlhalo),stat=merr)
      merr_all = merr_all + merr
      allocate(vz_me(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr)
      merr_all = merr_all + merr
      allocate(vy_me(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr)
      merr_all = merr_all + merr
      allocate(vx_me(1:n1,1:n2,kstart-lvlhalo:kend+lvlhalo), stat=merr)
      merr_all = merr_all + merr

      !-------------------------------------------------
      ! Arrays for linear system
      !-------------------------------------------------
      allocate(rhs(1:n1m,1:n2m,kstart:kend), stat=merr)
      merr_all = merr_all + merr
      allocate(dq(1:n1,1:n2,kstart:kend), stat=merr)
      merr_all = merr_all + merr                    
      allocate(ru1(1:n1,1:n2,kstart:kend), stat=merr)
      merr_all = merr_all + merr
      allocate(ru2(1:n1,1:n2,kstart:kend), stat=merr)
      merr_all = merr_all + merr
      allocate(ru3(1:n1,1:n2,kstart:kend), stat=merr)
      merr_all = merr_all + merr
      allocate(qcap(1:n1,1:n2,kstart:kend), stat=merr)
      merr_all = merr_all + merr
      allocate(hro(1:n1,1:n2,kstart:kend), stat=merr)
      merr_all = merr_all + merr
      allocate(ruro(1:n1,1:n2,kstart:kend), stat=merr)
      merr_all = merr_all + merr

      !---------------------------------------------------
      ! Arrays for salinity
      !---------------------------------------------------
      allocate(dsal(1:n1r,1:n2r,kstartr-lvlhalo:kendr+lvlhalo),stat=merr)
      merr_all = merr_all + merr
      allocate(qvap(1:n1r,1:n2r,kstartr-lvlhalo:kendr+lvlhalo),stat=merr)
      merr_all = merr_all + merr
      allocate(hsal(1:n1r,1:n2r,kstartr:kendr), stat=merr)
      merr_all = merr_all + merr      
      allocate(hqvap(1:n1r,1:n2r,kstartr:kendr), stat=merr)
      merr_all = merr_all + merr
      allocate(rusal(1:n1r,1:n2r,kstartr:kendr), stat=merr)
      merr_all = merr_all + merr
      allocate(ruqvap(1:n1r,1:n2r,kstartr:kendr), stat=merr)
      merr_all = merr_all + merr
      allocate(rhsr(1:n1mr,1:n2mr,kstartr:kendr), stat=merr)
      merr_all = merr_all + merr
      allocate(kpa(1:n1mr,1:n2mr,kstartr-lvlhalo:kendr+lvlhalo), stat=merr)
      merr_all = merr_all + merr
      allocate(T(1:n1r,1:n2r,kstartr-lvlhalo:kendr+lvlhalo), stat=merr)
      merr_all = merr_all + merr
      allocate(qs(1:n1r,1:n2r,kstartr-lvlhalo:kendr+lvlhalo), stat=merr)
      merr_all = merr_all + merr
      allocate(cond_term(1:n1r,1:n2r,kstartr-lvlhalo:kendr+lvlhalo), stat=merr)
      merr_all = merr_all + merr
      allocate(temp_me(1:n1r,1:n2r,kstartr-lvlhalo:kendr+lvlhalo),stat=merr)
      merr_all = merr_all + merr

      if(merr_all.ne.0)then
        write(*,*)myid, ' memory alloc error'
        write(*,*)merr_all
      endif

      return
      end subroutine mem_alloc

!==================================================      
      
      subroutine mem_dealloc
      use local_arrays
      use mpi_param
      use mpih
      use stat_arrays

      implicit none
      
      if(allocated(q1)) deallocate(q1)
      if(allocated(q2)) deallocate(q2)
      if(allocated(q3)) deallocate(q3)
      if(allocated(dens)) deallocate(dens)
      if(allocated(dsal)) deallocate(dsal)
      if(allocated(qvap)) deallocate(qvap)
      if(allocated(kpa)) deallocate(kpa)

      if(allocated(temp_me)) deallocate(temp_me)
      if(allocated(vx_me)) deallocate(vx_me)
      if(allocated(vy_me)) deallocate(vy_me)
      if(allocated(vz_me)) deallocate(vz_me)

      if(allocated(pr)) deallocate(pr)
      if(allocated(dph)) deallocate(dph)
      if(allocated(dq)) deallocate(dq)
      if(allocated(qcap)) deallocate(qcap)

      if(allocated(hro)) deallocate(hro)
      if(allocated(rhs)) deallocate(rhs)
      if(allocated(ru1)) deallocate(ru1)
      if(allocated(ru2)) deallocate(ru2)
      if(allocated(ru3)) deallocate(ru3)
      if(allocated(ruro)) deallocate(ruro)

      if(allocated(hsal)) deallocate(hsal)
      if(allocated(hqvap)) deallocate(hqvap)
      if(allocated(rusal)) deallocate(rusal)
      if(allocated(ruqvap)) deallocate(ruqvap)
      if(allocated(rhsr)) deallocate(rhsr)

      if(allocated(T)) deallocate(T)
      if(allocated(qs)) deallocate(qs)
      if(allocated(cond_term)) deallocate(cond_term)

      if(allocated(countj)) deallocate(countj)
      if(allocated(countk)) deallocate(countk)
      if(allocated(countjp)) deallocate(countjp)

      if(allocated(offsetj)) deallocate(offsetj)
      if(allocated(offsetk)) deallocate(offsetk)
      if(allocated(offsetjp)) deallocate(offsetjp)

      if(allocated(countjr)) deallocate(countjr)
      if(allocated(countkr)) deallocate(countkr)

      if(allocated(offsetjr)) deallocate(offsetjr)
      if(allocated(offsetkr)) deallocate(offsetkr)
      
      if(allocated(countf)) deallocate(countf)
      if(allocated(offsetf)) deallocate(offsetf)

      return    
      end subroutine mem_dealloc
!================================================
