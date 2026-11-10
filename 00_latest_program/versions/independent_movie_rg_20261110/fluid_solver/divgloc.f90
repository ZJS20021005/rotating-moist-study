!***********************************************************************
!EP      Output the location(s) of excessive divergence
!EP      i.e. i,j,k where dqcap(i,j,k) > resid
!***********************************************************************

      subroutine divgloc
      use param
      use local_arrays, only: q2,q3,q1
      use mpih
      use mpi_param, only: kstart,kend
      implicit none
      integer :: jc,kc,kp,jp,ic,ip
      real    :: dqcap
        
      if(myid.eq.0) write(*,*) "I   J   K   MYID"
      do kc=kstart,kend
        kp=kc+1
        do jc=1,n2m
          jp=jpv(jc)
            do ic=1,n1m
              ip=ipv(ic)
              dqcap= (q1(ip,jc,kc)-q1(ic,jc,kc))*dx1&
     &              +(q2(ic,jp,kc)-q2(ic,jc,kc))*dx2&
     &              +(q3(ic,jc,kp)-q3(ic,jc,kc))*udx3m(kc)
              if (dabs(dqcap).gt.resid) then
                 write(*,*) ic,jc,kc,myid
              endif
      enddo
      enddo
      enddo
      
      return     
      end         
