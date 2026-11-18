      subroutine mgrd_qvapc
      use param
      use local_arrays, only: qvap 
      use mgrd_arrays, only: qvapc
      use mpi_param
      use mpih
      implicit none
       
      integer ic,jc,kc,icr,jcr,kcr

      real ldzr, ldz, dmrefxy, qvalloc

      !m=========================================================
      !  single mesh

      IF(mref1.eq.1 .and. mref2.eq.1 .and. mref3.eq.1)then

        do kc=kstart-1,kend+1
         do jc=1,n2
         do ic=1,n1
           qvapc(ic,jc,kc) = qvap(ic,jc,kc)
         enddo
         enddo
        enddo

      !m=========================================================
      !  multiple mesh
      ELSE

        dmrefxy = 1.d0/dble(mref1*mref2)

        do kc=kstart,kend
         ldz = 1.d0/(zc(kc+1)-zc(kc))
         do jc=1,n2m
         do ic=1,n1m
           qvalloc = 0.d0
           do kcr=(kc-1)*mref3+1,kc*mref3
             ldzr = zcr(kcr+1)-zcr(kcr)
             do jcr=(jc-1)*mref2+1,jc*mref2
             do icr=(ic-1)*mref1+1,ic*mref1
               qvalloc=qvalloc+qvap(icr,jcr,kcr)*ldzr
             enddo
             enddo
           enddo
           qvapc(ic,jc,kc) = qvalloc*ldz*dmrefxy
         enddo
         enddo
        enddo

        !   periodic B.C.
        do kc=kstart,kend
         do jc=1,n2m
           qvapc(n1,jc,kc) = qvapc(1,jc,kc)
         enddo
         do ic=1,n1
           qvapc(ic,n2,kc) = qvapc(ic,1,kc)
         enddo
        enddo

        !  boundary value 
        if(kstart.eq.1)then
         do jc=1,n2
          do ic=1,n1
            qvapc(ic,jc,0) = qvapbot +  A_sbotmod*sin(2.0*pi*k_sbotmod*ym(jc))*sin(2.0*pi*k_sbotmod*xm(ic))
          enddo
         enddo
        endif

        if(kend.eq.n3m)then
         do jc=1,n2
          do ic=1,n1
           qvapc(ic,jc,n3) = qvaptop +  A_stopmod*sin(2.0*pi*k_stopmod*ym(jc))
          enddo
         enddo
        endif
        call update_both_ghosts(n1,n2,qvapc,kstart,kend)
      ENDIF

      return
      end  subroutine mgrd_qvapc
