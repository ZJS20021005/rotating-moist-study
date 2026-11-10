      subroutine initia
      use param
      use local_arrays
      implicit none
      
      q1=0.d0
      q2=0.d0
      q3=0.d0
      dens=0.d0
      pr=0.d0
      dph=0.d0

      dq=0.d0
      rhs=0.d0
      ru1=0.d0
      ru2=0.d0
      ru3=0.d0
      ruro=0.d0
      qcap=0.d0
      hro=0.d0

      dsal=0.d0
      qvap=0.d0
      hqvap=0.d0
      hsal=0.d0
      rusal=0.d0
      ruqvap=0.d0
      rhsr=0.d0

      !m==================================         
      xc = 0.d0
      xm = 0.d0
      yc = 0.d0
      ym = 0.d0
      ap3ck = 0.d0 
      ac3ck = 0.d0 
      am3ck = 0.d0
      ap3sk = 0.d0
      ac3sk = 0.d0
      am3sk = 0.d0
      am3ssk = 0.d0
      ac3ssk = 0.d0
      ap3ssk = 0.d0
      zc = 0.d0
      zm = 0.d0
      g3rc = 0.d0
      g3rm = 0.d0
           
      !m==================================         
      xcr = 0.d0
      xmr = 0.d0
      ycr = 0.d0
      ymr = 0.d0
      am3sskr = 0.d0
      ac3sskr = 0.d0
      ap3sskr = 0.d0
      zcr = 0.d0
      zmr = 0.d0
      g3rcr = 0.d0
      g3rmr = 0.d0
     
      return 
      end   
