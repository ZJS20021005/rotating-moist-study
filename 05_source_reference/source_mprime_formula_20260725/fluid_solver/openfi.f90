      subroutine openfi
      use param
      use mpih
      implicit none

      IF(myid.eq.0)then

        open(95, file='./data/vmax.out',status='unknown',access='sequential',position='append')
        open(96, file='./data/avgvar.out',status='unknown',access='sequential',position='append')
        open(97, file='./data/nusse_walls.out',status='unknown',access='sequential',position='append')
        open(98, file='./data/layer.out',status='unknown',access='sequential',position='append')
        open(99, file='./data/nu_profiles.out',status='unknown',access='sequential',position='append')
        open(100, file='./data/kh_energy.out',status='unknown',access='sequential',position='append')
        open(101, file='./data/mse_aggregation.out',status='unknown',access='sequential',position='append')
        open(102, file='./data/mse_variance_profile.out',status='unknown',access='sequential',position='append')
        open(103, file='./data/w_z075_spectral_length.out',status='unknown',access='sequential',position='append')
        open(104, file='./data/mse_aggregation_scales.out',status='unknown',access='sequential',position='append')
        open(105, file='./data/mprime_squared.out',status='unknown',access='sequential',position='append')

      ! reset the time history
      if(ireset.eq.1 .or. nread.eq.0)then
          rewind(95)
          rewind(96)
          rewind(97)
          rewind(98)
          rewind(99)
          rewind(100)
          rewind(101)
          rewind(102)
          rewind(103)
          rewind(104)
          rewind(105)
      endif

      ENDIF

      return
      end   
      
!==============================================

      subroutine closefi
      use mpih
      implicit none
      
      if(myid.eq.0)then

      close(95)
      close(96)
      close(97)
      close(98)
      close(99)
       close(100)
       close(101)
       close(102)
       close(103)
       close(104)
       close(105)
      endif
     
      return      
      end
