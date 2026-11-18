      subroutine meshes
      use param
      use mpih
      implicit none
 
      dx1=rext1/dble(n1m)
      dx2=rext2/dble(n2m)
      dx3=alx3/dble(n3m)
      dx1=1.d0/dx1
      dx2=1.d0/dx2
      dx3=1.d0/dx3
      dx1q=dx1*dx1                                                      
      dx2q=dx2*dx2                                                      
      dx3q=dx3*dx3   
      
      dx1r=rext1/dble(n1mr)
      dx2r=rext2/dble(n2mr)
      dx3r=alx3/dble(n3mr)
      dx1r=1.d0/dx1r
      dx2r=1.d0/dx2r
      dx3r=1.d0/dx3r
      dx1qr=dx1r*dx1r   
      dx2qr=dx2r*dx2r   
      dx3qr=dx3r*dx3r 
 
      return                                                            
      end                                                               
