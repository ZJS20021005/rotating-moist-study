! Online LMSE budget evaluated on the scalar grid at the end of a time step.
! This routine only reads the live state. It never calls stage operators
! (which overwrite the DNS history arrays) or saves volume-size buffers.
subroutine compute_lmse_budget(values, level_variances, gmse)
  use param
  use local_arrays, only: dsal, qvap, kpa
  use mgrd_arrays, only: q1lr, q2lr, q3lr
  use mpi_param, only: kstartr, kendr
  use mpih, only: MPI_COMM_WORLD, MPI_SUM, MDP, ierr, myid
  implicit none

  real, intent(out) :: values(5), level_variances(3), gmse(5)
  integer :: i, j, k, im, ip, jm, jp, km, kp, level
  integer :: first_plane, last_plane, target_levels(3)
  integer, parameter :: m_column=1, xflux_column=2, yflux_column=3, diffusion_column=4, wall_column=5
  real :: local(8), global(8), area_inv, height_inv, dz, half_dz, weight
  real :: center, east, west, north, south, below, above, anomaly
  real :: adv_h, flux_top, flux_bot, adv_v, grad_production
  real :: flux_west, flux_south, boundary_flux
  real :: pot, pot_e, pot_w, pot_n, pot_s, pot_t, pot_b, wall
  real :: coeff_e, coeff_w, coeff_n, coeff_s, coeff_t, coeff_b
  real :: diffusion_h, diffusion_v, a_top, a_bot, variance_plane
  real :: gamma_sm, column_anomaly, divergence, grad_x, grad_y
  real :: mbar(max(1,kstartr-1):min(n3mr,kendr+1))
  real, allocatable :: local_columns(:,:,:), global_columns(:,:,:)
  external :: MPI_ALLREDUCE, MPI_REDUCE

  ! Each z slab owns a complete periodic xy plane. The final stage has
  ! already exchanged scalar halos, so the neighboring averages are local.
  first_plane = max(1,kstartr-1)
  last_plane = min(n3mr,kendr+1)
  area_inv = 1.0d0 / (dble(n1mr)*dble(n2mr))
  height_inv = 1.0d0 / (zcr(n3r)-zcr(1))
  gamma_sm = gamma*Sm
  do k=first_plane,last_plane
    mbar(k) = sum(dsal(1:n1mr,1:n2mr,k) + gamma*qvap(1:n1mr,1:n2mr,k))*area_inv
  enddo

  target_levels(1) = minloc(abs(zmr(1:n3mr)-0.25d0*alx3),dim=1)
  target_levels(2) = minloc(abs(zmr(1:n3mr)-0.50d0*alx3),dim=1)
  target_levels(3) = minloc(abs(zmr(1:n3mr)-0.75d0*alx3),dim=1)
  local = 0.0d0
  allocate(local_columns(n1mr,n2mr,wall_column))
  local_columns = 0.0d0
  do k=kstartr,kendr
    dz = zcr(k+1)-zcr(k)
    half_dz = 0.5d0*dz
    weight = dz*area_inv*height_inv
    variance_plane = 0.0d0
    km = kmvr(k)
    kp = kpvr(k)
    a_top = 0.0d0
    a_bot = 0.0d0
    if (k.lt.n3mr) a_top = g3rmr(k+1)/(g3rmr(k)+g3rmr(k+1))
    if (k.gt.1) a_bot = g3rmr(k)/(g3rmr(k-1)+g3rmr(k))
    do j=1,n2mr
      jm = jmvr(j)
      jp = jpvr(j)
      do i=1,n1mr
        im = imvr(i)
        ip = ipvr(i)
        center = dsal(i,j,k)+gamma*qvap(i,j,k)
        east = dsal(ip,j,k)+gamma*qvap(ip,j,k)
        west = dsal(im,j,k)+gamma*qvap(im,j,k)
        north = dsal(i,jp,k)+gamma*qvap(i,jp,k)
        south = dsal(i,jm,k)+gamma*qvap(i,jm,k)
        anomaly = center-mbar(k)
        variance_plane = variance_plane+anomaly**2

        ! Same face fluxes as hdnlsa and hdnlqvap. The n3m test in
        ! that operator is kept here, including for a refined z grid.
        flux_west = q1lr(i,j,k)*(center+west)
        flux_south = q2lr(i,j,k)*(center+south)
        adv_h = 0.5d0*dx1r*(q1lr(i+1,j,k)*(east+center)-flux_west) &
              + 0.5d0*dx2r*(q2lr(i,j+1,k)*(north+center)-flux_south)
        flux_top = 0.0d0
        flux_bot = 0.0d0
        grad_production = 0.0d0
        if (k.lt.n3mr) then
          above = dsal(i,j,k+1)+gamma*qvap(i,j,k+1)
          if (k.ne.n3m .or. k.eq.1) then
            flux_top = q3lr(i,j,k+1)*(a_top*center+(1.0d0-a_top)*above)
          endif
          ! On the stretched grid the discrete product rule pairs the
          ! face mean gradient with the adjoint interpolation of m'.
          ! Each inter-layer face is counted once by its lower cell.
          grad_production = -q3lr(i,j,k+1)*(mbar(k+1)-mbar(k)) &
                            *((1.0d0-a_top)*anomaly+a_top*(above-mbar(k+1)))
        endif
        if (k.gt.1) then
          below = dsal(i,j,k-1)+gamma*qvap(i,j,k-1)
          flux_bot = q3lr(i,j,k)*(a_bot*below+(1.0d0-a_bot)*center)
        endif
        adv_v = (flux_top-flux_bot)*udx3mr(k)

        ! D_b and D_q use the same kpa face stencil, with Sm scaling
        ! the vapor operator. Linearity therefore gives D(b+gamma*Sm*q).
        pot = dsal(i,j,k)+gamma_sm*qvap(i,j,k)
        pot_e = dsal(ip,j,k)+gamma_sm*qvap(ip,j,k)
        pot_w = dsal(im,j,k)+gamma_sm*qvap(im,j,k)
        pot_n = dsal(i,jp,k)+gamma_sm*qvap(i,jp,k)
        pot_s = dsal(i,jm,k)+gamma_sm*qvap(i,jm,k)
        coeff_e = 0.5d0*(kpa(ip,j,k)+kpa(i,j,k))
        coeff_w = 0.5d0*(kpa(im,j,k)+kpa(i,j,k))
        coeff_n = 0.5d0*(kpa(i,jp,k)+kpa(i,j,k))
        coeff_s = 0.5d0*(kpa(i,jm,k)+kpa(i,j,k))
        diffusion_h = dx1qr*(coeff_e*(pot_e-pot)-coeff_w*(pot-pot_w)) &
                    + dx2qr*(coeff_n*(pot_n-pot)-coeff_s*(pot-pot_s))

        if (k.lt.n3mr) then
          coeff_t = (kpa(i,j,kp)*(zcr(kp)-zmr(k))+kpa(i,j,k)*(zmr(kp)-zcr(kp))) &
                    /(zmr(kp)-zmr(k))
          pot_t = dsal(i,j,kp)+gamma_sm*qvap(i,j,kp)
        endif
        if (k.gt.1) then
          coeff_b = (kpa(i,j,k)*(zcr(k)-zmr(km))+kpa(i,j,km)*(zmr(k)-zcr(k))) &
                    /(zmr(k)-zmr(km))
          pot_b = dsal(i,j,km)+gamma_sm*qvap(i,j,km)
        endif
        if (k.eq.1) then
          wall = dsalbot+gamma_sm*qvapbot &
               + A_sbotmod*sin(2.0d0*pi*k_sbotmod*ymr(j))*sin(2.0d0*pi*k_sbotmod*xmr(i))
          boundary_flux = kpa(i,j,k)*(pot-wall)/(zmr(k)-zcr(k))
          diffusion_v = (coeff_t*(pot_t-pot)/(zmr(kp)-zmr(k))-boundary_flux)/dz
          local_columns(i,j,wall_column) = local_columns(i,j,wall_column)-boundary_flux
        else if (k.eq.n3mr) then
          wall = dsaltop+gamma_sm*qvaptop+A_stopmod*sin(2.0d0*pi*k_stopmod*ymr(j))
          boundary_flux = kpa(i,j,k)*(wall-pot)/(zcr(k+1)-zmr(k))
          diffusion_v = (boundary_flux-coeff_b*(pot-pot_b)/(zmr(k)-zmr(km)))/dz
          local_columns(i,j,wall_column) = local_columns(i,j,wall_column)+boundary_flux
        else
          diffusion_v = (ap3sskr(k)*pot_t+0.5d0*ac3sskr(k)*pot)*coeff_t &
                      + (am3sskr(k)*pot_b+0.5d0*ac3sskr(k)*pot)*coeff_b
        endif

        local(2) = local(2)-weight*anomaly*adv_h
        local(3) = local(3)-weight*anomaly*adv_v
        local(4) = local(4)+area_inv*height_inv*grad_production
        local(5) = local(5)-weight*anomaly*(diffusion_h+diffusion_v)

        ! Integrate full columns without a height normalization. The west
        ! and south fluxes live on the scalar-cell faces of the DNS stencil.
        local_columns(i,j,m_column) = local_columns(i,j,m_column)+dz*anomaly
        local_columns(i,j,xflux_column) = local_columns(i,j,xflux_column)+half_dz*flux_west
        local_columns(i,j,yflux_column) = local_columns(i,j,yflux_column)+half_dz*flux_south
        local_columns(i,j,diffusion_column) = local_columns(i,j,diffusion_column)+dz*diffusion_h
      enddo
    enddo
    local(1) = local(1)+weight*variance_plane
    do level=1,3
      if (k.eq.target_levels(level)) local(level+5)=area_inv*variance_plane
    enddo
  enddo

  call MPI_ALLREDUCE(local,global,8,MDP,MPI_SUM,MPI_COMM_WORLD,ierr)
  values = global(1:5)
  level_variances = global(6:8)

  ! Only the output rank needs the assembled columns. MPI ignores recvbuf
  ! on the other ranks, which supply a valid one-element buffer.
  if (myid.eq.0) then
    allocate(global_columns(n1mr,n2mr,wall_column))
  else
    allocate(global_columns(1,1,1))
  endif
  call MPI_REDUCE(local_columns,global_columns,size(local_columns),MDP,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  gmse = 0.0d0
  if (myid.ne.0) return

  ! Project out horizontal means, including roundoff in the summed anomaly.
  global_columns(:,:,m_column) = global_columns(:,:,m_column)-sum(global_columns(:,:,m_column))*area_inv
  global_columns(:,:,wall_column) = global_columns(:,:,wall_column)-sum(global_columns(:,:,wall_column))*area_inv
  do j=1,n2mr
    jm = jmvr(j)
    jp = jpvr(j)
    do i=1,n1mr
      im = imvr(i)
      ip = ipvr(i)
      column_anomaly = global_columns(i,j,m_column)
      divergence = dx1r*(global_columns(ip,j,xflux_column)-global_columns(i,j,xflux_column)) &
                 + dx2r*(global_columns(i,jp,yflux_column)-global_columns(i,j,yflux_column))
      grad_x = dx1r*(column_anomaly-global_columns(im,j,m_column))
      grad_y = dx2r*(column_anomaly-global_columns(i,jm,m_column))
      gmse(1) = gmse(1)+column_anomaly**2
      gmse(2) = gmse(2)-column_anomaly*divergence
      ! Periodicity pairs the horizontal diffusion operator with minus the
      ! diffusivity-weighted gradient product, also when Sm differs from 1.
      gmse(3) = gmse(3)-column_anomaly*global_columns(i,j,diffusion_column)
      gmse(4) = gmse(4)+column_anomaly*global_columns(i,j,wall_column)
      gmse(5) = gmse(5)+global_columns(i,j,xflux_column)*grad_x+global_columns(i,j,yflux_column)*grad_y
    enddo
  enddo
  gmse = gmse*area_inv
end subroutine compute_lmse_budget

