module online_budget_ke
  use param
  use local_arrays, only: q1, q2, q3
  use mgrd_arrays, only: dsalc, q1lr, q2lr, q3lr
  use mpi_param, only: kstart, kend, kstartr, kendr
  use mpih, only: MPI_COMM_WORLD, MPI_SUM, MDP
  implicit none
  private
  public :: compute_kinetic_budgets
contains
  subroutine compute_kinetic_budgets(do_ke, do_ke2d, ke, ke2d)
    logical, intent(in) :: do_ke, do_ke2d
    ! Signed drag work is reported separately, and is zero when alpha=0.
    ! ke   = K, P_B, epsilon, drag.
    ! ke2d = K_2D, P_3Dto2D, epsilon_2D, W_wall, drag,
    !        legacy K_2D, legacy K_total, legacy fraction.
    real, intent(out) :: ke(4), ke2d(8)

    ke = 0.0d0
    ke2d = 0.0d0
    if (do_ke) call compute_total_ke(ke)
    if (do_ke2d) then
      call compute_depth_mean_ke(ke2d(1:5))
      call compute_legacy_ke(ke2d(6:8))
    endif
  end subroutine compute_kinetic_budgets

  subroutine compute_total_ke(ke)
    real, intent(out) :: ke(4)
    integer :: i, j, k, ip, jp, ierr_budget
    real :: local(3), totals(3), height, inv_area_count
    real :: cell_weight, face_weight, u, v, w, w_top, b_face
    real :: grad_h2, grad_z2, grad_w_h2

    height = zc(n3) - zc(1)
    inv_area_count = 1.0d0 / (dble(n1m) * dble(n2m))
    local = 0.0d0
    do k = kstart, kend
      cell_weight = g3rm(k) / (dx3 * height)
      face_weight = g3rc(k) / (dx3 * height)
      do j = 1, n2m
        jp = jpv(j)
        do i = 1, n1m
          ip = ipv(i)
          u = q1(i,j,k)
          v = q2(i,j,k)
          w = q3(i,j,k)

          ! u and v live at layer centres in z; w is on the dual cells.
          local(1) = local(1) + 0.5d0 * cell_weight * (u*u + v*v)
          if (k > 1) then
            local(1) = local(1) + 0.5d0 * face_weight * w*w
            ! Same w-face buoyancy as hdnlq3; dsalc is the base-grid buoyancy.
            b_face = 0.5d0 * (dsalc(i,j,k) + dsalc(i,j,k-1))
            local(2) = local(2) + face_weight * w * b_face
          endif

          ! Summation by parts of invtrq[123]'s native Laplacians.
          grad_h2 = dx1q * ((q1(ip,j,k)-u)**2 + (q2(ip,j,k)-v)**2) &
                  + dx2q * ((q1(i,jp,k)-u)**2 + (q2(i,jp,k)-v)**2)
          grad_w_h2 = 0.0d0
          if (k > 1) then
            grad_w_h2 = dx1q * (q3(ip,j,k)-w)**2 + dx2q * (q3(i,jp,k)-w)**2
            grad_z2 = (dx3/g3rc(k)) * ((u-q1(i,j,k-1))**2 + (v-q2(i,j,k-1))**2)
          else
            grad_z2 = (2.0d0*dx3/g3rc(1)) * dble(ubcbot) * (u*u+v*v)
            w = 0.0d0
          endif
          if (k < n3m) then
            w_top = q3(i,j,k+1)
          else
            w_top = 0.0d0
            grad_z2 = grad_z2 + (2.0d0*dx3/g3rc(n3)) * dble(ubctop) * (u*u+v*v)
          endif
          grad_z2 = grad_z2 + (dx3/g3rm(k)) * (w_top-w)**2
          local(3) = local(3) + cell_weight*grad_h2 + face_weight*grad_w_h2 + grad_z2/height
        enddo
      enddo
    enddo
    call MPI_ALLREDUCE(local, totals, 3, MDP, MPI_SUM, MPI_COMM_WORLD, ierr_budget)
    ke(1:3) = totals * inv_area_count
    ke(3) = nu * ke(3)
    ke(4) = -2.0d0 * alpha * ke(1)
  end subroutine compute_total_ke

  subroutine compute_depth_mean_ke(ke2d)
    real, intent(out) :: ke2d(5)
    integer :: i, j, k, ip, im, jp, jm, ierr_budget
    real :: partial(n1m,n2m,2), mean_velocity(n1m,n2m,2)
    real :: height, inv_area_count, dz_weight, shear_factor
    real :: ux, uy, up, vp, u_center, v_center, u_corner, v_corner
    real :: local(2), totals(2)

    height = zc(n3) - zc(1)
    inv_area_count = 1.0d0 / (dble(n1m) * dble(n2m))
    partial = 0.0d0
    do k = kstart, kend
      dz_weight = g3rm(k) / (dx3 * height)
      do j = 1, n2m
        do i = 1, n1m
          partial(i,j,1) = partial(i,j,1) + dz_weight * q1(i,j,k)
          partial(i,j,2) = partial(i,j,2) + dz_weight * q2(i,j,k)
        enddo
      enddo
    enddo
    call MPI_ALLREDUCE(partial, mean_velocity, size(partial), MDP, MPI_SUM, MPI_COMM_WORLD, ierr_budget)

    ke2d = 0.0d0
    do j = 1, n2m
      jp = jpv(j)
      do i = 1, n1m
        ip = ipv(i)
        ux = mean_velocity(i,j,1)
        uy = mean_velocity(i,j,2)
        ke2d(1) = ke2d(1) + 0.5d0*(ux*ux + uy*uy)
        ke2d(3) = ke2d(3) &
                  + dx1q*((mean_velocity(ip,j,1)-ux)**2 + (mean_velocity(ip,j,2)-uy)**2) &
                  + dx2q*((mean_velocity(i,jp,1)-ux)**2 + (mean_velocity(i,jp,2)-uy)**2)
      enddo
    enddo

    local = 0.0d0
    do k = kstart, kend
      dz_weight = g3rm(k) / (dx3 * height)
      shear_factor = 0.0d0
      ! Top derivative minus bottom derivative.  At a stationary no-slip
      ! wall the sign is negative in either boundary contribution.
      if (k == 1) shear_factor = shear_factor - dble(ubcbot)*2.0d0*dx3/g3rc(1)
      if (k == n3m) shear_factor = shear_factor - dble(ubctop)*2.0d0*dx3/g3rc(n3)
      do j = 1, n2m
        jp = jpv(j)
        jm = jmv(j)
        do i = 1, n1m
          ip = ipv(i)
          im = imv(i)
          ux = mean_velocity(i,j,1)
          uy = mean_velocity(i,j,2)
          up = q1(i,j,k) - ux
          vp = q2(i,j,k) - uy
          ! Normal stresses at scalar cell centers and shear stress at
          ! the common corner of the control volumes, as in hdnlq1/2.
          u_center = 0.5d0 * (up + q1(ip,j,k) - mean_velocity(ip,j,1))
          v_center = 0.5d0 * (vp + q2(i,jp,k) - mean_velocity(i,jp,2))
          u_corner = 0.5d0 * (up + q1(i,jm,k) - mean_velocity(i,jm,1))
          v_corner = 0.5d0 * (vp + q2(im,j,k) - mean_velocity(im,j,2))
          local(1) = local(1) + dz_weight * ( &
             u_center**2 * (mean_velocity(ip,j,1)-ux)*dx1 &
           + v_center**2 * (mean_velocity(i,jp,2)-uy)*dx2 &
           + u_corner*v_corner * ((ux-mean_velocity(i,jm,1))*dx2 + (uy-mean_velocity(im,j,2))*dx1))
          local(2) = local(2) + shear_factor * (ux*q1(i,j,k) + uy*q2(i,j,k))
        enddo
      enddo
    enddo
    call MPI_ALLREDUCE(local, totals, 2, MDP, MPI_SUM, MPI_COMM_WORLD, ierr_budget)
    ke2d(1) = ke2d(1) * inv_area_count
    ke2d(2) = totals(1) * inv_area_count
    ke2d(3) = ke2d(3) * nu * inv_area_count
    ke2d(4) = totals(2) * (nu/height) * inv_area_count
    ke2d(5) = -2.0d0 * alpha * ke2d(1)
  end subroutine compute_depth_mean_ke

  subroutine compute_legacy_ke(legacy)
    real, intent(out) :: legacy(3)
    integer :: i, j, k, ip, jp, ierr_budget
    real :: partial(n1mr,n2mr,2), mean_velocity(n1mr,n2mr,2)
    real :: height, inv_area_count, dz_weight, u, v, w
    real :: local_total, global_total

    ! Cell-centre interpolation and layer thickness of the existing movie
    ! diagnostic (mkmov_field and the offline reduce_fraction.py).
    height = zcr(n3r) - zcr(1)
    inv_area_count = 1.0d0 / (dble(n1mr) * dble(n2mr))
    partial = 0.0d0
    local_total = 0.0d0
    do k = kstartr, kendr
      dz_weight = (zcr(k+1) - zcr(k)) / height
      do j = 1, n2mr
        jp = jpvr(j)
        do i = 1, n1mr
          ip = ipvr(i)
          u = 0.5d0 * (q1lr(i,j,k) + q1lr(ip,j,k))
          v = 0.5d0 * (q2lr(i,j,k) + q2lr(i,jp,k))
          w = 0.5d0 * (q3lr(i,j,k) + q3lr(i,j,k+1))
          partial(i,j,1) = partial(i,j,1) + dz_weight*u
          partial(i,j,2) = partial(i,j,2) + dz_weight*v
          local_total = local_total + 0.5d0*dz_weight*(u*u+v*v+w*w)
        enddo
      enddo
    enddo
    call MPI_ALLREDUCE(partial, mean_velocity, size(partial), MDP, MPI_SUM, MPI_COMM_WORLD, ierr_budget)
    call MPI_ALLREDUCE(local_total, global_total, 1, MDP, MPI_SUM, MPI_COMM_WORLD, ierr_budget)
    legacy(1) = 0.5d0 * sum(mean_velocity**2) * inv_area_count
    legacy(2) = global_total * inv_area_count
    legacy(3) = 0.0d0
    if (legacy(2) > 0.0d0) legacy(3) = legacy(1) / legacy(2)
  end subroutine compute_legacy_ke
end module online_budget_ke
