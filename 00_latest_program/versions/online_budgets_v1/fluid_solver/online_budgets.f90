module online_budgets
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use param, only: time, nread, ireset
  use mpih
  use online_budget_ke, only: compute_kinetic_budgets
  implicit none
  private
  public :: initialize_budgets, sample_budgets, finish_budgets
  integer :: units(4), samples = 0
  double precision :: elapsed = 0.0d0
  real :: previous_time, previous_values(4), previous_rhs(4)
  logical :: has_previous = .false.

contains
  subroutine initialize_budgets
    integer :: status, ierr

    status = 0
    if (myid == 0) then
      call execute_command_line('mkdir -p diagnostics/budget', exitstat=status)
      if (status == 0) then
        call open_budget(1, 'lmse', 'time,sigma_m2,P_H,P_V,P_grad,D_m,LMSE_rhs,' // &
          'd_sigma_m2_dt,rhs_interval,closure,production_residual', status)
        call open_budget(2, 'ke', 'time,K,P_B,epsilon,KE_rhs,dK_dt,rhs_interval,closure,drag_work', status)
        call open_budget(3, 'ke2d', 'time,K_2D,P_3Dto2D,epsilon_2D,W_wall,KE2D_rhs,dK_2D_dt,' // &
          'rhs_interval,closure,drag_work,K_2D_legacy,K_legacy,E2D_over_Etot_legacy', status)
        call open_budget(4, 'gmse', 'time,GMSE_variance,GMSE_PH,GMSE_DH,GMSE_SB,GMSE_rhs,' // &
          'd_GMSE_variance_dt,rhs_interval,closure,closure_abs,closure_rel,GMSE_PH_ibp,GMSE_PH_residual', status)
      end if
      if (status == 0) write(*,'(a)') 'DNS online-budgets-v2: LMSE, KE, KE2D and GMSE budgets enabled; vortex tracking disabled'
    end if
    call MPI_BCAST(status, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (status /= 0) call MPI_ABORT(MPI_COMM_WORLD, status, ierr)
  end subroutine

  subroutine open_budget(index, name, header, status)
    integer, intent(in) :: index
    character(*), intent(in) :: name, header
    integer, intent(inout) :: status
    character(:), allocatable :: path
    logical :: new_file
    integer :: file_size

    if (status /= 0) return
    path = 'diagnostics/budget/' // name // '_budget.csv'
    inquire(file=path, size=file_size, iostat=status)
    if (status /= 0) return
    new_file = file_size <= 0 .or. nread == 0 .or. ireset == 1
    if (new_file) then
      open(newunit=units(index), file=path, status='replace', action='write', iostat=status)
      if (status == 0) write(units(index), '(a)', iostat=status) header
    else
      open(newunit=units(index), file=path, status='old', position='append', action='write', iostat=status)
    end if
    if (status /= 0) write(*,'(a)') 'Unable to open budget series: ' // path
  end subroutine

  subroutine sample_budgets(legacy_variance)
    real, intent(out) :: legacy_variance(4)
    real :: mse(5), ke(4), ke2d(8), gmse(5)
    real :: state(4), rhs(4), rate(4), interval_rhs(4), residual(4), gmse_relative
    double precision :: started
    external :: compute_lmse_budget

    started = MPI_WTIME()
    call compute_lmse_budget(mse, legacy_variance(1:3), gmse)
    legacy_variance(4) = mse(1)
    call compute_kinetic_budgets(.true., .true., ke, ke2d)
    samples = samples + 1
    if (myid /= 0) then
      elapsed = elapsed + MPI_WTIME() - started
      return
    end if

    state = [mse(1), ke(1), ke2d(1), gmse(1)]
    rhs = [2*(mse(2)+mse(3)-mse(5)), ke(2)-ke(3)+ke(4), &
      ke2d(2)-ke2d(3)+ke2d(4)+ke2d(5), 2*(gmse(2)-gmse(3)+gmse(4))]
    rate = ieee_value(0.0, ieee_quiet_nan)
    interval_rhs = rate
    gmse_relative = rate(4)
    if (has_previous) then
      if (time > previous_time) then
        rate = (state-previous_values)/(time-previous_time)
        interval_rhs = 0.5*(previous_rhs+rhs)
        gmse_relative = abs(rate(4)-interval_rhs(4))/max(abs(rate(4)), abs(interval_rhs(4)), tiny(rate(4)))
      end if
    end if
    residual = rate - interval_rhs
    write(units(1), '(*(es25.16e3,:,","))') time, mse, rhs(1), &
      rate(1), interval_rhs(1), residual(1), mse(2)+mse(3)-mse(4)
    write(units(2), '(*(es25.16e3,:,","))') time, ke(1:3), rhs(2), &
      rate(2), interval_rhs(2), residual(2), ke(4)
    write(units(3), '(*(es25.16e3,:,","))') time, ke2d(1:4), rhs(3), &
      rate(3), interval_rhs(3), residual(3), ke2d(5:8)
    write(units(4), '(*(es25.16e3,:,","))') time, gmse(1:4), rhs(4), &
      rate(4), interval_rhs(4), residual(4), abs(residual(4)), gmse_relative, gmse(5), gmse(2)-gmse(5)
    previous_time = time
    previous_values = state
    previous_rhs = rhs
    has_previous = .true.
    elapsed = elapsed + MPI_WTIME() - started
  end subroutine

  subroutine finish_budgets
    double precision :: max_elapsed
    integer :: i, status

    call MPI_REDUCE(elapsed, max_elapsed, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, status)
    if (myid /= 0) return
    do i = 1, size(units)
      close(units(i))
    end do
    write(*,'(a,i0,a,es12.4,a)') 'Online budgets: ', samples, ' samples; max rank time ', max_elapsed, ' s'
  end subroutine
end module
