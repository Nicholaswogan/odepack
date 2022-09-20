program test_lsoda
  use odepack_mod
  implicit none
  call test_simple()
  call test_rootfinding()
  call test_dense_jacobian()
contains
  subroutine test_simple()
    type(lsoda_class) :: ls
    integer :: neq, itask, istate
    real(dp) :: y(2), t, tout, rtol, atol(1)
    real(dp) :: correct_y(2) = [0.38583246250193476_dp, 4.602012234037773_dp]

    neq = 2
    call ls%initialize(rhs, neq, istate=istate)
    if (istate < 0) then
      print*,istate
      error stop '"test_simple" failed'
    endif

    y(:) = [5.0_dp, 0.8_dp]
    t = 0.0_dp
    tout = 10.0_dp
    rtol = 1.0e-8_dp
    atol = 1.0e-8_dp
    itask = 1
    istate = 1
    call ls%integrate(y, t, tout, rtol, atol, itask, istate)
    if (istate < 0) then
      print*,istate
      error stop '"test_simple" failed'
    endif

    if (.not.all(is_close(correct_y, y))) then
      print*,correct_y
      print*,y
      error stop '"test_simple" failed'
    endif

    print*,'"test_simple" passed'
  end subroutine

  subroutine test_rootfinding()
    type(lsoda_class) :: ls
    integer :: neq, itask, istate
    real(dp) :: y(2), t, tout, rtol, atol(1)

    neq = 2
    call ls%initialize(rhs, neq, g=root_fcn, ng=1, istate=istate)
    if (istate < 0) then
      print*,istate
      error stop '"test_rootfinding" failed'
    endif

    y(:) = [5.0_dp, 0.8_dp]
    t = 0.0_dp
    tout = 10.0_dp
    rtol = 1.0e-8_dp
    atol = 1.0e-8_dp
    itask = 1
    istate = 1
    call ls%integrate(y, t, tout, rtol, atol, itask, istate)
    if (istate < 0) then
      print*,istate
      error stop '"test_rootfinding" failed'
    endif

    if (istate /= 3 .or. .not.is_close(y(1),1.0_dp) .or. ls%jroot(1) /= -1) then
      error stop '"test_rootfinding" failed'
    endif

    print*,'"test_rootfinding" passed'
  end subroutine

  subroutine test_dense_jacobian()
    type(lsoda_class) :: ls
    integer :: neq, itask, istate, nje
    real(dp) :: y(3), t, tout, rtol, atol(1)

    neq = 3
    rtol = 1.0e-8_dp
    atol = 1.0e-8_dp
    itask = 1
    
    ! with finite differenced jacobian
    call ls%initialize(rhs_rober, neq, jt=2, istate=istate)
    if (istate < 0) then
      print*,istate
      error stop '"test_dense_jacobian" failed'
    endif

    y(:) = [1.0_dp,0.0_dp,0.0_dp]
    t = 0.0_dp
    tout = 1.0e5_dp
    istate = 1
    call ls%integrate(y, t, tout, rtol, atol, itask, istate)
    if (istate < 0) then
      print*,istate
      error stop '"test_dense_jacobian" failed'
    endif

    ! with analytical jacobian
    call ls%initialize(rhs_rober, neq, jt=1, jac=jac_rober, istate=istate)
    if (istate < 0) then
      print*,istate
      error stop '"test_dense_jacobian" failed'
    endif

    y(:) = [1.0_dp,0.0_dp,0.0_dp]
    t = 0.0_dp
    tout = 1.0e5_dp
    istate = 1
    call ls%integrate(y, t, tout, rtol, atol, itask, istate)
    if (istate < 0) then
      print*,istate
      error stop '"test_dense_jacobian" failed'
    endif

    print*,'"test_dense_jacobian" passed'
  end subroutine

  subroutine jac_rober(self, neq, t, u, ml, mu, pd, nrpd, ierr)
    class(lsoda_class), intent(inout) :: self
    integer, intent(in) :: neq
    real(dp), intent(in) :: t
    real(dp), intent(inout) :: u(neq)
    integer, intent(in) :: ml, mu
    real(dp), intent(out) :: pd(nrpd,neq)
    integer, intent(in) :: nrpd
    integer, intent(out) :: ierr
    
    real(dp), parameter :: k1 = 0.04_dp, &
                           k2 = 3.0e7_dp, &
                           k3 = 1.0e4_dp
    
    pd(:,1) = [-k1, k1, 0.0_dp]
    pd(:,2) = [k3*u(3), -2.0_dp*k2*u(2) - k3*u(3), 2.0_dp*k2*u(2)]
    pd(:,3) = [k3*u(2), -k3*u(2), 0.0_dp]
    
    ierr = 0
  end subroutine

  subroutine rhs_rober(self, neq, t, u, du, ierr)
    class(lsoda_class), intent(inout) :: self
    integer, intent(in) :: neq
    real(dp), intent(in) :: t
    real(dp), intent(in) :: u(neq)
    real(dp), intent(out) :: du(neq)
    integer, intent(out) :: ierr

    real(dp), parameter :: k1 = 0.04_dp, &
                           k2 = 3.0e7_dp, &
                           k3 = 1.0e4_dp
    
    du(1) = -k1*u(1) + k3*u(2)*u(3)
    du(2) =  k1*u(1) - k2*u(2)**2.0_dp - k3*u(2)*u(3)
    du(3) =  k2*u(2)**2.0_dp

    ierr = 0
  end subroutine

  subroutine rhs(self, neq, t, y, ydot, ierr)
    class(lsoda_class), intent(inout) :: self
    integer, intent(in) :: neq
    real(dp), intent(in) :: t
    real(dp), intent(in) :: y(neq)
    real(dp), intent(out) :: ydot(neq)
    integer, intent(out) :: ierr
    ydot(1) = y(1)-y(1)*y(2)
    ydot(2) = y(1)*y(2)-y(2)
    ierr = 0
  end subroutine

  subroutine root_fcn(self, neq, t, y, ng, gout, ierr)
    class(lsoda_class), intent(inout) :: self
    integer, intent(in) :: neq
    real(dp), intent(in) :: t
    real(dp), intent(in) :: y(neq)
    integer, intent(in) :: ng
    real(dp), intent(out) :: gout(ng)
    integer, intent(out) :: ierr
    gout(1) = y(1) - 1.0_dp
    ierr = 0
  end subroutine

  !> coppied from fortran stdlib v0.2.0
  elemental function is_close(a, b, tol, abs_tol, equal_nan) result(close)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    real(dp), intent(in) :: a, b
    real(dp), intent(in), optional :: tol, abs_tol
    logical, intent(in), optional :: equal_nan
    logical :: close

    real(dp) :: rel_tol_, abs_tol_
    logical :: equal_nan_

    if (present(tol)) then
      rel_tol_ = tol
    else
      rel_tol_ = 1.0e-5_dp
    endif

    if (present(abs_tol)) then
      abs_tol_ = abs_tol
    else
      abs_tol_ = 0.0_dp
    endif

    if (present(equal_nan)) then
      equal_nan_ = equal_nan
    else
      equal_nan_ = .false.
    endif

    if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
        close = merge(.true., .false., equal_nan_ .and. ieee_is_nan(a) .and. ieee_is_nan(b))
    else
        close = abs(a - b) <= max(abs(rel_tol_*max(abs(a), abs(b))), abs(abs_tol_))
    end if     

  end function
end program