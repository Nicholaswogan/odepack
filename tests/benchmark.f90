program main
  use odepack_mod
  implicit none
  call benchmark_lorenz
  call benchmark_rober()
contains

  subroutine benchmark_lorenz()
    type(lsoda_class) :: ls
    integer :: neq, itask, istate
    real(dp) :: y(3), t, tout, rtol, atol(1)
    real(dp) :: t_eval(1001)
    real(dp) :: ysol(3,1001)
    integer :: i, j
    integer, parameter :: nruns = 500
    real(dp) :: time(2)

    neq = 3
    call ls%initialize(rhs_lorenz, neq, istate=istate)
    if (istate < 0) then
      print*,istate
      error stop
    endif

    call linspace(0.0_dp, 100.0_dp, t_eval)
    rtol = 1.0e-8_dp
    atol = 1.0e-8_dp
    itask = 1

    call cpu_time(time(1))
    do j = 1,nruns
      istate = 1
      y(:) = [1.0_dp,0.0_dp,0.0_dp]
      t = 0.0_dp
      do i = 1,size(t_eval)
        tout = t_eval(i)
        call ls%integrate(y, t, tout, rtol, atol, itask, istate)
        if (istate < 0) then
          print*,istate
          error stop
        endif
        ysol(:,i) = y(:)
      enddo
    enddo
    call cpu_time(time(2))

    print*,'lorenz: ',1000.0_dp*(time(2) - time(1))/real(nruns,dp),'ms'

  end subroutine

  subroutine benchmark_rober()
    type(lsoda_class) :: ls
    integer :: neq, itask, istate
    real(dp) :: y(3), t, tout, rtol, atol(1)
    real(dp) :: t_eval(100)
    real(dp) :: ysol(3,100)
    integer :: i, j
    integer, parameter :: nruns = 5000
    real(dp) :: time(2)

    neq = 3
    call ls%initialize(rhs_rober, neq, istate=istate)
    if (istate < 0) then
      print*,istate
      error stop
    endif

    call linspace(0.0_dp, 1.0e5_dp, t_eval)
    rtol = 1.0e-8_dp
    atol = 1.0e-8_dp
    itask = 1

    call cpu_time(time(1))
    do j = 1,nruns
      istate = 1
      y(:) = [1.0_dp,0.0_dp,0.0_dp]
      t = 0.0_dp
      do i = 1,size(t_eval)
        tout = t_eval(i)
        call ls%integrate(y, t, tout, rtol, atol, itask, istate)
        if (istate < 0) then
          print*,istate
          error stop
        endif
        ysol(:,i) = y(:)
      enddo
      stop
    enddo
    call cpu_time(time(2))

    print*,'rober: ',1000.0_dp*(time(2) - time(1))/real(nruns,dp),'ms'
  end subroutine

  subroutine rhs_lorenz(self, neq, t, u, du)
    class(lsoda_class), intent(inout) :: self
    integer, intent(in) :: neq
    real(dp), intent(in) :: t
    real(dp), intent(in) :: u(neq)
    real(dp), intent(out) :: du(neq)
    
    real(dp) :: x, y, z
    real(dp), parameter :: sigma = 10.0_dp, &
                           rho = 28.0_dp, &
                           beta = 8.0_dp/3.0_dp
    
    x = u(1)
    y = u(2)
    z = u(3)
    du(1) = sigma * (y - x)
    du(2) = x * (rho - z) - y
    du(3) = x * y - beta * z

  end subroutine

  subroutine rhs_rober(self, neq, t, u, du)
    class(lsoda_class), intent(inout) :: self
    integer, intent(in) :: neq
    real(dp), intent(in) :: t
    real(dp), intent(in) :: u(neq)
    real(dp), intent(out) :: du(neq)
    
    real(dp), parameter :: k1 = 0.04_dp, &
                           k2 = 3.0e7_dp, &
                           k3 = 1.0e4_dp
    
    du(1) = -k1*u(1) + k3*u(2)*u(3)
    du(2) =  k1*u(1) - k2*u(2)**2.0_dp - k3*u(2)*u(3)
    du(3) =  k2*u(2)**2.0_dp

  end subroutine

  subroutine linspace(from, to, array)
    real(dp), intent(in) :: from, to
    real(dp), intent(out) :: array(:)
    real(dp) :: range
    integer :: n, i
    n = size(array)
    range = to - from

    if (n == 0) return

    if (n == 1) then
      array(1) = from
      return
    end if

    do i=1, n
      array(i) = from + range * (i - 1) / (n - 1)
    end do
  end subroutine

end program