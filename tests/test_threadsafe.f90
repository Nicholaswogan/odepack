program main
  use odepack_mod
  implicit none
  call test_threadsafe()
contains

  subroutine test_threadsafe()
    integer, parameter :: n = 100
    real(dp) :: y1(2,n)
    real(dp) :: y2(2,n)
    integer :: i

    !$omp parallel private(i)
    !$omp do
    do i = 1,n
      call test_threadsafe_helper(y1(:,i))
    enddo
    !$omp enddo
    !$omp end parallel

    do i = 1,n
      call test_threadsafe_helper(y2(:,i))
    enddo

    do i = 1,n
      if (any(y1(:,i) /= y2(:,i))) then
        print*,y1(:,i)
        print*,y2(:,i)
        error stop '"test_threadsafe" failed'
      endif
    enddo

    print*,'"test_threadsafe" passed'
  end subroutine

  subroutine test_threadsafe_helper(y)
    real(dp), intent(out) :: y(2)

    type(lsoda_class) :: ls
    integer :: neq, itask, istate
    real(dp) :: t, tout, rtol, atol(1)

    neq = 2
    call ls%initialize(rhs, neq, istate=istate)
    if (istate < 0) then
      print*,istate
      error stop '"test_threadsafe" failed'
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
      error stop '"test_threadsafe" failed'
    endif

  end subroutine

  subroutine rhs(self, neq, t, y, ydot)
    class(lsoda_class), intent(inout) :: self
    integer, intent(in) :: neq
    real(dp), intent(in) :: t
    real(dp), intent(in) :: y(neq)
    real(dp), intent(out) :: ydot(neq)
    ydot(1) = y(1)-y(1)*y(2)
    ydot(2) = y(1)*y(2)-y(2)
  end subroutine

end program