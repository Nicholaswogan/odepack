module odepack_mod
  use iso_fortran_env, only: dp => real64
  use odepack_common, only: odepack_common_data
  use odepack_interface, only: dlsoda, dlsodar
  implicit none
  private

  public :: dp, lsoda_class

  type :: lsoda_class
    integer :: neq
    procedure(lsoda_rhs_fcn), pointer :: f => NULL()
    integer :: jt
    procedure(lsoda_jac_fcn), pointer :: jac => NULL()
    integer :: ng
    procedure(lsoda_root_fcn), pointer :: g => NULL()
    integer, allocatable :: jroot(:)
    
    ! work memory
    integer :: lrw
    real(dp), allocatable :: rwork(:)
    integer :: liw
    integer, allocatable :: iwork(:)
    integer :: istate = 1

    ! common block data
    type(odepack_common_data) :: common_data
    
  contains
    procedure :: initialize => lsoda_initialize
    procedure :: integrate => lsoda_integrate
    ! procedure :: info => lsoda_get_info
  end type

  abstract interface
    subroutine lsoda_rhs_fcn(self, neq, t, y, ydot)
      import :: dp, lsoda_class
      implicit none
      class(lsoda_class), intent(inout) :: self
      integer, intent(in) :: neq
      real(dp), intent(in) :: t
      real(dp), intent(in) :: y(neq)
      real(dp), intent(out) :: ydot(neq)
    end subroutine

    subroutine lsoda_jac_fcn(self, neq, t, y, ml, mu, pd, nrowpd)
      import :: dp, lsoda_class
      implicit none
      class(lsoda_class), intent(inout) :: self
      integer, intent(in) :: neq
      real(dp), intent(in) :: t
      real(dp), intent(in) :: y(neq)
      integer, intent(in) :: ml
      integer, intent(in) :: mu
      real(dp), intent(out) :: pd(nrowpd,neq)
      integer, intent(in) :: nrowpd
    end subroutine

    subroutine lsoda_root_fcn(self, neq, t, y, ng, gout)
      import :: dp, lsoda_class
      implicit none
      class(lsoda_class), intent(inout) :: self
      integer, intent(in) :: neq
      real(dp), intent(in) :: t
      real(dp), intent(in) :: y(neq)
      integer, intent(in) :: ng
      real(dp), intent(out) :: gout(ng)
    end subroutine
  end interface

contains

  subroutine lsoda_initialize(self, f, neq, &
                              h0, hmax, hmin, ixpr, mxstep, mxhnil, mxordn, mxords, jac, jt, g, ng, &
                              istate)
    class(lsoda_class), intent(inout) :: self
    procedure(lsoda_rhs_fcn) :: f
    integer, intent(in) :: neq
    ! optional
    real(dp), optional, intent(in) :: h0
    real(dp), optional, intent(in) :: hmax
    real(dp), optional, intent(in) :: hmin
    integer, optional, intent(in) :: ixpr
    integer, optional, intent(in) :: mxstep
    integer, optional, intent(in) :: mxhnil
    integer, optional, intent(in) :: mxordn
    integer, optional, intent(in) :: mxords
    procedure(lsoda_jac_fcn), optional :: jac
    integer, optional, intent(in) :: jt
    procedure(lsoda_root_fcn), optional :: g
    integer, optional, intent(in) :: ng

    integer, intent(out) :: istate

    istate = 1

    self%f => f
    self%neq = neq

    ! allocate memory
    self%lrw = 22 + neq * max(16, neq + 9)
    if (allocated(self%rwork)) deallocate(self%rwork)
    allocate(self%rwork(self%lrw))
    self%liw = 20 + neq
    if (allocated(self%iwork)) deallocate(self%iwork)
    allocate(self%iwork(self%liw))

    ! optional inputs
    if (present(h0)) then
      self%rwork(5) = h0
    else
      self%rwork(5) = 0.0_dp
    endif
    if (present(hmax)) then
      self%rwork(6) = hmax
    else
      self%rwork(6) = huge(1.0_dp)
    endif
    if (present(hmin)) then
      self%rwork(7) = hmin
    else
      self%rwork(7) = 0.0_dp
    endif
    if (present(ixpr)) then
      self%iwork(5) = ixpr
    else
      self%iwork(5) = 0
    endif
    if (present(mxstep)) then
      self%iwork(6) = mxstep
    else
      self%iwork(6) = 10000
    endif
    if (present(mxhnil)) then
      self%iwork(7) = mxhnil
    else
      self%iwork(7) = 10
    endif
    if (present(mxordn)) then
      self%iwork(8) = mxordn
    else
      self%iwork(8) = 12
    endif
    if (present(mxords)) then
      self%iwork(9) = mxords
    else
      self%iwork(9) = 5
    endif

    self%istate = 1
    
    ! jacobian stuff
    if (present(jt)) then
      if (jt == 1 .or. jt == 4) then
        if (.not.present(jac)) then
          ! err = 'if jt is 1 or 4, then jac must always be present'
          istate = -3
          return
        endif
      endif
      self%jt = jt
    else
      self%jt = 2
    endif
    if (present(jac)) then
      if (.not.present(jt)) then
        ! err = 'if jac is present, then jt must always be present'
        istate = -3
        return
      endif
      self%jac => jac
    endif

    ! root function
    if (present(g)) then
      if (.not.present(ng)) then
        istate = -3
        return
      endif
      self%g => g
      self%ng = ng
      allocate(self%jroot(ng))
    else
      self%ng = 0
    endif

  end subroutine

  subroutine lsoda_integrate(self, y, t, tout, rtol, atol, itask, istate)
    class(lsoda_class), intent(inout) :: self
    real(dp), intent(inout) :: y(:)
    real(dp), intent(inout) :: t
    real(dp), intent(in) :: tout
    real(dp), intent(in) :: rtol
    real(dp), intent(in) :: atol(:)
    integer, intent(in) :: itask
    integer, intent(out) :: istate

    integer :: itol
    integer, parameter :: iopt = 1

    ! check dimensions
    if (size(y) /= self%neq) then
      ! err = 'lsoda_integrate: "y" has the wrong dimension.'
      istate = -3
      return
    endif

    if (size(atol) == 1) then
      itol = 1
    elseif (size(atol) == self%neq) then
      itol = 2
    else
      ! err = 'lsoda_integrate: "atol" can be size 1 or size neq.'
      istate = -3
      return
    endif

    if (associated(self%g)) then
      call dlsodar(f, self%neq, y, t, tout, itol, rtol, atol, itask, &
                   self%istate, iopt, self%rwork, self%lrw, self%iwork, self%liw, jac, self%jt, &
                   g, self%ng, self%jroot, self%common_data)
    else
      call dlsoda(f, self%neq, y, t, tout, itol, rtol, atol, itask, &
                  self%istate, iopt, self%rwork, self%lrw, self%iwork, self%liw, jac, self%jt, &
                  self%common_data)
    endif
    istate = self%istate

  contains
    subroutine f(neq_, t_, y_, ydot_)
      integer, intent(in) :: neq_
      real(dp), intent(in) :: t_
      real(dp), intent(in) :: y_(neq_)
      real(dp), intent(out) :: ydot_(neq_)
      call self%f(neq_, t_, y_, ydot_)
    end subroutine
    
    subroutine jac(neq_, t_, y_, ml_, mu_, pd_, nrowpd_)
      integer, intent(in) :: neq_
      real(dp), intent(in) :: t_
      real(dp), intent(in) :: y_(neq_)
      integer, intent(in) :: ml_
      integer, intent(in) :: mu_
      real(dp), intent(out) :: pd_(nrowpd_,neq_)
      integer, intent(in) :: nrowpd_
      call self%jac(neq_, t_, y_, ml_, mu_, pd_, nrowpd_)
    end subroutine

    subroutine g(neq_, t_, y_, ng_, gout_)
      integer, intent(in) :: neq_
      real(dp), intent(in) :: t_
      real(dp), intent(in) :: y_(neq_)
      integer, intent(in) :: ng_
      real(dp), intent(out) :: gout_(ng_)
      call self%g(neq_, t_, y_, ng_, gout_)
    end subroutine
  end subroutine

end module