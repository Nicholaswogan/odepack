module odepack_mod
  use iso_fortran_env, only: dp => real64
  use odepack_common, only: odepack_common_data
  use odepack_interface, only: dlsoda, dlsodar
  implicit none
  private

  public :: dp, lsoda_class

  type :: lsoda_class

    ! functions
    integer :: neq !! number of ODEs
    procedure(lsoda_rhs_fcn), pointer :: f => NULL() !! right-hand-side of ODEs
    integer :: jt !! Jacobian type indicator
                  !! 1 means a user-supplied full (NEQ by NEQ) Jacobian.
                  !! 2 means an internally generated (difference quotient) full
                  !!   Jacobian (using NEQ extra calls to F per df/dy value).
                  !! 4 means a user-supplied banded Jacobian.
                  !! 5 means an internally generated banded Jacobian (using
                  !!   ML+MU+1 extra calls to F per df/dy evaluation).
    procedure(lsoda_jac_fcn), pointer :: jac => NULL() !! jacobian of ODEs
    integer :: ng !! number of roots in g
    procedure(lsoda_root_fcn), pointer :: g => NULL() !! Root subroutine.
    integer, allocatable :: jroot(:) !! On a return with ISTATE = 3 (one or more roots found),
                                     !! JROOT(i) = 1 if g(i) has a root at T, or JROOT(i) = 0 if not.
    
    ! work memory
    integer :: lrw
    real(dp), allocatable :: rwork(:)
    integer :: liw
    integer, allocatable :: iwork(:)

    ! common block data
    type(odepack_common_data), private :: common_data
    
  contains
    procedure :: initialize => lsoda_initialize
    procedure :: integrate => lsoda_integrate
    procedure :: info => lsoda_info
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
                              ml, mu, &
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
    integer, optional, intent(in) :: ml
    integer, optional, intent(in) :: mu
    integer, intent(out) :: istate

    integer :: ml_, mu_

    istate = 1

    self%f => f
    self%neq = neq

    ! jacobian stuff
    if (present(ml) .and. present(mu)) then
      ml_ = ml
      mu_ = mu
    elseif (.not.present(ml) .and. .not.present(mu)) then
      ! nothing
    else
      ! err = '"ml" and "mu" most both be inputs'
      istate = -3
      return
    endif
    if (present(jt)) then
      if (jt == 1 .or. jt == 4) then
        if (.not.present(jac)) then
          ! err = 'if jt is 1 or 4, then jac must always be present'
          istate = -3
          return
        endif
      endif
      if (jt == 4 .or. jt == 5) then
        if (.not.present(ml)) then
          ! err = '"jt" is 4 or 5, therefore, "ml" and "mu" must be inputs.'
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
    if (present(ng)) then
      if (.not.present(g)) then
        istate = -3
        return
      endif
    endif
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

    ! allocate memory
    self%lrw = 22 + neq * max(16, neq + 9) + 3*self%ng
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
    if (present(mu)) then
      self%iwork(1) = ml_
      self%iwork(2) = mu_
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
    integer, intent(inout) :: istate

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
                   istate, iopt, self%rwork, self%lrw, self%iwork, self%liw, jac, self%jt, &
                   g, self%ng, self%jroot, self%common_data)
    else
      call dlsoda(f, self%neq, y, t, tout, itol, rtol, atol, itask, &
                  istate, iopt, self%rwork, self%lrw, self%iwork, self%liw, jac, self%jt, &
                  self%common_data)
    endif

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

  subroutine lsoda_info(self, hu, hcur, tcur, tolsf, tsw, nst, &
                        nfe, nje, nqu, nqcur, imxer, mused, mcur)
    class(lsoda_class), intent(inout) :: self
    real(dp), optional, intent(out) :: hu
    real(dp), optional, intent(out) :: hcur
    real(dp), optional, intent(out) :: tcur
    real(dp), optional, intent(out) :: tolsf
    real(dp), optional, intent(out) :: tsw
    integer, optional, intent(out) :: nst
    integer, optional, intent(out) :: nfe
    integer, optional, intent(out) :: nje
    integer, optional, intent(out) :: nqu
    integer, optional, intent(out) :: nqcur
    integer, optional, intent(out) :: imxer
    integer, optional, intent(out) :: mused
    integer, optional, intent(out) :: mcur

    if (present(hu)) then
      hu = self%rwork(11)
    endif
    if (present(hcur)) then
      hcur = self%rwork(12)
    endif
    if (present(tcur)) then
      tcur = self%rwork(13)
    endif
    if (present(tolsf)) then
      tolsf = self%rwork(14)
    endif
    if (present(tsw)) then
      tsw = self%rwork(15)
    endif

    if (present(nst)) then
      nst = self%iwork(11)
    endif
    if (present(nfe)) then
      nfe = self%iwork(12)
    endif
    if (present(nje)) then
      nje = self%iwork(13)
    endif
    if (present(nqu)) then
      nqu= self%iwork(14)
    endif
    if (present(nqcur)) then
      nqcur = self%iwork(15)
    endif
    if (present(imxer)) then
      imxer = self%iwork(16)
    endif
    if (present(mused)) then
      mused = self%iwork(19)
    endif
    if (present(mcur)) then
      mcur = self%iwork(20)
    endif

  end subroutine

end module