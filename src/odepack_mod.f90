module odepack_mod
  use iso_fortran_env, only: dp => real64
  use odepack_common, only: odepack_common_data
  use odepack_interface, only: dlsoda, dlsodar
  implicit none
  private

  public :: dp, lsoda_class

  !> Wrapper to the `DLSODA` and `DLSODAR` subroutines in
  !> ODEPACK. For additional details see documentation in 
  !> `odepack.f` 
  type :: lsoda_class

    ! functions
    integer :: neq !! number of ODEs
    procedure(lsoda_rhs_fcn), pointer :: f => NULL() !! right-hand-side of ODEs
    integer :: jt !! Jacobian type indicator
    procedure(lsoda_jac_fcn), pointer :: jac => NULL() !! jacobian of ODEs
    integer :: ng !! number of roots in g
    procedure(lsoda_root_fcn), pointer :: g => NULL() !! Root subroutine.
    integer, allocatable :: jroot(:) !! On a return with ISTATE = 3 (one or more roots found),
                                     !! JROOT(i) = 1 or -1 if g(i) has a root at T, or JROOT(i) = 0 if not.
                                     !! JROOT(i) = 1 means the value of gout(i) was increasing, and 
                                     !! JROOT(i) = -1 means the value of gout(i) was decreasing.
    real(dp), allocatable :: gout_input(:) !! a work array for root finding

    ! work memory
    integer :: lrw !! length of `rwork` array
    real(dp), allocatable :: rwork(:) !! real work array
    integer :: liw !! length of `iwork` array
    integer, allocatable :: iwork(:) !! integer work array

    !> Derived type that replaces the common blocks in 
    !> the original ODEPACK
    type(odepack_common_data), private :: common_data 
    
  contains
    procedure :: initialize => lsoda_initialize
    procedure :: integrate => lsoda_integrate
    procedure :: info => lsoda_info
  end type

  abstract interface

    !> Interface for the right-hand-side function defining the system
    !> of ODEs.
    subroutine lsoda_rhs_fcn(self, neq, t, y, ydot)
      import :: dp, lsoda_class
      implicit none
      class(lsoda_class), intent(inout) :: self
      integer, intent(in) :: neq !! number of ODEs
      real(dp), intent(in) :: t !! current time
      real(dp), intent(in) :: y(neq) !! state vector
      real(dp), intent(out) :: ydot(neq) !! derivative vector
    end subroutine

    !> Interface for the user-supplied jacobian function
    subroutine lsoda_jac_fcn(self, neq, t, y, ml, mu, pd, nrowpd)
      import :: dp, lsoda_class
      implicit none
      class(lsoda_class), intent(inout) :: self
      integer, intent(in) :: neq !! number of ODEs
      real(dp), intent(in) :: t !! current time
      real(dp), intent(inout) :: y(neq) !! state vector
      integer, intent(in) :: ml !! If the jacobian is banded (jt == 4)
                                !! then ml is the lower half-bandwidth. If the 
                                !! the jacobian is dense, then ignore this input.
      integer, intent(in) :: mu !! If the jacobian is banded (jt == 4)
                                !! then mu is the upper half-bandwidth. If the 
                                !! the jacobian is dense, then ignore this input.
      real(dp), intent(out) :: pd(nrowpd,neq)
      !! In the full matrix case (JT = 1), ML and MU are
      !! ignored, and the Jacobian is to be loaded into PD in
      !! columnwise manner, with df(i)/dy(j) loaded into PD(i,j).
      !!       In the band matrix case (JT = 4), the elements
      !! within the band are to be loaded into PD in columnwise
      !! manner, with diagonal lines of df/dy loaded into the rows
      !! of PD.  Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).
      !! ML and MU are the half-bandwidth parameters (see IWORK).
      !! The locations in PD in the two triangular areas which
      !! correspond to nonexistent matrix elements can be ignored
      !! or loaded arbitrarily, as they are overwritten by DLSODA.
      !!      JAC need not provide df/dy exactly.  A crude
      !! approximation (possibly with a smaller bandwidth) will do.
      !!      In either case, PD is preset to zero by the solver,
      !! so that only the nonzero elements need be loaded by JAC.
      !! Each call to JAC is preceded by a call to F with the same
      !! arguments NEQ, T, and Y.  Thus to gain some efficiency,
      !! intermediate quantities shared by both calculations may be
      !! saved in a user Common block by F and not recomputed by JAC,
      !! if desired.  Also, JAC may alter the Y array, if desired.
      integer, intent(in) :: nrowpd !! leading dimension of jacobian matrix
    end subroutine

    !> Interface for the root finding function
    subroutine lsoda_root_fcn(self, neq, t, y, ng, gout)
      import :: dp, lsoda_class
      implicit none
      class(lsoda_class), intent(inout) :: self
      integer, intent(in) :: neq !! number of ODEs
      real(dp), intent(in) :: t !! current time
      real(dp), intent(in) :: y(neq) !! state vector
      integer, intent(in) :: ng !! number of roots
      real(dp), intent(out) :: gout(ng) !! Roots. The program will find
                                        !! scenarios where gout(i) changes sign
    end subroutine
  end interface

contains

  !> Initializes the ODE integrator
  subroutine lsoda_initialize(self, f, neq, &
                              h0, hmax, hmin, ixpr, mxstep, mxhnil, mxordn, mxords, jac, jt, g, ng, &
                              ml, mu, &
                              istate)
    class(lsoda_class), intent(inout) :: self
    procedure(lsoda_rhs_fcn) :: f !! right-hand-side function defining the system of ODEs.
                                  !! See the `lsoda_rhs_fcn` interface for more information.
    integer, intent(in) :: neq !! number of ODEs
    ! optional
    real(dp), optional, intent(in) :: h0 !! The step size to be attempted on the first step.
                                         !! The default value is determined by the solver.
    real(dp), optional, intent(in) :: hmax !! The maximum absolute step size allowed.
                                           !! The default value is infinite.
    real(dp), optional, intent(in) :: hmin !! the minimum absolute step size allowed.
                                           !! The default value is 0.  (This lower bound is not
                                           !! enforced on the final step before reaching TCRIT
                                           !! when ITASK = 4 or 5.)
    integer, optional, intent(in) :: ixpr !! flag to generate extra printing at method switches.
                                          !! IXPR = 0 means no extra printing (the default).
                                          !! IXPR = 1 means print data on each switch.
                                          !! T, H, and NST will be printed on the same logical
                                          !! unit as used for error messages.
    integer, optional, intent(in) :: mxstep !! mxstep maximum number of (internally defined) steps
                                            !! allowed during one call to the solver.
                                            !! The default value is 10000.
    integer, optional, intent(in) :: mxhnil !! maximum number of messages printed (per problem)
                                            !! warning that T + H = T on a step (H = step size).
                                            !! This must be positive to result in a non-default
                                            !! value.  The default value is 10.
    integer, optional, intent(in) :: mxordn !! the maximum order to be allowed for the nonstiff
                                            !! (Adams) method.  The default value is 12.
                                            !! If MXORDN exceeds the default value, it will
                                            !! be reduced to the default value.
                                            !! MXORDN is held constant during the problem.
    integer, optional, intent(in) :: mxords !! the maximum order to be allowed for the stiff
                                            !! (BDF) method.  The default value is 5.
                                            !! If MXORDS exceeds the default value, it will
                                            !! be reduced to the default value.
                                            !! MXORDS is held constant during the problem.
    procedure(lsoda_jac_fcn), optional :: jac !! The user-supplied jaocobian matrix. See the
                                              !! `lsoda_jac_fcn` for more information.
    integer, optional, intent(in) :: jt !! Jacobian type indicator
                                        !!
                                        !! * 1 means a user-supplied full (NEQ by NEQ) Jacobian.
                                        !! * 2 means an internally generated (difference quotient) full
                                        !!     Jacobian (using NEQ extra calls to F per df/dy value).
                                        !! * 4 means a user-supplied banded Jacobian.
                                        !! * 5 means an internally generated banded Jacobian (using
                                        !!     ML+MU+1 extra calls to F per df/dy evaluation).
    procedure(lsoda_root_fcn), optional :: g !! The user-supplied root-finding subroutine. See
                                             !! the `lsoda_root_fcn` interface for more information.
    integer, optional, intent(in) :: ng !! Number of roots to be searched for with subroutine `g`
    integer, optional, intent(in) :: ml !! If the jacobian is banded (jt == 4 or 5), then ml is the lower half-bandwidth.
    integer, optional, intent(in) :: mu !! If the jacobian is banded (jt == 4 or 5), then ml is the lower upper-bandwidth.
    integer, intent(out) :: istate !! an index to specify success/failure. istate > 0 is success.

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
      if (allocated(self%jroot)) deallocate(self%jroot)
      allocate(self%jroot(ng))
      if (allocated(self%gout_input)) deallocate(self%gout_input)
      allocate(self%gout_input(ng))
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

  !> Integrates the ODEs forward in time from `t` until `tout`
  subroutine lsoda_integrate(self, y, t, tout, rtol, atol, itask, istate)
    class(lsoda_class), intent(inout) :: self
    real(dp), intent(inout) :: y(:)
    !! a real array for the vector of dependent variables, of
    !! length NEQ or more.  Used for both input and output on the
    !! first call (ISTATE = 1), and only for output on other calls.
    !! On the first call, Y must contain the vector of initial
    !! values.  On output, Y contains the computed solution vector,
    !! evaluated at T.  If desired, the Y array may be used
    !! for other purposes between calls to the solver.
    real(dp), intent(inout) :: t
    !! the independent variable.  On input, T is used only on the
    !! first call, as the initial point of the integration.
    !! On output, after each call, T is the value at which a
    !! computed solution y is evaluated (usually the same as TOUT).
    !! If a root was found, T is the computed location of the
    !! root reached first, on output.
    !! On an error return, T is the farthest point reached.
    real(dp), intent(in) :: tout
    !! the next value of t at which a computed solution is desired.
    !! Used only for input.
    !! 
    !! When starting the problem (ISTATE = 1), TOUT may be equal
    !! to T for one call, then should .ne. T for the next call.
    !! For the initial T, an input value of TOUT .ne. T is used
    !! in order to determine the direction of the integration
    !! (i.e. the algebraic sign of the step sizes) and the rough
    !! scale of the problem.  Integration in either direction
    !! (forward or backward in t) is permitted.
    !! 
    !! If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
    !! the first call (i.e. the first call with TOUT .ne. T).
    !! Otherwise, TOUT is required on every call.
    !! 
    !! If ITASK = 1, 3, or 4, the values of TOUT need not be
    !! monotone, but a value of TOUT which backs up is limited
    !! to the current internal T interval, whose endpoints are
    !! TCUR - HU and TCUR (see optional outputs, below, for
    !! TCUR and HU).
    real(dp), intent(in) :: rtol !! The relative error tolerance.
    real(dp), intent(in) :: atol(:) 
    !! The abolute error tolerance. Can be a length 1 or length `neq`.
    !! If it is length `neq`, then a separate tolerance is applied to
    !! each evolving variable `y(i)`. I length 1, then the same
    !! abosoluate tolerance is applied to every variable
    integer, intent(in) :: itask
    !! an index specifying the task to be performed.
    !! input only.  ITASK has the following values and meanings.
    !!
    !! * 1  means normal computation of output values of y(t) at
    !!      t = TOUT (by overshooting and interpolating).
    !! * 2  means take one step only and return.
    !! * 3  means stop at the first internal mesh point at or
    !!      beyond t = TOUT and return.
    !! * 4  means normal computation of output values of y(t) at
    !!      t = TOUT but without overshooting t = TCRIT.
    !!      TCRIT must be input as RWORK(1).  TCRIT may be equal to
    !!      or beyond TOUT, but not behind it in the direction of
    !!      integration.  This option is useful if the problem
    !!      has a singularity at or beyond t = TCRIT.
    !! * 5  means take one step, without passing TCRIT, and return.
    !!      TCRIT must be input as RWORK(1).
    !! 
    !! Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
    !! (within roundoff), it will return T = TCRIT (exactly) to
    !! indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
    !! in which case answers at t = TOUT are returned first).
    integer, intent(inout) :: istate
    !! an index used for input and output to specify the
    !! the state of the calculation.
    !!
    !! On input, the values of ISTATE are as follows.
    !!
    !! * 1  means this is the first call for the problem
    !!      (initializations will be done).  See note below.
    !! * 2  means this is not the first call, and the calculation
    !!      is to continue normally, with no change in any input
    !!      parameters except possibly TOUT and ITASK.
    !!      (If ITOL, RTOL, and/or ATOL are changed between calls
    !!      with ISTATE = 2, the new values will be used but not
    !!      tested for legality.)
    !! * 3  means this is not the first call, and the
    !!      calculation is to continue normally, but with
    !!      a change in input parameters other than
    !!      TOUT and ITASK.  Changes are allowed in
    !!      NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, JT, ML, MU,
    !!      and any optional inputs except H0, MXORDN, and MXORDS.
    !!      (See IWORK description for ML and MU.)
    !!      In addition, immediately following a return with
    !!      ISTATE = 3 (root found), NG and G may be changed.
    !!      (But changing NG from 0 to .gt. 0 is not allowed.)
    !!
    !! Note:  A preliminary call with TOUT = T is not counted
    !! as a first call here, as no initialization or checking of
    !! input is done.  (Such a call is sometimes useful for the
    !! purpose of outputting the initial conditions.)
    !! Thus the first call for which TOUT .ne. T requires
    !! ISTATE = 1 on input.
    !! 
    !! On output, ISTATE has the following values and meanings.
    !!
    !! * 1  means nothing was done; TOUT = t and ISTATE = 1 on input.
    !! * 2  means the integration was performed successfully, and
    !!      no roots were found.
    !! * 3  means the integration was successful, and one or more
    !!      roots were found before satisfying the stop condition
    !!      specified by ITASK.  See JROOT.
    !! * -1  means an excessive amount of work (more than MXSTEP
    !!       steps) was done on this call, before completing the
    !!       requested task, but the integration was otherwise
    !!       successful as far as T.  (MXSTEP is an optional input
    !!       and is normally 500.)  To continue, the user may
    !!       simply reset ISTATE to a value .gt. 1 and call again
    !!       (the excess work step counter will be reset to 0).
    !!       In addition, the user may increase MXSTEP to avoid
    !!       this error return (see below on optional inputs).
    !! * -2  means too much accuracy was requested for the precision
    !!       of the machine being used.  This was detected before
    !!       completing the requested task, but the integration
    !!       was successful as far as T.  To continue, the tolerance
    !!       parameters must be reset, and ISTATE must be set
    !!       to 3.  The optional output TOLSF may be used for this
    !!       purpose.  (Note: If this condition is detected before
    !!       taking any steps, then an illegal input return
    !!       (ISTATE = -3) occurs instead.)
    !! * -3  means illegal input was detected, before taking any
    !!       integration steps.  See written message for details.
    !!       Note:  If the solver detects an infinite loop of calls
    !!       to the solver with illegal input, it will cause
    !!       the run to stop.
    !! * -4  means there were repeated error test failures on
    !!       one attempted step, before completing the requested
    !!       task, but the integration was successful as far as T.
    !!       The problem may have a singularity, or the input
    !!       may be inappropriate.
    !! * -5  means there were repeated convergence test failures on
    !!       one attempted step, before completing the requested
    !!       task, but the integration was successful as far as T.
    !!       This may be caused by an inaccurate Jacobian matrix,
    !!       if one is being used.
    !! * -6  means EWT(i) became zero for some i during the
    !!       integration.  Pure relative error control (ATOL(i)=0.0)
    !!       was requested on a variable which has now vanished.
    !!       The integration was successful as far as T.
    !! * -7  means the length of RWORK and/or IWORK was too small to
    !!       proceed, but the integration was successful as far as T.
    !!       This happens when DLSODAR chooses to switch methods
    !!       but LRW and/or LIW is too small for the new method.
    !! 
    !! Note:  Since the normal output value of ISTATE is 2,
    !! it does not need to be reset for normal continuation.
    !! Also, since a negative input value of ISTATE will be
    !! regarded as illegal, a negative output value requires the
    !! user to change it, and possibly other inputs, before
    !! calling the solver again.

    integer :: i
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
      ! we compute the roots at input
      call self%g(self%neq, t, y, self%ng, self%gout_input)
      call dlsodar(f, self%neq, y, t, tout, itol, rtol, atol, itask, &
                   istate, iopt, self%rwork, self%lrw, self%iwork, self%liw, jac, self%jt, &
                   g, self%ng, self%jroot, self%common_data)
      if (istate == 3) then
        ! a root was found
        do i = 1,self%ng
          if (self%jroot(i) == 1) then
            ! if the root was descending, then replace with -1
            if (self%gout_input(i) > 0.0_dp) then
              self%jroot(i) = -1
            endif
          endif
        enddo
      endif
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
      real(dp), intent(inout) :: y_(neq_)
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

  !> Get information about the integration.
  subroutine lsoda_info(self, hu, hcur, tcur, tolsf, tsw, nst, &
                        nfe, nje, nqu, nqcur, imxer, mused, mcur)
    class(lsoda_class), intent(inout) :: self
    real(dp), optional, intent(out) :: hu !! the step size in t last used (successfully).
    real(dp), optional, intent(out) :: hcur !! the step size to be attempted on the next step.
    real(dp), optional, intent(out) :: tcur !! the current value of the independent variable
                                            !! which the solver has actually reached, i.e. the
                                            !! current internal mesh point in t.  On output, TCUR
                                            !! will always be at least as far as the argument
                                            !! T, but may be farther (if interpolation was done).
    real(dp), optional, intent(out) :: tolsf !! a tolerance scale factor, greater than 1.0,
                                             !! computed when a request for too much accuracy was
                                             !! detected (ISTATE = -3 if detected at the start of
                                             !! the problem, ISTATE = -2 otherwise).  If ITOL is
                                             !! left unaltered but RTOL and ATOL are uniformly
                                             !! scaled up by a factor of TOLSF for the next call,
                                             !! then the solver is deemed likely to succeed.
                                             !! (The user may also ignore TOLSF and alter the
                                             !! tolerance parameters in any other way appropriate.)
    real(dp), optional, intent(out) :: tsw !! the value of t at the time of the last method 
                                           !! switch, if any.
    integer, optional, intent(out) :: nst !! the number of steps taken for the problem so far.
    integer, optional, intent(out) :: nfe !! the number of f evaluations for the problem so far.
    integer, optional, intent(out) :: nje !! the number of Jacobian evaluations (and of matrix
                                          !! LU decompositions) for the problem so far.
    integer, optional, intent(out) :: nqu !! the method order last used (successfully).
    integer, optional, intent(out) :: nqcur !! the order to be attempted on the next step.
    integer, optional, intent(out) :: imxer !! the index of the component of largest magnitude in
                                            !! the weighted local error vector ( E(i)/EWT(i) ),
                                            !! on an error return with ISTATE = -4 or -5.
    integer, optional, intent(out) :: mused !! the method indicator for the last successful step:
                                            !! 1 means Adams (nonstiff), 2 means BDF (stiff).
    integer, optional, intent(out) :: mcur !! the current method indicator:
                                           !! 1 means Adams (nonstiff), 2 means BDF (stiff).
                                           !! This is the method to be attempted
                                           !! on the next step.  Thus it differs from MUSED
                                           !! only if a method switch has just been made.

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