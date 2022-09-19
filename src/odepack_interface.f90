module odepack_interface
  use iso_fortran_env, only: dp => real64
  use odepack_common, only: odepack_common_data
  implicit none
  private

  public :: dlsoda, dlsodar
  public :: DINTDY, DROOTS

  abstract interface
    subroutine odepack_f(neq, t, y, ydot)
      import :: dp
      implicit none
      integer, intent(in) :: neq
      real(dp), intent(in) :: t
      real(dp), intent(in) :: y(neq)
      real(dp), intent(out) :: ydot(neq)
    end subroutine
    
    subroutine odepack_jac(neq, t, y, ml, mu, pd, nrowpd)
      import :: dp
      implicit none
      integer, intent(in) :: neq
      real(dp), intent(in) :: t
      real(dp), intent(inout) :: y(neq)
      integer, intent(in) :: ml
      integer, intent(in) :: mu
      real(dp), intent(out) :: pd(nrowpd,neq)
      integer, intent(in) :: nrowpd
    end subroutine

    subroutine odepack_g(neq, t, y, ng, gout)
      import :: dp
      implicit none
      integer, intent(in) :: neq
      real(dp), intent(in) :: t
      real(dp), intent(in) :: y(neq)
      integer, intent(in) :: ng
      real(dp), intent(out) :: gout(ng)
    end subroutine
  end interface
  
  interface
    subroutine dlsoda(f, neq, y, t, tout, itol, rtol, atol, itask, &
                      istate, iopt, rwork, lrw, iwork, liw, jac, jt, &
                      common_data)
      import :: dp, odepack_common_data
      implicit none
      procedure(odepack_f) :: f
      integer, intent(in) :: neq
      real(dp), intent(inout) :: y(neq)
      real(dp), intent(inout) :: t
      real(dp), intent(in) :: tout
      integer, intent(in) :: itol
      real(dp), intent(in) :: rtol
      real(dp), intent(in) :: atol(*)
      integer, intent(in) :: itask
      integer, intent(inout) :: istate
      integer, intent(in) :: iopt
      real(dp), intent(inout) :: rwork(lrw)
      integer, intent(in) :: lrw
      integer, intent(inout) :: iwork(liw)
      integer, intent(in) :: liw
      procedure(odepack_jac) :: jac
      integer, intent(in) :: jt
      type(odepack_common_data), target, intent(inout) :: common_data
    end subroutine

    subroutine dlsodar(f, neq, y, t, tout, itol, rtol, atol, itask, &
                       istate, iopt, rwork, lrw, iwork, liw, jac, jt, g, ng, jroot, &
                       common_data)
      import :: dp, odepack_common_data
      implicit none
      procedure(odepack_f) :: f
      integer, intent(in) :: neq
      real(dp), intent(inout) :: y(neq)
      real(dp), intent(inout) :: t
      real(dp), intent(in) :: tout
      integer, intent(in) :: itol
      real(dp), intent(in) :: rtol
      real(dp), intent(in) :: atol(*)
      integer, intent(in) :: itask
      integer, intent(inout) :: istate
      integer, intent(in) :: iopt
      real(dp), intent(inout) :: rwork(lrw)
      integer, intent(in) :: lrw
      integer, intent(inout) :: iwork(liw)
      integer, intent(in) :: liw
      procedure(odepack_jac) :: jac
      integer, intent(in) :: jt
      procedure(odepack_g) :: g
      integer, intent(in) :: ng
      integer, intent(out) :: jroot(ng)
      type(odepack_common_data), target, intent(inout) :: common_data
    end subroutine

  end interface

  interface
    SUBROUTINE DINTDY (T, K, YH, NYH, DKY, IFLAG, common_data)
      import :: odepack_common_data
      type(odepack_common_data), target, intent(inout) :: common_data
      INTEGER K, NYH, IFLAG
      DOUBLE PRECISION T, YH, DKY
      DIMENSION YH(NYH,*), DKY(*)
    end subroutine

    SUBROUTINE DROOTS (NG, HMIN, JFLAG, X0, X1, G0, G1, GX, X, JROOT, common_data)
      import :: odepack_common_data
      type(odepack_common_data), target, intent(inout) :: common_data
      INTEGER NG, JFLAG, JROOT
      DOUBLE PRECISION HMIN, X0, X1, G0, G1, GX, X
      DIMENSION G0(NG), G1(NG), GX(NG), JROOT(NG)
    end subroutine
  end interface

end module