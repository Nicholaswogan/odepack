module odepack_interface
  use odepack_common, only: odepack_common_data
  implicit none
  private

  public :: DINTDY, DROOTS

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