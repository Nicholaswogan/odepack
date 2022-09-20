module odepack_common
  implicit none
  private

  public :: DLS001_type, DLSA01_type, DLSR01_type
  public :: odepack_common_data

  type :: DLS001_type
    double precision :: reals(218) = 0.0d0
    ! double precision :: ROWNS(209)
    ! double precision :: CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
    integer :: ints(37) = 0
    ! integer :: INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS(6)
    ! integer :: ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L
    ! integer :: LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER
    ! integer :: MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
  end type

  ! DOUBLE PRECISION ROWNS,
  ! 1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
  ! INTEGER INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS,
  ! 1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
  ! 2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
  ! 3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
  ! COMMON /DLS001/ ROWNS(209),
  ! 1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
  ! 2   INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS(6),
  ! 3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
  ! 4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
  ! 5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
  
  type :: DLSA01_type
    double precision :: reals(22) = 0.0d0
    ! double precision :: TSW, ROWNS2(20), PDNORM
    integer :: ints(9) = 0
    ! integer :: INSUFR, INSUFI, IXPR, IOWNS2(2), JTYP, MUSED, MXORDN, MXORDS
  end type

  ! DOUBLE PRECISION TSW, ROWNS2, PDNORM
  ! INTEGER INSUFR, INSUFI, IXPR, IOWNS2, JTYP, MUSED, MXORDN, MXORDS
  ! COMMON /DLSA01/ TSW, ROWNS2(20), PDNORM,
  ! 1   INSUFR, INSUFI, IXPR, IOWNS2(2), JTYP, MUSED, MXORDN, MXORDS

  type :: DLSR01_type
    double precision :: reals(5) = 0.0d0
    ! double precision :: ROWNR3(2), T0, TLAST, TOUTC
    integer :: ints(9) = 0
    ! integer :: LG0, LG1, LGX, IOWNR3(2), IRFND, ITASKC, NGC, NGE
  end type

  ! DOUBLE PRECISION ROWNR3, T0, TLAST, TOUTC
  ! INTEGER LG0, LG1, LGX, IOWNR3, IRFND, ITASKC, NGC, NGE
  ! COMMON /DLSR01/ ROWNR3(2), T0, TLAST, TOUTC,
  ! 1   LG0, LG1, LGX, IOWNR3(2), IRFND, ITASKC, NGC, NGE

  type :: odepack_common_data
    integer :: ierr
    type(DLS001_type) :: DLS001
    type(DLSA01_type) :: DLSA01
    type(DLSR01_type) :: DLSR01
  end type

end module