*DECK DUMACH
      DOUBLE PRECISION FUNCTION DUMACH ()
C***BEGIN PROLOGUE  DUMACH
C***PURPOSE  Compute the unit roundoff of the machine.
C***CATEGORY  R1
C***TYPE      DOUBLE PRECISION (RUMACH-S, DUMACH-D)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C *Usage:
C        DOUBLE PRECISION  A, DUMACH
C        A = DUMACH()
C
C *Function Return Values:
C     A : the unit roundoff of the machine.
C
C *Description:
C     The unit roundoff is defined as the smallest positive machine
C     number u such that  1.0 + u .ne. 1.0.  This is computed by DUMACH
C     in a machine-independent manner.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DUMSUM
C***REVISION HISTORY  (YYYYMMDD)
C   19930216  DATE WRITTEN
C   19930818  Added SLATEC-format prologue.  (FNF)
C   20030707  Added DUMSUM to force normal storage of COMP.  (ACH)
C***END PROLOGUE  DUMACH
C
      DOUBLE PRECISION U, COMP
C***FIRST EXECUTABLE STATEMENT  DUMACH
      U = 1.0D0
 10   U = U*0.5D0
      CALL DUMSUM(1.0D0, U, COMP)
      IF (COMP .NE. 1.0D0) GO TO 10
      DUMACH = U*2.0D0
      RETURN
C----------------------- End of Function DUMACH ------------------------
      END
      SUBROUTINE DUMSUM(A,B,C)
C     Routine to force normal storing of A + B, for DUMACH.
      DOUBLE PRECISION A, B, C
      C = A + B
      RETURN
      END
*DECK DCFODE
      SUBROUTINE DCFODE (METH, ELCO, TESCO)
C***BEGIN PROLOGUE  DCFODE
C***SUBSIDIARY
C***PURPOSE  Set ODE integrator coefficients.
C***TYPE      DOUBLE PRECISION (SCFODE-S, DCFODE-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  DCFODE is called by the integrator routine to set coefficients
C  needed there.  The coefficients for the current method, as
C  given by the value of METH, are set for all orders and saved.
C  The maximum order assumed here is 12 if METH = 1 and 5 if METH = 2.
C  (A smaller value of the maximum order is also allowed.)
C  DCFODE is called once at the beginning of the problem,
C  and is not called again unless and until METH is changed.
C
C  The ELCO array contains the basic method coefficients.
C  The coefficients el(i), 1 .le. i .le. nq+1, for the method of
C  order nq are stored in ELCO(i,nq).  They are given by a genetrating
C  polynomial, i.e.,
C      l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
C  For the implicit Adams methods, l(x) is given by
C      dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
C  For the BDF methods, l(x) is given by
C      l(x) = (x+1)*(x+2)* ... *(x+nq)/K,
C  where         K = factorial(nq)*(1 + 1/2 + ... + 1/nq).
C
C  The TESCO array contains test constants used for the
C  local error test and the selection of step size and/or order.
C  At order nq, TESCO(k,nq) is used for the selection of step
C  size at order nq - 1 if k = 1, at order nq if k = 2, and at order
C  nq + 1 if k = 3.
C
C***SEE ALSO  DLSODE
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/double precision versions. (ACH)
C***END PROLOGUE  DCFODE
C**End
      INTEGER METH
      INTEGER I, IB, NQ, NQM1, NQP1
      DOUBLE PRECISION ELCO, TESCO
      DOUBLE PRECISION AGAMQ, FNQ, FNQM1, PC, PINT, RAGQ,
     1   RQFAC, RQ1FAC, TSIGN, XPIN
      DIMENSION ELCO(13,12), TESCO(3,12)
      DIMENSION PC(12)
C
C***FIRST EXECUTABLE STATEMENT  DCFODE
      GO TO (100, 200), METH
C
 100  ELCO(1,1) = 1.0D0
      ELCO(2,1) = 1.0D0
      TESCO(1,1) = 0.0D0
      TESCO(2,1) = 2.0D0
      TESCO(1,2) = 1.0D0
      TESCO(3,12) = 0.0D0
      PC(1) = 1.0D0
      RQFAC = 1.0D0
      DO 140 NQ = 2,12
C-----------------------------------------------------------------------
C The PC array will contain the coefficients of the polynomial
C     p(x) = (x+1)*(x+2)*...*(x+nq-1).
C Initially, p(x) = 1.
C-----------------------------------------------------------------------
        RQ1FAC = RQFAC
        RQFAC = RQFAC/NQ
        NQM1 = NQ - 1
        FNQM1 = NQM1
        NQP1 = NQ + 1
C Form coefficients of p(x)*(x+nq-1). ----------------------------------
        PC(NQ) = 0.0D0
        DO 110 IB = 1,NQM1
          I = NQP1 - IB
 110      PC(I) = PC(I-1) + FNQM1*PC(I)
        PC(1) = FNQM1*PC(1)
C Compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
        PINT = PC(1)
        XPIN = PC(1)/2.0D0
        TSIGN = 1.0D0
        DO 120 I = 2,NQ
          TSIGN = -TSIGN
          PINT = PINT + TSIGN*PC(I)/I
 120      XPIN = XPIN + TSIGN*PC(I)/(I+1)
C Store coefficients in ELCO and TESCO. --------------------------------
        ELCO(1,NQ) = PINT*RQ1FAC
        ELCO(2,NQ) = 1.0D0
        DO 130 I = 2,NQ
 130      ELCO(I+1,NQ) = RQ1FAC*PC(I)/I
        AGAMQ = RQFAC*XPIN
        RAGQ = 1.0D0/AGAMQ
        TESCO(2,NQ) = RAGQ
        IF (NQ .LT. 12) TESCO(1,NQP1) = RAGQ*RQFAC/NQP1
        TESCO(3,NQM1) = RAGQ
 140    CONTINUE
      RETURN
C
 200  PC(1) = 1.0D0
      RQ1FAC = 1.0D0
      DO 230 NQ = 1,5
C-----------------------------------------------------------------------
C The PC array will contain the coefficients of the polynomial
C     p(x) = (x+1)*(x+2)*...*(x+nq).
C Initially, p(x) = 1.
C-----------------------------------------------------------------------
        FNQ = NQ
        NQP1 = NQ + 1
C Form coefficients of p(x)*(x+nq). ------------------------------------
        PC(NQP1) = 0.0D0
        DO 210 IB = 1,NQ
          I = NQ + 2 - IB
 210      PC(I) = PC(I-1) + FNQ*PC(I)
        PC(1) = FNQ*PC(1)
C Store coefficients in ELCO and TESCO. --------------------------------
        DO 220 I = 1,NQP1
 220      ELCO(I,NQ) = PC(I)/PC(2)
        ELCO(2,NQ) = 1.0D0
        TESCO(1,NQ) = RQ1FAC
        TESCO(2,NQ) = NQP1/ELCO(1,NQ)
        TESCO(3,NQ) = (NQ+2)/ELCO(1,NQ)
        RQ1FAC = RQ1FAC/FNQ
 230    CONTINUE
      RETURN
C----------------------- END OF SUBROUTINE DCFODE ----------------------
      END
*DECK DINTDY
      SUBROUTINE DINTDY (T, K, YH, NYH, DKY, IFLAG, common_data)
      use iso_c_binding, only: c_ptr, c_f_pointer, c_loc
      use odepack_common
      type(odepack_common_data), target, intent(inout) :: common_data
C***BEGIN PROLOGUE  DINTDY
C***SUBSIDIARY
C***PURPOSE  Interpolate solution derivatives.
C***TYPE      DOUBLE PRECISION (SINTDY-S, DINTDY-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  DINTDY computes interpolated values of the K-th derivative of the
C  dependent variable vector y, and stores it in DKY.  This routine
C  is called within the package with K = 0 and T = TOUT, but may
C  also be called by the user for any K up to the current order.
C  (See detailed instructions in the usage documentation.)
C
C  The computed values in DKY are gotten by interpolation using the
C  Nordsieck history array YH.  This array corresponds uniquely to a
C  vector-valued polynomial of degree NQCUR or less, and DKY is set
C  to the K-th derivative of this polynomial at T.
C  The formula for DKY is:
C               q
C   DKY(i)  =  sum  c(j,K) * (T - tn)**(j-K) * h**(-j) * YH(i,j+1)
C              j=K
C  where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, tn = TCUR, h = HCUR.
C  The quantities  nq = NQCUR, l = nq+1, N = NEQ, tn, and h are
C  communicated by COMMON.  The above sum is done in reverse order.
C  IFLAG is returned negative if either K or T is out of bounds.
C
C***SEE ALSO  DLSODE
C***ROUTINES CALLED  XERRWD
C***COMMON BLOCKS    DLS001
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/double precision versions. (ACH)
C   010418  Reduced size of Common block /DLS001/. (ACH)
C   031105  Restored 'own' variables to Common block /DLS001/, to
C           enable interrupt/restart feature. (ACH)
C   050427  Corrected roundoff decrement in TP. (ACH)
C***END PROLOGUE  DINTDY
C**End
      INTEGER K, NYH, IFLAG
      DOUBLE PRECISION T, YH, DKY
      DIMENSION YH(NYH,*), DKY(*)
      INTEGER, pointer :: IOWND(:), IOWNS(:),
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      DOUBLE PRECISION, pointer :: ROWNS(:),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
C      COMMON /DLS001/ ROWNS(209),
C     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
C     2   IOWND(6), IOWNS(6),
C     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
C     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
C     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER I, IC, J, JB, JB2, JJ, JJ1, JP1
      DOUBLE PRECISION C, R, S, TP
      CHARACTER*80 MSG
C
C     Common block pointers
      type(DLS001_type), pointer :: DLS001
      type(c_ptr) :: tmp_ptr
C-----------------------------------------------------------------------
C This code associates variables with common data
C-----------------------------------------------------------------------
      DLS001 => common_data%DLS001

      tmp_ptr = c_loc(DLS001%reals(1))
      call c_f_pointer(tmp_ptr,ROWNS,[209])
      tmp_ptr = c_loc(DLS001%reals(210))
      call c_f_pointer(tmp_ptr,CCMAX)
      tmp_ptr = c_loc(DLS001%reals(211))
      call c_f_pointer(tmp_ptr,EL0)
      tmp_ptr = c_loc(DLS001%reals(212))
      call c_f_pointer(tmp_ptr,H)
      tmp_ptr = c_loc(DLS001%reals(213))
      call c_f_pointer(tmp_ptr,HMIN)
      tmp_ptr = c_loc(DLS001%reals(214))
      call c_f_pointer(tmp_ptr,HMXI)
      tmp_ptr = c_loc(DLS001%reals(215))
      call c_f_pointer(tmp_ptr,HU)
      tmp_ptr = c_loc(DLS001%reals(216))
      call c_f_pointer(tmp_ptr,RC)
      tmp_ptr = c_loc(DLS001%reals(217))
      call c_f_pointer(tmp_ptr,TN)
      tmp_ptr = c_loc(DLS001%reals(218))
      call c_f_pointer(tmp_ptr,UROUND)

      tmp_ptr = c_loc(DLS001%ints(1))
      call c_f_pointer(tmp_ptr,IOWND,[6])
      tmp_ptr = c_loc(DLS001%ints(7))
      call c_f_pointer(tmp_ptr,IOWNS,[6])
      tmp_ptr = c_loc(DLS001%ints(13))
      call c_f_pointer(tmp_ptr,ICF)
      tmp_ptr = c_loc(DLS001%ints(14))
      call c_f_pointer(tmp_ptr,IERPJ)
      tmp_ptr = c_loc(DLS001%ints(15))
      call c_f_pointer(tmp_ptr,IERSL)
      tmp_ptr = c_loc(DLS001%ints(16))
      call c_f_pointer(tmp_ptr,JCUR)
      tmp_ptr = c_loc(DLS001%ints(17))
      call c_f_pointer(tmp_ptr,JSTART)
      tmp_ptr = c_loc(DLS001%ints(18))
      call c_f_pointer(tmp_ptr,KFLAG)
      tmp_ptr = c_loc(DLS001%ints(19))
      call c_f_pointer(tmp_ptr,L)
      tmp_ptr = c_loc(DLS001%ints(20))
      call c_f_pointer(tmp_ptr,LYH)
      tmp_ptr = c_loc(DLS001%ints(21))
      call c_f_pointer(tmp_ptr,LEWT)
      tmp_ptr = c_loc(DLS001%ints(22))
      call c_f_pointer(tmp_ptr,LACOR)
      tmp_ptr = c_loc(DLS001%ints(23))
      call c_f_pointer(tmp_ptr,LSAVF)
      tmp_ptr = c_loc(DLS001%ints(24))
      call c_f_pointer(tmp_ptr,LWM)
      tmp_ptr = c_loc(DLS001%ints(25))
      call c_f_pointer(tmp_ptr,LIWM)
      tmp_ptr = c_loc(DLS001%ints(26))
      call c_f_pointer(tmp_ptr,METH)
      tmp_ptr = c_loc(DLS001%ints(27))
      call c_f_pointer(tmp_ptr,MITER)
      tmp_ptr = c_loc(DLS001%ints(28))
      call c_f_pointer(tmp_ptr,MAXORD)
      tmp_ptr = c_loc(DLS001%ints(29))
      call c_f_pointer(tmp_ptr,MAXCOR)
      tmp_ptr = c_loc(DLS001%ints(30))
      call c_f_pointer(tmp_ptr,MSBP)
      tmp_ptr = c_loc(DLS001%ints(31))
      call c_f_pointer(tmp_ptr,MXNCF)
      tmp_ptr = c_loc(DLS001%ints(32))
      call c_f_pointer(tmp_ptr,N)
      tmp_ptr = c_loc(DLS001%ints(33))
      call c_f_pointer(tmp_ptr,NQ)
      tmp_ptr = c_loc(DLS001%ints(34))
      call c_f_pointer(tmp_ptr,NST)
      tmp_ptr = c_loc(DLS001%ints(35))
      call c_f_pointer(tmp_ptr,NFE)
      tmp_ptr = c_loc(DLS001%ints(36))
      call c_f_pointer(tmp_ptr,NJE)
      tmp_ptr = c_loc(DLS001%ints(37))
      call c_f_pointer(tmp_ptr,NQU)
C
C***FIRST EXECUTABLE STATEMENT  DINTDY
      IFLAG = 0
      IF (K .LT. 0 .OR. K .GT. NQ) GO TO 80
      TP = TN - HU -  100.0D0*UROUND*SIGN(ABS(TN) + ABS(HU), HU)
      IF ((T-TP)*(T-TN) .GT. 0.0D0) GO TO 90
C
      S = (T - TN)/H
      IC = 1
      IF (K .EQ. 0) GO TO 15
      JJ1 = L - K
      DO 10 JJ = JJ1,NQ
 10     IC = IC*JJ
 15   C = IC
      DO 20 I = 1,N
 20     DKY(I) = C*YH(I,L)
      IF (K .EQ. NQ) GO TO 55
      JB2 = NQ - K
      DO 50 JB = 1,JB2
        J = NQ - JB
        JP1 = J + 1
        IC = 1
        IF (K .EQ. 0) GO TO 35
        JJ1 = JP1 - K
        DO 30 JJ = JJ1,J
 30       IC = IC*JJ
 35     C = IC
        DO 40 I = 1,N
 40       DKY(I) = C*YH(I,JP1) + S*DKY(I)
 50     CONTINUE
      IF (K .EQ. 0) RETURN
 55   R = H**(-K)
      DO 60 I = 1,N
 60     DKY(I) = R*DKY(I)
      RETURN
C
 80   MSG = 'DINTDY-  K (=I1) illegal      '
      if (common_data%iprint == 1) then
      CALL XERRWD (MSG, 30, 51, 0, 1, K, 0, 0, 0.0D0, 0.0D0)
      endif
      IFLAG = -1
      RETURN
 90   MSG = 'DINTDY-  T (=R1) illegal      '
      if (common_data%iprint == 1) then
      CALL XERRWD (MSG, 30, 52, 0, 0, 0, 0, 1, T, 0.0D0)
      MSG='      T not in interval TCUR - HU (= R1) to TCUR (=R2)      '
      CALL XERRWD (MSG, 60, 52, 0, 0, 0, 0, 2, TP, TN)
      endif
      IFLAG = -2
      RETURN
C----------------------- END OF SUBROUTINE DINTDY ----------------------
      END
*DECK DSOLSY
      SUBROUTINE DSOLSY (WM, IWM, X, TEM, common_data)
      use iso_c_binding, only: c_ptr, c_f_pointer, c_loc
      use odepack_common
      type(odepack_common_data), target, intent(inout) :: common_data
C***BEGIN PROLOGUE  DSOLSY
C***SUBSIDIARY
C***PURPOSE  ODEPACK linear system solver.
C***TYPE      DOUBLE PRECISION (SSOLSY-S, DSOLSY-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  This routine manages the solution of the linear system arising from
C  a chord iteration.  It is called if MITER .ne. 0.
C  If MITER is 1 or 2, it calls DGESL to accomplish this.
C  If MITER = 3 it updates the coefficient h*EL0 in the diagonal
C  matrix, and then computes the solution.
C  If MITER is 4 or 5, it calls DGBSL.
C  Communication with DSOLSY uses the following variables:
C  WM    = real work space containing the inverse diagonal matrix if
C          MITER = 3 and the LU decomposition of the matrix otherwise.
C          Storage of matrix elements starts at WM(3).
C          WM also contains the following matrix-related data:
C          WM(1) = SQRT(UROUND) (not used here),
C          WM(2) = HL0, the previous value of h*EL0, used if MITER = 3.
C  IWM   = integer work space containing pivot information, starting at
C          IWM(21), if MITER is 1, 2, 4, or 5.  IWM also contains band
C          parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
C  X     = the right-hand side vector on input, and the solution vector
C          on output, of length N.
C  TEM   = vector of work space of length N, not used in this version.
C  IERSL = output flag (in COMMON).  IERSL = 0 if no trouble occurred.
C          IERSL = 1 if a singular matrix arose with MITER = 3.
C  This routine also uses the COMMON variables EL0, H, MITER, and N.
C
C***SEE ALSO  DLSODE
C***ROUTINES CALLED  DGBSL, DGESL
C***COMMON BLOCKS    DLS001
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/double precision versions. (ACH)
C   010418  Reduced size of Common block /DLS001/. (ACH)
C   031105  Restored 'own' variables to Common block /DLS001/, to
C           enable interrupt/restart feature. (ACH)
C***END PROLOGUE  DSOLSY
C**End
      INTEGER IWM
      DOUBLE PRECISION WM, X, TEM
      DIMENSION WM(*), IWM(*), X(*), TEM(*)
      INTEGER, pointer :: IOWND(:), IOWNS(:),
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      DOUBLE PRECISION, pointer :: ROWNS(:),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
C      COMMON /DLS001/ ROWNS(209),
C     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
C     2   IOWND(6), IOWNS(6),
C     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
C     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
C     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER I, MEBAND, ML, MU
      DOUBLE PRECISION DI, HL0, PHL0, R
C     Common block pointers
      type(DLS001_type), pointer :: DLS001
      type(c_ptr) :: tmp_ptr
C-----------------------------------------------------------------------
C This code associates variables with common data
C-----------------------------------------------------------------------
      DLS001 => common_data%DLS001

      tmp_ptr = c_loc(DLS001%reals(1))
      call c_f_pointer(tmp_ptr,ROWNS,[209])
      tmp_ptr = c_loc(DLS001%reals(210))
      call c_f_pointer(tmp_ptr,CCMAX)
      tmp_ptr = c_loc(DLS001%reals(211))
      call c_f_pointer(tmp_ptr,EL0)
      tmp_ptr = c_loc(DLS001%reals(212))
      call c_f_pointer(tmp_ptr,H)
      tmp_ptr = c_loc(DLS001%reals(213))
      call c_f_pointer(tmp_ptr,HMIN)
      tmp_ptr = c_loc(DLS001%reals(214))
      call c_f_pointer(tmp_ptr,HMXI)
      tmp_ptr = c_loc(DLS001%reals(215))
      call c_f_pointer(tmp_ptr,HU)
      tmp_ptr = c_loc(DLS001%reals(216))
      call c_f_pointer(tmp_ptr,RC)
      tmp_ptr = c_loc(DLS001%reals(217))
      call c_f_pointer(tmp_ptr,TN)
      tmp_ptr = c_loc(DLS001%reals(218))
      call c_f_pointer(tmp_ptr,UROUND)

      tmp_ptr = c_loc(DLS001%ints(1))
      call c_f_pointer(tmp_ptr,IOWND,[6])
      tmp_ptr = c_loc(DLS001%ints(7))
      call c_f_pointer(tmp_ptr,IOWNS,[6])
      tmp_ptr = c_loc(DLS001%ints(13))
      call c_f_pointer(tmp_ptr,ICF)
      tmp_ptr = c_loc(DLS001%ints(14))
      call c_f_pointer(tmp_ptr,IERPJ)
      tmp_ptr = c_loc(DLS001%ints(15))
      call c_f_pointer(tmp_ptr,IERSL)
      tmp_ptr = c_loc(DLS001%ints(16))
      call c_f_pointer(tmp_ptr,JCUR)
      tmp_ptr = c_loc(DLS001%ints(17))
      call c_f_pointer(tmp_ptr,JSTART)
      tmp_ptr = c_loc(DLS001%ints(18))
      call c_f_pointer(tmp_ptr,KFLAG)
      tmp_ptr = c_loc(DLS001%ints(19))
      call c_f_pointer(tmp_ptr,L)
      tmp_ptr = c_loc(DLS001%ints(20))
      call c_f_pointer(tmp_ptr,LYH)
      tmp_ptr = c_loc(DLS001%ints(21))
      call c_f_pointer(tmp_ptr,LEWT)
      tmp_ptr = c_loc(DLS001%ints(22))
      call c_f_pointer(tmp_ptr,LACOR)
      tmp_ptr = c_loc(DLS001%ints(23))
      call c_f_pointer(tmp_ptr,LSAVF)
      tmp_ptr = c_loc(DLS001%ints(24))
      call c_f_pointer(tmp_ptr,LWM)
      tmp_ptr = c_loc(DLS001%ints(25))
      call c_f_pointer(tmp_ptr,LIWM)
      tmp_ptr = c_loc(DLS001%ints(26))
      call c_f_pointer(tmp_ptr,METH)
      tmp_ptr = c_loc(DLS001%ints(27))
      call c_f_pointer(tmp_ptr,MITER)
      tmp_ptr = c_loc(DLS001%ints(28))
      call c_f_pointer(tmp_ptr,MAXORD)
      tmp_ptr = c_loc(DLS001%ints(29))
      call c_f_pointer(tmp_ptr,MAXCOR)
      tmp_ptr = c_loc(DLS001%ints(30))
      call c_f_pointer(tmp_ptr,MSBP)
      tmp_ptr = c_loc(DLS001%ints(31))
      call c_f_pointer(tmp_ptr,MXNCF)
      tmp_ptr = c_loc(DLS001%ints(32))
      call c_f_pointer(tmp_ptr,N)
      tmp_ptr = c_loc(DLS001%ints(33))
      call c_f_pointer(tmp_ptr,NQ)
      tmp_ptr = c_loc(DLS001%ints(34))
      call c_f_pointer(tmp_ptr,NST)
      tmp_ptr = c_loc(DLS001%ints(35))
      call c_f_pointer(tmp_ptr,NFE)
      tmp_ptr = c_loc(DLS001%ints(36))
      call c_f_pointer(tmp_ptr,NJE)
      tmp_ptr = c_loc(DLS001%ints(37))
      call c_f_pointer(tmp_ptr,NQU)
C
C***FIRST EXECUTABLE STATEMENT  DSOLSY
      IERSL = 0
      GO TO (100, 100, 300, 400, 400), MITER
C 100  CALL DGESL (WM(3), N, N, IWM(21), X, 0)
 100  call dgetrs ('N', n, 1, wm(3), n, iwm(21), x, n, ier)
      RETURN
C
 300  PHL0 = WM(2)
      HL0 = H*EL0
      WM(2) = HL0
      IF (HL0 .EQ. PHL0) GO TO 330
      R = HL0/PHL0
      DO 320 I = 1,N
        DI = 1.0D0 - R*(1.0D0 - 1.0D0/WM(I+2))
        IF (ABS(DI) .EQ. 0.0D0) GO TO 390
 320    WM(I+2) = 1.0D0/DI
 330  DO 340 I = 1,N
 340    X(I) = WM(I+2)*X(I)
      RETURN
 390  IERSL = 1
      RETURN
C
 400  ML = IWM(1)
      MU = IWM(2)
      MEBAND = 2*ML + MU + 1
C      CALL DGBSL (WM(3), MEBAND, N, ML, MU, IWM(21), X, 0)
      call dgbtrs ('N', n, ml, mu, 1, wm(3), meband, iwm(21), x, n, ier)
      RETURN
C----------------------- END OF SUBROUTINE DSOLSY ----------------------
      END
*DECK DEWSET
      SUBROUTINE DEWSET (N, ITOL, RTOL, ATOL, YCUR, EWT)
C***BEGIN PROLOGUE  DEWSET
C***SUBSIDIARY
C***PURPOSE  Set error weight vector.
C***TYPE      DOUBLE PRECISION (SEWSET-S, DEWSET-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  This subroutine sets the error weight vector EWT according to
C      EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N,
C  with the subscript on RTOL and/or ATOL possibly replaced by 1 above,
C  depending on the value of ITOL.
C
C***SEE ALSO  DLSODE
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/double precision versions. (ACH)
C***END PROLOGUE  DEWSET
C**End
      INTEGER N, ITOL
      INTEGER I
      DOUBLE PRECISION RTOL, ATOL, YCUR, EWT
      DIMENSION RTOL(*), ATOL(*), YCUR(N), EWT(N)
C
C***FIRST EXECUTABLE STATEMENT  DEWSET
      GO TO (10, 20, 30, 40), ITOL
 10   CONTINUE
      DO 15 I = 1,N
 15     EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(1)
      RETURN
 20   CONTINUE
      DO 25 I = 1,N
 25     EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(I)
      RETURN
 30   CONTINUE
      DO 35 I = 1,N
 35     EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(1)
      RETURN
 40   CONTINUE
      DO 45 I = 1,N
 45     EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(I)
      RETURN
C----------------------- END OF SUBROUTINE DEWSET ----------------------
      END
*DECK DSTODA
      SUBROUTINE DSTODA (NEQ, Y, YH, NYH, YH1, EWT, SAVF, ACOR,
     1   WM, IWM, F, JAC, PJAC, SLVS, common_data)
      use iso_c_binding, only: c_ptr, c_f_pointer, c_loc
      use odepack_common
      type(odepack_common_data), target, intent(inout) :: common_data
      EXTERNAL F, JAC, PJAC, SLVS
      INTEGER NEQ, NYH, IWM
      DOUBLE PRECISION Y, YH, YH1, EWT, SAVF, ACOR, WM
      DIMENSION NEQ(*), Y(*), YH(NYH,*), YH1(*), EWT(*), SAVF(*),
     1   ACOR(*), WM(*), IWM(*)
      INTEGER, pointer ::
     1   IOWND(:), IALTH, IPUP, LMAX, MEO, NQNYH, NSLP,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER, pointer ::
     1   IOWND2(:), ICOUNT, IRFLAG, JTYP, MUSED, MXORDN, MXORDS
      DOUBLE PRECISION, pointer ::
     1   CONIT, CRATE, EL(:), ELCO(:,:), HOLD, RMAX, TESCO(:,:),
     2   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION, pointer ::
     1   ROWND2, CM1(:), CM2(:), PDEST, PDLAST, RATIO,
     1   PDNORM
C      COMMON /DLS001/ CONIT, CRATE, EL(13), ELCO(13,12),
C     1   HOLD, RMAX, TESCO(3,12),
C     2   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
C     3   IOWND(6), IALTH, IPUP, LMAX, MEO, NQNYH, NSLP,
C     4   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
C     5   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
C     6   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
C      COMMON /DLSA01/ ROWND2, CM1(12), CM2(5), PDEST, PDLAST, RATIO,
C     1   PDNORM,
C     2   IOWND2(3), ICOUNT, IRFLAG, JTYP, MUSED, MXORDN, MXORDS
      INTEGER I, I1, IREDO, IRET, J, JB, M, NCF, NEWQ
      INTEGER LM1, LM1P1, LM2, LM2P1, NQM1, NQM2
      DOUBLE PRECISION DCON, DDN, DEL, DELP, DSM, DUP, EXDN, EXSM, EXUP,
     1   R, RH, RHDN, RHSM, RHUP, TOLD, DMNORM
      DOUBLE PRECISION ALPHA, DM1,DM2, EXM1,EXM2,
     1   PDH, PNORM, RATE, RH1, RH1IT, RH2, RM, SM1(12)
      SAVE SM1
      DATA SM1/0.5D0, 0.575D0, 0.55D0, 0.45D0, 0.35D0, 0.25D0,
     1   0.20D0, 0.15D0, 0.10D0, 0.075D0, 0.050D0, 0.025D0/
C-----------------------------------------------------------------------
C DSTODA performs one step of the integration of an initial value
C problem for a system of ordinary differential equations.
C Note: DSTODA is independent of the value of the iteration method
C indicator MITER, when this is .ne. 0, and hence is independent
C of the type of chord method used, or the Jacobian structure.
C Communication with DSTODA is done with the following variables:
C
C Y      = an array of length .ge. N used as the Y argument in
C          all calls to F and JAC.
C NEQ    = integer array containing problem size in NEQ(1), and
C          passed as the NEQ argument in all calls to F and JAC.
C YH     = an NYH by LMAX array containing the dependent variables
C          and their approximate scaled derivatives, where
C          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
C          j-th derivative of y(i), scaled by H**j/factorial(j)
C          (j = 0,1,...,NQ).  On entry for the first step, the first
C          two columns of YH must be set from the initial values.
C NYH    = a constant integer .ge. N, the first dimension of YH.
C YH1    = a one-dimensional array occupying the same space as YH.
C EWT    = an array of length N containing multiplicative weights
C          for local error measurements.  Local errors in y(i) are
C          compared to 1.0/EWT(i) in various error tests.
C SAVF   = an array of working storage, of length N.
C ACOR   = a work array of length N, used for the accumulated
C          corrections.  On a successful return, ACOR(i) contains
C          the estimated one-step local error in y(i).
C WM,IWM = real and integer work arrays associated with matrix
C          operations in chord iteration (MITER .ne. 0).
C PJAC   = name of routine to evaluate and preprocess Jacobian matrix
C          and P = I - H*EL0*Jac, if a chord method is being used.
C          It also returns an estimate of norm(Jac) in PDNORM.
C SLVS   = name of routine to solve linear system in chord iteration.
C CCMAX  = maximum relative change in H*EL0 before PJAC is called.
C H      = the step size to be attempted on the next step.
C          H is altered by the error control algorithm during the
C          problem.  H can be either positive or negative, but its
C          sign must remain constant throughout the problem.
C HMIN   = the minimum absolute value of the step size H to be used.
C HMXI   = inverse of the maximum absolute value of H to be used.
C          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
C          HMIN and HMXI may be changed at any time, but will not
C          take effect until the next change of H is considered.
C TN     = the independent variable. TN is updated on each step taken.
C JSTART = an integer used for input only, with the following
C          values and meanings:
C               0  perform the first step.
C           .gt.0  take a new step continuing from the last.
C              -1  take the next step with a new value of H,
C                    N, METH, MITER, and/or matrix parameters.
C              -2  take the next step with a new value of H,
C                    but with other inputs unchanged.
C          On return, JSTART is set to 1 to facilitate continuation.
C KFLAG  = a completion code with the following meanings:
C               0  the step was succesful.
C              -1  the requested error could not be achieved.
C              -2  corrector convergence could not be achieved.
C              -3  fatal error in PJAC or SLVS.
C          A return with KFLAG = -1 or -2 means either
C          ABS(H) = HMIN or 10 consecutive failures occurred.
C          On a return with KFLAG negative, the values of TN and
C          the YH array are as of the beginning of the last
C          step, and H is the last step size attempted.
C MAXORD = the maximum order of integration method to be allowed.
C MAXCOR = the maximum number of corrector iterations allowed.
C MSBP   = maximum number of steps between PJAC calls (MITER .gt. 0).
C MXNCF  = maximum number of convergence failures allowed.
C METH   = current method.
C          METH = 1 means Adams method (nonstiff)
C          METH = 2 means BDF method (stiff)
C          METH may be reset by DSTODA.
C MITER  = corrector iteration method.
C          MITER = 0 means functional iteration.
C          MITER = JT .gt. 0 means a chord iteration corresponding
C          to Jacobian type JT.  (The DLSODA/DLSODAR argument JT is
C          communicated here as JTYP, but is not used in DSTODA
C          except to load MITER following a method switch.)
C          MITER may be reset by DSTODA.
C N      = the number of first-order differential equations.
C-----------------------------------------------------------------------
C     Common block pointers
      type(DLS001_type), pointer :: DLS001
      type(DLSA01_type), pointer :: DLSA01
      type(c_ptr) :: tmp_ptr
C-----------------------------------------------------------------------
C This code associates variables with common data
C-----------------------------------------------------------------------
      DLS001 => common_data%DLS001

      tmp_ptr = c_loc(DLS001%reals(1))
      call c_f_pointer(tmp_ptr,CONIT)
      tmp_ptr = c_loc(DLS001%reals(2))
      call c_f_pointer(tmp_ptr,CRATE)
      tmp_ptr = c_loc(DLS001%reals(3))
      call c_f_pointer(tmp_ptr,EL,[13])
      tmp_ptr = c_loc(DLS001%reals(16))
      call c_f_pointer(tmp_ptr,ELCO,[13,12])
      tmp_ptr = c_loc(DLS001%reals(172))
      call c_f_pointer(tmp_ptr,HOLD)
      tmp_ptr = c_loc(DLS001%reals(173))
      call c_f_pointer(tmp_ptr,RMAX)
      tmp_ptr = c_loc(DLS001%reals(174))
      call c_f_pointer(tmp_ptr,TESCO,[3,12])
      tmp_ptr = c_loc(DLS001%reals(210))
      call c_f_pointer(tmp_ptr,CCMAX)
      tmp_ptr = c_loc(DLS001%reals(211))
      call c_f_pointer(tmp_ptr,EL0)
      tmp_ptr = c_loc(DLS001%reals(212))
      call c_f_pointer(tmp_ptr,H)
      tmp_ptr = c_loc(DLS001%reals(213))
      call c_f_pointer(tmp_ptr,HMIN)
      tmp_ptr = c_loc(DLS001%reals(214))
      call c_f_pointer(tmp_ptr,HMXI)
      tmp_ptr = c_loc(DLS001%reals(215))
      call c_f_pointer(tmp_ptr,HU)
      tmp_ptr = c_loc(DLS001%reals(216))
      call c_f_pointer(tmp_ptr,RC)
      tmp_ptr = c_loc(DLS001%reals(217))
      call c_f_pointer(tmp_ptr,TN)
      tmp_ptr = c_loc(DLS001%reals(218))
      call c_f_pointer(tmp_ptr,UROUND)

      tmp_ptr = c_loc(DLS001%ints(1))
      call c_f_pointer(tmp_ptr,IOWND,[6])
      tmp_ptr = c_loc(DLS001%ints(7))
      call c_f_pointer(tmp_ptr,IALTH)
      tmp_ptr = c_loc(DLS001%ints(8))
      call c_f_pointer(tmp_ptr,IPUP)
      tmp_ptr = c_loc(DLS001%ints(9))
      call c_f_pointer(tmp_ptr,LMAX)
      tmp_ptr = c_loc(DLS001%ints(10))
      call c_f_pointer(tmp_ptr,MEO)
      tmp_ptr = c_loc(DLS001%ints(11))
      call c_f_pointer(tmp_ptr,NQNYH)
      tmp_ptr = c_loc(DLS001%ints(12))
      call c_f_pointer(tmp_ptr,NSLP)
      tmp_ptr = c_loc(DLS001%ints(13))
      call c_f_pointer(tmp_ptr,ICF)
      tmp_ptr = c_loc(DLS001%ints(14))
      call c_f_pointer(tmp_ptr,IERPJ)
      tmp_ptr = c_loc(DLS001%ints(15))
      call c_f_pointer(tmp_ptr,IERSL)
      tmp_ptr = c_loc(DLS001%ints(16))
      call c_f_pointer(tmp_ptr,JCUR)
      tmp_ptr = c_loc(DLS001%ints(17))
      call c_f_pointer(tmp_ptr,JSTART)
      tmp_ptr = c_loc(DLS001%ints(18))
      call c_f_pointer(tmp_ptr,KFLAG)
      tmp_ptr = c_loc(DLS001%ints(19))
      call c_f_pointer(tmp_ptr,L)
      tmp_ptr = c_loc(DLS001%ints(20))
      call c_f_pointer(tmp_ptr,LYH)
      tmp_ptr = c_loc(DLS001%ints(21))
      call c_f_pointer(tmp_ptr,LEWT)
      tmp_ptr = c_loc(DLS001%ints(22))
      call c_f_pointer(tmp_ptr,LACOR)
      tmp_ptr = c_loc(DLS001%ints(23))
      call c_f_pointer(tmp_ptr,LSAVF)
      tmp_ptr = c_loc(DLS001%ints(24))
      call c_f_pointer(tmp_ptr,LWM)
      tmp_ptr = c_loc(DLS001%ints(25))
      call c_f_pointer(tmp_ptr,LIWM)
      tmp_ptr = c_loc(DLS001%ints(26))
      call c_f_pointer(tmp_ptr,METH)
      tmp_ptr = c_loc(DLS001%ints(27))
      call c_f_pointer(tmp_ptr,MITER)
      tmp_ptr = c_loc(DLS001%ints(28))
      call c_f_pointer(tmp_ptr,MAXORD)
      tmp_ptr = c_loc(DLS001%ints(29))
      call c_f_pointer(tmp_ptr,MAXCOR)
      tmp_ptr = c_loc(DLS001%ints(30))
      call c_f_pointer(tmp_ptr,MSBP)
      tmp_ptr = c_loc(DLS001%ints(31))
      call c_f_pointer(tmp_ptr,MXNCF)
      tmp_ptr = c_loc(DLS001%ints(32))
      call c_f_pointer(tmp_ptr,N)
      tmp_ptr = c_loc(DLS001%ints(33))
      call c_f_pointer(tmp_ptr,NQ)
      tmp_ptr = c_loc(DLS001%ints(34))
      call c_f_pointer(tmp_ptr,NST)
      tmp_ptr = c_loc(DLS001%ints(35))
      call c_f_pointer(tmp_ptr,NFE)
      tmp_ptr = c_loc(DLS001%ints(36))
      call c_f_pointer(tmp_ptr,NJE)
      tmp_ptr = c_loc(DLS001%ints(37))
      call c_f_pointer(tmp_ptr,NQU)
      
      DLSA01 => common_data%DLSA01

      tmp_ptr = c_loc(DLSA01%reals(1))
      call c_f_pointer(tmp_ptr,ROWND2)
      tmp_ptr = c_loc(DLSA01%reals(2))
      call c_f_pointer(tmp_ptr,CM1,[12])

      tmp_ptr = c_loc(DLSA01%reals(14))
      call c_f_pointer(tmp_ptr,CM2,[5])
      tmp_ptr = c_loc(DLSA01%reals(19))
      call c_f_pointer(tmp_ptr,PDEST)
      tmp_ptr = c_loc(DLSA01%reals(20))
      call c_f_pointer(tmp_ptr,PDLAST)
      tmp_ptr = c_loc(DLSA01%reals(21))
      call c_f_pointer(tmp_ptr,RATIO)
      tmp_ptr = c_loc(DLSA01%reals(22))
      call c_f_pointer(tmp_ptr,PDNORM)

      tmp_ptr = c_loc(DLSA01%ints(1))
      call c_f_pointer(tmp_ptr,IOWND2,[3])
      tmp_ptr = c_loc(DLSA01%ints(4))
      call c_f_pointer(tmp_ptr,ICOUNT)
      tmp_ptr = c_loc(DLSA01%ints(5))
      call c_f_pointer(tmp_ptr,IRFLAG)
      tmp_ptr = c_loc(DLSA01%ints(6))
      call c_f_pointer(tmp_ptr,JTYP)
      tmp_ptr = c_loc(DLSA01%ints(7))
      call c_f_pointer(tmp_ptr,MUSED)
      tmp_ptr = c_loc(DLSA01%ints(8))
      call c_f_pointer(tmp_ptr,MXORDN)
      tmp_ptr = c_loc(DLSA01%ints(9))
      call c_f_pointer(tmp_ptr,MXORDS)

C
      KFLAG = 0
      TOLD = TN
      NCF = 0
      IERPJ = 0
      IERSL = 0
      JCUR = 0
      ICF = 0
      DELP = 0.0D0
      IF (JSTART .GT. 0) GO TO 200
      IF (JSTART .EQ. -1) GO TO 100
      IF (JSTART .EQ. -2) GO TO 160
C-----------------------------------------------------------------------
C On the first call, the order is set to 1, and other variables are
C initialized.  RMAX is the maximum ratio by which H can be increased
C in a single step.  It is initially 1.E4 to compensate for the small
C initial H, but then is normally equal to 10.  If a failure
C occurs (in corrector convergence or error test), RMAX is set at 2
C for the next increase.
C DCFODE is called to get the needed coefficients for both methods.
C-----------------------------------------------------------------------
      LMAX = MAXORD + 1
      NQ = 1
      L = 2
      IALTH = 2
      RMAX = 10000.0D0
      RC = 0.0D0
      EL0 = 1.0D0
      CRATE = 0.7D0
      HOLD = H
      NSLP = 0
      IPUP = MITER
      IRET = 3
C Initialize switching parameters.  METH = 1 is assumed initially. -----
      ICOUNT = 20
      IRFLAG = 0
      PDEST = 0.0D0
      PDLAST = 0.0D0
      RATIO = 5.0D0
      CALL DCFODE (2, ELCO, TESCO)
      DO 10 I = 1,5
 10     CM2(I) = TESCO(2,I)*ELCO(I+1,I)
      CALL DCFODE (1, ELCO, TESCO)
      DO 20 I = 1,12
 20     CM1(I) = TESCO(2,I)*ELCO(I+1,I)
      GO TO 150
C-----------------------------------------------------------------------
C The following block handles preliminaries needed when JSTART = -1.
C IPUP is set to MITER to force a matrix update.
C If an order increase is about to be considered (IALTH = 1),
C IALTH is reset to 2 to postpone consideration one more step.
C If the caller has changed METH, DCFODE is called to reset
C the coefficients of the method.
C If H is to be changed, YH must be rescaled.
C If H or METH is being changed, IALTH is reset to L = NQ + 1
C to prevent further changes in H for that many steps.
C-----------------------------------------------------------------------
 100  IPUP = MITER
      LMAX = MAXORD + 1
      IF (IALTH .EQ. 1) IALTH = 2
      IF (METH .EQ. MUSED) GO TO 160
      CALL DCFODE (METH, ELCO, TESCO)
      IALTH = L
      IRET = 1
C-----------------------------------------------------------------------
C The el vector and related constants are reset
C whenever the order NQ is changed, or at the start of the problem.
C-----------------------------------------------------------------------
 150  DO 155 I = 1,L
 155    EL(I) = ELCO(I,NQ)
      NQNYH = NQ*NYH
      RC = RC*EL(1)/EL0
      EL0 = EL(1)
      CONIT = 0.5D0/(NQ+2)
      GO TO (160, 170, 200), IRET
C-----------------------------------------------------------------------
C If H is being changed, the H ratio RH is checked against
C RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
C L = NQ + 1 to prevent a change of H for that many steps, unless
C forced by a convergence or error test failure.
C-----------------------------------------------------------------------
 160  IF (H .EQ. HOLD) GO TO 200
      RH = H/HOLD
      H = HOLD
      IREDO = 3
      GO TO 175
 170  RH = MAX(RH,HMIN/ABS(H))
 175  RH = MIN(RH,RMAX)
      RH = RH/MAX(1.0D0,ABS(H)*HMXI*RH)
C-----------------------------------------------------------------------
C If METH = 1, also restrict the new step size by the stability region.
C If this reduces H, set IRFLAG to 1 so that if there are roundoff
C problems later, we can assume that is the cause of the trouble.
C-----------------------------------------------------------------------
      IF (METH .EQ. 2) GO TO 178
      IRFLAG = 0
      PDH = MAX(ABS(H)*PDLAST,0.000001D0)
      IF (RH*PDH*1.00001D0 .LT. SM1(NQ)) GO TO 178
      RH = SM1(NQ)/PDH
      IRFLAG = 1
 178  CONTINUE
      R = 1.0D0
      DO 180 J = 2,L
        R = R*RH
        DO 180 I = 1,N
 180      YH(I,J) = YH(I,J)*R
      H = H*RH
      RC = RC*RH
      IALTH = L
      IF (IREDO .EQ. 0) GO TO 690
C-----------------------------------------------------------------------
C This section computes the predicted values by effectively
C multiplying the YH array by the Pascal triangle matrix.
C RC is the ratio of new to old values of the coefficient  H*EL(1).
C When RC differs from 1 by more than CCMAX, IPUP is set to MITER
C to force PJAC to be called, if a Jacobian is involved.
C In any case, PJAC is called at least every MSBP steps.
C-----------------------------------------------------------------------
 200  IF (ABS(RC-1.0D0) .GT. CCMAX) IPUP = MITER
      IF (NST .GE. NSLP+MSBP) IPUP = MITER
      TN = TN + H
      I1 = NQNYH + 1
      DO 215 JB = 1,NQ
        I1 = I1 - NYH
CDIR$ IVDEP
        DO 210 I = I1,NQNYH
 210      YH1(I) = YH1(I) + YH1(I+NYH)
 215    CONTINUE
      PNORM = DMNORM (N, YH1, EWT)
C-----------------------------------------------------------------------
C Up to MAXCOR corrector iterations are taken.  A convergence test is
C made on the RMS-norm of each correction, weighted by the error
C weight vector EWT.  The sum of the corrections is accumulated in the
C vector ACOR(i).  The YH array is not altered in the corrector loop.
C-----------------------------------------------------------------------
 220  M = 0
      RATE = 0.0D0
      DEL = 0.0D0
      DO 230 I = 1,N
 230    Y(I) = YH(I,1)
      CALL F (NEQ, TN, Y, SAVF, common_data%ierr)
      if (common_data%ierr < 0) return
      NFE = NFE + 1
      IF (IPUP .LE. 0) GO TO 250
C-----------------------------------------------------------------------
C If indicated, the matrix P = I - H*EL(1)*J is reevaluated and
C preprocessed before starting the corrector iteration.  IPUP is set
C to 0 as an indicator that this has been done.
C-----------------------------------------------------------------------
      CALL PJAC (NEQ, Y, YH, NYH, EWT, ACOR, SAVF, WM, IWM, F, JAC,
     1           common_data)
      if (common_data%ierr < 0) return
      IPUP = 0
      RC = 1.0D0
      NSLP = NST
      CRATE = 0.7D0
      IF (IERPJ .NE. 0) GO TO 430
 250  DO 260 I = 1,N
 260    ACOR(I) = 0.0D0
 270  IF (MITER .NE. 0) GO TO 350
C-----------------------------------------------------------------------
C In the case of functional iteration, update Y directly from
C the result of the last function evaluation.
C-----------------------------------------------------------------------
      DO 290 I = 1,N
        SAVF(I) = H*SAVF(I) - YH(I,2)
 290    Y(I) = SAVF(I) - ACOR(I)
      DEL = DMNORM (N, Y, EWT)
      DO 300 I = 1,N
        Y(I) = YH(I,1) + EL(1)*SAVF(I)
 300    ACOR(I) = SAVF(I)
      GO TO 400
C-----------------------------------------------------------------------
C In the case of the chord method, compute the corrector error,
C and solve the linear system with that as right-hand side and
C P as coefficient matrix.
C-----------------------------------------------------------------------
 350  DO 360 I = 1,N
 360    Y(I) = H*SAVF(I) - (YH(I,2) + ACOR(I))
      CALL SLVS (WM, IWM, Y, SAVF, common_data)
      IF (IERSL .LT. 0) GO TO 430
      IF (IERSL .GT. 0) GO TO 410
      DEL = DMNORM (N, Y, EWT)
      DO 380 I = 1,N
        ACOR(I) = ACOR(I) + Y(I)
 380    Y(I) = YH(I,1) + EL(1)*ACOR(I)
C-----------------------------------------------------------------------
C Test for convergence.  If M .gt. 0, an estimate of the convergence
C rate constant is stored in CRATE, and this is used in the test.
C
C We first check for a change of iterates that is the size of
C roundoff error.  If this occurs, the iteration has converged, and a
C new rate estimate is not formed.
C In all other cases, force at least two iterations to estimate a
C local Lipschitz constant estimate for Adams methods.
C On convergence, form PDEST = local maximum Lipschitz constant
C estimate.  PDLAST is the most recent nonzero estimate.
C-----------------------------------------------------------------------
 400  CONTINUE
      IF (DEL .LE. 100.0D0*PNORM*UROUND) GO TO 450
      IF (M .EQ. 0 .AND. METH .EQ. 1) GO TO 405
      IF (M .EQ. 0) GO TO 402
      RM = 1024.0D0
      IF (DEL .LE. 1024.0D0*DELP) RM = DEL/DELP
      RATE = MAX(RATE,RM)
      CRATE = MAX(0.2D0*CRATE,RM)
 402  DCON = DEL*MIN(1.0D0,1.5D0*CRATE)/(TESCO(2,NQ)*CONIT)
      IF (DCON .GT. 1.0D0) GO TO 405
      PDEST = MAX(PDEST,RATE/ABS(H*EL(1)))
      IF (PDEST .NE. 0.0D0) PDLAST = PDEST
      GO TO 450
 405  CONTINUE
      M = M + 1
      IF (M .EQ. MAXCOR) GO TO 410
      IF (M .GE. 2 .AND. DEL .GT. 2.0D0*DELP) GO TO 410
      DELP = DEL
      CALL F (NEQ, TN, Y, SAVF, common_data%ierr)
      if (common_data%ierr < 0) return
      NFE = NFE + 1
      GO TO 270
C-----------------------------------------------------------------------
C The corrector iteration failed to converge.
C If MITER .ne. 0 and the Jacobian is out of date, PJAC is called for
C the next try.  Otherwise the YH array is retracted to its values
C before prediction, and H is reduced, if possible.  If H cannot be
C reduced or MXNCF failures have occurred, exit with KFLAG = -2.
C-----------------------------------------------------------------------
 410  IF (MITER .EQ. 0 .OR. JCUR .EQ. 1) GO TO 430
      ICF = 1
      IPUP = MITER
      GO TO 220
 430  ICF = 2
      NCF = NCF + 1
      RMAX = 2.0D0
      TN = TOLD
      I1 = NQNYH + 1
      DO 445 JB = 1,NQ
        I1 = I1 - NYH
CDIR$ IVDEP
        DO 440 I = I1,NQNYH
 440      YH1(I) = YH1(I) - YH1(I+NYH)
 445    CONTINUE
      IF (IERPJ .LT. 0 .OR. IERSL .LT. 0) GO TO 680
      IF (ABS(H) .LE. HMIN*1.00001D0) GO TO 670
      IF (NCF .EQ. MXNCF) GO TO 670
      RH = 0.25D0
      IPUP = MITER
      IREDO = 1
      GO TO 170
C-----------------------------------------------------------------------
C The corrector has converged.  JCUR is set to 0
C to signal that the Jacobian involved may need updating later.
C The local error test is made and control passes to statement 500
C if it fails.
C-----------------------------------------------------------------------
 450  JCUR = 0
      IF (M .EQ. 0) DSM = DEL/TESCO(2,NQ)
      IF (M .GT. 0) DSM = DMNORM (N, ACOR, EWT)/TESCO(2,NQ)
      IF (DSM .GT. 1.0D0) GO TO 500
C-----------------------------------------------------------------------
C After a successful step, update the YH array.
C Decrease ICOUNT by 1, and if it is -1, consider switching methods.
C If a method switch is made, reset various parameters,
C rescale the YH array, and exit.  If there is no switch,
C consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
C If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
C use in a possible order increase on the next step.
C If a change in H is considered, an increase or decrease in order
C by one is considered also.  A change in H is made only if it is by a
C factor of at least 1.1.  If not, IALTH is set to 3 to prevent
C testing for that many steps.
C-----------------------------------------------------------------------
      KFLAG = 0
      IREDO = 0
      NST = NST + 1
      HU = H
      NQU = NQ
      MUSED = METH
      DO 460 J = 1,L
        DO 460 I = 1,N
 460      YH(I,J) = YH(I,J) + EL(J)*ACOR(I)
      ICOUNT = ICOUNT - 1
      IF (ICOUNT .GE. 0) GO TO 488
      IF (METH .EQ. 2) GO TO 480
C-----------------------------------------------------------------------
C We are currently using an Adams method.  Consider switching to BDF.
C If the current order is greater than 5, assume the problem is
C not stiff, and skip this section.
C If the Lipschitz constant and error estimate are not polluted
C by roundoff, go to 470 and perform the usual test.
C Otherwise, switch to the BDF methods if the last step was
C restricted to insure stability (irflag = 1), and stay with Adams
C method if not.  When switching to BDF with polluted error estimates,
C in the absence of other information, double the step size.
C
C When the estimates are OK, we make the usual test by computing
C the step size we could have (ideally) used on this step,
C with the current (Adams) method, and also that for the BDF.
C If NQ .gt. MXORDS, we consider changing to order MXORDS on switching.
C Compare the two step sizes to decide whether to switch.
C The step size advantage must be at least RATIO = 5 to switch.
C-----------------------------------------------------------------------
      IF (NQ .GT. 5) GO TO 488
      IF (DSM .GT. 100.0D0*PNORM*UROUND .AND. PDEST .NE. 0.0D0)
     1   GO TO 470
      IF (IRFLAG .EQ. 0) GO TO 488
      RH2 = 2.0D0
      NQM2 = MIN(NQ,MXORDS)
      GO TO 478
 470  CONTINUE
      EXSM = 1.0D0/L
      RH1 = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
      RH1IT = 2.0D0*RH1
      PDH = PDLAST*ABS(H)
      IF (PDH*RH1 .GT. 0.00001D0) RH1IT = SM1(NQ)/PDH
      RH1 = MIN(RH1,RH1IT)
      IF (NQ .LE. MXORDS) GO TO 474
         NQM2 = MXORDS
         LM2 = MXORDS + 1
         EXM2 = 1.0D0/LM2
         LM2P1 = LM2 + 1
         DM2 = DMNORM (N, YH(1,LM2P1), EWT)/CM2(MXORDS)
         RH2 = 1.0D0/(1.2D0*DM2**EXM2 + 0.0000012D0)
         GO TO 476
 474  DM2 = DSM*(CM1(NQ)/CM2(NQ))
      RH2 = 1.0D0/(1.2D0*DM2**EXSM + 0.0000012D0)
      NQM2 = NQ
 476  CONTINUE
      IF (RH2 .LT. RATIO*RH1) GO TO 488
C THE SWITCH TEST PASSED.  RESET RELEVANT QUANTITIES FOR BDF. ----------
 478  RH = RH2
      ICOUNT = 20
      METH = 2
      MITER = JTYP
      PDLAST = 0.0D0
      NQ = NQM2
      L = NQ + 1
      GO TO 170
C-----------------------------------------------------------------------
C We are currently using a BDF method.  Consider switching to Adams.
C Compute the step size we could have (ideally) used on this step,
C with the current (BDF) method, and also that for the Adams.
C If NQ .gt. MXORDN, we consider changing to order MXORDN on switching.
C Compare the two step sizes to decide whether to switch.
C The step size advantage must be at least 5/RATIO = 1 to switch.
C If the step size for Adams would be so small as to cause
C roundoff pollution, we stay with BDF.
C-----------------------------------------------------------------------
 480  CONTINUE
      EXSM = 1.0D0/L
      IF (MXORDN .GE. NQ) GO TO 484
         NQM1 = MXORDN
         LM1 = MXORDN + 1
         EXM1 = 1.0D0/LM1
         LM1P1 = LM1 + 1
         DM1 = DMNORM (N, YH(1,LM1P1), EWT)/CM1(MXORDN)
         RH1 = 1.0D0/(1.2D0*DM1**EXM1 + 0.0000012D0)
         GO TO 486
 484  DM1 = DSM*(CM2(NQ)/CM1(NQ))
      RH1 = 1.0D0/(1.2D0*DM1**EXSM + 0.0000012D0)
      NQM1 = NQ
      EXM1 = EXSM
 486  RH1IT = 2.0D0*RH1
      PDH = PDNORM*ABS(H)
      IF (PDH*RH1 .GT. 0.00001D0) RH1IT = SM1(NQM1)/PDH
      RH1 = MIN(RH1,RH1IT)
      RH2 = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
      IF (RH1*RATIO .LT. 5.0D0*RH2) GO TO 488
      ALPHA = MAX(0.001D0,RH1)
      DM1 = (ALPHA**EXM1)*DM1
      IF (DM1 .LE. 1000.0D0*UROUND*PNORM) GO TO 488
C The switch test passed.  Reset relevant quantities for Adams. --------
      RH = RH1
      ICOUNT = 20
      METH = 1
      MITER = 0
      PDLAST = 0.0D0
      NQ = NQM1
      L = NQ + 1
      GO TO 170
C
C No method switch is being made.  Do the usual step/order selection. --
 488  CONTINUE
      IALTH = IALTH - 1
      IF (IALTH .EQ. 0) GO TO 520
      IF (IALTH .GT. 1) GO TO 700
      IF (L .EQ. LMAX) GO TO 700
      DO 490 I = 1,N
 490    YH(I,LMAX) = ACOR(I)
      GO TO 700
C-----------------------------------------------------------------------
C The error test failed.  KFLAG keeps track of multiple failures.
C Restore TN and the YH array to their previous values, and prepare
C to try the step again.  Compute the optimum step size for this or
C one lower order.  After 2 or more failures, H is forced to decrease
C by a factor of 0.2 or less.
C-----------------------------------------------------------------------
 500  KFLAG = KFLAG - 1
      TN = TOLD
      I1 = NQNYH + 1
      DO 515 JB = 1,NQ
        I1 = I1 - NYH
CDIR$ IVDEP
        DO 510 I = I1,NQNYH
 510      YH1(I) = YH1(I) - YH1(I+NYH)
 515    CONTINUE
      RMAX = 2.0D0
      IF (ABS(H) .LE. HMIN*1.00001D0) GO TO 660
      IF (KFLAG .LE. -3) GO TO 640
      IREDO = 2
      RHUP = 0.0D0
      GO TO 540
C-----------------------------------------------------------------------
C Regardless of the success or failure of the step, factors
C RHDN, RHSM, and RHUP are computed, by which H could be multiplied
C at order NQ - 1, order NQ, or order NQ + 1, respectively.
C In the case of failure, RHUP = 0.0 to avoid an order increase.
C The largest of these is determined and the new order chosen
C accordingly.  If the order is to be increased, we compute one
C additional scaled derivative.
C-----------------------------------------------------------------------
 520  RHUP = 0.0D0
      IF (L .EQ. LMAX) GO TO 540
      DO 530 I = 1,N
 530    SAVF(I) = ACOR(I) - YH(I,LMAX)
      DUP = DMNORM (N, SAVF, EWT)/TESCO(3,NQ)
      EXUP = 1.0D0/(L+1)
      RHUP = 1.0D0/(1.4D0*DUP**EXUP + 0.0000014D0)
 540  EXSM = 1.0D0/L
      RHSM = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
      RHDN = 0.0D0
      IF (NQ .EQ. 1) GO TO 550
      DDN = DMNORM (N, YH(1,L), EWT)/TESCO(1,NQ)
      EXDN = 1.0D0/NQ
      RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
C If METH = 1, limit RH according to the stability region also. --------
 550  IF (METH .EQ. 2) GO TO 560
      PDH = MAX(ABS(H)*PDLAST,0.000001D0)
      IF (L .LT. LMAX) RHUP = MIN(RHUP,SM1(L)/PDH)
      RHSM = MIN(RHSM,SM1(NQ)/PDH)
      IF (NQ .GT. 1) RHDN = MIN(RHDN,SM1(NQ-1)/PDH)
      PDEST = 0.0D0
 560  IF (RHSM .GE. RHUP) GO TO 570
      IF (RHUP .GT. RHDN) GO TO 590
      GO TO 580
 570  IF (RHSM .LT. RHDN) GO TO 580
      NEWQ = NQ
      RH = RHSM
      GO TO 620
 580  NEWQ = NQ - 1
      RH = RHDN
      IF (KFLAG .LT. 0 .AND. RH .GT. 1.0D0) RH = 1.0D0
      GO TO 620
 590  NEWQ = L
      RH = RHUP
      IF (RH .LT. 1.1D0) GO TO 610
      R = EL(L)/L
      DO 600 I = 1,N
 600    YH(I,NEWQ+1) = ACOR(I)*R
      GO TO 630
 610  IALTH = 3
      GO TO 700
C If METH = 1 and H is restricted by stability, bypass 10 percent test.
 620  IF (METH .EQ. 2) GO TO 622
      IF (RH*PDH*1.00001D0 .GE. SM1(NEWQ)) GO TO 625
 622  IF (KFLAG .EQ. 0 .AND. RH .LT. 1.1D0) GO TO 610
 625  IF (KFLAG .LE. -2) RH = MIN(RH,0.2D0)
C-----------------------------------------------------------------------
C If there is a change of order, reset NQ, L, and the coefficients.
C In any case H is reset according to RH and the YH array is rescaled.
C Then exit from 690 if the step was OK, or redo the step otherwise.
C-----------------------------------------------------------------------
      IF (NEWQ .EQ. NQ) GO TO 170
 630  NQ = NEWQ
      L = NQ + 1
      IRET = 2
      GO TO 150
C-----------------------------------------------------------------------
C Control reaches this section if 3 or more failures have occured.
C If 10 failures have occurred, exit with KFLAG = -1.
C It is assumed that the derivatives that have accumulated in the
C YH array have errors of the wrong order.  Hence the first
C derivative is recomputed, and the order is set to 1.  Then
C H is reduced by a factor of 10, and the step is retried,
C until it succeeds or H reaches HMIN.
C-----------------------------------------------------------------------
 640  IF (KFLAG .EQ. -10) GO TO 660
      RH = 0.1D0
      RH = MAX(HMIN/ABS(H),RH)
      H = H*RH
      DO 645 I = 1,N
 645    Y(I) = YH(I,1)
      CALL F (NEQ, TN, Y, SAVF, common_data%ierr)
      if (common_data%ierr < 0) return
      NFE = NFE + 1
      DO 650 I = 1,N
 650    YH(I,2) = H*SAVF(I)
      IPUP = MITER
      IALTH = 5
      IF (NQ .EQ. 1) GO TO 200
      NQ = 1
      L = 2
      IRET = 3
      GO TO 150
C-----------------------------------------------------------------------
C All returns are made through this section.  H is saved in HOLD
C to allow the caller to change H on the next step.
C-----------------------------------------------------------------------
 660  KFLAG = -1
      GO TO 720
 670  KFLAG = -2
      GO TO 720
 680  KFLAG = -3
      GO TO 720
 690  RMAX = 10.0D0
 700  R = 1.0D0/TESCO(2,NQU)
      DO 710 I = 1,N
 710    ACOR(I) = ACOR(I)*R
 720  HOLD = H
      JSTART = 1
      RETURN
C----------------------- End of Subroutine DSTODA ----------------------
      END
*DECK DPRJA
      SUBROUTINE DPRJA (NEQ, Y, YH, NYH, EWT, FTEM, SAVF, WM, IWM,
     1   F, JAC, common_data)
      use iso_c_binding, only: c_ptr, c_f_pointer, c_loc
      use odepack_common
      type(odepack_common_data), target, intent(inout) :: common_data
      EXTERNAL F, JAC
      INTEGER NEQ, NYH, IWM
      DOUBLE PRECISION Y, YH, EWT, FTEM, SAVF, WM
      DIMENSION NEQ(*), Y(*), YH(NYH,*), EWT(*), FTEM(*), SAVF(*),
     1   WM(*), IWM(*)
      INTEGER, pointer :: IOWND(:), IOWNS(:),
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER, pointer :: 
     1   IOWND2(:), IOWNS2(:), JTYP, MUSED, MXORDN, MXORDS
      DOUBLE PRECISION, pointer :: ROWNS(:),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION, pointer :: ROWND2, ROWNS2(:), PDNORM
C      COMMON /DLS001/ ROWNS(209),
C     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
C     2   IOWND(6), IOWNS(6),
C     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
C     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
C     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
C      COMMON /DLSA01/ ROWND2, ROWNS2(20), PDNORM,
C     1   IOWND2(3), IOWNS2(2), JTYP, MUSED, MXORDN, MXORDS
      INTEGER I, I1, I2, IER, II, J, J1, JJ, LENP,
     1   MBA, MBAND, MEB1, MEBAND, ML, ML3, MU, NP1
      DOUBLE PRECISION CON, FAC, HL0, R, R0, SRUR, YI, YJ, YJJ,
     1   DMNORM, DFNORM, DBNORM
C-----------------------------------------------------------------------
C DPRJA is called by DSTODA to compute and process the matrix
C P = I - H*EL(1)*J , where J is an approximation to the Jacobian.
C Here J is computed by the user-supplied routine JAC if
C MITER = 1 or 4 or by finite differencing if MITER = 2 or 5.
C J, scaled by -H*EL(1), is stored in WM.  Then the norm of J (the
C matrix norm consistent with the weighted max-norm on vectors given
C by DMNORM) is computed, and J is overwritten by P.  P is then
C subjected to LU decomposition in preparation for later solution
C of linear systems with P as coefficient matrix.  This is done
C by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.
C
C In addition to variables described previously, communication
C with DPRJA uses the following:
C Y     = array containing predicted values on entry.
C FTEM  = work array of length N (ACOR in DSTODA).
C SAVF  = array containing f evaluated at predicted y.
C WM    = real work space for matrices.  On output it contains the
C         LU decomposition of P.
C         Storage of matrix elements starts at WM(3).
C         WM also contains the following matrix-related data:
C         WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
C IWM   = integer work space containing pivot information, starting at
C         IWM(21).   IWM also contains the band parameters
C         ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
C EL0   = EL(1) (input).
C PDNORM= norm of Jacobian matrix. (Output).
C IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
C         P matrix found to be singular.
C JCUR  = output flag = 1 to indicate that the Jacobian matrix
C         (or approximation) is now current.
C This routine also uses the Common variables EL0, H, TN, UROUND,
C MITER, N, NFE, and NJE.
C-----------------------------------------------------------------------
C     Common block pointers
      type(DLS001_type), pointer :: DLS001
      type(DLSA01_type), pointer :: DLSA01
      type(c_ptr) :: tmp_ptr
C-----------------------------------------------------------------------
C This code associates variables with common data
C-----------------------------------------------------------------------
      DLS001 => common_data%DLS001

      tmp_ptr = c_loc(DLS001%reals(1))
      call c_f_pointer(tmp_ptr,ROWNS,[209])
      tmp_ptr = c_loc(DLS001%reals(210))
      call c_f_pointer(tmp_ptr,CCMAX)
      tmp_ptr = c_loc(DLS001%reals(211))
      call c_f_pointer(tmp_ptr,EL0)
      tmp_ptr = c_loc(DLS001%reals(212))
      call c_f_pointer(tmp_ptr,H)
      tmp_ptr = c_loc(DLS001%reals(213))
      call c_f_pointer(tmp_ptr,HMIN)
      tmp_ptr = c_loc(DLS001%reals(214))
      call c_f_pointer(tmp_ptr,HMXI)
      tmp_ptr = c_loc(DLS001%reals(215))
      call c_f_pointer(tmp_ptr,HU)
      tmp_ptr = c_loc(DLS001%reals(216))
      call c_f_pointer(tmp_ptr,RC)
      tmp_ptr = c_loc(DLS001%reals(217))
      call c_f_pointer(tmp_ptr,TN)
      tmp_ptr = c_loc(DLS001%reals(218))
      call c_f_pointer(tmp_ptr,UROUND)

      tmp_ptr = c_loc(DLS001%ints(1))
      call c_f_pointer(tmp_ptr,IOWND,[6])
      tmp_ptr = c_loc(DLS001%ints(7))
      call c_f_pointer(tmp_ptr,IOWNS,[6])
      tmp_ptr = c_loc(DLS001%ints(13))
      call c_f_pointer(tmp_ptr,ICF)
      tmp_ptr = c_loc(DLS001%ints(14))
      call c_f_pointer(tmp_ptr,IERPJ)
      tmp_ptr = c_loc(DLS001%ints(15))
      call c_f_pointer(tmp_ptr,IERSL)
      tmp_ptr = c_loc(DLS001%ints(16))
      call c_f_pointer(tmp_ptr,JCUR)
      tmp_ptr = c_loc(DLS001%ints(17))
      call c_f_pointer(tmp_ptr,JSTART)
      tmp_ptr = c_loc(DLS001%ints(18))
      call c_f_pointer(tmp_ptr,KFLAG)
      tmp_ptr = c_loc(DLS001%ints(19))
      call c_f_pointer(tmp_ptr,L)
      tmp_ptr = c_loc(DLS001%ints(20))
      call c_f_pointer(tmp_ptr,LYH)
      tmp_ptr = c_loc(DLS001%ints(21))
      call c_f_pointer(tmp_ptr,LEWT)
      tmp_ptr = c_loc(DLS001%ints(22))
      call c_f_pointer(tmp_ptr,LACOR)
      tmp_ptr = c_loc(DLS001%ints(23))
      call c_f_pointer(tmp_ptr,LSAVF)
      tmp_ptr = c_loc(DLS001%ints(24))
      call c_f_pointer(tmp_ptr,LWM)
      tmp_ptr = c_loc(DLS001%ints(25))
      call c_f_pointer(tmp_ptr,LIWM)
      tmp_ptr = c_loc(DLS001%ints(26))
      call c_f_pointer(tmp_ptr,METH)
      tmp_ptr = c_loc(DLS001%ints(27))
      call c_f_pointer(tmp_ptr,MITER)
      tmp_ptr = c_loc(DLS001%ints(28))
      call c_f_pointer(tmp_ptr,MAXORD)
      tmp_ptr = c_loc(DLS001%ints(29))
      call c_f_pointer(tmp_ptr,MAXCOR)
      tmp_ptr = c_loc(DLS001%ints(30))
      call c_f_pointer(tmp_ptr,MSBP)
      tmp_ptr = c_loc(DLS001%ints(31))
      call c_f_pointer(tmp_ptr,MXNCF)
      tmp_ptr = c_loc(DLS001%ints(32))
      call c_f_pointer(tmp_ptr,N)
      tmp_ptr = c_loc(DLS001%ints(33))
      call c_f_pointer(tmp_ptr,NQ)
      tmp_ptr = c_loc(DLS001%ints(34))
      call c_f_pointer(tmp_ptr,NST)
      tmp_ptr = c_loc(DLS001%ints(35))
      call c_f_pointer(tmp_ptr,NFE)
      tmp_ptr = c_loc(DLS001%ints(36))
      call c_f_pointer(tmp_ptr,NJE)
      tmp_ptr = c_loc(DLS001%ints(37))
      call c_f_pointer(tmp_ptr,NQU)

      DLSA01 => common_data%DLSA01

      tmp_ptr = c_loc(DLSA01%reals(1))
      call c_f_pointer(tmp_ptr,ROWND2)
      tmp_ptr = c_loc(DLSA01%reals(2))
      call c_f_pointer(tmp_ptr,ROWNS2,[20])
      tmp_ptr = c_loc(DLSA01%reals(22))
      call c_f_pointer(tmp_ptr,PDNORM)

      tmp_ptr = c_loc(DLSA01%ints(1))
      call c_f_pointer(tmp_ptr,IOWND2,[3])
      tmp_ptr = c_loc(DLSA01%ints(4))
      call c_f_pointer(tmp_ptr,IOWNS2,[2])
      tmp_ptr = c_loc(DLSA01%ints(6))
      call c_f_pointer(tmp_ptr,JTYP)
      tmp_ptr = c_loc(DLSA01%ints(7))
      call c_f_pointer(tmp_ptr,MUSED)
      tmp_ptr = c_loc(DLSA01%ints(8))
      call c_f_pointer(tmp_ptr,MXORDN)
      tmp_ptr = c_loc(DLSA01%ints(9))
      call c_f_pointer(tmp_ptr,MXORDS)
C
      NJE = NJE + 1
      IERPJ = 0
      JCUR = 1
      HL0 = H*EL0
      GO TO (100, 200, 300, 400, 500), MITER
C If MITER = 1, call JAC and multiply by scalar. -----------------------
 100  LENP = N*N
      DO 110 I = 1,LENP
 110    WM(I+2) = 0.0D0
      CALL JAC (NEQ, TN, Y, 0, 0, WM(3), N, common_data%ierr)
      if (common_data%ierr < 0) return
      CON = -HL0
      DO 120 I = 1,LENP
 120    WM(I+2) = WM(I+2)*CON
      GO TO 240
C If MITER = 2, make N calls to F to approximate J. --------------------
 200  FAC = DMNORM (N, SAVF, EWT)
      R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
      IF (R0 .EQ. 0.0D0) R0 = 1.0D0
      SRUR = WM(1)
      J1 = 2
      DO 230 J = 1,N
        YJ = Y(J)
        R = MAX(SRUR*ABS(YJ),R0/EWT(J))
        Y(J) = Y(J) + R
        FAC = -HL0/R
        CALL F (NEQ, TN, Y, FTEM, common_data%ierr)
        if (common_data%ierr < 0) return
        DO 220 I = 1,N
 220      WM(I+J1) = (FTEM(I) - SAVF(I))*FAC
        Y(J) = YJ
        J1 = J1 + N
 230    CONTINUE
      NFE = NFE + N
 240  CONTINUE
C Compute norm of Jacobian. --------------------------------------------
      PDNORM = DFNORM (N, WM(3), EWT)/ABS(HL0)
C Add identity matrix. -------------------------------------------------
      J = 3
      NP1 = N + 1
      DO 250 I = 1,N
        WM(J) = WM(J) + 1.0D0
 250    J = J + NP1
C Do LU decomposition on P. --------------------------------------------
C      CALL DGEFA (WM(3), N, N, IWM(21), IER)
      call dgetrf (n, n, wm(3), n, iwm(21), ier)
      IF (IER .NE. 0) IERPJ = 1
      RETURN
C Dummy block only, since MITER is never 3 in this routine. ------------
 300  RETURN
C If MITER = 4, call JAC and multiply by scalar. -----------------------
 400  ML = IWM(1)
      MU = IWM(2)
      ML3 = ML + 3
      MBAND = ML + MU + 1
      MEBAND = MBAND + ML
      LENP = MEBAND*N
      DO 410 I = 1,LENP
 410    WM(I+2) = 0.0D0
      CALL JAC (NEQ, TN, Y, ML, MU, WM(ML3), MEBAND, common_data%ierr)
      if (common_data%ierr < 0) return
      CON = -HL0
      DO 420 I = 1,LENP
 420    WM(I+2) = WM(I+2)*CON
      GO TO 570
C If MITER = 5, make MBAND calls to F to approximate J. ----------------
 500  ML = IWM(1)
      MU = IWM(2)
      MBAND = ML + MU + 1
      MBA = MIN(MBAND,N)
      MEBAND = MBAND + ML
      MEB1 = MEBAND - 1
      SRUR = WM(1)
      FAC = DMNORM (N, SAVF, EWT)
      R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
      IF (R0 .EQ. 0.0D0) R0 = 1.0D0
      DO 560 J = 1,MBA
        DO 530 I = J,N,MBAND
          YI = Y(I)
          R = MAX(SRUR*ABS(YI),R0/EWT(I))
 530      Y(I) = Y(I) + R
        CALL F (NEQ, TN, Y, FTEM, common_data%ierr)
        if (common_data%ierr < 0) return
        DO 550 JJ = J,N,MBAND
          Y(JJ) = YH(JJ,1)
          YJJ = Y(JJ)
          R = MAX(SRUR*ABS(YJJ),R0/EWT(JJ))
          FAC = -HL0/R
          I1 = MAX(JJ-MU,1)
          I2 = MIN(JJ+ML,N)
          II = JJ*MEB1 - ML + 2
          DO 540 I = I1,I2
 540        WM(II+I) = (FTEM(I) - SAVF(I))*FAC
 550      CONTINUE
 560    CONTINUE
      NFE = NFE + MBA
 570  CONTINUE
C Compute norm of Jacobian. --------------------------------------------
      PDNORM = DBNORM (N, WM(ML+3), MEBAND, ML, MU, EWT)/ABS(HL0)
C Add identity matrix. -------------------------------------------------
      II = MBAND + 2
      DO 580 I = 1,N
        WM(II) = WM(II) + 1.0D0
 580    II = II + MEBAND
C Do LU decomposition of P. --------------------------------------------
C      CALL DGBFA (WM(3), MEBAND, N, ML, MU, IWM(21), IER)
      call dgbtrf (n, n, ml, mu, wm(3), meband, iwm(21), ier)
      IF (IER .NE. 0) IERPJ = 1
      RETURN
C----------------------- End of Subroutine DPRJA -----------------------
      END
*DECK DMNORM
      DOUBLE PRECISION FUNCTION DMNORM (N, V, W)
C-----------------------------------------------------------------------
C This function routine computes the weighted max-norm
C of the vector of length N contained in the array V, with weights
C contained in the array w of length N:
C   DMNORM = MAX(i=1,...,N) ABS(V(i))*W(i)
C-----------------------------------------------------------------------
      INTEGER N,   I
      DOUBLE PRECISION V, W,   VM
      DIMENSION V(N), W(N)
      VM = 0.0D0
      DO 10 I = 1,N
 10     VM = MAX(VM,ABS(V(I))*W(I))
      DMNORM = VM
      RETURN
C----------------------- End of Function DMNORM ------------------------
      END
*DECK DFNORM
      DOUBLE PRECISION FUNCTION DFNORM (N, A, W)
C-----------------------------------------------------------------------
C This function computes the norm of a full N by N matrix,
C stored in the array A, that is consistent with the weighted max-norm
C on vectors, with weights stored in the array W:
C   DFNORM = MAX(i=1,...,N) ( W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j) )
C-----------------------------------------------------------------------
      INTEGER N,   I, J
      DOUBLE PRECISION A,   W, AN, SUM
      DIMENSION A(N,N), W(N)
      AN = 0.0D0
      DO 20 I = 1,N
        SUM = 0.0D0
        DO 10 J = 1,N
 10       SUM = SUM + ABS(A(I,J))/W(J)
        AN = MAX(AN,SUM*W(I))
 20     CONTINUE
      DFNORM = AN
      RETURN
C----------------------- End of Function DFNORM ------------------------
      END
*DECK DBNORM
      DOUBLE PRECISION FUNCTION DBNORM (N, A, NRA, ML, MU, W)
C-----------------------------------------------------------------------
C This function computes the norm of a banded N by N matrix,
C stored in the array A, that is consistent with the weighted max-norm
C on vectors, with weights stored in the array W.
C ML and MU are the lower and upper half-bandwidths of the matrix.
C NRA is the first dimension of the A array, NRA .ge. ML+MU+1.
C In terms of the matrix elements a(i,j), the norm is given by:
C   DBNORM = MAX(i=1,...,N) ( W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j) )
C-----------------------------------------------------------------------
      INTEGER N, NRA, ML, MU
      INTEGER I, I1, JLO, JHI, J
      DOUBLE PRECISION A, W
      DOUBLE PRECISION AN, SUM
      DIMENSION A(NRA,N), W(N)
      AN = 0.0D0
      DO 20 I = 1,N
        SUM = 0.0D0
        I1 = I + MU + 1
        JLO = MAX(I-ML,1)
        JHI = MIN(I+MU,N)
        DO 10 J = JLO,JHI
 10       SUM = SUM + ABS(A(I1-J,J))/W(J)
        AN = MAX(AN,SUM*W(I))
 20     CONTINUE
      DBNORM = AN
      RETURN
C----------------------- End of Function DBNORM ------------------------
      END
*DECK DSRCMA
      SUBROUTINE DSRCMA (RSAV, ISAV, JOB, common_data)
      use iso_c_binding, only: c_ptr, c_f_pointer, c_loc
      use odepack_common
      type(odepack_common_data), target, intent(inout) :: common_data
C-----------------------------------------------------------------------
C This routine saves or restores (depending on JOB) the contents of
C the Common blocks DLS001, DLSA01, which are used
C internally by one or more ODEPACK solvers.
C
C RSAV = real array of length 240 or more.
C ISAV = integer array of length 46 or more.
C JOB  = flag indicating to save or restore the Common blocks:
C        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
C        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
C        A call with JOB = 2 presumes a prior call with JOB = 1.
C-----------------------------------------------------------------------
      INTEGER ISAV, JOB
      INTEGER, pointer :: ILS(:), ILSA(:)
      INTEGER I, LENRLS, LENILS, LENRLA, LENILA
      DOUBLE PRECISION RSAV
      DOUBLE PRECISION, pointer :: RLS(:), RLSA(:)
      DIMENSION RSAV(*), ISAV(*)
      SAVE LENRLS, LENILS, LENRLA, LENILA
C      COMMON /DLS001/ RLS(218), ILS(37)
C      COMMON /DLSA01/ RLSA(22), ILSA(9)
      DATA LENRLS/218/, LENILS/37/, LENRLA/22/, LENILA/9/
C     Common block pointers
      type(DLS001_type), pointer :: DLS001
      type(DLSA01_type), pointer :: DLSA01
      type(c_ptr) :: tmp_ptr
C-----------------------------------------------------------------------
C This code associates variables with common data
C-----------------------------------------------------------------------
      DLS001 => common_data%DLS001

      tmp_ptr = c_loc(DLS001%reals(1))
      call c_f_pointer(tmp_ptr,RLS,[218])

      tmp_ptr = c_loc(DLS001%ints(1))
      call c_f_pointer(tmp_ptr,ILS,[37])

      DLSA01 => common_data%DLSA01

      tmp_ptr = c_loc(DLSA01%reals(1))
      call c_f_pointer(tmp_ptr,RLSA,[22])

      tmp_ptr = c_loc(DLSA01%ints(1))
      call c_f_pointer(tmp_ptr,ILSA,[9])
      
C
      IF (JOB .EQ. 2) GO TO 100
      DO 10 I = 1,LENRLS
 10     RSAV(I) = RLS(I)
      DO 15 I = 1,LENRLA
 15     RSAV(LENRLS+I) = RLSA(I)
C
      DO 20 I = 1,LENILS
 20     ISAV(I) = ILS(I)
      DO 25 I = 1,LENILA
 25     ISAV(LENILS+I) = ILSA(I)
C
      RETURN
C
 100  CONTINUE
      DO 110 I = 1,LENRLS
 110     RLS(I) = RSAV(I)
      DO 115 I = 1,LENRLA
 115     RLSA(I) = RSAV(LENRLS+I)
C
      DO 120 I = 1,LENILS
 120     ILS(I) = ISAV(I)
      DO 125 I = 1,LENILA
 125     ILSA(I) = ISAV(LENILS+I)
C
      RETURN
C----------------------- End of Subroutine DSRCMA ----------------------
      END
*DECK DRCHEK
      SUBROUTINE DRCHEK (JOB, G, NEQ, Y, YH,NYH, G0, G1, GX, JROOT, IRT,
     1                   common_data)
      use iso_c_binding, only: c_ptr, c_f_pointer, c_loc
      use odepack_interface, only: DINTDY, DROOTS
      use odepack_common
      type(odepack_common_data), target, intent(inout) :: common_data
      EXTERNAL G
      INTEGER JOB, NEQ, NYH, JROOT, IRT
      DOUBLE PRECISION Y, YH, G0, G1, GX
      DIMENSION NEQ(*), Y(*), YH(NYH,*), G0(*), G1(*), GX(*), JROOT(*)
      INTEGER, pointer :: IOWND(:), IOWNS(:),
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER, pointer :: IOWND3(:), IOWNR3(:), IRFND, ITASKC, NGC,
     1   NGE
      DOUBLE PRECISION, pointer :: ROWNS(:),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION, pointer :: ROWNR3(:), T0, TLAST, TOUTC
C      COMMON /DLS001/ ROWNS(209),
C     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
C     2   IOWND(6), IOWNS(6),
C     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
C     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
C     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
C      COMMON /DLSR01/ ROWNR3(2), T0, TLAST, TOUTC,
C     1   IOWND3(3), IOWNR3(2), IRFND, ITASKC, NGC, NGE
      INTEGER I, IFLAG, JFLAG
      DOUBLE PRECISION HMING, T1, TEMP1, TEMP2, X
      LOGICAL ZROOT
C-----------------------------------------------------------------------
C This routine checks for the presence of a root in the vicinity of
C the current T, in a manner depending on the input flag JOB.  It calls
C Subroutine DROOTS to locate the root as precisely as possible.
C
C In addition to variables described previously, DRCHEK
C uses the following for communication:
C JOB    = integer flag indicating type of call:
C          JOB = 1 means the problem is being initialized, and DRCHEK
C                  is to look for a root at or very near the initial T.
C          JOB = 2 means a continuation call to the solver was just
C                  made, and DRCHEK is to check for a root in the
C                  relevant part of the step last taken.
C          JOB = 3 means a successful step was just taken, and DRCHEK
C                  is to look for a root in the interval of the step.
C G0     = array of length NG, containing the value of g at T = T0.
C          G0 is input for JOB .ge. 2, and output in all cases.
C G1,GX  = arrays of length NG for work space.
C IRT    = completion flag:
C          IRT = 0  means no root was found.
C          IRT = -1 means JOB = 1 and a root was found too near to T.
C          IRT = 1  means a legitimate root was found (JOB = 2 or 3).
C                   On return, T0 is the root location, and Y is the
C                   corresponding solution vector.
C T0     = value of T at one endpoint of interval of interest.  Only
C          roots beyond T0 in the direction of integration are sought.
C          T0 is input if JOB .ge. 2, and output in all cases.
C          T0 is updated by DRCHEK, whether a root is found or not.
C TLAST  = last value of T returned by the solver (input only).
C TOUTC  = copy of TOUT (input only).
C IRFND  = input flag showing whether the last step taken had a root.
C          IRFND = 1 if it did, = 0 if not.
C ITASKC = copy of ITASK (input only).
C NGC    = copy of NG (input only).
C-----------------------------------------------------------------------
C
C     Common block pointers
      type(DLS001_type), pointer :: DLS001
      type(DLSR01_type), pointer :: DLSR01
      type(c_ptr) :: tmp_ptr
C-----------------------------------------------------------------------
C This code associates variables with common data
C-----------------------------------------------------------------------
      DLS001 => common_data%DLS001

      tmp_ptr = c_loc(DLS001%reals(1))
      call c_f_pointer(tmp_ptr,ROWNS,[209])
      tmp_ptr = c_loc(DLS001%reals(210))
      call c_f_pointer(tmp_ptr,CCMAX)
      tmp_ptr = c_loc(DLS001%reals(211))
      call c_f_pointer(tmp_ptr,EL0)
      tmp_ptr = c_loc(DLS001%reals(212))
      call c_f_pointer(tmp_ptr,H)
      tmp_ptr = c_loc(DLS001%reals(213))
      call c_f_pointer(tmp_ptr,HMIN)
      tmp_ptr = c_loc(DLS001%reals(214))
      call c_f_pointer(tmp_ptr,HMXI)
      tmp_ptr = c_loc(DLS001%reals(215))
      call c_f_pointer(tmp_ptr,HU)
      tmp_ptr = c_loc(DLS001%reals(216))
      call c_f_pointer(tmp_ptr,RC)
      tmp_ptr = c_loc(DLS001%reals(217))
      call c_f_pointer(tmp_ptr,TN)
      tmp_ptr = c_loc(DLS001%reals(218))
      call c_f_pointer(tmp_ptr,UROUND)

      tmp_ptr = c_loc(DLS001%ints(1))
      call c_f_pointer(tmp_ptr,IOWND,[6])
      tmp_ptr = c_loc(DLS001%ints(7))
      call c_f_pointer(tmp_ptr,IOWNS,[6])
      tmp_ptr = c_loc(DLS001%ints(13))
      call c_f_pointer(tmp_ptr,ICF)
      tmp_ptr = c_loc(DLS001%ints(14))
      call c_f_pointer(tmp_ptr,IERPJ)
      tmp_ptr = c_loc(DLS001%ints(15))
      call c_f_pointer(tmp_ptr,IERSL)
      tmp_ptr = c_loc(DLS001%ints(16))
      call c_f_pointer(tmp_ptr,JCUR)
      tmp_ptr = c_loc(DLS001%ints(17))
      call c_f_pointer(tmp_ptr,JSTART)
      tmp_ptr = c_loc(DLS001%ints(18))
      call c_f_pointer(tmp_ptr,KFLAG)
      tmp_ptr = c_loc(DLS001%ints(19))
      call c_f_pointer(tmp_ptr,L)
      tmp_ptr = c_loc(DLS001%ints(20))
      call c_f_pointer(tmp_ptr,LYH)
      tmp_ptr = c_loc(DLS001%ints(21))
      call c_f_pointer(tmp_ptr,LEWT)
      tmp_ptr = c_loc(DLS001%ints(22))
      call c_f_pointer(tmp_ptr,LACOR)
      tmp_ptr = c_loc(DLS001%ints(23))
      call c_f_pointer(tmp_ptr,LSAVF)
      tmp_ptr = c_loc(DLS001%ints(24))
      call c_f_pointer(tmp_ptr,LWM)
      tmp_ptr = c_loc(DLS001%ints(25))
      call c_f_pointer(tmp_ptr,LIWM)
      tmp_ptr = c_loc(DLS001%ints(26))
      call c_f_pointer(tmp_ptr,METH)
      tmp_ptr = c_loc(DLS001%ints(27))
      call c_f_pointer(tmp_ptr,MITER)
      tmp_ptr = c_loc(DLS001%ints(28))
      call c_f_pointer(tmp_ptr,MAXORD)
      tmp_ptr = c_loc(DLS001%ints(29))
      call c_f_pointer(tmp_ptr,MAXCOR)
      tmp_ptr = c_loc(DLS001%ints(30))
      call c_f_pointer(tmp_ptr,MSBP)
      tmp_ptr = c_loc(DLS001%ints(31))
      call c_f_pointer(tmp_ptr,MXNCF)
      tmp_ptr = c_loc(DLS001%ints(32))
      call c_f_pointer(tmp_ptr,N)
      tmp_ptr = c_loc(DLS001%ints(33))
      call c_f_pointer(tmp_ptr,NQ)
      tmp_ptr = c_loc(DLS001%ints(34))
      call c_f_pointer(tmp_ptr,NST)
      tmp_ptr = c_loc(DLS001%ints(35))
      call c_f_pointer(tmp_ptr,NFE)
      tmp_ptr = c_loc(DLS001%ints(36))
      call c_f_pointer(tmp_ptr,NJE)
      tmp_ptr = c_loc(DLS001%ints(37))
      call c_f_pointer(tmp_ptr,NQU)

      DLSR01 => common_data%DLSR01

      tmp_ptr = c_loc(DLSR01%reals(1))
      call c_f_pointer(tmp_ptr,ROWNR3,[2])
      tmp_ptr = c_loc(DLSR01%reals(3))
      call c_f_pointer(tmp_ptr,T0)
      tmp_ptr = c_loc(DLSR01%reals(4))
      call c_f_pointer(tmp_ptr,TLAST)
      tmp_ptr = c_loc(DLSR01%reals(5))
      call c_f_pointer(tmp_ptr,TOUTC)

      tmp_ptr = c_loc(DLSR01%ints(1))
      call c_f_pointer(tmp_ptr,IOWND3,[3])
      tmp_ptr = c_loc(DLSR01%ints(4))
      call c_f_pointer(tmp_ptr,IOWNR3,[2])
      tmp_ptr = c_loc(DLSR01%ints(6))
      call c_f_pointer(tmp_ptr,IRFND)
      tmp_ptr = c_loc(DLSR01%ints(7))
      call c_f_pointer(tmp_ptr,ITASKC)
      tmp_ptr = c_loc(DLSR01%ints(8))
      call c_f_pointer(tmp_ptr,NGC)
      tmp_ptr = c_loc(DLSR01%ints(9))
      call c_f_pointer(tmp_ptr,NGE)

C
C
      IRT = 0
      DO 10 I = 1,NGC
 10     JROOT(I) = 0
      HMING = (ABS(TN) + ABS(H))*UROUND*100.0D0
C
      GO TO (100, 200, 300), JOB
C
C Evaluate g at initial T, and check for zero values. ------------------
 100  CONTINUE
      T0 = TN
      CALL G (NEQ, T0, Y, NGC, G0, common_data%ierr)
      if (common_data%ierr < 0) return
      NGE = 1
      ZROOT = .FALSE.
      DO 110 I = 1,NGC
 110    IF (ABS(G0(I)) .LE. 0.0D0) ZROOT = .TRUE.
      IF (.NOT. ZROOT) GO TO 190
C g has a zero at T.  Look at g at T + (small increment). --------------
      TEMP2 = MAX(HMING/ABS(H), 0.1D0)
      TEMP1 = TEMP2*H
      T0 = T0 + TEMP1
      DO 120 I = 1,N
 120    Y(I) = Y(I) + TEMP2*YH(I,2)
      CALL G (NEQ, T0, Y, NGC, G0, common_data%ierr)
      if (common_data%ierr < 0) return
      NGE = NGE + 1
      ZROOT = .FALSE.
      DO 130 I = 1,NGC
 130    IF (ABS(G0(I)) .LE. 0.0D0) ZROOT = .TRUE.
      IF (.NOT. ZROOT) GO TO 190
C g has a zero at T and also close to T.  Take error return. -----------
      IRT = -1
      RETURN
C
 190  CONTINUE
      RETURN
C
C
 200  CONTINUE
      IF (IRFND .EQ. 0) GO TO 260
C If a root was found on the previous step, evaluate G0 = g(T0). -------
      CALL DINTDY (T0, 0, YH, NYH, Y, IFLAG, common_data)
      CALL G (NEQ, T0, Y, NGC, G0, common_data%ierr)
      if (common_data%ierr < 0) return
      NGE = NGE + 1
      ZROOT = .FALSE.
      DO 210 I = 1,NGC
 210    IF (ABS(G0(I)) .LE. 0.0D0) ZROOT = .TRUE.
      IF (.NOT. ZROOT) GO TO 260
C g has a zero at T0.  Look at g at T + (small increment). -------------
      TEMP1 = SIGN(HMING,H)
      T0 = T0 + TEMP1
      IF ((T0 - TN)*H .LT. 0.0D0) GO TO 230
      TEMP2 = TEMP1/H
      DO 220 I = 1,N
 220    Y(I) = Y(I) + TEMP2*YH(I,2)
      GO TO 240
 230  CALL DINTDY (T0, 0, YH, NYH, Y, IFLAG, common_data)
 240  CALL G (NEQ, T0, Y, NGC, G0, common_data%ierr)
      if (common_data%ierr < 0) return
      NGE = NGE + 1
      ZROOT = .FALSE.
      DO 250 I = 1,NGC
        IF (ABS(G0(I)) .GT. 0.0D0) GO TO 250
        JROOT(I) = 1
        ZROOT = .TRUE.
 250    CONTINUE
      IF (.NOT. ZROOT) GO TO 260
C g has a zero at T0 and also close to T0.  Return root. ---------------
      IRT = 1
      RETURN
C G0 has no zero components.  Proceed to check relevant interval. ------
 260  IF (TN .EQ. TLAST) GO TO 390
C
 300  CONTINUE
C Set T1 to TN or TOUTC, whichever comes first, and get g at T1. -------
      IF (ITASKC.EQ.2 .OR. ITASKC.EQ.3 .OR. ITASKC.EQ.5) GO TO 310
      IF ((TOUTC - TN)*H .GE. 0.0D0) GO TO 310
      T1 = TOUTC
      IF ((T1 - T0)*H .LE. 0.0D0) GO TO 390
      CALL DINTDY (T1, 0, YH, NYH, Y, IFLAG, common_data)
      GO TO 330
 310  T1 = TN
      DO 320 I = 1,N
 320    Y(I) = YH(I,1)
 330  CALL G (NEQ, T1, Y, NGC, G1, common_data%ierr)
      if (common_data%ierr < 0) return
      NGE = NGE + 1
C Call DROOTS to search for root in interval from T0 to T1. ------------
      JFLAG = 0
 350  CONTINUE
      CALL DROOTS (NGC, HMING, JFLAG, T0, T1, G0, G1, GX, X, JROOT,
     1             common_data)
      IF (JFLAG .GT. 1) GO TO 360
      CALL DINTDY (X, 0, YH, NYH, Y, IFLAG, common_data)
      CALL G (NEQ, X, Y, NGC, GX, common_data%ierr)
      if (common_data%ierr < 0) return
      NGE = NGE + 1
      GO TO 350
 360  T0 = X
      CALL DCOPY (NGC, GX, 1, G0, 1)
      IF (JFLAG .EQ. 4) GO TO 390
C Found a root.  Interpolate to X and return. --------------------------
      CALL DINTDY (X, 0, YH, NYH, Y, IFLAG, common_data)
      IRT = 1
      RETURN
C
 390  CONTINUE
      RETURN
C----------------------- End of Subroutine DRCHEK ----------------------
      END
*DECK DROOTS
      SUBROUTINE DROOTS (NG, HMIN, JFLAG, X0, X1, G0, G1, GX, X, JROOT,
     1                   common_data)
      use iso_c_binding, only: c_ptr, c_f_pointer, c_loc
      use odepack_common
      type(odepack_common_data), target, intent(inout) :: common_data
      INTEGER NG, JFLAG, JROOT
      DOUBLE PRECISION HMIN, X0, X1, G0, G1, GX, X
      DIMENSION G0(NG), G1(NG), GX(NG), JROOT(NG)
      INTEGER, pointer :: IOWND3(:), IMAX, LAST, IDUM3(:)
      DOUBLE PRECISION, pointer :: ALPHA, X2, RDUM3(:)
C      COMMON /DLSR01/ ALPHA, X2, RDUM3(3),
C     1   IOWND3(3), IMAX, LAST, IDUM3(4)
C-----------------------------------------------------------------------
C This subroutine finds the leftmost root of a set of arbitrary
C functions gi(x) (i = 1,...,NG) in an interval (X0,X1).  Only roots
C of odd multiplicity (i.e. changes of sign of the gi) are found.
C Here the sign of X1 - X0 is arbitrary, but is constant for a given
C problem, and -leftmost- means nearest to X0.
C The values of the vector-valued function g(x) = (gi, i=1...NG)
C are communicated through the call sequence of DROOTS.
C The method used is the Illinois algorithm.
C
C Reference:
C Kathie L. Hiebert and Lawrence F. Shampine, Implicitly Defined
C Output Points for Solutions of ODEs, Sandia Report SAND80-0180,
C February 1980.
C
C Description of parameters.
C
C NG     = number of functions gi, or the number of components of
C          the vector valued function g(x).  Input only.
C
C HMIN   = resolution parameter in X.  Input only.  When a root is
C          found, it is located only to within an error of HMIN in X.
C          Typically, HMIN should be set to something on the order of
C               100 * UROUND * MAX(ABS(X0),ABS(X1)),
C          where UROUND is the unit roundoff of the machine.
C
C JFLAG  = integer flag for input and output communication.
C
C          On input, set JFLAG = 0 on the first call for the problem,
C          and leave it unchanged until the problem is completed.
C          (The problem is completed when JFLAG .ge. 2 on return.)
C
C          On output, JFLAG has the following values and meanings:
C          JFLAG = 1 means DROOTS needs a value of g(x).  Set GX = g(X)
C                    and call DROOTS again.
C          JFLAG = 2 means a root has been found.  The root is
C                    at X, and GX contains g(X).  (Actually, X is the
C                    rightmost approximation to the root on an interval
C                    (X0,X1) of size HMIN or less.)
C          JFLAG = 3 means X = X1 is a root, with one or more of the gi
C                    being zero at X1 and no sign changes in (X0,X1).
C                    GX contains g(X) on output.
C          JFLAG = 4 means no roots (of odd multiplicity) were
C                    found in (X0,X1) (no sign changes).
C
C X0,X1  = endpoints of the interval where roots are sought.
C          X1 and X0 are input when JFLAG = 0 (first call), and
C          must be left unchanged between calls until the problem is
C          completed.  X0 and X1 must be distinct, but X1 - X0 may be
C          of either sign.  However, the notion of -left- and -right-
C          will be used to mean nearer to X0 or X1, respectively.
C          When JFLAG .ge. 2 on return, X0 and X1 are output, and
C          are the endpoints of the relevant interval.
C
C G0,G1  = arrays of length NG containing the vectors g(X0) and g(X1),
C          respectively.  When JFLAG = 0, G0 and G1 are input and
C          none of the G0(i) should be zero.
C          When JFLAG .ge. 2 on return, G0 and G1 are output.
C
C GX     = array of length NG containing g(X).  GX is input
C          when JFLAG = 1, and output when JFLAG .ge. 2.
C
C X      = independent variable value.  Output only.
C          When JFLAG = 1 on output, X is the point at which g(x)
C          is to be evaluated and loaded into GX.
C          When JFLAG = 2 or 3, X is the root.
C          When JFLAG = 4, X is the right endpoint of the interval, X1.
C
C JROOT  = integer array of length NG.  Output only.
C          When JFLAG = 2 or 3, JROOT indicates which components
C          of g(x) have a root at X.  JROOT(i) is 1 if the i-th
C          component has a root, and JROOT(i) = 0 otherwise.
C-----------------------------------------------------------------------
      INTEGER I, IMXOLD, NXLAST
      DOUBLE PRECISION T2, TMAX, FRACINT, FRACSUB, ZERO,HALF,TENTH,FIVE
      LOGICAL ZROOT, SGNCHG, XROOT
      SAVE ZERO, HALF, TENTH, FIVE
      DATA ZERO/0.0D0/, HALF/0.5D0/, TENTH/0.1D0/, FIVE/5.0D0/
C     Common block pointers
      type(DLSR01_type), pointer :: DLSR01
      type(c_ptr) :: tmp_ptr
C-----------------------------------------------------------------------
C This code associates variables with common data
C-----------------------------------------------------------------------
      DLSR01 => common_data%DLSR01

      tmp_ptr = c_loc(DLSR01%reals(1))
      call c_f_pointer(tmp_ptr,ALPHA)
      tmp_ptr = c_loc(DLSR01%reals(2))
      call c_f_pointer(tmp_ptr,X2)
      tmp_ptr = c_loc(DLSR01%reals(3))
      call c_f_pointer(tmp_ptr,RDUM3,[3])

      tmp_ptr = c_loc(DLSR01%ints(1))
      call c_f_pointer(tmp_ptr,IOWND3,[3])
      tmp_ptr = c_loc(DLSR01%ints(4))
      call c_f_pointer(tmp_ptr,IMAX)
      tmp_ptr = c_loc(DLSR01%ints(5))
      call c_f_pointer(tmp_ptr,LAST)
      tmp_ptr = c_loc(DLSR01%ints(6))
      call c_f_pointer(tmp_ptr,IDUM3,[4])

C
      IF (JFLAG .EQ. 1) GO TO 200
C JFLAG .ne. 1.  Check for change in sign of g or zero at X1. ----------
      IMAX = 0
      TMAX = ZERO
      ZROOT = .FALSE.
      DO 120 I = 1,NG
        IF (ABS(G1(I)) .GT. ZERO) GO TO 110
        ZROOT = .TRUE.
        GO TO 120
C At this point, G0(i) has been checked and cannot be zero. ------------
 110    IF (SIGN(1.0D0,G0(I)) .EQ. SIGN(1.0D0,G1(I))) GO TO 120
          T2 = ABS(G1(I)/(G1(I)-G0(I)))
          IF (T2 .LE. TMAX) GO TO 120
            TMAX = T2
            IMAX = I
 120    CONTINUE
      IF (IMAX .GT. 0) GO TO 130
      SGNCHG = .FALSE.
      GO TO 140
 130  SGNCHG = .TRUE.
 140  IF (.NOT. SGNCHG) GO TO 400
C There is a sign change.  Find the first root in the interval. --------
      XROOT = .FALSE.
      NXLAST = 0
      LAST = 1
C
C Repeat until the first root in the interval is found.  Loop point. ---
 150  CONTINUE
      IF (XROOT) GO TO 300
      IF (NXLAST .EQ. LAST) GO TO 160
      ALPHA = 1.0D0
      GO TO 180
 160  IF (LAST .EQ. 0) GO TO 170
      ALPHA = 0.5D0*ALPHA
      GO TO 180
 170  ALPHA = 2.0D0*ALPHA
 180  X2 = X1 - (X1 - X0)*G1(IMAX) / (G1(IMAX) - ALPHA*G0(IMAX))
C If X2 is too close to X0 or X1, adjust it inward, by a fractional ----
C distance that is between 0.1 and 0.5. --------------------------------
      IF (ABS(X2 - X0) < HALF*HMIN) THEN
        FRACINT = ABS(X1 - X0)/HMIN
        FRACSUB = TENTH
        IF (FRACINT .LE. FIVE) FRACSUB = HALF/FRACINT
        X2 = X0 + FRACSUB*(X1 - X0)
      ENDIF
      IF (ABS(X1 - X2) < HALF*HMIN) THEN
        FRACINT = ABS(X1 - X0)/HMIN
        FRACSUB = TENTH
        IF (FRACINT .LE. FIVE) FRACSUB = HALF/FRACINT
        X2 = X1 - FRACSUB*(X1 - X0)
      ENDIF
      JFLAG = 1
      X = X2
C Return to the calling routine to get a value of GX = g(X). -----------
      RETURN
C Check to see in which interval g changes sign. -----------------------
 200  IMXOLD = IMAX
      IMAX = 0
      TMAX = ZERO
      ZROOT = .FALSE.
      DO 220 I = 1,NG
        IF (ABS(GX(I)) .GT. ZERO) GO TO 210
        ZROOT = .TRUE.
        GO TO 220
C Neither G0(i) nor GX(i) can be zero at this point. -------------------
 210    IF (SIGN(1.0D0,G0(I)) .EQ. SIGN(1.0D0,GX(I))) GO TO 220
          T2 = ABS(GX(I)/(GX(I) - G0(I)))
          IF (T2 .LE. TMAX) GO TO 220
            TMAX = T2
            IMAX = I
 220    CONTINUE
      IF (IMAX .GT. 0) GO TO 230
      SGNCHG = .FALSE.
      IMAX = IMXOLD
      GO TO 240
 230  SGNCHG = .TRUE.
 240  NXLAST = LAST
      IF (.NOT. SGNCHG) GO TO 250
C Sign change between X0 and X2, so replace X1 with X2. ----------------
      X1 = X2
      CALL DCOPY (NG, GX, 1, G1, 1)
      LAST = 1
      XROOT = .FALSE.
      GO TO 270
 250  IF (.NOT. ZROOT) GO TO 260
C Zero value at X2 and no sign change in (X0,X2), so X2 is a root. -----
      X1 = X2
      CALL DCOPY (NG, GX, 1, G1, 1)
      XROOT = .TRUE.
      GO TO 270
C No sign change between X0 and X2.  Replace X0 with X2. ---------------
 260  CONTINUE
      CALL DCOPY (NG, GX, 1, G0, 1)
      X0 = X2
      LAST = 0
      XROOT = .FALSE.
 270  IF (ABS(X1-X0) .LE. HMIN) XROOT = .TRUE.
      GO TO 150
C
C Return with X1 as the root.  Set JROOT.  Set X = X1 and GX = G1. -----
 300  JFLAG = 2
      X = X1
      CALL DCOPY (NG, G1, 1, GX, 1)
      DO 320 I = 1,NG
        JROOT(I) = 0
        IF (ABS(G1(I)) .GT. ZERO) GO TO 310
          JROOT(I) = 1
          GO TO 320
 310    IF (SIGN(1.0D0,G0(I)) .NE. SIGN(1.0D0,G1(I))) JROOT(I) = 1
 320    CONTINUE
      RETURN
C
C No sign change in the interval.  Check for zero at right endpoint. ---
 400  IF (.NOT. ZROOT) GO TO 420
C
C Zero value at X1 and no sign change in (X0,X1).  Return JFLAG = 3. ---
      X = X1
      CALL DCOPY (NG, G1, 1, GX, 1)
      DO 410 I = 1,NG
        JROOT(I) = 0
        IF (ABS(G1(I)) .LE. ZERO) JROOT (I) = 1
 410  CONTINUE
      JFLAG = 3
      RETURN
C
C No sign changes in this interval.  Set X = X1, return JFLAG = 4. -----
 420  CALL DCOPY (NG, G1, 1, GX, 1)
      X = X1
      JFLAG = 4
      RETURN
C----------------------- End of Subroutine DROOTS ----------------------
      END
*DECK DSRCAR
      SUBROUTINE DSRCAR (RSAV, ISAV, JOB, common_data)
      use iso_c_binding, only: c_ptr, c_f_pointer, c_loc
      use odepack_common
      type(odepack_common_data), target, intent(inout) :: common_data
C-----------------------------------------------------------------------
C This routine saves or restores (depending on JOB) the contents of
C the Common blocks DLS001, DLSA01, DLSR01, which are used
C internally by one or more ODEPACK solvers.
C
C RSAV = real array of length 245 or more.
C ISAV = integer array of length 55 or more.
C JOB  = flag indicating to save or restore the Common blocks:
C        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
C        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
C        A call with JOB = 2 presumes a prior call with JOB = 1.
C-----------------------------------------------------------------------
      INTEGER ISAV, JOB
      INTEGER, pointer :: ILS(:), ILSA(:), ILSR(:)
      INTEGER I, IOFF, LENRLS, LENILS, LENRLA, LENILA, LENRLR, LENILR
      DOUBLE PRECISION RSAV
      DOUBLE PRECISION, pointer :: RLS(:), RLSA(:), RLSR(:)
      DIMENSION RSAV(*), ISAV(*)
      SAVE LENRLS, LENILS, LENRLA, LENILA, LENRLR, LENILR
C      COMMON /DLS001/ RLS(218), ILS(37)
C      COMMON /DLSA01/ RLSA(22), ILSA(9)
C      COMMON /DLSR01/ RLSR(5), ILSR(9)
      DATA LENRLS/218/, LENILS/37/, LENRLA/22/, LENILA/9/
      DATA LENRLR/5/, LENILR/9/
C     Common block pointers
      type(DLS001_type), pointer :: DLS001
      type(DLSA01_type), pointer :: DLSA01
      type(DLSR01_type), pointer :: DLSR01
      type(c_ptr) :: tmp_ptr
C-----------------------------------------------------------------------
C This code associates variables with common data
C-----------------------------------------------------------------------
      DLS001 => common_data%DLS001

      tmp_ptr = c_loc(DLS001%reals(1))
      call c_f_pointer(tmp_ptr,RLS,[218])

      tmp_ptr = c_loc(DLS001%ints(1))
      call c_f_pointer(tmp_ptr,ILS,[37])

      DLSA01 => common_data%DLSA01

      tmp_ptr = c_loc(DLSA01%reals(1))
      call c_f_pointer(tmp_ptr,RLSA,[22])

      tmp_ptr = c_loc(DLSA01%ints(1))
      call c_f_pointer(tmp_ptr,ILSA,[9])

      DLSR01 => common_data%DLSR01

      tmp_ptr = c_loc(DLSR01%reals(1))
      call c_f_pointer(tmp_ptr,RLSR,[5])

      tmp_ptr = c_loc(DLSR01%ints(1))
      call c_f_pointer(tmp_ptr,ILSR,[9])

C
      IF (JOB .EQ. 2) GO TO 100
      DO 10 I = 1,LENRLS
 10     RSAV(I) = RLS(I)
       DO 15 I = 1,LENRLA
 15     RSAV(LENRLS+I) = RLSA(I)
      IOFF = LENRLS + LENRLA
      DO 20 I = 1,LENRLR
 20     RSAV(IOFF+I) = RLSR(I)
C
      DO 30 I = 1,LENILS
 30     ISAV(I) = ILS(I)
      DO 35 I = 1,LENILA
 35     ISAV(LENILS+I) = ILSA(I)
      IOFF = LENILS + LENILA
      DO 40 I = 1,LENILR
 40     ISAV(IOFF+I) = ILSR(I)
C
      RETURN
C
 100  CONTINUE
      DO 110 I = 1,LENRLS
 110     RLS(I) = RSAV(I)
      DO 115 I = 1,LENRLA
 115     RLSA(I) = RSAV(LENRLS+I)
      IOFF = LENRLS + LENRLA
      DO 120 I = 1,LENRLR
 120     RLSR(I) = RSAV(IOFF+I)
C
      DO 130 I = 1,LENILS
 130     ILS(I) = ISAV(I)
      DO 135 I = 1,LENILA
 135     ILSA(I) = ISAV(LENILS+I)
      IOFF = LENILS + LENILA
      DO 140 I = 1,LENILR
 140     ILSR(I) = ISAV(IOFF+I)
C
      RETURN
C----------------------- End of Subroutine DSRCAR ----------------------
      END