import numpy as np
from numba import njit, float64, int64, void, types
from numba.experimental import jitclass

from lapack import dgetrs, dgbtrs

@njit(float64())
def dumach():
    """
    Compute the unit roundoff (smallest u such that 1.0 + u != 1.0)
    using the same halving loop as the Fortran DUMACH.
    """
    u = 1.0
    while True:
        u *= 0.5
        if 1.0 + u == 1.0:
            break
    return u * 2.0


@njit(void(int64, float64[:, ::1], float64[:, ::1]))
def dcfode(meth, elco, tesco):
    """
    In-place version of DCFODE: fill ELCO and TESCO for the given method.

    Parameters
    ----------
    meth : int
        1 for Adams (max order 12), 2 for BDF (max order 5).
    elco : (13, 12) float64 ndarray
    tesco: (3, 12) float64 ndarray
    """
    if meth == 1:
        pc = np.zeros(12, dtype=np.float64)
        elco[0, 0] = 1.0
        elco[1, 0] = 1.0
        tesco[0, 0] = 0.0
        tesco[1, 0] = 2.0
        tesco[0, 1] = 1.0
        tesco[2, 11] = 0.0
        pc[0] = 1.0
        rqfac = 1.0

        for nq in range(2, 13):  # orders 2..12
            rq1fac = rqfac
            rqfac = rqfac / nq
            nqm1 = nq - 1
            fnqm1 = float(nqm1)
            nqp1 = nq + 1

            pc[nq - 1] = 0.0
            for ib in range(1, nqm1 + 1):
                i = nqp1 - ib  # Fortran 1-based
                idx = i - 1
                pc[idx] = pc[idx - 1] + fnqm1 * pc[idx]
            pc[0] = fnqm1 * pc[0]

            pint = pc[0]
            xpin = pc[0] * 0.5
            tsign = 1.0
            for i in range(2, nq + 1):
                tsign = -tsign
                pint += tsign * pc[i - 1] / i
                xpin += tsign * pc[i - 1] / (i + 1)

            elco[0, nq - 1] = pint * rq1fac
            elco[1, nq - 1] = 1.0
            for i in range(2, nq + 1):
                elco[i, nq - 1] = rq1fac * pc[i - 1] / i

            agamq = rqfac * xpin
            ragq = 1.0 / agamq
            tesco[1, nq - 1] = ragq
            if nq < 12:
                tesco[0, nqp1 - 1] = ragq * rqfac / nqp1
            tesco[2, nqm1 - 1] = ragq

    elif meth == 2:
        pc = np.zeros(12, dtype=np.float64)
        pc[0] = 1.0
        rq1fac = 1.0

        for nq in range(1, 6):  # orders 1..5
            fnq = float(nq)
            nqp1 = nq + 1

            pc[nqp1 - 1] = 0.0
            for ib in range(1, nq + 1):
                i = nq + 2 - ib
                idx = i - 1
                pc[idx] = pc[idx - 1] + fnq * pc[idx]
            pc[0] = fnq * pc[0]

            for i in range(1, nqp1 + 1):
                elco[i - 1, nq - 1] = pc[i - 1] / pc[1]
            elco[1, nq - 1] = 1.0

            tesco[0, nq - 1] = rq1fac
            tesco[1, nq - 1] = nqp1 / elco[0, nq - 1]
            tesco[2, nq - 1] = (nq + 2) / elco[0, nq - 1]
            rq1fac = rq1fac / fnq

    else:
        raise ValueError("meth must be 1 (Adams) or 2 (BDF)")


# Common-data structure mirroring odepack_common_data derived type.
common_spec = [
    ("iprint", int64),
    ("error_message", types.unicode_type),
    ("ierr", int64),
    ("DLS001_reals", float64[::1]),
    ("DLS001_ints", int64[::1]),
    ("DLSA01_reals", float64[::1]),
    ("DLSA01_ints", int64[::1]),
    ("DLSR01_reals", float64[::1]),
    ("DLSR01_ints", int64[::1]),
]


@jitclass(common_spec)
class CommonData:
    def __init__(self, iprint, ierr, error_message):
        self.iprint = iprint
        self.error_message = error_message
        self.ierr = ierr
        self.DLS001_reals = np.zeros(218, dtype=np.float64)
        self.DLS001_ints = np.zeros(37, dtype=np.int64)
        self.DLSA01_reals = np.zeros(22, dtype=np.float64)
        self.DLSA01_ints = np.zeros(9, dtype=np.int64)
        self.DLSR01_reals = np.zeros(5, dtype=np.float64)
        self.DLSR01_ints = np.zeros(9, dtype=np.int64)


@njit(void(float64, int64, float64[:, ::1], int64, float64[:], int64[:], CommonData.class_type.instance_type))
def dintdy(t, k, yh, nyh, dky, iflag, common):
    """
    Numba version of DINTDY.
    iflag is a length-1 int64 array for output status.
    """
    # map Fortran common layout offsets (1-based in Fortran, 0-based here)
    dls_reals = common.DLS001_reals
    dls_ints = common.DLS001_ints

    h = dls_reals[211]    # reals(212)
    hu = dls_reals[214]   # reals(215)
    tn = dls_reals[216]   # reals(217)
    uround = dls_reals[217]  # reals(218)

    l = dls_ints[18]   # ints(19)
    n = dls_ints[31]   # ints(32)
    nq = dls_ints[32]  # ints(33)

    iflag[0] = 0
    if k < 0 or k > nq:
        iflag[0] = -1
        return

    sign_hu = 1.0 if hu >= 0.0 else -1.0
    tp = tn - hu - 100.0 * uround * sign_hu * (abs(tn) + abs(hu))
    if (t - tp) * (t - tn) > 0.0:
        iflag[0] = -2
        return

    s = (t - tn) / h
    ic = 1
    if k != 0:
        for jj in range(l - k, nq + 1):
            ic *= jj
    c = ic

    for i in range(n):
        dky[i] = c * yh[i, l - 1]

    if k != nq:
        jb2 = nq - k
        for jb in range(1, jb2 + 1):
            j = nq - jb
            jp1 = j + 1
            ic = 1
            if k != 0:
                for jj in range(jp1 - k, j + 1):
                    ic *= jj
            c = ic
            for i in range(n):
                dky[i] = c * yh[i, jp1 - 1] + s * dky[i]

    if k != 0:
        r = h ** (-k)
        for i in range(n):
            dky[i] = r * dky[i]


@njit
def dsolsy(wm, iwm, x, tem, common):
    """
    Numba version of DSOLSY.

    Parameters mirror Fortran:
    wm : 1D float64 work array
    iwm: 1D int64 work array
    x  : RHS on input, solution on output
    tem: workspace (unused)
    common: CommonData jitclass
    """
    el0 = common.DLS001_reals[210]
    h = common.DLS001_reals[211]
    miter = common.DLS001_ints[26]
    n = common.DLS001_ints[31]

    ml = iwm[0]
    mu = iwm[1]
    ldab = 2 * ml + mu + 1

    if miter == 3:
        phl0 = wm[1]
        hl0 = h * el0
        wm[1] = hl0
        if hl0 != phl0:
            r = hl0 / phl0
            for i in range(n):
                di = 1.0 - r * (1.0 - 1.0 / wm[2 + i])
                wm[2 + i] = 1.0 / di
        for i in range(n):
            x[i] = wm[2 + i] * x[i]
        common.DLS001_ints[14] = 0  # IERSL
        return

    if miter in (1, 2):
        # wm stores the LU factors in column-major order starting at offset 2.
        a = np.empty((n, n), dtype=np.float64)
        idx = 0
        for j in range(n):
            for i in range(n):
                a[i, j] = wm[2 + idx]
                idx += 1
        piv = iwm[20 : 20 + n]
        dgetrs(a, piv, x)
        common.DLS001_ints[14] = 0
        return

    if miter in (4, 5):
        # Band LU is stored column-major in wm.
        ab = np.empty((ldab, n), dtype=np.float64)
        idx = 0
        for j in range(n):
            for i in range(ldab):
                ab[i, j] = wm[2 + idx]
                idx += 1
        piv = iwm[20 : 20 + n]
        dgbtrs(ab, n, n, ml, mu, piv, x, ldab)
        common.DLS001_ints[14] = 0
        return

    common.DLS001_ints[14] = -1
