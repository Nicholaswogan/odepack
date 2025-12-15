import numpy as np
from numba import njit, float64, int64, void


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
