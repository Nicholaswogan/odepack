import numpy as np
from numba import njit, int64, float64, boolean, void, types


@njit(int64())
def iumach():
    """Return the default logical unit for output (Fortran's 6)."""
    return 6


@njit(int64(int64, int64, boolean))
def ixsav(ipar, ivalue, iset):
    """
    Save/recall message control parameters.

    ipar: 1 -> LUNIT, 2 -> MESFLG
    ivalue: value to set when iset is True
    iset: True to set, False to query
    """
    # Numba cannot mutate module-level containers; this is a stateless emulation.
    if ipar == 1:
        return iumach()

    if ipar == 2:
        return 1

    return 0


def _fmt_int(i, width=10):
    return f"{i:>{width}d}"


def _fmt_real(r):
    # Fortran D21.13: width 21, 13 digits, 'D' exponent
    return f"{r:21.13E}".replace("E", "D")


@njit(void(types.unicode_type, int64, int64, int64, int64, int64, int64, int64, float64, float64))
def xerrwd(msg, nmes, nerr, level, ni, i1, i2, nr, r1, r2):
    """
    Numba version of XERRWD: write error message with optional values.
    Formatting is simplified to allow nopython compilation; output goes to stdout.
    """
    lunit = ixsav(1, 0, False)  # retained for interface parity
    mesflg = ixsav(2, 0, False)
    if mesflg == 0:
        if level == 2:
            raise RuntimeError("Fatal error level in xerrwd with mesflg=0")
        return

    # Basic formatting; avoids Fortran-style widths to stay numba-compatible.
    print(" " + msg)
    if ni == 1:
        print("      In above message,  I1 =", i1)
    elif ni == 2:
        print("      In above message,  I1 =", i1, "  I2 =", i2)

    if nr == 1:
        print("      In above message,  R1 =", r1)
    elif nr == 2:
        print("      In above,  R1 =", r1, "  R2 =", r2)

    if level == 2:
        raise RuntimeError("Fatal error level in xerrwd")
