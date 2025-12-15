import numpy as np
from numba import njit, int64, float64


@njit
def dgetrf(a, ipiv):
    """
    In-place LU factorization with partial pivoting (Numba).

    Parameters
    ----------
    a : 2D float64 array (Fortran-order preferred), shape (m, n)
    ipiv : 1D int64 array, length min(m, n)

    Returns
    -------
    info : int64 (0 on success; k if U[k,k] was zero)
    """
    m, n = a.shape
    info = 0
    kmax = min(m, n)
    for k in range(kmax):
        # Find pivot in column k
        piv = k
        pivmax = abs(a[k, k])
        for i in range(k + 1, m):
            if abs(a[i, k]) > pivmax:
                pivmax = abs(a[i, k])
                piv = i
        ipiv[k] = piv + 1  # 1-based
        if a[piv, k] == 0.0:
            if info == 0:
                info = k + 1
            continue
        # Swap rows if needed
        if piv != k:
            for j in range(n):
                tmp = a[k, j]
                a[k, j] = a[piv, j]
                a[piv, j] = tmp
        # Compute multipliers and update trailing submatrix
        for i in range(k + 1, m):
            a[i, k] = a[i, k] / a[k, k]
            for j in range(k + 1, n):
                a[i, j] -= a[i, k] * a[k, j]
    return info


@njit
def dgetrs(a, ipiv, b):
    """
    Solve A*x = b using LU factors from dgetrf (no transpose).

    Parameters
    ----------
    a : LU factors from dgetrf (2D float64 array, shape (n, n))
    ipiv : pivot array (1D int64, length n)
    b : RHS on input, overwritten with solution (1D float64 length n)
    """
    n = a.shape[0]
    # Apply pivots to RHS
    for i in range(n):
        piv = ipiv[i] - 1
        if piv != i:
            tmp = b[i]
            b[i] = b[piv]
            b[piv] = tmp
    # Forward solve Ly = Pb (L unit diagonal)
    for i in range(n):
        for j in range(i):
            b[i] -= a[i, j] * b[j]
    # Backward solve Ux = y
    for i in range(n - 1, -1, -1):
        for j in range(i + 1, n):
            b[i] -= a[i, j] * b[j]
        b[i] /= a[i, i]


@njit
def dgbtrf(ab, m, n, kl, ku, ipiv, ldab):
    """
    LU factorization of a band matrix in LAPACK storage.

    To keep implementation simple and robust, this converts to dense,
    factors with dgetrf, then packs the LU back into band storage.

    Parameters
    ----------
    ab : 2D float64 array, shape (ldab, n), band storage (modified in-place)
    m, n : matrix dimensions
    kl, ku : lower/upper bandwidths
    ipiv : 1D int64 pivot array, length min(m, n)
    ldab : leading dimension (>= 2*kl+ku+1)

    Returns
    -------
    info : int64
    """
    # Form dense copy
    Ad = np.zeros((m, n), dtype=np.float64)
    for j in range(n):
        i_start = max(0, j - ku)
        i_end = min(m - 1, j + kl)
        for i in range(i_start, i_end + 1):
            Ad[i, j] = ab[ku + i - j, j]
    info = dgetrf(Ad, ipiv)
    # Pack LU back into band storage (overwrite ab)
    ab[:] = 0.0
    for j in range(n):
        i_start = max(0, j - ku)
        i_end = min(m - 1, j + kl)
        for i in range(i_start, i_end + 1):
            ab[ku + i - j, j] = Ad[i, j]
    return info


@njit
def dgbtrs(ab, m, n, kl, ku, ipiv, b, ldab):
    """
    Solve A*x = b for band matrix using LU from dgbtrf (no transpose).

    Parameters
    ----------
    ab : LU factors in band storage (ldab x n)
    m, n : dimensions (assume m==n)
    kl, ku : bandwidths
    ipiv : pivot array from dgbtrf
    b : RHS on input, overwritten with solution
    ldab : leading dimension
    """
    # Reconstruct dense LU from band storage
    Ad = np.zeros((m, n), dtype=np.float64)
    for j in range(n):
        i_start = max(0, j - ku)
        i_end = min(m - 1, j + kl)
        for i in range(i_start, i_end + 1):
            Ad[i, j] = ab[ku + i - j, j]
    # Use dense solve
    dgetrs(Ad, ipiv, b)
