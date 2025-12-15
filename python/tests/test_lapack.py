import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "python"))
from lapack import dgbtrs, dgbtrf, dgetrf, dgetrs  # noqa: E402


def test_dgetrs_matches_fortran(run_dir_lapack):
    bin_path = run_dir_lapack / "dgetrs_ref.bin"
    data = memoryview(bin_path.read_bytes())
    pos = 0
    info = int(np.frombuffer(data[pos : pos + 4], dtype=np.int32)[0]); pos += 4
    assert info == 0
    ipiv = np.frombuffer(data[pos : pos + 3 * 4], dtype=np.int32).astype(np.int64); pos += 3 * 4
    b_ref = np.frombuffer(data[pos : pos + 3 * 8], dtype=np.float64).copy()

    a = np.array([[4.0, 2.0, 0.0], [1.0, 5.0, 2.0], [0.0, 1.0, 3.0]], dtype=np.float64)
    ipiv_py = np.zeros_like(ipiv)
    a_lu = a.copy(order="F")
    info_py = dgetrf(a_lu, ipiv_py)
    assert info_py == 0
    x = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    dgetrs(a_lu, ipiv_py, x)
    np.testing.assert_allclose(x, b_ref, rtol=0.0, atol=1e-12)


def test_dgbtrs_matches_fortran(run_dir_lapack):
    bin_path = run_dir_lapack / "dgbtrs_ref.bin"
    data = memoryview(bin_path.read_bytes())
    pos = 0
    info = int(np.frombuffer(data[pos : pos + 4], dtype=np.int32)[0]); pos += 4
    assert info == 0
    b_ref = np.frombuffer(data[pos : pos + 4 * 8], dtype=np.float64).copy()

    n = 4
    kl = 1
    ku = 1
    ldab = 2 * kl + ku + 1
    A = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        A[i, i] = 2.0
        if i > 0:
            A[i, i - 1] = -1.0
            A[i - 1, i] = -1.0
    ab_py = np.zeros((ldab, n), dtype=np.float64, order="F")
    for j in range(n):
        i_start = max(0, j - ku)
        i_end = min(n - 1, j + kl)
        for i in range(i_start, i_end + 1):
            ab_py[ku + i - j, j] = A[i, j]
    ipiv_py = np.zeros(n, dtype=np.int64)
    info_py = dgbtrf(ab_py, n, n, kl, ku, ipiv_py, ldab)
    assert info_py == 0
    x = np.array([1.0, 0.0, 0.0, 1.0], dtype=np.float64)
    dgbtrs(ab_py, n, n, kl, ku, ipiv_py, x, ldab)
    np.testing.assert_allclose(x, b_ref, rtol=0.0, atol=1e-12)
