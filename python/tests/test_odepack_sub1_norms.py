import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "python"))
from odepack_sub1 import dbnorm, dfnorm, dmnorm  # noqa: E402


def test_norms_match_fortran(run_dir_norms):
    bin_path = run_dir_norms / "norms.bin"
    data = memoryview(bin_path.read_bytes())
    pos = 0
    dm_ref = float(np.frombuffer(data[pos : pos + 8], dtype=np.float64)[0]); pos += 8
    df_ref = float(np.frombuffer(data[pos : pos + 8], dtype=np.float64)[0]); pos += 8
    db_ref = float(np.frombuffer(data[pos : pos + 8], dtype=np.float64)[0]); pos += 8

    n = 4
    ml = 1
    mu = 1
    nra = ml + mu + 1
    v = np.array([1.0, -2.0, 0.5, -1.5], dtype=np.float64)
    w = np.array([1.0, 0.5, 2.0, 1.5], dtype=np.float64)
    A = np.array(
        [
            [2.0, -1.0, 0.0, 0.0],
            [-1.0, 2.0, -1.0, 0.0],
            [0.0, -1.0, 2.0, -1.0],
            [0.0, 0.0, -1.0, 2.0],
        ],
        dtype=np.float64,
    )
    ab = np.zeros((nra, n), dtype=np.float64)
    for j in range(n):
        for i in range(max(0, j - mu), min(n, j + ml + 1)):
            ab[mu + i - j, j] = A[i, j]

    np.testing.assert_allclose(dmnorm(n, v, w), dm_ref, rtol=0.0, atol=1e-14)
    np.testing.assert_allclose(dfnorm(n, A, w), df_ref, rtol=0.0, atol=1e-14)
    np.testing.assert_allclose(dbnorm(n, ab, nra, ml, mu, w), db_ref, rtol=0.0, atol=1e-14)
