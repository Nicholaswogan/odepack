import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "python"))
from odepack_sub1 import CommonData, dsolsy  # noqa: E402


def test_dsolsy_matches_fortran(run_dir_dsolsy):
    bin_path = run_dir_dsolsy / "dsolsy_dense.bin"
    data = memoryview(bin_path.read_bytes())
    pos = 0
    info = int(np.frombuffer(data[pos : pos + 4], dtype=np.int32)[0]); pos += 4
    assert info == 0
    ipiv = np.frombuffer(data[pos : pos + 3 * 4], dtype=np.int32).astype(np.int64); pos += 3 * 4
    lu = np.frombuffer(data[pos : pos + 3 * 3 * 8], dtype=np.float64).reshape((3, 3), order="F"); pos += 3 * 3 * 8
    b_ref = np.frombuffer(data[pos : pos + 3 * 8], dtype=np.float64).copy()

    n = 3
    wm = np.zeros(2 + n * n, dtype=np.float64)
    iwm = np.zeros(20 + n, dtype=np.int64)
    wm[2:] = lu.ravel(order="F")
    iwm[20 : 20 + n] = ipiv
    common = CommonData(0, 0, "")
    common.DLS001_reals[210] = 1.0  # EL0
    common.DLS001_reals[211] = 1.0  # H
    common.DLS001_ints[26] = 1      # MITER dense
    common.DLS001_ints[31] = n
    x = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    tem = np.zeros_like(x)
    dsolsy(wm, iwm, x, tem, common)
    np.testing.assert_allclose(x, b_ref, rtol=0.0, atol=1e-12)
