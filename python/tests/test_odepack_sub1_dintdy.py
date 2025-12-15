import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "python"))
from odepack_sub1 import CommonData, dintdy  # noqa: E402


@pytest.fixture(scope="session")
def fortran_reference(run_dir_dintdy):
    bin_path = run_dir_dintdy / "dintdy_ref.bin"
    data = memoryview(bin_path.read_bytes())
    ref = []
    pos = 0
    for _ in range(3):
        k = int(np.frombuffer(data[pos : pos + 4], dtype=np.int32)[0]); pos += 4
        iflag = int(np.frombuffer(data[pos : pos + 4], dtype=np.int32)[0]); pos += 4
        dky = np.frombuffer(data[pos : pos + 16], dtype=np.float64).copy(); pos += 16
        ref.append((k, iflag, dky))
    return ref


def test_dintdy_matches_fortran(fortran_reference):
    common = CommonData(0, 0, "")
    common.DLS001_ints[18] = 3   # L
    common.DLS001_ints[31] = 2   # N
    common.DLS001_ints[32] = 2   # NQ
    common.DLS001_reals[211] = 0.5
    common.DLS001_reals[214] = 0.5
    common.DLS001_reals[216] = 1.0
    common.DLS001_reals[217] = 1e-16

    yh = np.zeros((2, 3), dtype=np.float64)
    yh[:, 0] = np.array([1.0, 2.0])
    yh[:, 1] = np.array([0.1, 0.2])
    yh[:, 2] = np.array([0.01, 0.02])

    for k_ref, iflag_ref, dky_ref in fortran_reference:
        dky = np.zeros(2, dtype=np.float64)
        iflag = np.zeros(1, dtype=np.int64)
        dintdy(1.0, k_ref, yh, 2, dky, iflag, common)
        assert iflag[0] == iflag_ref
        np.testing.assert_allclose(dky, dky_ref, rtol=0.0, atol=1e-12)
