import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "python"))
from odepack_sub1 import dewset  # noqa: E402


def test_dewset_matches_fortran(run_dir_dewset):
    bin_path = run_dir_dewset / "dewset.bin"
    data = memoryview(bin_path.read_bytes())
    pos = 0

    n = 4
    rtol = np.array([1.0e-3, 2.0e-3, 3.0e-3, 4.0e-3], dtype=np.float64)
    atol = np.array([1.0e-6, 2.0e-6, 3.0e-6, 4.0e-6], dtype=np.float64)
    ycur = np.array([1.0, -2.0, 0.5, -1.5], dtype=np.float64)

    for itol in (1, 2, 3, 4):
        ewt_ref = np.frombuffer(data[pos : pos + n * 8], dtype=np.float64).copy()
        pos += n * 8
        ewt = np.zeros(n, dtype=np.float64)
        dewset(n, itol, rtol, atol, ycur, ewt)
        np.testing.assert_allclose(ewt, ewt_ref, rtol=0.0, atol=1e-14)
