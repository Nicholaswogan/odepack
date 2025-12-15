import sys
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

# Ensure we can import the Python implementation without installing a package.
REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "python"))
from odepack_sub1 import dcfode, dumach  # noqa: E402


@pytest.fixture(scope="session")
def fortran_reference(run_dir_sub1):
    bin_path = Path(run_dir_sub1) / "odepack_sub1.bin"
    raw = np.fromfile(bin_path, dtype=np.float64)
    expected_len = 1 + 156 + 36 + 156 + 36  # dumach + elco1 + tesco1 + elco2 + tesco2
    assert raw.size == expected_len, f"Unexpected reference size {raw.size}, expected {expected_len}"

    offset = 0
    f_dumach = raw[offset]
    offset += 1
    elco1 = raw[offset : offset + 156].reshape((13, 12), order="F")
    offset += 156
    tesco1 = raw[offset : offset + 36].reshape((3, 12), order="F")
    offset += 36
    elco2 = raw[offset : offset + 156].reshape((13, 12), order="F")
    offset += 156
    tesco2 = raw[offset : offset + 36].reshape((3, 12), order="F")

    return {
        "dumach": f_dumach,
        "elco1": elco1,
        "tesco1": tesco1,
        "elco2": elco2,
        "tesco2": tesco2,
    }


def test_dumach_matches_fortran(fortran_reference):
    py_val = dumach()
    f_val = fortran_reference["dumach"]
    assert_allclose(py_val, f_val, rtol=0.0, atol=0.0)


@pytest.mark.parametrize("meth,max_order", [(1, 12), (2, 5)])
def test_dcfode_matches_fortran(fortran_reference, meth, max_order):
    elco_py = np.full((13, 12), 7.7, dtype=np.float64, order="C")
    tesco_py = np.full((3, 12), 8.8, dtype=np.float64, order="C")
    dcfode(meth, elco_py, tesco_py)

    if meth == 1:
        elco_ref, tesco_ref = fortran_reference["elco1"], fortran_reference["tesco1"]
    else:
        elco_ref, tesco_ref = fortran_reference["elco2"], fortran_reference["tesco2"]

    col_limit = max_order
    assert_allclose(elco_py[:, :col_limit], elco_ref[:, :col_limit], rtol=0.0, atol=1e-15)
    assert_allclose(tesco_py[:, :col_limit], tesco_ref[:, :col_limit], rtol=0.0, atol=1e-15)
