import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

# Ensure we can import the Python implementation without installing a package.
REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "python"))
from odepack_sub1 import dcfode, dumach


def _check_tools():
    if shutil.which("cmake") is None:
        pytest.skip("cmake not available; skipping Fortran reference build")
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran not available; skipping Fortran reference build")


@pytest.fixture(scope="session")
def fortran_reference(tmp_path_factory):
    _check_tools()
    build_dir = tmp_path_factory.mktemp("cmake_build")
    src_dir = REPO_ROOT

    # Configure and build the dedicated driver target using the project's CMake.
    subprocess.run(
        ["cmake", "-S", str(src_dir), "-B", str(build_dir), "-DENABLE_PYTHON_TESTS=ON"],
        check=True,
    )
    subprocess.run(["cmake", "--build", str(build_dir), "--target", "odepack_sub1_driver"], check=True)

    candidates = [
        Path(build_dir) / "python" / "tests" / "odepack_sub1_driver",
        Path(build_dir) / "test" / "odepack_sub1_driver",
        Path(build_dir) / "odepack_sub1_driver",
    ]
    exe_path = None
    for cand in candidates:
        if sys.platform.startswith("win"):
            cand_exe = cand.with_suffix(".exe")
            if cand_exe.exists():
                exe_path = cand_exe
                break
        if cand.exists():
            exe_path = cand
            break
    if exe_path is None:
        pytest.skip("odepack_sub1_driver not built or not found")

    # Run the driver; it writes odepack_sub1.bin in its working directory.
    subprocess.run([str(exe_path)], check=True, cwd=exe_path.parent)

    bin_path = exe_path.parent / "odepack_sub1.bin"
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
