import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
_build_note = None


def _check_tools():
    missing = []
    if shutil.which("cmake") is None:
        missing.append("cmake")
    if shutil.which("gfortran") is None:
        missing.append("gfortran")
    if missing:
        raise RuntimeError(f"Required build tools missing: {', '.join(missing)}")


@pytest.fixture(scope="session")
def build_dir(tmp_path_factory, request):
    _check_tools()
    build_dir = tmp_path_factory.mktemp("cmake_build")
    subprocess.run(
        ["cmake", "-S", str(REPO_ROOT), "-B", str(build_dir), "-DENABLE_PYTHON_TESTS=ON"],
        check=True,
    )
    subprocess.run(
        [
            "cmake",
            "--build",
            str(build_dir),
            "--target",
            "odepack_sub1_driver",
            "odepack_sub2_driver",
            "odepack_sub1_dintdy_driver",
            "lapack_drivers",
        ],
        check=True,
    )
    return Path(build_dir)


def _find_exe(build_dir: Path, target: str) -> Path:
    candidates = [
        build_dir / "python" / "tests" / target,
        build_dir / "test" / target,
        build_dir / target,
    ]
    for cand in candidates:
        cand_exe = cand.with_suffix(".exe") if sys.platform.startswith("win") else cand
        if cand_exe.exists():
            return cand_exe
    raise FileNotFoundError(f"{target} not found in build dir {build_dir}")


@pytest.fixture(scope="session")
def run_dir_sub2(build_dir):
    exe_path = _find_exe(build_dir, "odepack_sub2_driver")
    subprocess.run([str(exe_path)], check=True, cwd=exe_path.parent)
    return exe_path.parent


@pytest.fixture(scope="session")
def run_dir_dintdy(build_dir):
    exe_path = _find_exe(build_dir, "odepack_sub1_dintdy_driver")
    subprocess.run([str(exe_path)], check=True, cwd=exe_path.parent)
    return exe_path.parent


@pytest.fixture(scope="session")
def run_dir_lapack(build_dir):
    exe_path = _find_exe(build_dir, "lapack_drivers")
    subprocess.run([str(exe_path)], check=True, cwd=exe_path.parent)
    return exe_path.parent


@pytest.fixture(scope="session")
def run_dir_sub1(build_dir):
    exe_path = _find_exe(build_dir, "odepack_sub1_driver")
    subprocess.run([str(exe_path)], check=True, cwd=exe_path.parent)
    return exe_path.parent


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    if _build_note:
        terminalreporter.write_line(_build_note)
