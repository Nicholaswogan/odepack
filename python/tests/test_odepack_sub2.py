import sys
from io import StringIO
from pathlib import Path
from contextlib import redirect_stdout

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "python"))
from odepack_sub2 import ixsav, xerrwd  # noqa: E402


def _normalize_line(line):
    toks = line.strip().replace(",", " ").split()
    norm = []
    for t in toks:
        t_repl = t.replace("D", "E")
        try:
            val = float(t_repl)
            if val.is_integer():
                norm.append(int(val))
            else:
                norm.append(val)
            continue
        except ValueError:
            pass
        norm.append(t.upper())
    return norm


def _normalize_output(lines):
    return [_normalize_line(line) for line in lines if line.strip() != ""]


@pytest.fixture(scope="session")
def fortran_output(run_dir_sub2):
    txt_path = run_dir_sub2 / "odepack_sub2.txt"
    return txt_path.read_text().splitlines()


def test_xerrwd_matches_fortran(fortran_output):
    buf = StringIO()
    with redirect_stdout(buf):
        xerrwd("TEST MSG", 8, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
        xerrwd("WITH I", 6, 0, 0, 1, 123, 0, 0, 0.0, 0.0)
        xerrwd("WITH I12", 7, 0, 0, 2, 12, 34, 0, 0.0, 0.0)
        xerrwd("WITH R", 6, 0, 0, 0, 0, 0, 1, 1.23456789, 0.0)
        xerrwd("WITH R12", 7, 0, 0, 0, 0, 0, 2, 1.0, -2.5)

    py_norm = _normalize_output(buf.getvalue().splitlines())
    f_norm = _normalize_output(fortran_output)

    assert py_norm == f_norm
