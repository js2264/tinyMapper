"""End-to-end integration tests for all tinyMapper modes.

Each test invokes ``tinymapper`` via ``micromamba run -n autotinymapper_tm``
and verifies the command exits cleanly.  The entire suite is skipped when the
``autotinymapper_tm`` micromamba environment is not available.

Run from the repo root so that relative paths like ``tests/testChIP`` resolve
correctly against the bundled test FASTQ files.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).parent.parent
TESTS_DIR = REPO_ROOT / "tests"
ENV_NAME = "autotinymapper_tm"


# ---------------------------------------------------------------------------
# Skip marker
# ---------------------------------------------------------------------------


def _env_available() -> bool:
    """Return True if the autotinymapper_tm micromamba env exists and has tinymapper."""
    micromamba = shutil.which("micromamba")
    if not micromamba:
        return False
    result = subprocess.run(
        [micromamba, "run", "-n", ENV_NAME, "tinymapper", "--version"],
        capture_output=True,
        timeout=30,
    )
    return result.returncode == 0


requires_env = pytest.mark.skipif(
    not _env_available(),
    reason=f"micromamba env '{ENV_NAME}' not available",
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _run(args: list[str], output_dir: Path) -> subprocess.CompletedProcess[str]:
    """Run tinymapper with *args* via the autotinymapper_tm env."""
    cmd = ["micromamba", "run", "-n", ENV_NAME, "tinymapper"] + args + ["--output", str(output_dir)]
    print(f"Running command: {' '.join(cmd)}")
    return subprocess.run(
        cmd,
        capture_output=False,
        cwd=str(REPO_ROOT),
        timeout=600,
    )


def _sample(name: str) -> str:
    return str(TESTS_DIR / name)


def _genome(name: str) -> str:
    return str(TESTS_DIR / name)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@requires_env
def test_chip_sample_only(tmp_path):
    result = _run(
        ["-m", "ChIP", "-s", _sample("testChIP"), "-g", _genome("R64-1-1/R64-1-1")],
        tmp_path,
    )
    assert result.returncode == 0, "ChIP (sample only) failed"


@requires_env
def test_chip_with_input(tmp_path):
    result = _run(
        [
            "-m",
            "ChIP",
            "--sample",
            _sample("testChIP"),
            "--input",
            _sample("testChIP.input"),
            "--genome",
            _genome("R64-1-1/R64-1-1"),
            "-k",
        ],
        tmp_path,
    )
    assert result.returncode == 0, "ChIP (sample + input) failed"


@requires_env
def test_chip_with_input_and_calibration(tmp_path):
    result = _run(
        [
            "-m",
            "ChIP",
            "--sample",
            _sample("testChIP"),
            "--input",
            _sample("testChIP.input"),
            "--genome",
            _genome("R64-1-1/R64-1-1"),
            "--calibration",
            _genome("CBS138/CBS138"),
            "-k",
        ],
        tmp_path,
    )
    assert result.returncode == 0, "ChIP (sample + input + calibration) failed"


@requires_env
def test_rna(tmp_path):
    result = _run(
        ["-m", "RNA", "-s", _sample("testRNA"), "-g", _genome("R64-1-1/R64-1-1")],
        tmp_path,
    )
    assert result.returncode == 0, "RNA failed"


@requires_env
def test_mnase(tmp_path):
    result = _run(
        [
            "-m",
            "MNase",
            "-s",
            _sample("testMNase"),
            "-g",
            _genome("R64-1-1/R64-1-1"),
            "--MNaseSizes",
            "70,250",
        ],
        tmp_path,
    )
    assert result.returncode == 0, "MNase failed"


@requires_env
def test_atac_paired(tmp_path):
    result = _run(
        ["-m", "ATAC", "-s", _sample("testATAC"), "-g", _genome("R64-1-1/R64-1-1")],
        tmp_path,
    )
    assert result.returncode == 0, "ATAC (paired-end) failed"


@requires_env
def test_atac_single_end(tmp_path):
    result = _run(
        ["-m", "ATAC", "-s", _sample("testATAC_se"), "-g", _genome("R64-1-1/R64-1-1")],
        tmp_path,
    )
    assert result.returncode == 0, "ATAC (single-end) failed"


@requires_env
def test_hic(tmp_path):
    result = _run(
        [
            "-m",
            "HiC",
            "-s",
            _sample("testHiC"),
            "-g",
            _genome("R64-1-1/R64-1-1"),
            "--binning",
            "2000",
            "--restriction",
            "DpnII,HinfI",
        ],
        tmp_path,
    )
    assert result.returncode == 0, "HiC failed"


@requires_env
def test_shotgun(tmp_path):
    result = _run(
        ["-m", "shotgun", "-s", _sample("testShotgun"), "-g", _genome("R64-1-1/R64-1-1")],
        tmp_path,
    )
    assert result.returncode == 0, "shotgun failed"
