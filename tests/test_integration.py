"""End-to-end integration tests for all tinyMapper modes.

Each test invokes ``tinymapper`` and verifies the command exits cleanly.
The suite runs if *either* of the following conditions is met:

1. The ``autotinymapper_tm`` micromamba environment is available and contains
   a working ``tinymapper`` installation.
2. All CLI dependencies listed in ``env/tinymapper.yaml`` are present on
   ``PATH`` (i.e. the current environment already has them installed).

The entire suite is skipped when neither condition is satisfied.

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
TESTS_RESULTS_DIR = REPO_ROOT / "tests_results"
ENV_NAME = "autotinymapper_tm"


# ---------------------------------------------------------------------------
# Purge results directory before running tests
# ---------------------------------------------------------------------------
@pytest.fixture(scope="session", autouse=True)
def purge_results():
    """Delete the tests_results directory before running any tests."""
    if TESTS_RESULTS_DIR.exists():
        shutil.rmtree(TESTS_RESULTS_DIR)


# ---------------------------------------------------------------------------
# Skip marker
# ---------------------------------------------------------------------------

# Conda package name → executable(s) to check on PATH.
# At least one executable per entry must be found.
_CONDA_TO_BINS: dict[str, list[str]] = {
    "bwa-mem2": ["bwa-mem2"],
    "bowtie2": ["bowtie2"],
    "star": ["STAR", "star"],
    "samtools": ["samtools"],
    "deeptools": ["bamCoverage"],
    "macs3": ["macs3"],
    "hicstuff": ["hicstuff"],
    "cooler": ["cooler"],
    "bedtools": ["bedtools"],
    "ucsc-bedgraphtobigwig": ["bedGraphToBigWig"],
    "mawk": ["mawk"],
    "tree": ["tree"],
    "java-jdk": ["java"],
}


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


def _deps_available() -> bool:
    """Return True if all CLI tools from env/tinymapper.yaml are on PATH."""
    import re

    yaml_path = REPO_ROOT / "env" / "tinymapper.yaml"
    if not yaml_path.exists():
        return False
    # Parse dependency names (strip version specifiers and leading dashes)
    with yaml_path.open() as fh:
        lines = fh.readlines()
    in_deps = False
    declared: list[str] = []
    for line in lines:
        stripped = line.strip()
        if stripped == "dependencies:":
            in_deps = True
            continue
        if in_deps:
            if stripped.startswith("- pip:"):
                break  # skip pip sub-section
            if stripped.startswith("- ") and not stripped.startswith("- pip"):
                pkg = re.split(r"[>=<!\s]", stripped[2:])[0].lower()
                declared.append(pkg)
    # Check executables
    for pkg in declared:
        bins = _CONDA_TO_BINS.get(pkg)
        if bins is None:
            continue  # no executable to check (e.g. pip)
        if not any(shutil.which(b) for b in bins):
            return False
    # Also require tinymapper itself
    return shutil.which("tinymapper") is not None


_USE_ENV = _env_available()
_USE_DIRECT = not _USE_ENV and _deps_available()

requires_env = pytest.mark.skipif(
    not (_USE_ENV or _USE_DIRECT),
    reason=(
        f"micromamba env '{ENV_NAME}' not available and "
        "not all dependencies from env/tinymapper.yaml are installed"
    ),
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _run(args: list[str], output_dir: Path) -> subprocess.CompletedProcess[str]:
    """Run tinymapper with *args*, via the micromamba env or directly."""
    output_dir.mkdir(parents=True, exist_ok=True)
    if _USE_ENV:
        cmd = (
            ["micromamba", "run", "-n", ENV_NAME, "tinymapper"]
            + args
            + ["--output", str(output_dir)]
        )
    else:
        cmd = ["tinymapper"] + args + ["--output", str(output_dir)]
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


def _blacklist(name: str) -> str:
    return str(TESTS_DIR / name)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@requires_env
def test_chip_sample_only():
    result = _run(
        ["-m", "ChIP", "-s", _sample("testChIP"), "-g", _genome("R64-1-1/R64-1-1"), "-k"],
        TESTS_RESULTS_DIR / "chip_sample_only",
    )
    assert result.returncode == 0, "ChIP (sample only) failed"


@requires_env
def test_chip_with_input():
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
            "--blacklist",
            _blacklist("blacklist.bed"),
            "-k",
        ],
        TESTS_RESULTS_DIR / "chip_with_input",
    )
    assert result.returncode == 0, "ChIP (sample + input) failed"


@requires_env
def test_chip_with_input_and_calibration():
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
        TESTS_RESULTS_DIR / "chip_with_input_and_calibration",
    )
    assert result.returncode == 0, "ChIP (sample + input + calibration) failed"


@requires_env
def test_rna():
    result = _run(
        ["-m", "RNA", "-s", _sample("testRNA"), "-g", _genome("R64-1-1/R64-1-1"), "-k"],
        TESTS_RESULTS_DIR / "rna",
    )
    assert result.returncode == 0, "RNA failed"


@requires_env
def test_mnase():
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
            "-k",
        ],
        TESTS_RESULTS_DIR / "mnase",
    )
    assert result.returncode == 0, "MNase failed"


@requires_env
def test_atac_paired():
    result = _run(
        [
            "-m",
            "ATAC",
            "-s",
            _sample("testATAC"),
            "-g",
            _genome("R64-1-1/R64-1-1"),
            "--blacklist",
            _blacklist("blacklist.bed"),
            "-k",
        ],
        TESTS_RESULTS_DIR / "atac_paired",
    )
    assert result.returncode == 0, "ATAC (paired-end) failed"


@requires_env
def test_atac_single_end():
    result = _run(
        ["-m", "ATAC", "-s", _sample("testATAC_se"), "-g", _genome("R64-1-1/R64-1-1"), "-k"],
        TESTS_RESULTS_DIR / "atac_single_end",
    )
    assert result.returncode == 1, "ATAC (single-end w/o extend) failed"


@requires_env
def test_atac_single_end_with_extend():
    result = _run(
        [
            "-m",
            "ATAC",
            "-s",
            _sample("testATAC_se"),
            "-g",
            _genome("R64-1-1/R64-1-1"),
            "--extend-reads",
            "500",
            "-k",
        ],
        TESTS_RESULTS_DIR / "atac_single_end",
    )
    assert result.returncode == 0, "ATAC (single-end w/ extend) failed"


@requires_env
def test_hic():
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
            "-k",
        ],
        TESTS_RESULTS_DIR / "hic",
    )
    assert result.returncode == 0, "HiC failed"


@requires_env
def test_shotgun():
    result = _run(
        ["-m", "shotgun", "-s", _sample("testShotgun"), "-g", _genome("R64-1-1/R64-1-1"), "-k"],
        TESTS_RESULTS_DIR / "shotgun",
    )
    assert result.returncode == 0, "shotgun failed"


@requires_env
def test_shotgun_as_se_noextend():
    result = _run(
        [
            "-m",
            "shotgun",
            "-s",
            _sample("testShotgun"),
            "-g",
            _genome("R64-1-1/R64-1-1"),
            "--as-single-end",
            "-k",
        ],
        TESTS_RESULTS_DIR / "shotgun",
    )
    assert result.returncode == 1, "shotgun (as se, no extend) failed"


@requires_env
def test_shotgun_as_se():
    result = _run(
        [
            "-m",
            "shotgun",
            "-s",
            _sample("testShotgun"),
            "-g",
            _genome("R64-1-1/R64-1-1"),
            "--as-single-end",
            "--extend-reads",
            "2000",
            "-k",
        ],
        TESTS_RESULTS_DIR / "shotgun",
    )
    assert result.returncode == 0, "shotgun (as se) failed"


@requires_env
def test_shotgun_se():
    result = _run(
        [
            "-m",
            "shotgun",
            "-s",
            _sample("testShotgun_se"),
            "-g",
            _genome("R64-1-1/R64-1-1"),
            "--blacklist",
            _blacklist("blacklist.bed"),
            "--extend-reads",
            "200",
            "-k",
        ],
        TESTS_RESULTS_DIR / "shotgun",
    )
    assert result.returncode == 0, "shotgun failed"
