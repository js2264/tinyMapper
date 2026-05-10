"""Tests for JobSpec model validation and to_cli_args round-trip."""

from __future__ import annotations

import pytest

from tinymapper.models import JobSpec, TinyMapperMode, _generate_hash

# ------------------------------------------------------------------ #
#  _generate_hash                                                      #
# ------------------------------------------------------------------ #


def test_generate_hash_format():
    h = _generate_hash()
    assert len(h) == 6
    assert h.isalnum()
    assert h == h.upper()


def test_generate_hash_is_random():
    hashes = {_generate_hash() for _ in range(50)}
    assert len(hashes) > 1


# ------------------------------------------------------------------ #
#  TinyMapperMode                                                      #
# ------------------------------------------------------------------ #


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("chip", TinyMapperMode.chip),
        ("ChIP", TinyMapperMode.chip),
        ("CHIP", TinyMapperMode.chip),
        ("rna", TinyMapperMode.rna),
        ("RNA", TinyMapperMode.rna),
        ("atac", TinyMapperMode.atac),
        ("ATAC", TinyMapperMode.atac),
        ("mnase", TinyMapperMode.mnase),
        ("MNase", TinyMapperMode.mnase),
        ("hic", TinyMapperMode.hic),
        ("HiC", TinyMapperMode.hic),
        ("shotgun", TinyMapperMode.shotgun),
    ],
)
def test_mode_normalisation(raw, expected):
    spec = JobSpec(mode=raw, sample="/reads/JS001", genome="/genome/W303/W303")
    assert spec.mode == expected


# ------------------------------------------------------------------ #
#  JobSpec defaults                                                    #
# ------------------------------------------------------------------ #


def test_jobspec_defaults():
    spec = JobSpec(mode="chip", sample="/reads/JS001", genome="/genome/W303/W303")
    assert spec.threads == 8
    assert spec.alignment == "--maxins 1000"
    assert spec.filter_opts == "-f 0x001 -f 0x002 -F 0x004 -F 0x008 -q 10"
    assert spec.gsize == "13000000"
    assert not spec.keep_duplicates
    assert not spec.keep_intermediate
    assert not spec.do_input
    assert not spec.do_calibration
    assert spec.do_peaks  # chip → peaks


def test_jobspec_hash_auto_populated():
    spec = JobSpec(mode="chip", sample="/reads/JS001", genome="/genome/W303/W303")
    assert len(spec.hash) == 6
    assert spec.hash.isalnum()


def test_jobspec_hash_preserved_when_set():
    spec = JobSpec(mode="chip", sample="/reads/JS001", genome="/genome/W303/W303", hash="ABC123")
    assert spec.hash == "ABC123"


# ------------------------------------------------------------------ #
#  Derived properties                                                  #
# ------------------------------------------------------------------ #


def test_sample_base():
    spec = JobSpec(mode="chip", sample="/reads/JS001", genome="/genome/W303/W303")
    assert spec.sample_base == "JS001"


def test_genome_base():
    spec = JobSpec(mode="chip", sample="/reads/JS001", genome="/genome/W303/W303")
    assert spec.genome_base == "W303"


def test_input_base_empty():
    spec = JobSpec(mode="chip", sample="/reads/JS001", genome="/g/W303/W303")
    assert spec.input_base == ""


def test_input_base_set():
    spec = JobSpec(mode="chip", sample="/reads/JS001", genome="/g/W303/W303", input="/reads/JS002")
    assert spec.input_base == "JS002"


def test_do_peaks_chip():
    spec = JobSpec(mode="chip", sample="/reads/JS001", genome="/g/W303/W303")
    assert spec.do_peaks is True


def test_do_peaks_rna():
    spec = JobSpec(mode="rna", sample="/reads/JS001", genome="/g/W303/W303")
    assert spec.do_peaks is False


def test_do_peaks_atac():
    spec = JobSpec(mode="atac", sample="/reads/JS001", genome="/g/W303/W303")
    assert spec.do_peaks is True


def test_do_peaks_mnase():
    spec = JobSpec(mode="mnase", sample="/reads/JS001", genome="/g/W303/W303")
    assert spec.do_peaks is False


def test_do_peaks_hic():
    spec = JobSpec(mode="hic", sample="/reads/JS001", genome="/g/W303/W303")
    assert spec.do_peaks is False


def test_remove_duplicates_flag_default():
    spec = JobSpec(mode="chip", sample="/reads/JS001", genome="/g/W303/W303")
    assert spec.remove_duplicates_flag == "-r"


def test_remove_duplicates_flag_keep():
    spec = JobSpec(mode="chip", sample="/reads/JS001", genome="/g/W303/W303", keep_duplicates=True)
    assert spec.remove_duplicates_flag == ""


def test_ignore_duplicates_flag_default():
    spec = JobSpec(mode="chip", sample="/reads/JS001", genome="/g/W303/W303")
    assert spec.ignore_duplicates_flag == "--ignoreDuplicates"


def test_blacklist_options_empty():
    spec = JobSpec(mode="chip", sample="/reads/JS001", genome="/g/W303/W303")
    assert spec.blacklist_options == ""


def test_blacklist_options_set():
    spec = JobSpec(
        mode="chip",
        sample="/reads/JS001",
        genome="/g/W303/W303",
        blacklist="/path/to/blacklist.bed",
    )
    assert spec.blacklist_options == "--blackListFileName /path/to/blacklist.bed"


def test_mnase_sizes_parsing():
    spec = JobSpec(
        mode="mnase", sample="/reads/JS001", genome="/g/W303/W303", mnase_sizes="150,200"
    )
    assert spec.mnase_min_size == "150"
    assert spec.mnase_max_size == "200"


def test_hic_base_rez_single():
    spec = JobSpec(mode="hic", sample="/reads/JS001", genome="/g/W303/W303", binning="1000")
    assert spec.hic_base_rez == "1000"


def test_hic_base_rez_multi():
    spec = JobSpec(
        mode="hic", sample="/reads/JS001", genome="/g/W303/W303", binning="1000,2000,4000"
    )
    assert spec.hic_base_rez == "1000"


# ------------------------------------------------------------------ #
#  Validation                                                          #
# ------------------------------------------------------------------ #


def test_calibration_without_input_raises():
    with pytest.raises(ValueError, match="--calibration requires --input"):
        JobSpec(
            mode="chip",
            sample="/reads/JS001",
            genome="/g/W303/W303",
            calibration="/g/CBS138/CBS138",
        )


def test_calibration_with_input_ok():
    spec = JobSpec(
        mode="chip",
        sample="/reads/JS001",
        genome="/g/W303/W303",
        input="/reads/JS002",
        calibration="/g/CBS138/CBS138",
    )
    assert spec.do_calibration


# ------------------------------------------------------------------ #
#  to_cli_args round-trip                                              #
# ------------------------------------------------------------------ #


def test_to_cli_args_minimal():
    spec = JobSpec(mode="chip", sample="/reads/JS001", genome="/g/W303/W303", hash="TSTRUN")
    args = spec.to_cli_args()
    assert "--mode" in args
    assert "ChIP" in args
    assert "--sample" in args
    assert "/reads/JS001" in args
    assert "--genome" in args
    assert "/g/W303/W303" in args
    # Default flags should NOT be in args
    assert "--duplicates" not in args
    assert "--keepIntermediate" not in args
    assert "--blacklist" not in args


def test_to_cli_args_with_input():
    spec = JobSpec(
        mode="chip",
        sample="/reads/JS001",
        genome="/g/W303/W303",
        input="/reads/JS002",
        hash="TSTRUN",
    )
    args = spec.to_cli_args()
    assert "--input" in args
    assert "/reads/JS002" in args


def test_to_cli_args_with_flags():
    spec = JobSpec(
        mode="chip",
        sample="/reads/JS001",
        genome="/g/W303/W303",
        keep_duplicates=True,
        keep_intermediate=True,
        hash="TSTRUN",
    )
    args = spec.to_cli_args()
    assert "--duplicates" in args
    assert "--keepIntermediate" in args


def test_to_cli_args_hic_includes_hic_opts():
    spec = JobSpec(
        mode="hic",
        sample="/reads/JS001",
        genome="/g/W303/W303",
        binning="1000",
        hash="TSTRUN",
    )
    args = spec.to_cli_args()
    assert "--restriction" in args
    assert "--binning" in args
    assert "1000" in args
    assert "--hicstuff" in args


def test_to_cli_args_mnase_includes_sizes():
    spec = JobSpec(
        mode="mnase",
        sample="/reads/JS001",
        genome="/g/W303/W303",
        mnase_sizes="150,200",
        hash="TSTRUN",
    )
    args = spec.to_cli_args()
    assert "--MNaseSizes" in args
    assert "150,200" in args


def test_to_cli_args_non_default_filter():
    spec = JobSpec(
        mode="chip",
        sample="/reads/JS001",
        genome="/g/W303/W303",
        filter_opts="-F 4 -q 5",
        hash="TSTRUN",
    )
    args = spec.to_cli_args()
    assert "--filter" in args
    idx = args.index("--filter")
    assert args[idx + 1] == "-F 4 -q 5"


def test_to_cli_args_non_default_alignment():
    spec = JobSpec(
        mode="chip",
        sample="/reads/JS001",
        genome="/g/W303/W303",
        alignment="--very-sensitive",
        hash="TSTRUN",
    )
    args = spec.to_cli_args()
    assert "--alignment" in args


def test_to_cli_args_default_alignment_omitted():
    """Default alignment option should be omitted from CLI args."""
    spec = JobSpec(mode="chip", sample="/reads/JS001", genome="/g/W303/W303", hash="TSTRUN")
    args = spec.to_cli_args()
    assert "--alignment" not in args
