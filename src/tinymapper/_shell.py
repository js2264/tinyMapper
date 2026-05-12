"""Thin subprocess wrapper for external CLI tool invocations.

All mode modules use ``run_cmd`` so that:
- every command is logged before execution
- dry-run mode is honoured globally
- errors propagate as exceptions with a clear message
"""

from __future__ import annotations

import logging
import subprocess
import sys
from pathlib import Path

logger = logging.getLogger(__name__)


_err_file: Path | None = None


def configure_err_file(path: Path | None) -> None:
    """Set the global err file for all subsequent ``run_cmd`` calls."""
    global _err_file
    _err_file = path


def run_cmd(
    cmd: str,
    log_file: Path | None = None,
    dry_run: bool = False,
    err_file: Path | None = None,
) -> None:
    """Execute *cmd* as a shell string, logging it first.

    Parameters
    ----------
    cmd:        Shell command string (passed to ``subprocess.run`` with
                ``shell=True``).
    log_file:   Optional path to append the ``[EXEC]`` record to.
    dry_run:    When ``True`` the command is logged but not executed.
    err_file:   Optional path to append stderr output to.  Falls back to
                the module-level err file set via ``configure_err_file``.

    Raises
    ------
    subprocess.CalledProcessError
        Propagated on non-zero exit when not in dry-run mode.
    """
    logger.info("[EXEC] %s", cmd)
    if log_file is not None:
        with open(log_file, "a") as fh:
            fh.write(f"[EXEC] {cmd}\n")
    if dry_run:
        return
    _ef = err_file if err_file is not None else _err_file
    result = subprocess.run(cmd, shell=True, check=False, capture_output=True, text=True)
    if result.stdout:
        print(result.stdout, end="")
        if log_file is not None:
            with open(log_file, "a") as fh:
                fh.write(result.stdout)
    if result.stderr:
        print(result.stderr, end="", file=sys.stderr)
        if _ef is not None:
            with open(_ef, "a") as fh:
                fh.write(result.stderr)
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd)
