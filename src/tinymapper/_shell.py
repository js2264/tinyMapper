"""Thin subprocess wrapper for external CLI tool invocations.

All mode modules use ``run_cmd`` so that:
- every command is logged before execution
- dry-run mode is honoured globally
- errors propagate as exceptions with a clear message
"""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)


def run_cmd(
    cmd: str,
    log_file: Path | None = None,
    dry_run: bool = False,
) -> None:
    """Execute *cmd* as a shell string, logging it first.

    Parameters
    ----------
    cmd:        Shell command string (passed to ``subprocess.run`` with
                ``shell=True``).
    log_file:   Optional path to append the ``[EXEC]`` record to.
    dry_run:    When ``True`` the command is logged but not executed.

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
    subprocess.run(cmd, shell=True, check=True)
