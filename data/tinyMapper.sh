#!/usr/bin/env bash
# tinyMapper.sh — backward-compatible wrapper for the tinymapper Python CLI.
#
# This script is a thin shim: all arguments are forwarded verbatim to the
# `tinymapper` command installed in the same environment.  The Python package
# accepts exactly the same flags as this script used to.
#
# Kept for compatibility with autotinymapper and any existing Slurm scripts
# that invoke `tinyMapper.sh`.  See:
#   src/tinymapper/cli.py   — Python CLI implementation
#   README.md               — updated documentation
exec tinymapper "$@"
