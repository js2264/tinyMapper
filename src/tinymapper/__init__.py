"""tinyMapper — map and process ChIP / RNA / ATAC / MNase / HiC / shotgun reads."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("tinymapper")
except PackageNotFoundError:
    __version__ = "unknown"
