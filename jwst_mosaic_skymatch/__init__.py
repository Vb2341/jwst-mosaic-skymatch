# Licensed under a 3-clause BSD style license - see LICENSE.rst

from importlib.metadata import version as _version, PackageNotFoundError
try:
    __version__ = _version(__name__)
except PackageNotFoundError:
    pass


# Expose subpackage API at package level.
from .pixel_skymatch_step import PixelSkyMatchStep
