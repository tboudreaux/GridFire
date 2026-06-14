try:
    import fourdst as fst
except ImportError as e:
    raise ImportError(
        "gridfire requires the fourdst package (its C++ types and shared "
        "libraries come from there). pip install fourdst."
    ) from e

from ._gridfire import *

from ._gridfire import *
import sys

from ._gridfire import type, utils, engine, solver, exceptions, partition, reaction, screening, io, policy, config

sys.modules['gridfire.type'] = type
sys.modules['gridfire.utils'] = utils
sys.modules['gridfire.engine'] = engine
sys.modules['gridfire.solver'] = solver
sys.modules['gridfire.exceptions'] = exceptions
sys.modules['gridfire.partition'] = partition
sys.modules['gridfire.reaction'] = reaction
sys.modules['gridfire.screening'] = screening
sys.modules['gridfire.policy'] = policy
sys.modules['gridfire.io'] = io
sys.modules['gridfire.config'] = config

__all__ = ['type', 'utils', 'engine', 'solver', 'exceptions', 'partition', 'reaction', 'screening', 'io', 'policy', 'config']

import importlib.metadata

try:
    _meta = importlib.metadata.metadata('gridfire')
    __version__ = _meta['Version']
    __author__ = _meta['Authors']
    __license__ = _meta['License']
    __email__ = _meta['Author-email']
    __url__ = _meta['Home-page'] or _meta.get('Project-URL', '').split(',')[0].split(' ')[-1].strip()
    __description__ = _meta['Summary']
except importlib.metadata.PackageNotFoundError :
    __version__ = 'unknown - Package not installed'
    __author__ = 'Emily M. Boudreaux'
    __license__ = 'GNU General Public License v3.0'
    __email__ = 'emily.boudreaux@dartmouth.edu'
    __url__ = 'https://github.com/4D-STAR/GridFire'

def gf_metadata():
    return {
        'version': __version__,
        'author': __author__,
        'license': __license__,
        'email': __email__,
        'url': __url__,
        'description': __description__
    }

def gf_version():
    return __version__

def gf_author():
    return __author__

def gf_license():
    return __license__

def gf_email():
    return __email__

def gf_url():
    return __url__

def gf_description():
    return __description__

def gf_collaboration():
    return "4D-STAR Collaboration"

def gf_credits():
    return [
        "Emily M. Boudreaux - Lead Developer",
        "Aaron Dotter - Co-Developer",
        "4D-STAR Collaboration - Contributors"
    ]

import os
from pathlib import Path
from typing import List

_PACKAGE_DIR = Path(__file__).resolve().parent
def gf_get_include_dirs():
    return [
        os.fspath(_PACKAGE_DIR / "include"),
        os.fspath(_PACKAGE_DIR / "include" / "gridfire" / "vendor"),
    ]

def gf_get_lib_dirs():
    return [
        os.fspath(_PACKAGE_DIR / "lib"),
    ]

def gf_get_rpath_flags() -> List[str]:
    return ["-Wl,-rpath," + os.fspath(_PACKAGE_DIR / "lib")]


def gf_get_lib_flags() -> List[str]:
    flags = ["-L" + d for d in gf_get_lib_dirs()]
    flags += ["-lgridfire"]
    return flags

def gf_get_include_flags() -> List[str]:
    return ["-I" + d for d in gf_get_include_dirs()]

def gf_get_extra_flags() -> List[str]:
    return ['--std=c++23', '-fPIC']
def gf_compiler_flags(just_gridfire=False):
    flags = []
    if not just_gridfire:
        flags.extend(fst.get_compiler_flags())
    flags.extend(gf_get_rpath_flags())
    flags.extend(gf_get_lib_flags())
    flags.extend(gf_get_include_flags())
    flags.extend(gf_get_extra_flags())
    return flags

def gf_get_compiler_flags_formatted(just_gridfire=False) -> int:
    flags = gf_compiler_flags(just_gridfire)
    print(" ".join(flags))
    return 0
