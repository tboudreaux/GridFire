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


# ---------------------------------------------------------------------------------------------------------------------
# gf-config: compiler and linker flags for programs that use GridFire from C++, C or Fortran
# ---------------------------------------------------------------------------------------------------------------------

_COMPONENTS = ("cpp", "c", "fortran")


def gf_get_fortran_source() -> str:
    """Path of gridfire_mod.f90, the Fortran interface; compile it with the same compiler as your program."""
    return os.fspath(_PACKAGE_DIR / "share" / "fortran" / "gridfire_mod.f90")


def gf_get_cmake_prefix() -> str:
    """Directory to add to CMAKE_PREFIX_PATH so that find_package(GridFire) finds this installation."""
    return os.fspath(_PACKAGE_DIR)


def _require_c_api():
    lib_dir = _PACKAGE_DIR / "lib"
    if not any(lib_dir.glob("libgridfire_extern.*")):
        raise RuntimeError(
            f"This gridfire installation ({_PACKAGE_DIR}) has no C API library (libgridfire_extern). "
            "Install a gridfire wheel built with the C API."
        )


def gf_cflags(component: str = "cpp") -> List[str]:
    """Compile flags for a program using GridFire through the C++ library, the C API, or the Fortran module."""
    if component == "cpp":
        return gf_get_include_flags() + fst.get_include_flags() + gf_get_extra_flags()
    if component in ("c", "fortran"):
        _require_c_api()
        return ["-I" + os.fspath(_PACKAGE_DIR / "include")]
    raise ValueError(f"Unknown component {component!r}; expected one of {', '.join(_COMPONENTS)}")


def gf_libs(component: str = "cpp") -> List[str]:
    """Link flags (library paths, libraries and run-time search paths) for the given component."""
    lib_dir = os.fspath(_PACKAGE_DIR / "lib")
    if component == "cpp":
        return ["-L" + lib_dir, "-lgridfire", "-Wl,-rpath," + lib_dir] + fst.get_lib_flags()
    if component in ("c", "fortran"):
        _require_c_api()
        # The C API library finds libgridfire and the fourdst libraries through its own run-time paths; the extra
        # -rpath entries let linkers that resolve indirect dependencies at link time (GNU ld) find them too.
        flags = ["-L" + lib_dir, "-lgridfire_extern", "-Wl,-rpath," + lib_dir]
        flags += ["-Wl,-rpath," + d for d in fst.get_lib_dirs()]
        return flags
    raise ValueError(f"Unknown component {component!r}; expected one of {', '.join(_COMPONENTS)}")


def gf_config_main(argv=None) -> int:
    """
    gf-config: print what a build needs to use this GridFire installation.

    Examples:
        gf-config --cflags --libs                      # C++ library
        gf-config --component c --cflags --libs        # C API (gridfire/extern/gridfire_extern.h)
        gf-config --component fortran --libs           # Fortran: link flags
        gf-config --fortran-source                     # Fortran: the module source to compile with your compiler
        PKG_CONFIG_PATH=$(gf-config --pkg-config-path) pkg-config --cflags --libs gridfire_extern
        cmake -DCMAKE_PREFIX_PATH="$(gf-config --cmake-prefix)" ...   # then find_package(GridFire)
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog="gf-config",
        description="Compiler and linker flags for programs that use GridFire.",
    )
    parser.add_argument("--component", choices=_COMPONENTS, default="cpp",
                        help="interface to build against: the C++ library (default), the C API, or the Fortran module")
    parser.add_argument("--cflags", action="store_true", help="compile flags")
    parser.add_argument("--libs", action="store_true", help="link flags")
    parser.add_argument("--fortran-source", action="store_true", help="path of gridfire_mod.f90")
    parser.add_argument("--include-dir", action="store_true", help="GridFire's header directory")
    parser.add_argument("--lib-dir", action="store_true", help="GridFire's library directory")
    parser.add_argument("--pkg-config-path", action="store_true",
                        help="directory of gridfire.pc, gridfire_extern.pc and gridfire_fortran.pc, for PKG_CONFIG_PATH")
    parser.add_argument("--cmake-prefix", action="store_true",
                        help="directory to add to CMAKE_PREFIX_PATH for find_package(GridFire)")
    parser.add_argument("--version", action="store_true", help="GridFire version")
    args = parser.parse_args(argv)

    out: List[str] = []
    try:
        if args.version:
            out.append(gf_version())
        if args.include_dir:
            out.append(os.fspath(_PACKAGE_DIR / "include"))
        if args.lib_dir:
            out.append(os.fspath(_PACKAGE_DIR / "lib"))
        if args.fortran_source:
            out.append(gf_get_fortran_source())
        if args.pkg_config_path:
            out.append(os.fspath(_PACKAGE_DIR / "lib" / "pkgconfig"))
        if args.cmake_prefix:
            out.append(gf_get_cmake_prefix())
        if args.cflags:
            out.append(" ".join(gf_cflags(args.component)))
        if args.libs:
            out.append(" ".join(gf_libs(args.component)))
    except RuntimeError as e:
        print(f"gf-config: {e}", file=sys.stderr)
        return 1
    if not out:
        parser.print_help()
        return 1
    print("\n".join(out))  # one item per line; unquoted $(gf-config --cflags --libs) joins them with spaces
    return 0
