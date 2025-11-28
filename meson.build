# ***********************************************************************
#
#   Copyright (C) 2025 -- The 4D-STAR Collaboration
#   File Author: Emily Boudreaux
#   Last Modified: June 21, 2025
#
#   GridFire is free software; you can use it and/or modify
#   it under the terms and restrictions the GNU General Library Public
#   License version 3 (GPLv3) as published by the Free Software Foundation.
#
#   GridFire is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#   See the GNU Library General Public License for more details.
#
#   You should have received a copy of the GNU Library General Public License
#   along with this software; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
# *********************************************************************** #
project('GridFire', ['c', 'cpp', 'fortran'], version: 'v0.7.1_rc2', default_options: ['cpp_std=c++23'], meson_version: '>=1.5.0')

# Add default visibility for all C++ targets
add_project_arguments('-fvisibility=default', language: 'cpp')

message('Found CXX compiler: ' + meson.get_compiler('cpp').get_id())
message('Found FORTRAN compiler: ' + meson.get_compiler('fortran').get_id())
message('C++ standard set to: ' + get_option('cpp_std'))
message('Fortran standard set to: ' + get_option('fortran_std'))


if meson.get_compiler('cpp').get_id() == 'clang'
    # We disable these because of CppAD
    message('disabling bitwise-instead-of-logical warnings for clang')
    add_project_arguments('-Wno-bitwise-instead-of-logical', language: 'cpp')
endif

if meson.get_compiler('cpp').get_id() == 'gcc'
    # We disable these because of boost notes about abi changes from C++10 -> C++17 make the build too noisey
    message('disabling psabi warnings for gcc')
    add_project_arguments('-Wno-psabi', language: 'cpp')
    
    if (meson.get_compiler('cpp').version().version_compare('<14.0'))
      error('g++ version must be at least 14.0, found ' + meson.get_compiler('cpp').version())
    endif
endif

build_fortran = get_option('build-fortran')
if (build_fortran)
  add_languages('fortran', native: true)
  message('Building fortran module (gridfire_mod.mod)')
  fc = meson.get_compiler('fortran')
  if not get_option('unsafe-fortran')
    if fc.get_id() != 'gcc'
      error('The only supported fortran compiler for GridFire is gfortran (version >= 14.0), found ' + fc  + '. GridFire has not been tested with any other compilers. You can disable this check with the -Dunsafe-fortran=true flag to try other compilers')
    endif
  endif
endif



if (meson.get_compiler('fortran').version().version_compare('<14.0'))
  error('gfortran version must be at least 14.0, found ' + meson.get_compiler('fortran').version())
endif

# For Eigen
add_project_arguments('-Wno-deprecated-declarations', language: 'cpp')

llevel = get_option('log-level')

logbase='QUILL_COMPILE_ACTIVE_LOG_LEVEL_'

if (llevel == 'traceL3')
    message('Setting log level to TRACE_L3')
    log_argument = logbase + 'TRACE_L3'
elif (llevel == 'traceL2')
    message('Setting log level to TRACE_L2')
    log_argument = logbase + 'TRACE_L2'
elif (llevel == 'traceL1')
    message('Setting log level to TRACE_L1')
    log_argument = logbase + 'TRACE_L1'
elif (llevel == 'debug')
    message('Setting log level to DEBUG')
    log_argument = logbase + 'DEBUG'
elif (llevel == 'info')
    message('Setting log level to INFO')
    log_argument = logbase + 'INFO'
elif (llevel == 'warning')
    message('Setting log level to WARNING')
    log_argument = logbase + 'WARNING'
elif (llevel == 'error')
    message('Setting log level to ERROR')
    log_argument = logbase + 'ERROR'
elif (llevel == 'critical')
    message('Setting log level to CRITICAL')
    log_argument = logbase + 'CRITICAL'
endif

log_argument = '-DQUILL_COMPILE_ACTIVE_LOG_LEVEL=' + log_argument
add_project_arguments(log_argument, language: 'cpp')

cpp = meson.get_compiler('cpp')

subdir('build-config')

subdir('src')

if get_option('build-python')
    message('Configuring Python bindings...')
    subdir('build-python')
else
    message('Skipping Python bindings...')
endif

if get_option('build-tests')
    message('Setting up tests for GridFire...')
    subdir('tests')
else
    message('Skipping tests for GridFire...')
endif


if get_option('pkg-config')
    message('Generating pkg-config file for GridFire...')
    pkg = import('pkgconfig')
    pkg.generate(
        name: 'gridfire',
        description: 'GridFire nuclear reaction network solver',
        version: meson.project_version(),
        libraries: [
            libgridfire,
            libcomposition,
            libconfig,
            libconst,
            liblogging
        ],
        subdirs: ['gridfire'],
        filebase: 'gridfire',
        install_dir: join_paths(get_option('libdir'), 'pkgconfig')
    )
endif


