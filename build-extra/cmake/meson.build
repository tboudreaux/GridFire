# GridFireConfig.cmake and GridFireConfigVersion.cmake, for CMake's find_package(GridFire). Installed next to the
# libraries (<libdir>/cmake/GridFire), both for a normal install and in the Python wheel, where
# `gf-config --cmake-prefix` prints the directory to add to CMAKE_PREFIX_PATH.
gridfire_cmake_version = meson.project_version().strip('v')
gridfire_cmakedir = gridfire_libdir / 'cmake' / 'GridFire'

cmake_conf = configuration_data()
cmake_conf.set('GF_VERSION', gridfire_cmake_version)
if get_option('build_python')
    # <site-packages>/gridfire/lib/cmake/GridFire; fourdst is the fourdst wheel in the same site-packages
    cmake_conf.set('GF_LAYOUT', 'wheel')
    cmake_conf.set('GF_PREFIX_FROM_CMAKEDIR', '../../..')
    cmake_conf.set('GF_INCLUDEDIR', 'include')
    cmake_conf.set('GF_LIBDIR', 'lib')
    cmake_conf.set('GF_FORTRAN_SOURCE', 'share/fortran/gridfire_mod.f90')
    cmake_conf.set('GF_FOURDST_PREFIX', '../fourdst')
    cmake_conf.set('GF_FOURDST_LIBDIR', 'lib')
    cmake_conf.set('GF_FOURDST_VENDOR_LIBDIR', 'lib/vendor')
    cmake_conf.set('GF_FOURDST_INCLUDEDIR', 'include')
else
    # Paths relative to the install prefix (absolute if an install directory was given outside it). fourdst is
    # installed into the same prefix when built as a subproject; otherwise CMake's usual search paths are tried too.
    libdir_opt = get_option('libdir')
    if libdir_opt.startswith('/')
        gf_prefix_from_cmakedir = get_option('prefix')
    else
        gf_prefix_from_cmakedir = '../..'
        foreach component : libdir_opt.split('/')
            gf_prefix_from_cmakedir += '/..'
        endforeach
    endif
    cmake_conf.set('GF_LAYOUT', 'install')
    cmake_conf.set('GF_PREFIX_FROM_CMAKEDIR', gf_prefix_from_cmakedir)
    cmake_conf.set('GF_INCLUDEDIR', get_option('includedir'))
    cmake_conf.set('GF_LIBDIR', libdir_opt)
    cmake_conf.set('GF_FORTRAN_SOURCE', gridfire_fortran_srcdir / 'gridfire_mod.f90')
    cmake_conf.set('GF_FOURDST_PREFIX', '.')
    cmake_conf.set('GF_FOURDST_LIBDIR', libdir_opt)
    cmake_conf.set('GF_FOURDST_VENDOR_LIBDIR', libdir_opt / 'fourdst' / 'vendor')
    cmake_conf.set('GF_FOURDST_INCLUDEDIR', get_option('includedir'))
endif

configure_file(
    input: 'GridFireConfig.cmake.in',
    output: 'GridFireConfig.cmake',
    configuration: cmake_conf,
    install: true,
    install_dir: gridfire_cmakedir,
)

# Meson's cmake.write_basic_package_version_file is not used: meson-python cannot map its install path into the wheel
configure_file(
    input: 'GridFireConfigVersion.cmake.in',
    output: 'GridFireConfigVersion.cmake',
    configuration: cmake_conf,
    install: true,
    install_dir: gridfire_cmakedir,
)
