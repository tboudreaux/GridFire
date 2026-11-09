subdir('include/gridfire/utils')

cpp = meson.get_compiler('cpp')

gridfire_sources = files(
    'lib/engine/engine_graph.cpp',
    'lib/engine/views/engine_adaptive.cpp',
    'lib/engine/views/engine_defined.cpp',
    'lib/engine/views/engine_multiscale.cpp',
    'lib/engine/views/engine_priming.cpp',
    'lib/engine/procedures/priming.cpp',
    'lib/engine/procedures/construction.cpp',
    'lib/engine/procedures/energy_channels.cpp',
    'lib/engine/diagnostics/dynamic_engine_diagnostics.cpp',
    'lib/engine/types/jacobian.cpp',
    'lib/engine/scratchpads/engine_multiscale_scratchpad.cpp',
    'lib/reaction/reaction.cpp',
    'lib/reaction/reaclib.cpp',
    'lib/reaction/weak/weak.cpp',
    'lib/reaction/weak/weak_interpolator.cpp',
    'lib/reaction/weak/neutrino_energies.cpp',
    'lib/io/network_file.cpp',
    'lib/io/generative/python.cpp',
    'lib/io/generative/mesa.cpp',
    'lib/solver/strategies/PointSolver.cpp',
    'lib/solver/strategies/GridSolver.cpp',
    'lib/solver/strategies/rates.cpp',
    'lib/solver/strategies/triggers/engine_partitioning_trigger.cpp',
    'lib/screening/screening_types.cpp',
    'lib/screening/screening_weak.cpp',
    'lib/screening/screening_bare.cpp',
    'lib/screening/screening_blended.cpp',
    'lib/partition/partition_rauscher_thielemann.cpp',
    'lib/partition/partition_ground.cpp',
    'lib/partition/composite/partition_composite.cpp',
    'lib/policy/chains.cpp',
    'lib/policy/policy_logical.cpp',
    'lib/policy/stellar_policy.cpp',
    'lib/utils/logging.cpp',
    'lib/utils/provenance.cpp',
    'lib/utils/gf_omp.cpp',
)


gridfire_build_dependencies = [
    const_dep,
    config_dep,
    composition_dep,
    cppad_dep,
    log_dep,
    xxhash_dep,
    eigen_dep,
    sundials_dep,
    json_dep,
]

if get_option('use_mimalloc')
    gridfire_build_dependencies += [mimalloc_dep]
endif

if get_option('plugin_support')
    gridfire_build_dependencies += [plugin_dep]
endif

if get_option('openmp_support')
    libomp_prefix = get_option('libomp_prefix')
    if libomp_prefix != ''
        omp_lib = cpp.find_library('omp', dirs: [libomp_prefix + '/lib'])
        openmp_dep = declare_dependency(
            compile_args: ['-Xpreprocessor', '-fopenmp', '-I' + libomp_prefix + '/include'],
            dependencies: omp_lib,
        )
    else
        openmp_dep = dependency('openmp', required: true)
    endif
    gridfire_build_dependencies += [openmp_dep]
endif

gridfire_link_whole = [libcvode_static, libkinsol_static, libcppad_static]

gridfire_link_args = cpp.get_supported_link_arguments(
    '-Wl,--exclude-libs,libcvode-static.a:libkinsol_static.a'
)

# Installed libraries find their dependencies (the fourdst libraries, mimalloc, each other) in their own directory
# (normal installs; the wheel sets its own run-time paths below).
# A program's RUNPATH covers only its direct dependencies on Linux, so without this a program linked against an
# install outside the system library path fails to load libgridfire's dependencies. The fourdst libraries carry
# the same $ORIGIN. macOS needs nothing: dependents are found through the program's @rpath.
if host_machine.system() == 'darwin'
    gridfire_install_rpath = ''
else
    gridfire_install_rpath = '$ORIGIN'
endif

if get_option('build_python')
    if get_option('plugin_support')
        error('plugin_support is not available when build_python=true: '
              + 'the fourdst wheel does not ship libplugin, and bundling it '
              + 'would break cross-package type compatibility.')
    endif

    # ---- rpaths ---------------------------------------------------------
    # Wheel layout (both packages share one site-packages):
    #   <site-packages>/gridfire/_gridfire.<tag>.so
    #   <site-packages>/gridfire/lib/libgridfire.{so,dylib}
    #   <site-packages>/fourdst/lib/lib{composition,logging,const}.*
    #   <site-packages>/fourdst/lib/vendor/libreflect_cpp.*
    #
    # libgridfire (in gridfire/lib/) needs:  itself-dir, ../../fourdst/lib{,/vendor}
    # _gridfire   (in gridfire/)     needs:  ./lib, ../fourdst/lib{,/vendor}
    #
    # Platform split:
    #   * macOS / Mach-O: LC_RPATH cannot hold colon-separated lists,
    #     meson carries a user-supplied build_rpath as a single string, and
    #     install_rpath is only applied by `meson install` (which
    #     meson-python never runs). So on darwin we emit one raw
    #     `-Wl,-rpath,<path>` link arg per path, and deliberately do NOT
    #     set build_rpath/install_rpath (a real `meson install` would then
    #     add install_rpath on top of the link args, producing duplicate
    #     LC_RPATH entries, which macOS 26+ dyld treats as a hard load
    #     error).
    #   * Linux / ELF: colon-joined DT_RUNPATH is the native format, and
    #     install_rpath is baked into the binary at link time, so the
    #     kwargs work for both meson-python wheels and `meson install`.
    #
    # NOTE: compose paths with '+', never the '/' join operator —
    # 'x' / '/y' silently discards 'x' in Meson.
    #
    # gridfire_ext_rpath_args / gridfire_ext_rpath are consumed later by
    # build-python/meson.build for the _gridfire extension module (this
    # subdir is processed before build-python in the root meson.build).
    if host_machine.system() == 'darwin'
        gridfire_lib_rpath_args = [
            '-Wl,-rpath,@loader_path',
            '-Wl,-rpath,@loader_path/../../fourdst/lib',
            '-Wl,-rpath,@loader_path/../../fourdst/lib/vendor',
        ]
        gridfire_lib_rpath = ''

        gridfire_ext_rpath_args = [
            '-Wl,-rpath,@loader_path/lib',
            '-Wl,-rpath,@loader_path/../fourdst/lib',
            '-Wl,-rpath,@loader_path/../fourdst/lib/vendor',
        ]
        gridfire_ext_rpath = ''
    else
        gridfire_lib_rpath_args = []
        gridfire_lib_rpath = '$ORIGIN:' + '$ORIGIN/../../fourdst/lib:' + '$ORIGIN/../../fourdst/lib/vendor'

        gridfire_ext_rpath_args = []
        gridfire_ext_rpath = '$ORIGIN/lib:' + '$ORIGIN/../fourdst/lib:' + '$ORIGIN/../fourdst/lib/vendor'
    endif

    libgridfire = shared_library('gridfire',
        gridfire_sources,
        include_directories: include_directories('include'),
        dependencies: gridfire_build_dependencies,
        cpp_args: gridfire_args,
        link_whole: gridfire_link_whole,
        link_args: gridfire_link_args + gridfire_lib_rpath_args,
        install: true,
        install_dir: gridfire_libdir,
        build_rpath: gridfire_lib_rpath,
        install_rpath: gridfire_lib_rpath,
    )
else
    libgridfire = library('gridfire',
        gridfire_sources,
        include_directories: include_directories('include'),
        dependencies: gridfire_build_dependencies,
        link_whole: gridfire_link_whole,
        link_args: gridfire_link_args,
        cpp_args: gridfire_args,
        install_rpath: gridfire_install_rpath,
        install : true
    )
endif

# Plain version (the project version carries a leading 'v') so parents can require e.g. '>=1.2.0'
gridfire_dep_version = meson.project_version().strip('v')

if get_option('build_python')
    gridfire_iface_dep = declare_dependency(
        dependencies: gridfire_build_dependencies,
    ).partial_dependency(compile_args: true, includes: true)

    # Consumers of libgridfire (the extension, the C API, the tools) take only its compile flags from the static
    # dependencies (SUNDIALS and CppAD live inside libgridfire and must not be linked twice), but they call the
    # shared fourdst libraries and the OpenMP runtime directly, so those are linked in full.
    gridfire_shared_link_deps = [composition_dep, log_dep, const_dep, config_dep]
    if get_option('openmp_support')
        gridfire_shared_link_deps += [openmp_dep]
    endif

    gridfire_dep = declare_dependency(
        version: gridfire_dep_version,
        include_directories: include_directories('include'),
        link_with: libgridfire,
        dependencies: [gridfire_iface_dep] + gridfire_shared_link_deps,
        compile_args: gridfire_args,
    )
else
    gridfire_dep = declare_dependency(
        version: gridfire_dep_version,
        include_directories: include_directories('include'),
        link_with: libgridfire,
        dependencies: gridfire_build_dependencies,
        compile_args: gridfire_args,
    )
endif

meson.override_dependency('gridfire', gridfire_dep)

install_subdir('include/gridfire',
    install_dir: gridfire_includedir,
    exclude_files: ['utils/config.h.in'],
)


if not get_option('build_c_api') and get_option('build_fortran')
    error('Cannot build fortran without C API. Set -Dbuild-c-api=true and -Dbuild-fortran=true')
endif
if get_option('build_c_api')
    message('Configuring C API...')
    subdir('extern')
endif
