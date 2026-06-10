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
    'lib/engine/diagnostics/dynamic_engine_diagnostics.cpp',
    'lib/engine/types/jacobian.cpp',
    'lib/reaction/reaction.cpp',
    'lib/reaction/reaclib.cpp',
    'lib/reaction/weak/weak.cpp',
    'lib/reaction/weak/weak_interpolator.cpp',
    'lib/io/network_file.cpp',
    'lib/io/generative/python.cpp',
    'lib/io/generative/mesa.cpp',
    'lib/solver/strategies/PointSolver.cpp',
    'lib/solver/strategies/GridSolver.cpp',
    'lib/solver/strategies/triggers/engine_partitioning_trigger.cpp',
    'lib/screening/screening_types.cpp',
    'lib/screening/screening_weak.cpp',
    'lib/screening/screening_bare.cpp',
    'lib/partition/partition_rauscher_thielemann.cpp',
    'lib/partition/partition_ground.cpp',
    'lib/partition/composite/partition_composite.cpp',
    'lib/policy/chains.cpp',
    'lib/policy/policy_logical.cpp',
    'lib/policy/stellar_policy.cpp',
    'lib/utils/logging.cpp',
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
    openmp_dep = dependency('openmp', required: true)
    gridfire_build_dependencies += [openmp_dep]
endif

gridfire_link_whole = [libcvode_static, libkinsol_static, libcppad_static]

gridfire_link_args = cpp.get_supported_link_arguments(
    '-Wl,--exclude-libs,libcvode-static.a:libkinsol_static.a'
)

if get_option('build_python')
    gridfire_link_whole += [libcomposition, libconst, liblogging]
    if get_option('plugin_support')
        gridfire_link_whole += [libplugin]
    endif
    libgridfire = shared_library('gridfire',
                                 gridfire_sources,
                                 include_directories: include_directories('include'),
                                 dependencies: gridfire_build_dependencies,
                                 cpp_args: gridfire_args,
                                 link_whole: gridfire_link_whole,
                                 link_args: gridfire_link_args,
                                 install: true,
                                 install_dir: gridfire_libdir)
else
    libgridfire = library('gridfire',
                          gridfire_sources,
                          include_directories: include_directories('include'),
                          dependencies: gridfire_build_dependencies,
                          link_whole: gridfire_link_whole,
                          link_args: gridfire_link_args,
                          cpp_args: gridfire_args,
                          install : true)
endif

if get_option('build_python')
    gridfire_iface_deps = []
    foreach d : gridfire_build_dependencies
        gridfire_iface_deps += d.partial_dependency(compile_args: true, includes: true)
    endforeach
    gridfire_dep = declare_dependency(
        include_directories: include_directories('include'),
        link_with: libgridfire,
        dependencies: gridfire_iface_deps,
        compile_args: gridfire_args,
    )
else
    gridfire_dep = declare_dependency(
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

