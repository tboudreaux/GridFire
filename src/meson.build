# Define the library
gridfire_sources = files(
    'lib/network.cpp',
    'lib/engine/engine_approx8.cpp',
    'lib/engine/engine_graph.cpp',
    'lib/engine/views/engine_adaptive.cpp',
    'lib/engine/views/engine_defined.cpp',
    'lib/engine/views/engine_multiscale.cpp',
    'lib/engine/views/engine_priming.cpp',
    'lib/engine/procedures/priming.cpp',
    'lib/engine/procedures/construction.cpp',
    'lib/reaction/reaction.cpp',
    'lib/reaction/reaclib.cpp',
    'lib/io/network_file.cpp',
    'lib/solver/solver.cpp',
    'lib/screening/screening_types.cpp',
    'lib/screening/screening_weak.cpp',
    'lib/screening/screening_bare.cpp',
    'lib/partition/partition_rauscher_thielemann.cpp',
    'lib/partition/partition_ground.cpp',
    'lib/partition/composite/partition_composite.cpp',
    'lib/utils/logging.cpp',
)


gridfire_build_dependencies = [
    boost_dep,
    const_dep,
    config_dep,
    composition_dep,
    cppad_dep,
    log_dep,
    xxhash_dep,
    eigen_dep,
]

# Define the libnetwork library so it can be linked against by other parts of the build system
libgridfire = library('gridfire',
    gridfire_sources,
    include_directories: include_directories('include'),
    dependencies: gridfire_build_dependencies,
    install : true)

gridfire_dep = declare_dependency(
    include_directories: include_directories('include'),
    link_with: libgridfire,
    sources: gridfire_sources,
    dependencies: gridfire_build_dependencies,
)

# Make headers accessible
gridfire_headers = files(
    'include/gridfire/network.h',
    'include/gridfire/engine/engine_abstract.h',
    'include/gridfire/engine/views/engine_view_abstract.h',
    'include/gridfire/engine/engine_approx8.h',
    'include/gridfire/engine/engine_graph.h',
    'include/gridfire/engine/views/engine_adaptive.h',
    'include/gridfire/engine/views/engine_defined.h',
    'include/gridfire/engine/views/engine_multiscale.h',
    'include/gridfire/engine/views/engine_priming.h',
    'include/gridfire/engine/procedures/priming.h',
    'include/gridfire/engine/procedures/construction.h',
    'include/gridfire/reaction/reaction.h',
    'include/gridfire/reaction/reaclib.h',
    'include/gridfire/io/network_file.h',
    'include/gridfire/solver/solver.h',
    'include/gridfire/screening/screening_abstract.h',
    'include/gridfire/screening/screening_bare.h',
    'include/gridfire/screening/screening_weak.h',
    'include/gridfire/screening/screening_types.h',
    'include/gridfire/partition/partition_abstract.h',
    'include/gridfire/partition/partition_rauscher_thielemann.h',
    'include/gridfire/partition/partition_ground.h',
    'include/gridfire/partition/composite/partition_composite.h',
    'include/gridfire/utils/logging.h',
)
install_headers(gridfire_headers, subdir : 'gridfire')

subdir('python')
