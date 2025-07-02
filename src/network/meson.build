# Define the library
network_sources = files(
    'lib/network.cpp',
    'lib/engine/engine_approx8.cpp',
    'lib/engine/engine_graph.cpp',
    'lib/engine/views/engine_adaptive.cpp',
    'lib/engine/views/engine_defined.cpp',
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


dependencies = [
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
libnetwork = library('network',
    network_sources, 
    include_directories: include_directories('include'),
    dependencies: dependencies,
    install : true)

network_dep = declare_dependency(
    include_directories: include_directories('include'),
    link_with: libnetwork,
    sources: network_sources,
    dependencies: dependencies,
)

# Make headers accessible
network_headers = files(
    'include/gridfire/network.h',
    'include/gridfire/engine/engine_abstract.h',
    'include/gridfire/engine/views/engine_view_abstract.h',
    'include/gridfire/engine/engine_approx8.h',
    'include/gridfire/engine/engine_graph.h',
    'include/gridfire/engine/views/engine_adaptive.h',
    'include/gridfire/engine/views/engine_defined.h',
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
install_headers(network_headers, subdir : 'gridfire')
