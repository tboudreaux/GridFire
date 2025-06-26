# Define the library
network_sources = files(
    'lib/network.cpp',
    'lib/engine/engine_approx8.cpp',
    'lib/engine/engine_graph.cpp',
    'lib/reaction/reaction.cpp',
    'lib/solver/solver.cpp',
)


dependencies = [
    boost_dep,
    const_dep,
    config_dep,
    composition_dep,
    reaclib_reactions_dep,
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
    'include/gridfire/engine/engine_approx8.h',
    'include/gridfire/engine/engine_graph.h',
    'include/gridfire/reaction/reaction.h',
    'include/gridfire/solver/solver.h',
)
install_headers(network_headers, subdir : 'gridfire')
