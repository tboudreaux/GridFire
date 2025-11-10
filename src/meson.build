# Define the library
gridfire_sources = files(
    'lib/network.cpp',
    'lib/engine/engine_graph.cpp',
    'lib/engine/views/engine_adaptive.cpp',
    'lib/engine/views/engine_defined.cpp',
    'lib/engine/views/engine_multiscale.cpp',
    'lib/engine/views/engine_priming.cpp',
    'lib/engine/procedures/priming.cpp',
    'lib/engine/procedures/construction.cpp',
    'lib/engine/diagnostics/dynamic_engine_diagnostics.cpp',
    'lib/reaction/reaction.cpp',
    'lib/reaction/reaclib.cpp',
    'lib/reaction/weak/weak.cpp',
    'lib/reaction/weak/weak_interpolator.cpp',
    'lib/io/network_file.cpp',
    'lib/io/generative/python.cpp',
    'lib/solver/strategies/CVODE_solver_strategy.cpp',
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
    boost_dep,
    const_dep,
    config_dep,
    composition_dep,
    cppad_dep,
    log_dep,
    xxhash_dep,
    eigen_dep,
    plugin_dep,
    sundials_dep,
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

install_subdir('include/gridfire', install_dir: get_option('includedir'))

if get_option('build-python')
    message('Configuring Python bindings...')
    subdir('python')
else
    message('Skipping Python bindings...')
endif