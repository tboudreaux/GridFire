# Define the library
network_sources = files(
    'private/network.cpp',
    'private/approx8.cpp',
    'private/netgraph.cpp'
)

network_headers = files(
    'public/network.h',
    'public/approx8.h',
    'public/netgraph.h',
)

dependencies = [
    boost_dep,
    const_dep,
    quill_dep,
    mfem_dep,
    config_dep,
    probe_dep,
    species_weight_dep,
    composition_dep,
    reaclib_reactions_dep,
]

# Define the libnetwork library so it can be linked against by other parts of the build system
libnetwork = library('network',
    network_sources, 
    include_directories: include_directories('public'),
    dependencies: dependencies,
    install : true)

network_dep = declare_dependency(
    include_directories: include_directories('public'),
    link_with: libnetwork,
    sources: network_sources,
    dependencies: dependencies,
)

# Make headers accessible
install_headers(network_headers, subdir : '4DSSE/network')