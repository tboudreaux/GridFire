gridfire_extern_sources = files(
    'lib/gridfire_context.cpp',
    'lib/gridfire_extern.cpp',
)

gridfire_extern_dependencies = [
    gridfire_dep
]

libgridfire_extern = library('gridfire_extern',
    gridfire_extern_sources,
    include_directories: include_directories('include'),
    dependencies: gridfire_extern_dependencies,
    install : true
)

gridfire_extern_dep = declare_dependency(
    include_directories: include_directories('include'),
    link_with: libgridfire_extern,
    sources: gridfire_extern_sources,
    dependencies: gridfire_extern_dependencies,
)

install_subdir('include/gridfire', install_dir: get_option('includedir'))

if get_option('build-fortran')
    message('Configuring Fortran bindings...')
    subdir('fortran')
endif
