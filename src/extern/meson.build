gridfire_extern_sources = files(
    'lib/gridfire_context.cpp',
    'lib/gridfire_extern.cpp',
)

gridfire_extern_dependencies = [
    gridfire_dep
]

if get_option('build_python')
    # Lives next to libgridfire in the wheel (<site-packages>/gridfire/lib), with the same rpaths. Always shared, like
    # libgridfire in the wheel (whose build uses default_library=static for its subprojects).
    gridfire_extern_link_args = gridfire_lib_rpath_args
    gridfire_extern_build_rpath = gridfire_lib_rpath
    gridfire_extern_rpath = gridfire_lib_rpath
    gridfire_extern_library_kind = 'shared_library'
else
    gridfire_extern_link_args = []
    gridfire_extern_build_rpath = ''
    gridfire_extern_rpath = gridfire_install_rpath  # $ORIGIN on Linux, as for libgridfire (src/meson.build)
    gridfire_extern_library_kind = 'library'
endif

libgridfire_extern = build_target('gridfire_extern',
    gridfire_extern_sources,
    target_type: gridfire_extern_library_kind,
    include_directories: include_directories('include'),
    dependencies: gridfire_extern_dependencies,
    link_args: gridfire_extern_link_args,
    build_rpath: gridfire_extern_build_rpath,
    install_rpath: gridfire_extern_rpath,
    install: true,
    install_dir: gridfire_libdir,
)

gridfire_extern_dep = declare_dependency(
    version: gridfire_dep_version,
    include_directories: include_directories('include'),
    link_with: libgridfire_extern,
    dependencies: gridfire_extern_dependencies,
)
# For Meson projects that use GridFire as a subproject: dependency('gridfire_extern')
meson.override_dependency('gridfire_extern', gridfire_extern_dep)

install_subdir('include/gridfire', install_dir: gridfire_includedir)

# Fortran users compile the module with their own compiler (see gf-config --fortran-source), so the source is
# installed whether or not this build compiles it
install_data('fortran/gridfire_mod.f90', install_dir: gridfire_fortran_srcdir)

if get_option('build_fortran')
    message('Configuring Fortran bindings...')
    subdir('fortran')
endif
