gridfire_fortran_lib = library('gridfire_fortran',
    'gridfire_mod.f90',
    link_with: libgridfire_extern,
    install_rpath: gridfire_install_rpath,
    install: true,
    install_dir: gridfire_libdir
)

gridfire_fortran_dep = declare_dependency(
    version: gridfire_dep_version,
    link_with: [gridfire_fortran_lib, libgridfire_extern],
    include_directories: include_directories('.')
)
# For Meson projects that use GridFire as a subproject (with build_fortran=true): dependency('gridfire_fortran').
# The module is compiled by the parent's own Fortran compiler, so its .mod files match.
meson.override_dependency('gridfire_fortran', gridfire_fortran_dep)
