gridfire_fortran_lib = library('gridfire_fortran',
    'gridfire_mod.f90',
    link_with: libgridfire_extern,
    install: true,
    install_dir: get_option('libdir')
)

gridfire_fortran_dep = declare_dependency(
    link_with: [gridfire_fortran_lib, libgridfire_extern],
    include_directories: include_directories('.')
)