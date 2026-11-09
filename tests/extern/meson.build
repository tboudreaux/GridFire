if get_option('build_c_api')
  subdir('C')
endif

if get_option('build_fortran')
  subdir('fortran')
endif

if get_option('build_fortran')
    # Every bind(C) interface in gridfire_mod.f90 must match the C header (argument count, pointer vs value, kinds)
    test('extern_interface_consistency',
        find_program('python3'),
        args: [
            files('check_interface_consistency.py'),
            meson.project_source_root() / 'src/extern/include/gridfire/extern/gridfire_extern.h',
            meson.project_source_root() / 'src/extern/fortran/gridfire_mod.f90',
        ],
        suite: 'extern',
    )
endif
