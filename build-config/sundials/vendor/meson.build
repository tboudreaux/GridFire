# The SUNDIALS subprojects are in GridFire's own subprojects/ or the parent project's (gridfire_subproject_roots);
# both their sources and CMake build directories exist by now (configured above), so take whichever is there.
sundials_subproject_dir = ''
foreach root : gridfire_subproject_roots
    if sundials_subproject_dir == '' and fs.is_dir(root[0] / 'cvode-7.5.0' / 'include')
        sundials_subproject_dir = root[0]
        sundials_subproject_build_dir = root[1]
    endif
endforeach
if sundials_subproject_dir == ''
    error('Cannot find the SUNDIALS (cvode-7.5.0) subproject sources to vendor their headers.')
endif

sundials_src_include = sundials_subproject_dir / 'cvode-7.5.0' / 'include'
kinsol_src_include   = sundials_subproject_dir / 'kinsol-7.5.0' / 'include'

sundials_vendor_dirs = [
    'sundials',           # core (includes sundials/priv/)
    'cvode',
    'nvector',
    'sunmatrix',
    'sunlinsol',
    'sunnonlinsol',
    'sunmemory',
    'sunadaptcontroller',
]

foreach d : sundials_vendor_dirs
    custom_target(
        'vendor_sundials_' + d,
        command: copytree_cmd + [sundials_src_include / d, '@OUTPUT@'],
        output: d,
        install: true,
        install_dir: gridfire_vendor_includedir,
    )
endforeach

custom_target(
    'vendor_sundials_kinsol',
    command: copytree_cmd + [kinsol_src_include / 'kinsol', '@OUTPUT@'],
    output: 'kinsol',
    install: true,
    install_dir: gridfire_vendor_includedir,
)

sundials_cmake_include = sundials_subproject_build_dir / 'cvode-7.5.0' / '__CMake_build' / 'include' / 'sundials'
if not fs.is_file(sundials_cmake_include / 'sundials_config.h')
    error('Cannot find the configured SUNDIALS headers in ' + sundials_cmake_include)
endif

custom_target(
    'vendor_sundials_generated',
    command: copyfiles_cmd + [
        sundials_cmake_include / 'sundials_config.h',
        sundials_cmake_include / 'sundials_export.h',
        '@OUTDIR@',
    ],
    output: ['sundials_config.h', 'sundials_export.h'],
    depends: libcvode_static,
    install: true,
    install_dir: gridfire_vendor_includedir / 'sundials',
)
