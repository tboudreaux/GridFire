sundials_src_include = meson.project_source_root() / 'subprojects' / 'cvode-7.5.0' / 'include'
kinsol_src_include   = meson.project_source_root() / 'subprojects' / 'kinsol-7.5.0' / 'include'

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

sundials_cmake_include = meson.global_build_root() / 'subprojects' / 'cvode-7.5.0' / '__CMake_build' / 'include' / 'sundials'

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
