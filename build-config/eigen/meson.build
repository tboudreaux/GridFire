eigen_sp = subproject('eigen')
eigen_dep = eigen_sp.get_variable('eigen_dep').as_system()

custom_target(
    'vendor_eigen_headers',
    command: copytree_cmd + [
        meson.project_source_root() / 'subprojects' / 'eigen-3.4.0' / 'Eigen',
        '@OUTPUT@',
    ],
    output: 'Eigen',
    install: true,
    install_dir: gridfire_vendor_includedir,
)
