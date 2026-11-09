eigen_sp = subproject('eigen')
eigen_dep = eigen_sp.get_variable('eigen_dep').as_system()

eigen_src_dir = ''
foreach root : gridfire_subproject_roots
    if eigen_src_dir == '' and fs.is_dir(root[0] / 'eigen-3.4.0' / 'Eigen')
        eigen_src_dir = root[0] / 'eigen-3.4.0' / 'Eigen'
    endif
endforeach
if eigen_src_dir == ''
    error('Cannot find the Eigen subproject sources to vendor their headers.')
endif

custom_target(
    'vendor_eigen_headers',
    command: copytree_cmd + [
        eigen_src_dir,
        '@OUTPUT@',
    ],
    output: 'Eigen',
    install: true,
    install_dir: gridfire_vendor_includedir,
)
