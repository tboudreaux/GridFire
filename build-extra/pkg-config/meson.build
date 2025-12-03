if get_option('pkg_config')
    message('Generating pkg-config file for GridFire...')
    pkg = import('pkgconfig')
    pkg.generate(
        name: 'gridfire',
        description: 'GridFire nuclear reaction network solver',
        version: meson.project_version(),
        libraries: [
            libgridfire,
            libcomposition,
            libconfig,
            libconst,
            liblogging
        ],
        subdirs: ['gridfire'],
        filebase: 'gridfire',
        install_dir: join_paths(get_option('libdir'), 'pkgconfig')
    )
endif
