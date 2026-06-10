if get_option('build_python')
    message('Generating wheel-relocatable pkg-config file for GridFire...')
    wheel_pc_conf = configuration_data()
    wheel_pc_conf.set('VERSION', meson.project_version())
    configure_file(
        input: 'gridfire-wheel.pc.in',
        output: 'gridfire.pc',
        configuration: wheel_pc_conf,
        install: true,
        install_dir: gridfire_pcdir,
    )
elif get_option('pkg_config')
    message('Generating pkg-config file for GridFire...')
    pkg = import('pkgconfig')
    pkg.generate(
        libgridfire,
        name: 'gridfire',
        description: 'GridFire nuclear reaction network solver',
        filebase: 'gridfire',
        subdirs: ['gridfire'],
        requires: [
            'fourdst_composition',
            'fourdst_config',
            'fourdst_constants',
            'fourdst_logging',
        ],
        extra_cflags: ['-I${includedir}' / 'gridfire' / 'vendor'],
    )
endif
