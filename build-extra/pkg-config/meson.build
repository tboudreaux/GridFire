gridfire_pc_version = meson.project_version().strip('v')

if get_option('build_python')
    # Relocatable .pc files for the wheel (<site-packages>/gridfire/lib/pkgconfig); `gf-config --pkg-config-path`
    # prints their directory for PKG_CONFIG_PATH.
    message('Generating wheel-relocatable pkg-config files for GridFire...')
    wheel_pc_conf = configuration_data()
    wheel_pc_conf.set('VERSION', gridfire_pc_version)
    wheel_pc_files = [['gridfire-wheel.pc.in', 'gridfire.pc']]
    if get_option('build_c_api')
        wheel_pc_files += [
            ['gridfire_extern-wheel.pc.in', 'gridfire_extern.pc'],
            ['gridfire_fortran-wheel.pc.in', 'gridfire_fortran.pc'],
        ]
    endif
    foreach f : wheel_pc_files
        configure_file(
            input: f[0],
            output: f[1],
            configuration: wheel_pc_conf,
            install: true,
            install_dir: gridfire_pcdir,
        )
    endforeach
elif get_option('pkg_config')
    message('Generating pkg-config files for GridFire...')
    pkg = import('pkgconfig')
    pkg.generate(
        name: 'gridfire',
        description: 'GridFire nuclear reaction network solver',
        version: gridfire_pc_version,
        filebase: 'gridfire',
        libraries: [libgridfire, '-Wl,-rpath,${libdir}'],
        requires: ['fourdst_composition', 'fourdst_config', 'fourdst_constants', 'fourdst_logging'],
        extra_cflags: ['-I${includedir}' / 'gridfire' / 'vendor'],
    )
    if get_option('build_c_api')
        # The C API is plain C: consumers need only its header directory and the library
        pkg.generate(
            name: 'gridfire_extern',
            description: 'GridFire nuclear reaction network solver (C API)',
            version: gridfire_pc_version,
            filebase: 'gridfire_extern',
            libraries: [libgridfire_extern, '-Wl,-rpath,${libdir}'],
            libraries_private: [libgridfire],
        )
        # The Fortran module is installed as source (variable fortran_source); compile it with your own compiler
        pkg.generate(
            name: 'gridfire_fortran',
            description: 'GridFire nuclear reaction network solver (Fortran interface over the C API)',
            version: gridfire_pc_version,
            filebase: 'gridfire_fortran',
            libraries: [libgridfire_extern, '-Wl,-rpath,${libdir}'],
            libraries_private: [libgridfire],
            variables: ['fortran_source=${prefix}' / gridfire_fortran_srcdir / 'gridfire_mod.f90'],
        )
    endif
endif
