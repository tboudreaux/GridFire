if get_option('build_python')
    fourdst_inc_probe = run_command(py_installation, '-c',
                                    'import fourdst; print("\\n".join(fourdst.get_include_dirs()))',
                                    check: false)
    if fourdst_inc_probe.returncode() != 0
        error('Could not interrogate the fourdst wheel:\n' + fourdst_inc_probe.stderr()
              + '\nIs fourdst installed in the build environment?')
    endif
    fourdst_inc_dirs = fourdst_inc_probe.stdout().strip().split('\n')

    fourdst_lib_probe = run_command(py_installation, '-c',
                                    'import fourdst; print("\\n".join(fourdst.get_lib_dirs()))',
                                    check: false)
    if fourdst_lib_probe.returncode() != 0
        error('Could not interrogate the fourdst wheel:\n' + fourdst_lib_probe.stderr())
    endif
    fourdst_lib_dirs = fourdst_lib_probe.stdout().strip().split('\n')
    fourdst_inc_args = []
    foreach d : fourdst_inc_dirs
        fourdst_inc_args += ['-I' + d]
    endforeach

    cpp = meson.get_compiler('cpp')
    # Check the libraries exist, but link them by name: linking the found files directly makes Meson add the build
    # environment's absolute fourdst/lib to every binary's rpath, which is meaningless in the installed wheel (the
    # relative rpaths set in src/meson.build are what find fourdst there).
    fourdst_L_args = []
    foreach d : fourdst_lib_dirs
        fourdst_L_args += ['-L' + d]
    endforeach
    foreach name : ['composition', 'logging', 'const', 'reflect_cpp']
        cpp.find_library(name, dirs: fourdst_lib_dirs)
    endforeach
    comp_lib  = declare_dependency(link_args: fourdst_L_args + ['-lcomposition'])
    log_lib   = declare_dependency(link_args: fourdst_L_args + ['-llogging'])
    const_lib = declare_dependency(link_args: fourdst_L_args + ['-lconst'])
    refl_lib  = declare_dependency(link_args: fourdst_L_args + ['-lreflect_cpp'])

    composition_dep = declare_dependency(compile_args: fourdst_inc_args, dependencies: [comp_lib, refl_lib])
    log_dep         = declare_dependency(compile_args: fourdst_inc_args, dependencies: log_lib)
    const_dep       = declare_dependency(compile_args: fourdst_inc_args, dependencies: const_lib)
    config_dep      = declare_dependency(compile_args: fourdst_inc_args)   # header-only libconfig
else
    fourdst_build_lib_all = true
    if not get_option('plugin_support')
        fourdst_build_lib_all=false
        message('Disabling fourdst plugin support as per user request.')
    endif

    fourdst_default_options = [
        'build_tests=' + get_option('build_tests').to_string(),
        'build_python=' + get_option('build_python').to_string(),
        'build_lib_all=' + fourdst_build_lib_all.to_string(),
        'build_lib_comp=true',
        'build_lib_config=true',
        'build_lib_log=true',
        'build_lib_const=true',
        'pkg_config=' + get_option('pkg_config').to_string(),
    ]

    if get_option('build_python')
        fourdst_default_options += ['default_library=static']
    endif

    fourdst_sp = subproject('fourdst', default_options: fourdst_default_options)

    composition_dep = fourdst_sp.get_variable('composition_dep')
    log_dep = fourdst_sp.get_variable('log_dep')
    const_dep = fourdst_sp.get_variable('const_dep')
    config_dep = fourdst_sp.get_variable('config_dep')
    if get_option('plugin_support')
        warning('Including plugin library from fourdst. Note this will bring in minizip-ng and openssl, which can cause build issues with cross compilation due to their complexity.')
        plugin_dep = fourdst_sp.get_variable('plugin_dep')
    endif

    libcomposition = fourdst_sp.get_variable('libcomposition')
    libconst = fourdst_sp.get_variable('libconst')
    liblogging = fourdst_sp.get_variable('liblogging')

    if get_option('plugin_support')
        warning('Including plugin library from fourdst. Note this will bring in minizip-ng and openssl, which can cause build issues with cross compilation due to their complexity.')
        libplugin = fourdst_sp.get_variable('libplugin')
    endif
endif
