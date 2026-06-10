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

if get_option('build_python')
    sp_root = meson.project_source_root() / 'subprojects'

    fourdst_header_trees = [

        ['config',
         sp_root / 'libconfig' / 'src' / 'config' / 'include' / 'fourdst' / 'config',
         gridfire_includedir / 'fourdst'],
        ['composition',
         sp_root / 'libcomposition' / 'src' / 'composition' / 'include' / 'fourdst' / 'composition',
         gridfire_includedir / 'fourdst'],
        ['atomic',
         sp_root / 'libcomposition' / 'src' / 'composition' / 'include' / 'fourdst' / 'atomic',
         gridfire_includedir / 'fourdst'],
        ['constants',
         sp_root / 'libconstants' / 'src' / 'constants' / 'include' / 'fourdst' / 'constants',
         gridfire_includedir / 'fourdst'],
        ['logging',
         sp_root / 'liblogging' / 'src' / 'logging' / 'include' / 'fourdst' / 'logging',
         gridfire_includedir / 'fourdst'],
        ['toml++',
         sp_root / 'libconfig' / 'build-config' / 'tomlpp' / 'vendor' / 'include' / 'toml++',
         gridfire_fourdst_vendor_includedir],
        ['quill',
         sp_root / 'quill' / 'include' / 'quill',
         gridfire_fourdst_vendor_includedir],
        ['CLI',
         sp_root / 'CLI11-2.6.1' / 'include' / 'CLI',
         gridfire_fourdst_vendor_includedir],
    ]

    foreach t : fourdst_header_trees
        custom_target(
            'wheel_headers_' + t[0].underscorify(),
            command: copytree_cmd + [t[1], '@OUTPUT@'],
            output: t[0],
            install: true,
            install_dir: t[2],
        )
    endforeach
endif
