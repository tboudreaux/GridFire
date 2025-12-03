# bring in all of the fourdst utility repositories

fourdst_build_lib_all = true
if not get_option('plugin_support')
    fourdst_build_lib_all=false
    message('Disabling fourdst plugin support as per user request.')
endif

fourdst_sp = subproject('fourdst',
                        default_options:
                            ['build_tests=' + get_option('build_tests').to_string(),
                             'build_python=' + get_option('build_python').to_string(),
                             'build_lib_all=' + fourdst_build_lib_all.to_string(),
                             'build_lib_comp=true',
                             'build_lib_config=true',
                             'build_lib_log=true',
                             'build_lib_const=true',
                             'pkg_config=' + get_option('pkg_config').to_string(),
                            ]
)

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
libconfig = fourdst_sp.get_variable('libconfig')
liblogging = fourdst_sp.get_variable('liblogging')

if get_option('plugin_support')
    warning('Including plugin library from fourdst. Note this will bring in minizip-ng and openssl, which can cause build issues with cross compilation due to their complexity.')
    libplugin = fourdst_sp.get_variable('libplugin')
endif
