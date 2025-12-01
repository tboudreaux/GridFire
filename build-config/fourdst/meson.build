# bring in all of the fourdst utility repositories

fourdst_build_lib_all = true
if get_option('unity-safe')
    fourdst_build_lib_all=false
endif

fourdst_sp = subproject('fourdst',
                        default_options:
                            ['build-tests=' + get_option('build-tests').to_string(),
                             'build-python=' + get_option('build-python').to_string(),
                             'build-lib-all=' + fourdst_build_lib_all.to_string(),
                             'pkg-config=' + get_option('pkg-config').to_string(),
                             'build-lib-log=true'
                            ]
)

composition_dep = fourdst_sp.get_variable('composition_dep')
log_dep = fourdst_sp.get_variable('log_dep')
const_dep = fourdst_sp.get_variable('const_dep')
config_dep = fourdst_sp.get_variable('config_dep')
if not get_option('unity-safe')
    plugin_dep = fourdst_sp.get_variable('plugin_dep')
endif

libcomposition = fourdst_sp.get_variable('libcomposition')
libconst = fourdst_sp.get_variable('libconst')
libconfig = fourdst_sp.get_variable('libconfig')
liblogging = fourdst_sp.get_variable('liblogging')

if not get_option('unity-safe')
    libplugin = fourdst_sp.get_variable('libplugin')
endif
