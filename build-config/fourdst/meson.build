# bring in all of the fourdst utility repositories

fourdst_sp = subproject('fourdst',
                        default_options:
                            ['build-tests=' + get_option('build-tests').to_string(),
                             'build-python=' + get_option('build-python').to_string()
                            ]
)

composition_dep = fourdst_sp.get_variable('composition_dep')
log_dep = fourdst_sp.get_variable('log_dep')
const_dep = fourdst_sp.get_variable('const_dep')
config_dep = fourdst_sp.get_variable('config_dep')
plugin_dep = fourdst_sp.get_variable('plugin_dep')

libcomposition = fourdst_sp.get_variable('libcomposition')
libconst = fourdst_sp.get_variable('libconst')
libconfig = fourdst_sp.get_variable('libconfig')
liblogging = fourdst_sp.get_variable('liblogging')
libplugin = fourdst_sp.get_variable('libplugin')