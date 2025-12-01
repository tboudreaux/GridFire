py_installation = import('python').find_installation('python3', pure: false)


if meson.is_cross_build() and host_machine.system() == 'darwin'
  py_ver = get_option('python-target-version')
  message('Cross build on Darwin, using python version ' + py_ver)
  py_inc_dir = include_directories('../../cross/python_includes/python-' + py_ver + '/include/python' + py_ver)
  py_dep = declare_dependency(include_directories: py_inc_dir)
  py_module_prefix = ''
  py_module_suffic = 'so'
  meson.override_dependency('python3', py_dep)
else
  py_dep = py_installation.dependency()
  py_module_prefix = ''
  py_module_suffic = 'so'
endif
