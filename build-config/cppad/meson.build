#cppad_inc = include_directories('include', is_system: true)
#cppad_dep = declare_dependency(
#  include_directories: cppad_inc,
#)
#
#message('Registering CppAD headers for installation...')
#install_subdir('include/cppad', install_dir: get_option('includedir'))
#message('Done registering CppAD headers for installation!')
#


cppad_cmake_options = cmake.subproject_options()

cppad_cmake_options.add_cmake_defines({
    'cppad_static_lib': 'true',
    'cpp_mas_num_threads': '10',
    'cppad_debug_and_release': 'false',
    'include_doc': 'false'
})

cppad_cmake_options.set_install(false)

cppad_sp = cmake.subproject(
    'cppad',
    options: cppad_cmake_options,
)

cppad_dep = cppad_sp.dependency('cppad_lib').as_system()