cppad_inc = include_directories('include', is_system: true)
cppad_dep = declare_dependency(
  include_directories: cppad_inc,
)

message('Registering CppAD headers for installation...')
install_subdir('include/cppad', install_dir: get_option('includedir'))
message('Done registering CppAD headers for installation!')


