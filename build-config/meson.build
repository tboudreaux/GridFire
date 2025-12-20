cmake = import('cmake')

if get_option('build_python')
  subdir('python')
  subdir('pybind')
endif
subdir('fourdst')
subdir('sundials')

subdir('cppad')
subdir('xxHash')
subdir('eigen')

subdir('json')

subdir('CLI11')
subdir('unordered_dense')

if get_option('use_mimalloc')
  subdir('mimalloc')
endif
