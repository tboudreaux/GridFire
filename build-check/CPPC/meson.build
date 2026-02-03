cppc = meson.get_compiler('cpp')

if cppc.get_id() == 'clang'
    message('disabling bitwise-instead-of-logical warnings for clang')
    add_project_arguments('-Wno-bitwise-instead-of-logical', language: 'cpp')
endif

if cppc.get_id() == 'gcc'
    message('disabling psabi warnings for gcc')
    add_project_arguments('-Wno-psabi', language: 'cpp')
    
    if (cppc.version().version_compare('<14.0'))
      error('g++ version must be at least 14.0, found ' + cppc.version())
    endif
endif

if not cppc.has_header('print')
    error('C++ standard library header <print> not found. Please ensure your compiler and standard library supports C++23. We have already validated your compiler version so this is likely an issue with your standard library installation.')
endif
if not cppc.has_header('format')
    error('C++ standard library header <format> not found. Please ensure your compiler and standard library supports C++23. We have already validated your compiler version so this is likely an issue with your standard library installation.')
endif

# For Eigen
add_project_arguments('-Wno-deprecated-declarations', language: 'cpp')

if get_option('build_python')
    message('enabling hidden visibility for C++ symbols when building Python extension. This reduces the size of the resulting shared library.')
    add_project_arguments('-fvisibility=hidden', language: 'cpp')
else
    message('enabling default visibility for C++ symbols')
    add_project_arguments('-fvisibility=default', language: 'cpp')
endif

if get_option('openmp_support')
    gridfire_args += ['-DGF_USE_OPENMP']
endif

if get_option('asan') and get_option('buildtype') != 'debug' and get_option('buildtype') != 'debugoptimized'
    error('AddressSanitizer (ASan) can only be enabled for debug or debugoptimized builds')
endif

if get_option('asan') and (get_option('buildtype') == 'debugoptimized' or get_option('buildtype') == 'debug')
    message('enabling AddressSanitizer (ASan) support')
    add_project_arguments('-fsanitize=address,undefined', language: 'cpp')
    add_project_arguments('-fno-omit-frame-pointer', language: 'cpp')

    add_project_link_arguments('-fsanitize=address,undefined', language: 'cpp')
    add_project_link_arguments('-fno-omit-frame-pointer', language: 'cpp')
endif
