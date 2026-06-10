cppad_cmake_options = cmake.subproject_options()

cppad_cmake_options.add_cmake_defines({
                                          'cppad_static_lib': 'true',
                                          'cpp_mas_num_threads': '10',
                                          'cppad_debug_and_release': 'false',
                                          'include_doc': 'false',
                                          'CMAKE_POSITION_INDEPENDENT_CODE': true
                                      })
cppad_cmake_options.set_install(false)

cppad_sp = cmake.subproject(
    'cppad',
    options: cppad_cmake_options,
)

cppad_target = cppad_sp.target('cppad_lib')
cppad_objs = [cppad_target.extract_all_objects(recursive: true)]
cppad_incs = cppad_sp.include_directories('cppad_lib')

empty_cppad_file = configure_file(output: 'cppad_dummy_ar.cpp', command: ['echo'], capture: true)

libcppad_static = static_library(
    'cppad-static',
    empty_cppad_file,
    objects: cppad_objs,
    include_directories: cppad_incs,
    pic: true,
    install: false
)

cppad_dep = declare_dependency(
    include_directories: cppad_incs
)

message('Staging vendored CppAD headers for ' + gridfire_vendor_includedir)

custom_target(
    'vendor_cppad_headers',
    command: copytree_cmd + [
        meson.project_source_root() / 'subprojects' / 'CppAD-20250000.2' / 'include' / 'cppad',
        '@OUTPUT@',
    ],
    output: 'cppad',
    depends: libcppad_static,
    install: true,
    install_dir: gridfire_vendor_includedir,
)