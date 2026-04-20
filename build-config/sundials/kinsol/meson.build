cmake = import('cmake')

kinsol_cmake_options = cmake.subproject_options()

kinsol_cmake_options.add_cmake_defines({
    'CMAKE_CXX_FLAGS'   : '-Wno-deprecated-declarations',
    'CMAKE_C_FLAGS'     : '-Wno-deprecated-declarations',
    'BUILD_SHARED_LIBS' : 'OFF',
    'BUILD_STATIC_LIBS' : 'ON',
    'EXAMPLES_ENABLE_C' : 'OFF',
    'CMAKE_POSITION_INDEPENDENT_CODE': true,
    'CMAKE_PLATFORM_NO_VERSIONED_SONAME': 'ON',
    'SUNDIALS_LOGGING_LEVEL': 1
})

kinsol_cmake_options.add_cmake_defines({
    'CMAKE_INSTALL_LIBDIR': get_option('libdir'),
    'CMAKE_INSTALL_INCLUDEDIR': get_option('includedir')
})

kinsol_cmake_options.set_install(false)

kinsol_sp = cmake.subproject(
    'kinsol',
    options: kinsol_cmake_options,
)

sundials_kinsol_static_tgt = kinsol_sp.target('sundials_kinsol_obj_static')
kinsol_includes = kinsol_sp.include_directories('sundials_kinsol_obj_static')

kinsol_objs = [sundials_kinsol_static_tgt.extract_all_objects(recursive: false)]

empty_kinsol_file = configure_file(
  output: 'kinsol_dummy_ar.cpp',
  command: ['echo'],
  capture: true
  )


libkinsol_static = static_library(
  'kinsol_static',
  empty_kinsol_file,
  objects: kinsol_objs,
  include_directories: kinsol_includes,
  pic: true,
  install: false
)


kinsol_dep = declare_dependency(
  link_with: libkinsol_static,
  include_directories: kinsol_includes
)


