cmake = import('cmake')
cvode_cmake_options = cmake.subproject_options()

cvode_cmake_options.add_cmake_defines({
    'CMAKE_CXX_FLAGS'   : '-Wno-deprecated-declarations',
    'CMAKE_C_FLAGS'     : '-Wno-deprecated-declarations',
    'BUILD_SHARED_LIBS' : 'OFF',
    'BUILD_STATIC_LIBS' : 'ON',
    'EXAMPLES_ENABLE_C' : 'OFF',
    'CMAKE_POSITION_INDEPENDENT_CODE': true

})

cvode_cmake_options.add_cmake_defines({
    'CMAKE_INSTALL_LIBDIR': get_option('libdir'),
    'CMAKE_INSTALL_INCLUDEDIR': get_option('includedir')
})

if meson.is_cross_build() and host_machine.system() == 'emscripten'
  cvode_cmake_options.add_cmake_defines({
      'CMAKE_C_FLAGS': '-s MEMORY64=1 -s ALLOW_MEMORY_GROWTH=1',
      'CMAKE_CXX_FLAGS': '-s MEMORY64=1 -s ALLOW_MEMORY_GROWTH=1',
      'CMAKE_SHARED_LINKER_FLAGS': '-s MEMORY64=1 -s ALLOW_MEMORY_GROWTH=1',
      'CMAKE_EXE_LINKER_FLAGS': '-s MEMORY64=1 -s ALLOW_MEMORY_GROWTH=1'
  })
endif

cvode_sp = cmake.subproject(
    'cvode',
    options: cvode_cmake_options,
)

sundials_core_tgt = cvode_sp.target('sundials_core_static')
sundials_cvode_tgt = cvode_sp.target('sundials_cvode_static')
sundials_nvecserial_tgt = cvode_sp.target('sundials_nvecserial_static')
sundials_sunmatrixdense_tgt = cvode_sp.target('sundials_sunmatrixdense_static')
sundials_sunlinsoldense_tgt = cvode_sp.target('sundials_sunlinsoldense_static')

cvode_objs = [
  sundials_core_tgt.extract_all_objects(recursive: true),
  sundials_cvode_tgt.extract_all_objects(recursive: true),
  sundials_nvecserial_tgt.extract_all_objects(recursive: true),
  sundials_sunmatrixdense_tgt.extract_all_objects(recursive: true),
  sundials_sunlinsoldense_tgt.extract_all_objects(recursive: true),
]

sundials_core_includes = cvode_sp.include_directories('sundials_core_static')
sundials_cvode_includes = cvode_sp.include_directories('sundials_cvode_static')
sundials_nvecserial_includes = cvode_sp.include_directories('sundials_nvecserial_static')
sundials_sunmatrixdense_includes = cvode_sp.include_directories('sundials_sunmatrixdense_static')
sundials_sunlinsoldense_includes = cvode_sp.include_directories('sundials_sunlinsoldense_static')

cvode_includes = [
  sundials_core_includes,
  sundials_cvode_includes,
  sundials_nvecserial_includes,
  sundials_sunmatrixdense_includes,
  sundials_sunlinsoldense_includes
]


empty_cvode_file = configure_file(
  output: 'cvode_dummy_ar.cpp',
  command: ['echo'],
  capture: true
  )



libcvode_static = static_library(
   'cvode-static',
   empty_cvode_file,
   objects: cvode_objs,
   include_directories: cvode_includes,
   pic: true,
   install: false
)


cvode_dep = declare_dependency(
  link_with: libcvode_static,
  include_directories: cvode_includes,
)


