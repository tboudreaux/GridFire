cmake = import('cmake')

kinsol_cmake_options = cmake.subproject_options()

kinsol_cmake_options.add_cmake_defines({
    'CMAKE_CXX_FLAGS'   : '-Wno-deprecated-declarations',
    'CMAKE_C_FLAGS'     : '-Wno-deprecated-declarations',
    'BUILD_SHARED_LIBS' : 'OFF',
    'BUILD_STATIC_LIBS' : 'ON',
    'CMAKE_BUILD_WITH_INSTALL_RPATH': 'ON',
    'EXAMPLES_ENABLE_C' : 'OFF',
    'CMAKE_POSITION_INDEPENDENT_CODE': 'ON'
})

kinsol_cmake_options.add_cmake_defines({
    'CMAKE_INSTALL_LIBDIR': get_option('libdir'),
    'CMAKE_INSTALL_INCLUDEDIR': get_option('includedir')
})

kinsol_sp = cmake.subproject(
    'kinsol',
    options: kinsol_cmake_options,
)

sundials_kinsol_shared = kinsol_sp.dependency('sundials_kinsol_static')

kinsol_dep = declare_dependency(
    dependencies: [
        sundials_kinsol_shared,
    ]
)

