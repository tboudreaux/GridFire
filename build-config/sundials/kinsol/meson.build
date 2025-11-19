cmake = import('cmake')

kinsol_cmake_options = cmake.subproject_options()

kinsol_cmake_options.add_cmake_defines({
    'CMAKE_CXX_FLAGS'   : '-Wno-deprecated-declarations',
    'CMAKE_C_FLAGS'     : '-Wno-deprecated-declarations',
    'BUILD_SHARED_LIBS' : 'ON',
    'BUILD_STATIC_LIBS' : 'OFF',
})

kinsol_cmake_options.add_cmake_defines({
    'CMAKE_INSTALL_LIBDIR': get_option('libdir'),
    'CMAKE_INSTALL_INCLUDEDIR': get_option('includedir')
})

kinsol_sp = cmake.subproject(
    'kinsol',
    options: kinsol_cmake_options,
)

sundials_kinsol_shared = kinsol_sp.dependency('sundials_kinsol_shared')

kinsol_dep = declare_dependency(
    dependencies: [
        sundials_kinsol_shared,
    ]
)

