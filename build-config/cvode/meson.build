cmake = import('cmake')
cvode_cmake_options = cmake.subproject_options()

cvode_cmake_options.add_cmake_defines({
    'CMAKE_CXX_FLAGS'   : '-Wno-deprecated-declarations',
    'CMAKE_C_FLAGS'     : '-Wno-deprecated-declarations',
    'BUILD_SHARED_LIBS' : 'ON',
    'BUILD_STATIC_LIBS' : 'OFF',

})

cvode_cmake_options.add_cmake_defines({
    'CMAKE_INSTALL_LIBDIR': get_option('libdir'),
    'CMAKE_INSTALL_INCLUDEDIR': get_option('includedir')
})

cvode_sp = cmake.subproject(
    'cvode',
    options: cvode_cmake_options,
)

# For the core SUNDIALS library (SUNContext, etc.)
sundials_core_dep = cvode_sp.dependency('sundials_core_shared')

# For the CVODE integrator library
sundials_cvode_dep = cvode_sp.dependency('sundials_cvode_shared')

# For the serial NVector library
sundials_nvecserial_dep = cvode_sp.dependency('sundials_nvecserial_shared')

# For the dense matrix library
sundials_sunmatrixdense_dep = cvode_sp.dependency('sundials_sunmatrixdense_shared')

# For the dense linear solver library
sundials_sunlinsoldense_dep = cvode_sp.dependency('sundials_sunlinsoldense_shared')

sundials_dep = declare_dependency(
    dependencies: [
        sundials_core_dep,
        sundials_cvode_dep,
        sundials_nvecserial_dep,
        sundials_sunmatrixdense_dep,
        sundials_sunlinsoldense_dep,
    ],
)


