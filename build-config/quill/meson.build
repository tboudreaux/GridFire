quill_cmake_options = cmake.subproject_options()
quill_cmake_options.add_cmake_defines({
                                     'BUILD_SHARED_LIBS': 'ON',
                                     'CMAKE_SKIP_INSTALL_RULES': 'ON'
})
quill_sp = cmake.subproject(
    'quill',
    options: quill_cmake_options,
)
quill_dep = quill_sp.dependency('quill')