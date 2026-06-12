# Google Test dependency
#gtest_dep = dependency('gtest', main: true, required : true)
#gtest_main = dependency('gtest_main', required: true)
#gtest_nomain_dep = dependency('gtest', main: false, required : true)

# Subdirectories for unit and integration tests
if get_option('build_tests') and not get_option('build_python')
    subdir('graphnet_sandbox')
    subdir('extern')
else
    message('Tests disabled by build options! To enable them build with build_tests = true and build_python=false')
endif
