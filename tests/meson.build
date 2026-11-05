# Google Test dependency. include_type 'non-system' turns gtest's -isystem into -I so that, listed first, its headers
# take precedence over any other googletest copy reachable through gridfire's include path (e.g. a conda prefix).
gtest_dep = dependency('gtest', main: false, required: true, include_type: 'non-system')
gtest_main_dep = dependency('gtest_main', required: true, include_type: 'non-system')

if get_option('build_tests')
    # Shared fixtures and helpers (tests/support)
    test_support_lib = static_library(
        'gf_test_support',
        'support/gf_test_support.cpp',
        dependencies: [gtest_dep, gridfire_dep], # gtest first so its headers win over any copy on gridfire's include path
    )
    test_support_dep = declare_dependency(
        link_with: test_support_lib,
        include_directories: include_directories('support'),
        dependencies: [gtest_dep, gtest_main_dep, gridfire_dep],
    )

    # Each entry: [source, suite, timeout in seconds]. Run one tier with `meson test --suite <suite>`.
    gridfire_tests = [
        ['unit/energy_channels_test.cpp', 'unit', 120],
        ['unit/network_data_test.cpp', 'unit', 120],
        ['unit/screening_test.cpp', 'unit', 120],
        ['unit/screening_cd09_test.cpp', 'unit', 300],
        ['unit/config_test.cpp', 'unit', 300],
        ['unit/policy_test.cpp', 'unit', 300],
        ['property/engine_properties_test.cpp', 'property', 900],
        ['integration/solver_integration_test.cpp', 'integration', 1800],
    ]

    foreach t : gridfire_tests
        test_name = t[0].split('/')[-1].split('.')[0]
        test_exe = executable(
            test_name,
            t[0],
            dependencies: [test_support_dep],
            install: false,
        )
        test(test_name, test_exe, suite: t[1], timeout: t[2], protocol: 'gtest')
    endforeach

    subdir('validation')
    subdir('graphnet_sandbox')
    subdir('extern')
else
    message('Tests disabled by build options! To enable them build with build_tests = true')
endif
