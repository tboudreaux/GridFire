# MESA (bbq) comparison. Reference data lives in reference/ and is produced by generate_reference.py on a machine with
# MESA + bbq; without it the comparison test skips. Run with `meson test --suite validation`.
mesa_validation_args = ['-DGF_MESA_VALIDATION_DIR="' + meson.current_source_dir() + '"']

mesa_case_lib = static_library(
    'gf_mesa_case',
    'mesa_case.cpp',
    cpp_args: mesa_validation_args,
    dependencies: [gtest_dep, gridfire_dep, json_dep],
)
mesa_case_dep = declare_dependency(
    link_with: mesa_case_lib,
    include_directories: include_directories('.'),
    dependencies: [json_dep],
)

# Exports a case's network and initial state for bbq (used by generate_reference.py)
gf_mesa_setup = executable(
    'gf_mesa_setup',
    'gf_mesa_setup.cpp',
    dependencies: [gridfire_dep, mesa_case_dep],
    install: false,
)

# Per-reaction energy generation at one state, the GridFire counterpart of mesa_rates/ (diagnostics only)
gf_mesa_rates = executable(
    'gf_mesa_rates',
    'gf_mesa_rates.cpp',
    dependencies: [gridfire_dep, mesa_case_dep],
    install: false,
)

mesa_comparison_test = executable(
    'mesa_comparison_test',
    'mesa_comparison_test.cpp',
    dependencies: [test_support_dep, mesa_case_dep],
    install: false,
)
test('mesa_comparison_test', mesa_comparison_test, suite: 'validation', timeout: 3600, protocol: 'gtest')
