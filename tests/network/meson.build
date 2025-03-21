# Test files for network
test_sources = [
    'approx8Test.cpp',
]

foreach test_file : test_sources
  exe_name = test_file.split('.')[0]
  message('Building test: ' + exe_name)

  # Create an executable target for each test
  test_exe = executable(
      exe_name,
      test_file,
      dependencies: [gtest_dep, network_dep, gtest_main],
      include_directories: include_directories('../../src/network/public'),
      link_with: libnetwork,  # Link the dobj library
      install_rpath: '@loader_path/../../src'  # Ensure runtime library path resolves correctly
  )

  # Add the executable as a test
  test(
    exe_name,
    test_exe,
    env: ['MESON_SOURCE_ROOT=' + meson.project_source_root(), 'MESON_BUILD_ROOT=' + meson.project_build_root()])
endforeach
