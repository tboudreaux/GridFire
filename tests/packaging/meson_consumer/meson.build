# A downstream Meson project that uses GridFire as a subproject: the C++ library, the C API, and the Fortran module
# compiled by this project's own Fortran compiler. See gridfire.wrap for the wrap file such a project ships.
project('gridfire_meson_consumer', ['c', 'cpp', 'fortran'],
    # b_ndebug=if-release: with assertions on, Meson 1.9 adds -D_LIBCPP_ENABLE_ASSERTIONS, which recent Apple libc++
    # rejects (for every C++ project, GridFire's dependencies included)
    default_options: ['cpp_std=c++23', 'b_ndebug=if-release'],
    meson_version: '>=1.5.0')

# Require a version: an older GridFire found on the system (e.g. a stale /usr/local install) is then skipped in
# favour of the subproject, instead of being mixed with it. The first lookup configures the subproject; its options
# go here.
gridfire_version = '>=1.2.0'
gridfire_dep = dependency('gridfire', version: gridfire_version, default_options: [
    'build_fortran=true',
    'build_tests=false',
    'build_tools=false',
    'build_examples=false',
    'build_benchmarks=false',
    'pkg_config=false',
])
gridfire_extern_dep = dependency('gridfire_extern', version: gridfire_version)
gridfire_fortran_dep = dependency('gridfire_fortran', version: gridfire_version)

test('cpp', executable('consumer_cpp', 'consumer.cpp', dependencies: gridfire_dep))
test('c', executable('consumer_c', 'consumer.c', dependencies: gridfire_extern_dep))
test('fortran', executable('consumer_f', 'consumer.f90', dependencies: gridfire_fortran_dep))
