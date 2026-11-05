executable('gf_fortran_single_zone_test', 'gridfire_evolve_single.f90',
  install: false,
  fortran_args: ['-Wall', '-Wextra'],
  dependencies: [gridfire_fortran_dep]
)

executable('gf_fortran_multi_zone_test', 'gridfire_evolve_multi.f90',
           install: false,
           fortran_args: ['-Wall', '-Wextra'],
           dependencies: [gridfire_fortran_dep]
)

gf_fortran_config_test = executable('gf_fortran_config_test', 'gridfire_config.f90',
           install: false,
           fortran_args: ['-Wall', '-Wextra'],
           dependencies: [gridfire_fortran_dep]
)
test('extern_fortran_config', gf_fortran_config_test, args: [meson.current_build_dir()], suite: 'extern', timeout: 600)

gf_fortran_evaluate_test = executable('gf_fortran_evaluate_test', 'gridfire_evaluate.f90',
           install: false,
           fortran_args: ['-Wall', '-Wextra'],
           dependencies: [gridfire_fortran_dep]
)
test('extern_fortran_evaluate', gf_fortran_evaluate_test, suite: 'extern', timeout: 600)

gf_fortran_zone_status_test = executable('gf_fortran_zone_status_test', 'gridfire_zone_status.f90',
           install: false,
           fortran_args: ['-Wall', '-Wextra'],
           dependencies: [gridfire_fortran_dep]
)
test('extern_fortran_zone_status', gf_fortran_zone_status_test, suite: 'extern', timeout: 600)
