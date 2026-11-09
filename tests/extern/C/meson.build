executable('gf_extern_c_single', 'gridfire_evolve_single.c',
  install: false,
  c_args: ['-Wall', '-Wextra'],
  dependencies: [gridfire_extern_dep]
)

executable('gf_extern_c_multi', 'gridfire_evolve_multi.c',
           install: false,
           c_args: ['-Wall', '-Wextra'],
           dependencies: [gridfire_extern_dep]
)

gf_extern_c_config = executable('gf_extern_c_config', 'gridfire_config.c',
           install: false,
           c_args: ['-Wall', '-Wextra'],
           dependencies: [gridfire_extern_dep, meson.get_compiler('c').find_library('m', required: false)]
)
test('extern_c_config', gf_extern_c_config, args: [meson.current_build_dir()], suite: 'extern', timeout: 600)

gf_extern_c_evaluate = executable('gf_extern_c_evaluate', 'gridfire_evaluate.c',
           install: false,
           c_args: ['-Wall', '-Wextra'],
           dependencies: [gridfire_extern_dep, meson.get_compiler('c').find_library('m', required: false)]
)
test('extern_c_evaluate', gf_extern_c_evaluate, suite: 'extern', timeout: 600)

gf_extern_c_zone_status = executable('gf_extern_c_zone_status', 'gridfire_zone_status.c',
           install: false,
           c_args: ['-Wall', '-Wextra'],
           dependencies: [gridfire_extern_dep]
)
test('extern_c_zone_status', gf_extern_c_zone_status, args: [meson.current_build_dir()], suite: 'extern', timeout: 600)

gf_extern_c_he_policy = executable('gf_extern_c_he_policy', 'gridfire_he_policy.c',
           install: false,
           c_args: ['-Wall', '-Wextra'],
           dependencies: [gridfire_extern_dep]
)
test('extern_c_he_policy', gf_extern_c_he_policy, suite: 'extern', timeout: 600)

gf_extern_c_shared_composition = executable('gf_extern_c_shared_composition', 'gridfire_shared_composition.c',
           install: false,
           c_args: ['-Wall', '-Wextra'],
           dependencies: [gridfire_extern_dep]
)
test('extern_c_shared_composition', gf_extern_c_shared_composition, suite: 'extern', timeout: 600)

gf_extern_c_dsep_composition = executable('gf_extern_c_dsep_composition', 'gridfire_dsep_composition.c',
           install: false,
           c_args: ['-Wall', '-Wextra'],
           dependencies: [gridfire_extern_dep, meson.get_compiler('c').find_library('m', required: false)]
)
test('extern_c_dsep_composition', gf_extern_c_dsep_composition, suite: 'extern', timeout: 600)

gf_extern_c_rate_multipliers = executable('gf_extern_c_rate_multipliers', 'gridfire_rate_multipliers.c',
           install: false,
           c_args: ['-Wall', '-Wextra'],
           dependencies: [gridfire_extern_dep, meson.get_compiler('c').find_library('m', required: false)]
)
test('extern_c_rate_multipliers', gf_extern_c_rate_multipliers, suite: 'extern', timeout: 600)

if get_option('openmp_support')
    # Separate contexts used from the caller's own OpenMP threads must not interfere
    gf_extern_c_threads = executable('gf_extern_c_threads', 'gridfire_threads.c',
               install: false,
               c_args: ['-Wall', '-Wextra'],
               dependencies: [gridfire_extern_dep, dependency('openmp', language: 'c'),
                              meson.get_compiler('c').find_library('m', required: false)]
    )
    test('extern_c_threads', gf_extern_c_threads, suite: 'extern', timeout: 900, env: ['OMP_NUM_THREADS=8'])
endif
