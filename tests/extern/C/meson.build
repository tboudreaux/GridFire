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
