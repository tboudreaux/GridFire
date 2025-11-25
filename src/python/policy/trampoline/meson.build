gf_policy_trampoline_sources = files('py_policy.cpp')

gf_policy_trapoline_dependencies = [
    gridfire_dep,
    pybind11_dep,
    python3_dep,
]

gf_policy_trampoline_lib = static_library(
    'policy_trampolines',
    gf_policy_trampoline_sources,
    include_directories: include_directories('.'),
    dependencies: gf_policy_trapoline_dependencies,
    install: false,
)

gr_policy_trampoline_dep = declare_dependency(
    link_with: gf_policy_trampoline_lib,
    include_directories: ('.'),
    dependencies: gf_policy_trapoline_dependencies,
)
