gf_solver_trampoline_sources = files('py_solver.cpp')

gf_solver_trapoline_dependencies = [
    gridfire_dep,
    pybind11_dep,
    python3_dep,
]

gf_solver_trampoline_lib = static_library(
    'solver_trampolines',
    gf_solver_trampoline_sources,
    include_directories: include_directories('.'),
    dependencies: gf_solver_trapoline_dependencies,
    install: false,
)

gr_solver_trampoline_dep = declare_dependency(
    link_with: gf_solver_trampoline_lib,
    include_directories: ('.'),
    dependencies: gf_solver_trapoline_dependencies,
)
