gf_screening_trampoline_sources = files('py_screening.cpp')

gf_screening_trapoline_dependencies = [
    gridfire_dep,
    pybind11_dep,
    python3_dep,
]

gf_screening_trampoline_lib = static_library(
    'screening_trampolines',
    gf_screening_trampoline_sources,
    include_directories: include_directories('.'),
    dependencies: gf_screening_trapoline_dependencies,
    install: false,
)

gr_screening_trampoline_dep = declare_dependency(
    link_with: gf_screening_trampoline_lib,
    include_directories: ('.'),
    dependencies: gf_screening_trapoline_dependencies,
)
