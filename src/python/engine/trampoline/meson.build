gf_engine_trampoline_sources = files('py_engine.cpp')

gf_engine_trapoline_dependencies = [
    gridfire_dep,
    pybind11_dep,
    python3_dep,
]

gf_engine_trampoline_lib = static_library(
    'engine_trampolines',
    gf_engine_trampoline_sources,
    include_directories: include_directories('.'),
    dependencies: gf_engine_trapoline_dependencies,
    install: false,
)

gr_engine_trampoline_dep = declare_dependency(
    link_with: gf_engine_trampoline_lib,
    include_directories: ('.'),
    dependencies: gf_engine_trapoline_dependencies,
)
