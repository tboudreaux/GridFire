gf_io_trampoline_sources = files('py_io.cpp')

gf_io_trapoline_dependencies = [
    gridfire_dep,
    pybind11_dep,
    python3_dep,
]

gf_io_trampoline_lib = static_library(
    'io_trampolines',
    gf_io_trampoline_sources,
    include_directories: include_directories('.'),
    dependencies: gf_io_trapoline_dependencies,
    install: false,
)

gr_io_trampoline_dep = declare_dependency(
    link_with: gf_io_trampoline_lib,
    include_directories: ('.'),
    dependencies: gf_io_trapoline_dependencies,
)
