gf_partition_trampoline_sources = files('py_partition.cpp')

gf_partition_trapoline_dependencies = [
    gridfire_dep,
    pybind11_dep,
    python3_dep,
]

gf_partition_trampoline_lib = static_library(
    'partition_trampolines',
    gf_partition_trampoline_sources,
    include_directories: include_directories('.'),
    dependencies: gf_partition_trapoline_dependencies,
    install: false,
)

gr_partition_trampoline_dep = declare_dependency(
    link_with: gf_partition_trampoline_lib,
    include_directories: ('.'),
    dependencies: gf_partition_trapoline_dependencies,
)
