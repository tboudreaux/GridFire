# --- Python Extension Setup ---
py_installation = import('python').find_installation('python3', pure: false)

gridfire_py_deps = [
    pybind11_dep,
    const_dep,
    config_dep,
    composition_dep,
    gridfire_dep
]

py_mod = py_installation.extension_module(
    '_gridfire', # Name of the generated .so/.pyd file (without extension)
    sources: [
        meson.project_source_root() + '/src/python/bindings.cpp',
        meson.project_source_root() + '/src/python/types/bindings.cpp',
        meson.project_source_root() + '/src/python/partition/bindings.cpp',
        meson.project_source_root() + '/src/python/partition/trampoline/py_partition.cpp',
        meson.project_source_root() + '/src/python/reaction/bindings.cpp',
        meson.project_source_root() + '/src/python/screening/bindings.cpp',
        meson.project_source_root() + '/src/python/screening/trampoline/py_screening.cpp',
        meson.project_source_root() + '/src/python/io/bindings.cpp',
        meson.project_source_root() + '/src/python/io/trampoline/py_io.cpp',
        meson.project_source_root() + '/src/python/exceptions/bindings.cpp',
        meson.project_source_root() + '/src/python/engine/bindings.cpp',
        meson.project_source_root() + '/src/python/engine/trampoline/py_engine.cpp',
        meson.project_source_root() + '/src/python/solver/bindings.cpp',
        meson.project_source_root() + '/src/python/solver/trampoline/py_solver.cpp',
        meson.project_source_root() + '/src/python/policy/bindings.cpp',
        meson.project_source_root() + '/src/python/policy/trampoline/py_policy.cpp',
        meson.project_source_root() + '/src/python/utils/bindings.cpp',
    ],
    dependencies : gridfire_py_deps,
    cpp_args : ['-UNDEBUG'], # Example: Ensure assertions are enabled if needed
    install : true,
    subdir: 'gridfire',
)


py_installation.install_sources(
    files(
        meson.project_source_root() + '/src/python/gridfire/__init__.py',
        meson.project_source_root() + '/stubs/gridfire/_gridfire/__init__.pyi',
        meson.project_source_root() + '/stubs/gridfire/_gridfire/exceptions.pyi',
        meson.project_source_root() + '/stubs/gridfire/_gridfire/partition.pyi',
        meson.project_source_root() + '/stubs/gridfire/_gridfire/reaction.pyi',
        meson.project_source_root() + '/stubs/gridfire/_gridfire/screening.pyi',
        meson.project_source_root() + '/stubs/gridfire/_gridfire/io.pyi',
        meson.project_source_root() + '/stubs/gridfire/_gridfire/solver.pyi',
        meson.project_source_root() + '/stubs/gridfire/_gridfire/policy.pyi',
        meson.project_source_root() + '/stubs/gridfire/_gridfire/type.pyi'
    ),
    subdir: 'gridfire',
)

py_installation.install_sources(
    files(
        meson.project_source_root() + '/stubs/gridfire/_gridfire/engine/__init__.pyi',
        meson.project_source_root() + '/stubs/gridfire/_gridfire/engine/diagnostics.pyi',
    ),
    subdir: 'gridfire/engine',
)

py_installation.install_sources(
    files(
        meson.project_source_root() + '/stubs/gridfire/_gridfire/utils/__init__.pyi',
    ),
    subdir: 'gridfire/utils',
)

py_installation.install_sources(
    files(
        meson.project_source_root() + '/stubs/gridfire/_gridfire/utils/hashing/__init__.pyi',
        meson.project_source_root() + '/stubs/gridfire/_gridfire/utils/hashing/reaction.pyi',
    ),
    subdir: 'gridfire/utils/hashing',
)
