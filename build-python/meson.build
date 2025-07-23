# --- Python Extension Setup ---
py_installation = import('python').find_installation('python3')

py_mod = py_installation.extension_module(
    'gridfire', # Name of the generated .so/.pyd file (without extension)
    sources: [
        meson.project_source_root() + '/src/python/bindings.cpp',
        meson.project_source_root() + '/src/python/types/bindings.cpp',
        meson.project_source_root() + '/src/python/expectations/bindings.cpp',
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
        meson.project_source_root() + '/src/python/utils/bindings.cpp',
    ],
    dependencies : [
        pybind11_dep,
        const_dep,
        config_dep,
        composition_dep,
        gridfire_dep,
    ],
    cpp_args : ['-UNDEBUG'], # Example: Ensure assertions are enabled if needed
    install : true,
)