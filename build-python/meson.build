if get_option('build_python')
    message('Building Python bindings...')

    gridfire_py_deps = [
        py_dep,
        pybind11_dep,
        gridfire_dep
    ]

    if host_machine.system() == 'darwin'
        gridfire_ext_rpath = '@loader_path/lib'
    else
        gridfire_ext_rpath = '$ORIGIN/lib'
    endif

    py_sources =  [
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
        meson.project_source_root() + '/src/python/config/bindings.cpp',
        meson.project_source_root() + '/src/python/engine/scratchpads/bindings.cpp',
    ]


    py_mod = py_installation.extension_module(
        '_gridfire',
        sources: py_sources,
        dependencies: gridfire_py_deps,
        install: true,
        link_args: gridfire_ext_rpath_args,
        build_rpath: gridfire_ext_rpath,
        install_rpath: gridfire_ext_rpath,
        subdir: 'gridfire',
    )

    py_installation.install_sources(
        files(
            meson.project_source_root() + '/src/python/gridfire/__init__.py',
            meson.project_source_root() + '/stubs/gridfire/_gridfire/__init__.pyi',
            meson.project_source_root() + '/stubs/gridfire/_gridfire/config.pyi',
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
            meson.project_source_root() + '/stubs/gridfire/_gridfire/engine/scratchpads.pyi'
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
else
    message('Python bindings disabled')
endif
