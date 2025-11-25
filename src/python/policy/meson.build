# Define the library
subdir('trampoline')

bindings_sources = files('bindings.cpp')
bindings_headers = files('bindings.h')

dependencies = [
    gridfire_dep,
    python3_dep,
    pybind11_dep,
]

shared_module('py_gf_policy',
              bindings_sources,
              cpp_args: ['-fvisibility=default'],
              install : true,
              dependencies: dependencies,
              include_directories: include_directories('.')
)
