subdir('trampoline')

# Define the library
bindings_sources = files('bindings.cpp')
bindings_headers = files('bindings.h')

dependencies = [
    gridfire_dep,
    python3_dep,
    pybind11_dep,
]

message('⏳ Python bindings for GridFire Engine are being registered...')
shared_module('py_gf_engine',
              bindings_sources,
              cpp_args: ['-fvisibility=default'],
              install : true,
              dependencies: dependencies,
              include_directories: include_directories('.')
)
message('✅ Python bindings for GridFire Engine registered successfully!')