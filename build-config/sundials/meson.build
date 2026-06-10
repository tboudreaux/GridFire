subdir('cvode')
subdir('kinsol')

sundials_dep = declare_dependency(
    dependencies: [
        cvode_dep,
        kinsol_dep,
    ],
)

# Vendor the SUNDIALS public headers (see vendor/meson.build). GridFire's
# installed public headers #include <sundials/...> etc., so consumers of
# gridfire.pc must compile against exactly the headers this build used --
# never whatever SUNDIALS happens to be in /usr/local/include (which may be
# a different version or configuration, e.g. an MPI-enabled build whose
# sundials_config.h pulls in <mpi.h>).
subdir('vendor')
