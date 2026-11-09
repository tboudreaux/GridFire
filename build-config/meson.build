cmake = import('cmake')

if get_option('build_python')
  subdir('python')
  subdir('pybind')
endif

if get_option('build_python')
  gridfire_pkg_dir    = py_installation.get_install_dir() / 'gridfire'
  gridfire_includedir = gridfire_pkg_dir / 'include'
  gridfire_libdir     = gridfire_pkg_dir / 'lib'
  gridfire_pcdir      = gridfire_libdir / 'pkgconfig'
  gridfire_datadir    = gridfire_pkg_dir / 'share'
else
  gridfire_includedir = get_option('includedir')
  gridfire_libdir     = get_option('libdir')
  gridfire_pcdir      = get_option('libdir') / 'pkgconfig'
  gridfire_datadir    = get_option('datadir') / 'gridfire'
endif
# The Fortran interface is installed as source: a compiled .mod only works with the gfortran version that wrote it
gridfire_fortran_srcdir = gridfire_datadir / 'fortran'

gridfire_vendor_includedir = gridfire_includedir / 'gridfire' / 'vendor'

# Where dependency subprojects are unpacked, as [sources, build] pairs to try in order. As the top-level project
# that is GridFire's own subprojects/ directory; as a subproject of another project, Meson may unpack them in the
# parent's top-level subprojects/ instead. Used to find headers to vendor (see gridfire_subproject_dir_of below).
fs = import('fs')
gridfire_subproject_roots = [
    [meson.project_source_root() / 'subprojects', meson.project_build_root() / 'subprojects'],
    [meson.global_source_root() / 'subprojects', meson.global_build_root() / 'subprojects'],
]

gridfire_fourdst_vendor_includedir = gridfire_includedir / 'fourdst' / 'vendor'

vendor_py = import('python').find_installation('python3')
copytree_cmd = [
  vendor_py, '-c',
  '''import sys, os, stat, shutil
src, dst = sys.argv[1], sys.argv[2]
shutil.rmtree(dst, ignore_errors=True)
shutil.copytree(src, dst)
for root, dirs, files in os.walk(dst):
    os.chmod(root, 0o755)
    for f in files:
        p = os.path.join(root, f)
        os.chmod(p, os.stat(p).st_mode | 0o644)
''',
]
copyfiles_cmd = [
  vendor_py, '-c',
  'import sys, shutil; [shutil.copy2(p, sys.argv[-1]) for p in sys.argv[1:-1]]',
]

subdir('fourdst')
subdir('sundials')

subdir('cppad')
subdir('xxHash')
subdir('eigen')

subdir('json')

subdir('CLI11')

if get_option('use_mimalloc')
  subdir('mimalloc')
endif
