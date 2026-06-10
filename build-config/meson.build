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
else
  gridfire_includedir = get_option('includedir')
  gridfire_libdir     = get_option('libdir')
  gridfire_pcdir      = get_option('libdir') / 'pkgconfig'
endif

gridfire_vendor_includedir = gridfire_includedir / 'gridfire' / 'vendor'

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
