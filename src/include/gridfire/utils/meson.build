python_exe = import('python').find_installation('python3')

version_parser = '''
import sys, re
ver = sys.argv[1]
if ver.startswith("v"): ver = ver[1:]
m = re.match(r"^(\d+)\.(\d+)\.(\d+)(.*)$", ver)
if m:
    print(f"{m.group(1)};{m.group(2)};{m.group(3)};{m.group(4)}")
else:
    print("0;0;0;unknown")
'''

ver_res = run_command(python_exe, '-c', version_parser, meson.project_version(), check: true)
ver_parts = ver_res.stdout().strip().split(';')

conf_data = configuration_data()
conf_data.set('GF_VERSION_MAJOR', ver_parts[0])
conf_data.set('GF_VERSION_MINOR', ver_parts[1])
conf_data.set('GF_VERSION_PATCH', ver_parts[2])
conf_data.set('GF_VERSION_TAG',   ver_parts[3])

message('Configuring include/utils/config.h with version ' + meson.project_version())
message('  Major: ' + ver_parts[0])
message('  Minor: ' + ver_parts[1])
message('  Patch: ' + ver_parts[2])
message('  Tag:   ' + ver_parts[3])

configure_file(
    input : 'config.h.in',
    output : 'config.h',
    configuration : conf_data,
    install : true,
    install_dir : gridfire_includedir / 'gridfire/utils'
)
