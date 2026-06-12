cppc = meson.get_compiler('cpp')
fs = import('fs')

# --- minimum supported versions ------------------------------------------
gridfire_gcc_min        = '14.0'   # g++
gridfire_clang_min      = '17.0'   # upstream LLVM clang++
gridfire_appleclang_min = '16.0'   # Apple clang (Xcode 16 toolchain)

# --- identify the active toolchain ----------------------------------------
compiler_id = cppc.get_id()
compiler_version = cppc.version()

is_apple_clang = false
if compiler_id == 'clang' and cppc.get_define('__apple_build_version__') != ''
    is_apple_clang = true
endif

if is_apple_clang
    toolchain_desc = 'Apple clang (Xcode) ' + compiler_version
    required_min = gridfire_appleclang_min
elif compiler_id == 'clang'
    toolchain_desc = 'LLVM clang ' + compiler_version
    required_min = gridfire_clang_min
elif compiler_id == 'gcc'
    toolchain_desc = 'GCC ' + compiler_version
    required_min = gridfire_gcc_min
else
    toolchain_desc = compiler_id + ' ' + compiler_version
    required_min = ''
    warning('GridFire is developed and tested with GCC and Clang; compiler "'
            + compiler_id + '" is untested. Proceeding, but the capability '
            + 'checks below must pass.')
endif

message('C++ toolchain: ' + toolchain_desc)
message('C++ standard set to: ' + get_option('cpp_std'))

# --- per-compiler warning tweaks (unchanged behavior) ----------------------
if compiler_id == 'clang'
    message('disabling bitwise-instead-of-logical warnings for clang')
    add_project_arguments('-Wno-bitwise-instead-of-logical', language: 'cpp')
endif
if compiler_id == 'gcc'
    message('disabling psabi warnings for gcc')
    add_project_arguments('-Wno-psabi', language: 'cpp')
endif

# --- functional capability probes (authoritative) --------------------------
have_print_hdr  = cppc.has_header('print')
have_format_hdr = cppc.has_header('format')
println_works = cppc.compiles(
    '''#include <print>
       #include <format>
       int main() { std::println("gridfire toolchain probe {}", 23); }''',
    args: ['-std=c++23'],
    name: 'std::println usable at -std=c++23',
)

toolchain_functional = have_print_hdr and have_format_hdr and println_works
version_sufficient = required_min == '' ? true : compiler_version.version_compare('>=' + required_min)

# --- failure analysis and reporting ----------------------------------------
if toolchain_functional and not version_sufficient
    warning(toolchain_desc + ' is below the minimum GridFire tests against ('
            + required_min + '), but all C++23 capability probes passed. '
            + 'Proceeding; if you hit compiler errors deep in the build, '
            + 'upgrade your toolchain before reporting a bug.')
endif

if not toolchain_functional
    failure_detail = ''
    if not have_print_hdr
        failure_detail += '\n  * C++ standard library header <print> not found.'
    endif
    if not have_format_hdr
        failure_detail += '\n  * C++ standard library header <format> not found.'
    endif
    if have_print_hdr and have_format_hdr and not println_works
        failure_detail += '\n  * Headers exist but std::println failed to compile at -std=c++23 (incomplete library support).'
    endif
    if not version_sufficient
        failure_detail += '\n  * ' + toolchain_desc + ' is below the required minimum (' + required_min + ').'
    elif compiler_id == 'clang' and not is_apple_clang
        failure_detail += ('\n  * clang itself is new enough; on Linux clang uses the system '
            + 'libstdc++, so the GNU C++ runtime is likely too old. Install GCC >= '
            + gridfire_gcc_min + ' (clang will pick up its libstdc++), or use -Dcpp_args=-stdlib=libc++ with libc++ >= 17 installed.')
    endif

    candidate_names = []
    if compiler_id == 'gcc'
        candidate_names += ['g++-16', 'g++-15', 'g++-14', 'clang++-21', 'clang++-20', 'clang++-19', 'clang++-18', 'clang++-17']
    else
        candidate_names += ['clang++-21', 'clang++-20', 'clang++-19', 'clang++-18', 'clang++-17', 'g++-16', 'g++-15', 'g++-14']
    endif
    if host_machine.system() == 'darwin'
        candidate_names += [
            '/opt/homebrew/opt/llvm/bin/clang++',
            '/usr/local/opt/llvm/bin/clang++',
        ]
    endif

    suitable_cxx = ''
    candidates_report = ''
    foreach cand : candidate_names
        p = find_program(cand, required: false)
        if p.found()
            cand_ver = p.version()
            cand_min = cand.contains('clang') ? gridfire_clang_min : gridfire_gcc_min
            if cand_ver != 'unknown' and cand_ver.version_compare('>=' + cand_min)
                candidates_report += '\n  [OK]       ' + p.full_path() + '  (version ' + cand_ver + ')'
                if suitable_cxx == ''
                    suitable_cxx = p.full_path()
                endif
            else
                candidates_report += '\n  [too old]  ' + p.full_path() + '  (version ' + cand_ver + ')'
            endif
        endif
    endforeach

    os_help = ''
    if host_machine.system() == 'darwin'
        os_help = '''
How to get a suitable compiler on macOS:
  * Apple clang: update Xcode / Command Line Tools to Xcode 16 or newer
      softwareupdate --list   (or install from developer.apple.com)
  * Homebrew LLVM:   brew install llvm     -> /opt/homebrew/opt/llvm/bin/clang++
  * Homebrew GCC:    brew install gcc      -> g++-14 (or newer) on your PATH'''
    elif host_machine.system() == 'linux'
        distro_hint = ''
        if fs.exists('/etc/os-release')
            os_release = fs.read('/etc/os-release')
            if os_release.contains('ubuntu') or os_release.contains('debian')
                distro_hint = '''
  Detected Debian/Ubuntu:
      sudo apt install g++-14            (or newer)
      sudo apt install clang-18 libstdc++-14-dev   (clang needs a modern libstdc++ too)'''
            elif os_release.contains('fedora') or os_release.contains('rhel') or os_release.contains('centos')
                distro_hint = '''
  Detected Fedora/RHEL family:
      sudo dnf install gcc-c++           (Fedora 40+ ships GCC 14)
      sudo dnf install clang'''
            elif os_release.contains('arch')
                distro_hint = '''
  Detected Arch:
      sudo pacman -S gcc clang           (both current in the repos)'''
            endif
        endif
        os_help = '''
How to get a suitable compiler on Linux:''' + distro_hint + '''
  Generic: install GCC >= ''' + gridfire_gcc_min + ''' or LLVM clang >= ''' + gridfire_clang_min + ''' from your
  package manager; many distros package versioned binaries (g++-14, clang++-18).'''
    else
        os_help = '\nInstall GCC >= ' + gridfire_gcc_min + ' or LLVM clang >= ' + gridfire_clang_min + ' for your platform.'
    endif

    if suitable_cxx != ''
        action = ('\nA suitable compiler IS already installed. Meson cannot switch compilers '
            + 'after configuration starts, so re-run setup pointing at it:\n\n'
            + '    CXX=' + suitable_cxx + ' meson setup --wipe <builddir>\n')
    else
        action = ('\nNo suitable alternate compiler was found on this system.\n' + os_help + '\n\n'
            + 'After installing one, configure with:\n\n'
            + '    CXX=<new compiler> meson setup --wipe <builddir>\n')
    endif

    scanned = candidates_report == '' ? '\n  (none of the common versioned compiler names were found on PATH)' : candidates_report

    error('GridFire requires a C++23 toolchain (GCC >= ' + gridfire_gcc_min
          + ', LLVM clang >= ' + gridfire_clang_min + ', or Apple clang >= ' + gridfire_appleclang_min + ').\n'
          + 'Active compiler: ' + toolchain_desc + '\n'
          + '\nProblems detected:' + failure_detail + '\n'
          + '\nAlternate compilers scanned:' + scanned + '\n'
          + action)
endif


add_project_arguments('-Wno-deprecated-declarations', language: 'cpp')
add_project_arguments('-fvisibility=default', language: 'cpp')

if get_option('openmp_support')
    gridfire_args += ['-DGF_USE_OPENMP']
endif

if get_option('asan') and get_option('buildtype') != 'debug' and get_option('buildtype') != 'debugoptimized'
    error('AddressSanitizer (ASan) can only be enabled for debug or debugoptimized builds')
endif

if get_option('asan') and (get_option('buildtype') == 'debugoptimized' or get_option('buildtype') == 'debug')
    message('enabling AddressSanitizer (ASan) support')
    add_project_arguments('-fsanitize=address,undefined', language: 'cpp')
    add_project_arguments('-fno-omit-frame-pointer', language: 'cpp')

    add_project_link_arguments('-fsanitize=address,undefined', language: 'cpp')
    add_project_link_arguments('-fno-omit-frame-pointer', language: 'cpp')
endif
