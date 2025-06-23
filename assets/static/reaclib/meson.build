required_headers = [
    'gridfire/reaclib.h',
    'gridfire/reactions.h',
]

foreach h : required_headers
    if not cpp.has_header(h, include_directories: include_directories('include'))
        error('GridFire requires the header file ' + h + ' to be present in the assets/static/reaclib/include/gridfire directory.')
    endif
endforeach
reaclib_reactions_dep = declare_dependency(
    include_directories: include_directories('include'),
)
message('✅ GridFire reaclib_reactions dependency declared')

to_install_headers = [
    'include/gridfire/reaclib.h',
    'include/gridfire/reactions.h',
]
install_headers(to_install_headers, subdir: 'gridfire/gridfire')
