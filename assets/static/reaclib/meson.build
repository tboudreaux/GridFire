required_headers = [
    'reaclib.h',
]

foreach h : required_headers
    if not cpp.has_header(h, include_directories: include_directories('include'))
        error('SERiF requires the header file ' + h + ' to be present in the assets/static/reaclib/include directory.')
    endif
endforeach
reaclib_reactions_dep = declare_dependency(
    include_directories: include_directories('include'),
)
message('✅ SERiF reaclib_reactions dependency declared')
