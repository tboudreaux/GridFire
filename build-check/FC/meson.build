if get_option('build_fortran')
  add_languages('fortran', native: true)
  message('Found FORTRAN compiler: ' + meson.get_compiler('fortran').get_id())
  message('Fortran standard set to: ' + get_option('fortran_std'))
  message('Building fortran module (gridfire_mod.mod)')
  fc = meson.get_compiler('fortran')
  if not get_option('unsafe_fortran')
    if fc.get_id() != 'gcc'
      error('The only supported fortran compiler for GridFire is gfortran (version >= 14.0), found ' + fc  + '. GridFire has not been tested with any other compilers. You can disable this check with the -Dunsafe-fortran=true flag to try other compilers')
    endif
  endif
  if (fc.version().version_compare('<14.0'))
      error('gfortran version must be at least 14.0, found ' + fc.version())
  endif
endif
