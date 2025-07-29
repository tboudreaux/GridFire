boost_dep = dependency('boost', version: '>=1.8.0')

if not boost_dep.found()
    message('❌ Boost library not found. Note that Boost is a system dependency and must be installed seperatly.')
    message('GridFire includes a Boost auto installation script, attempting to run it now (this will require some user input)...')
    run_res = run_command('./install.sh', check: true)
    boost_dep = dependency('boost', version: '>=1.8.0', required: true)
endif


