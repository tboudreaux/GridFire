llevel = get_option('log_level')

logbase='QUILL_COMPILE_ACTIVE_LOG_LEVEL_'

if (llevel == 'traceL3')
    message('Setting log level to TRACE_L3')
    log_argument = logbase + 'TRACE_L3'
elif (llevel == 'traceL2')
    message('Setting log level to TRACE_L2')
    log_argument = logbase + 'TRACE_L2'
elif (llevel == 'traceL1')
    message('Setting log level to TRACE_L1')
    log_argument = logbase + 'TRACE_L1'
elif (llevel == 'debug')
    message('Setting log level to DEBUG')
    log_argument = logbase + 'DEBUG'
elif (llevel == 'info')
    message('Setting log level to INFO')
    log_argument = logbase + 'INFO'
elif (llevel == 'warning')
    message('Setting log level to WARNING')
    log_argument = logbase + 'WARNING'
elif (llevel == 'error')
    message('Setting log level to ERROR')
    log_argument = logbase + 'ERROR'
elif (llevel == 'critical')
    message('Setting log level to CRITICAL')
    log_argument = logbase + 'CRITICAL'
endif

log_argument = '-DQUILL_COMPILE_ACTIVE_LOG_LEVEL=' + log_argument
gridfire_args += [log_argument]
