from .c8_tracer_ext import *  # import compiled lib

# by default, turn this off because the low-level logger is quite loud even on INFO
c8_tracer_ext.logging.logger_tracer.set_level(c8_tracer_ext.logging.LogLevel.ERROR)
