#pragma once
#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>

namespace c8_tracer
{

  enum class LogLevel
  {
    ALL,
    TRACE,
    DEBUG,
    INFO,
    WARNING,
    ERROR
  };

  class Logger
  {
  public:
    Logger(LogLevel level = LogLevel::INFO) : level_(level) {}

    void set_level(LogLevel level)
    {
      level_ = level;
    }

    void log(LogLevel msg_level, const std::string &message, const std::string &func_name)
    {
      if (msg_level < level_)
        return;

      std::cout << "[" << level_to_string(msg_level) << "] "
                << (func_name.empty() ? "" : "[" + func_name + "] ")
                << message << std::endl;
    }

    // These should only be used as pybdingings
    void all(const std::string &msg) { log(LogLevel::ALL, msg, ""); }
    void trace(const std::string &msg) { log(LogLevel::TRACE, msg, ""); }
    void debug(const std::string &msg) { log(LogLevel::DEBUG, msg, ""); }
    void info(const std::string &msg) { log(LogLevel::INFO, msg, ""); }
    void warning(const std::string &msg) { log(LogLevel::WARNING, msg, ""); }
    void error(const std::string &msg) { log(LogLevel::ERROR, msg, ""); }

  private:
    LogLevel level_;

    std::string level_to_string(LogLevel level) const
    {
      switch (level)
      {
      case LogLevel::ALL:
        return "ALL";
      case LogLevel::TRACE:
        return "TRACE";
      case LogLevel::DEBUG:
        return "DEBUG";
      case LogLevel::INFO:
        return "INFO";
      case LogLevel::WARNING:
        return "WARNING";
      case LogLevel::ERROR:
        return "ERROR";
      default:
        return "UNKNOWN";
      }
    }
  };

  // Global logger instance (optional)
  inline Logger logger;
  inline Logger logger_tracer;

#define LOG_ERROR(msg) c8_tracer::logger.log(LogLevel::ERROR, msg, __func__)
#define LOG_WARNING(msg) c8_tracer::logger.log(LogLevel::WARNING, msg, __func__)
#define LOG_INFO(msg) c8_tracer::logger.log(LogLevel::INFO, msg, __func__)
#define LOG_DEBUG(msg) c8_tracer::logger.log(LogLevel::DEBUG, msg, __func__)
#define LOG_TRACE(msg) c8_tracer::logger.log(LogLevel::TRACE, msg, __func__)
#define LOG_ALL(msg) c8_tracer::logger.log(LogLevel::ALL, msg, __func__)


#define TRACER_LOG_ERROR(msg) c8_tracer::logger_tracer.log(LogLevel::ERROR, msg, __func__)
#define TRACER_LOG_WARNING(msg) c8_tracer::logger_tracer.log(LogLevel::WARNING, msg, __func__)
#define TRACER_LOG_INFO(msg) c8_tracer::logger_tracer.log(LogLevel::INFO, msg, __func__)
#define TRACER_LOG_DEBUG(msg) c8_tracer::logger_tracer.log(LogLevel::DEBUG, msg, __func__)
#define TRACER_LOG_TRACE(msg) c8_tracer::logger_tracer.log(LogLevel::TRACE, msg, __func__)
#define TRACER_LOG_ALL(msg) c8_tracer::logger_tracer.log(LogLevel::ALL, msg, __func__)


}