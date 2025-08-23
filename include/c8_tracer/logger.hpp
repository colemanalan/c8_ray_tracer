#pragma once
#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>
#include <sstream>

namespace c8_tracer
{

    enum class LogLevel
    {
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

        void log(LogLevel msg_level, const std::string &message)
        {
            if (msg_level < level_)
                return;

            std::cout << "[" << timestamp() << "] "
                      << "[" << level_to_string(msg_level) << "] "
                      << message << std::endl;
        }

        void debug(const std::string &msg) { log(LogLevel::DEBUG, msg); }
        void info(const std::string &msg) { log(LogLevel::INFO, msg); }
        void warning(const std::string &msg) { log(LogLevel::WARNING, msg); }
        void error(const std::string &msg) { log(LogLevel::ERROR, msg); }

    private:
        LogLevel level_;

        std::string timestamp() const
        {
            auto now = std::chrono::system_clock::now();
            auto time = std::chrono::system_clock::to_time_t(now);
            std::stringstream ss;
            ss << std::put_time(std::localtime(&time), "%Y-%m-%d %H:%M:%S");
            return ss.str();
        }

        std::string level_to_string(LogLevel level) const
        {
            switch (level)
            {
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
}