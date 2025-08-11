#ifndef LOG_HPP
#define LOG_HPP

#include <iostream>
#include <fstream>
#include <mutex>

namespace stoat {

#define LOG_FATAL(msg)   Logger::instance().fatal((msg))
#define LOG_ERROR(msg)   Logger::instance().error((msg))
#define LOG_WARN(msg)    Logger::instance().warn((msg))
#define LOG_INFO(msg)    Logger::instance().info((msg))
#define LOG_DEBUG(msg)   Logger::instance().debug((msg))
#define LOG_TRACE(msg)   Logger::instance().trace((msg))

enum class LogLevel {
    Error = 0,
    Warning = 1,
    Info = 2,
    Debug = 3,
    Trace = 4
};

class Logger {
public:
    static Logger& instance();

    void setLevel(LogLevel level);
    void log(LogLevel level, const std::string& message);

    void debug(const std::string& msg);
    void info(const std::string& msg);
    void warn(const std::string& msg);
    void error(const std::string& msg);
    void fatal(const std::string& msg);  // logs error and exits
    void trace(const std::string& msg);

private:
    LogLevel logLevel = LogLevel::Info;
    std::mutex mutex;

    Logger() = default;

    std::string levelToString(LogLevel level) const {
        switch (level) {
            case LogLevel::Error: return "ERROR: ";
            case LogLevel::Warning: return "WARNING: ";
            case LogLevel::Info: return "";
            case LogLevel::Debug: return "DEBUG: ";
            case LogLevel::Trace: return "TRACE: ";
            default: return "UNKNOWN";
        }
    }
};

} //end stoat namespace

#endif