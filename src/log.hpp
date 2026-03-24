#ifndef LOG_HPP
#define LOG_HPP

#include <iostream>
#include <fstream>
#include <mutex>
#include <sstream>

namespace stoat {

#define LOG_ERROR(msg)   Logger::instance().error(msg)
#define LOG_WARN(msg, count_type_warning)    Logger::instance().warn(msg, count_type_warning)
#define LOG_INFO(msg)    Logger::instance().info(msg)
#define LOG_DEBUG(msg)   Logger::instance().debug(msg)
#define LOG_TRACE(msg)   Logger::instance().trace(msg)
#define LOG_SILENTE(msg) Logger::instance().silente(msg)

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

    // Warning log that check if too many warning of the same type are printing in the terminal
    // if it's the case just write them in the log without print them
    void log_warning(LogLevel level, const std::string& message, const size_t& count_type_warning);

    /// Is the logger at least at this level of verbosity?
    bool at_level(const LogLevel& level) const { return logLevel >= level; }

    /// Check an assertion and if it is false, print the message to the appropriate log level
    /// This should only really be used with at_level() so that the assertion check doesn't happen all the time
    void log_assert(LogLevel level, bool assertion, const std::string& message);

    void debug(const std::string& msg);
    void info(const std::string& msg);
    void warn(const std::string& msg, const size_t& count_warn_type);
    void error(const std::string& msg);
    void trace(const std::string& msg);
    void silente(const std::string& msg);

    void log(LogLevel level, const std::stringstream& message);
    void silente_log(const std::string& message);

    void debug(const std::stringstream& msg);
    void info(const std::stringstream& msg);
    void warn(const std::stringstream& msg, const size_t& count_warn_type);
    void error(const std::stringstream& msg);
    void trace(const std::stringstream& msg);
    void silente(const std::stringstream& msg);

    void setLogFile(const std::string& filename);

private:
    LogLevel logLevel = LogLevel::Info;
    std::mutex mutex;
    const size_t warning_count_threshold = 4; // 5 warning will be printed

    Logger() = default;

    std::ofstream logFile;
    bool fileLoggingEnabled = false;

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
