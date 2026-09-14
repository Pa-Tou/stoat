#ifndef LOG_HPP
#define LOG_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>

namespace stoat {

#define LOG_ERROR(msg)   Logger::instance().error(msg)
#define LOG_WARN(msg, warn_code)    Logger::instance().warn(msg, warn_code)
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

    // Warning log that checks if too many warning of the same type have been printed to the terminal
    // if it's the case just write them in the log without printing them
    void log_warning(LogLevel level, const std::string& message, const std::string& message_code);

    /// Is the logger at least at this level of verbosity?
    bool at_level(const LogLevel& level) const { return logLevel >= level; }

    /// Check an assertion and if it is false, print the message to the appropriate log level
    /// This should only really be used with at_level() so that the assertion check doesn't happen all the time
    void log_assert(LogLevel level, bool assertion, const std::string& message);

    void debug(const std::string& msg);
    void info(const std::string& msg);
    void warn(const std::string& msg, const std::string& warn_code);
    void error(const std::string& msg);
    void trace(const std::string& msg);
    void silente(const std::string& msg);

    void log(LogLevel level, const std::stringstream& message);
    void silente_log(const std::string& message);

    void debug(const std::stringstream& msg);
    void info(const std::stringstream& msg);
    void warn(const std::stringstream& msg, const std::string& warn_code);
    void error(const std::stringstream& msg);
    void trace(const std::stringstream& msg);
    void silente(const std::stringstream& msg);

    void setLogFile(const std::string& filename);

private:
    LogLevel logLevel = LogLevel::Info;
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

    // Map the warning to the count of times it has been printed.
    // TODO: I don't like this
    std::unordered_map<std::string, size_t> warning_to_count;
};

} //end stoat namespace

#endif
