#include "log.hpp"

// EXAMPLE :
// stoat::LOG_SILENT("something only present in log file");
// stoat::LOG_INFO("Program started");
// stoat::LOG_DEBUG("Loaded " + std::to_string(node_count) + " nodes");
// stoat::LOG_WARN("Using fallback parameter");
// stoat::LOG_ERROR("Cannot open file");
// stoat::LOG_TRACE("Detailed trace info...");

namespace stoat {

Logger& Logger::instance() {
    static Logger _instance;
    return _instance;
}

void Logger::setLevel(LogLevel level) {
    logLevel = level;
}

void Logger::log(LogLevel level, const std::string& message) {
    if (level <= logLevel) {
        std::lock_guard<std::mutex> lock(mutex);
        const std::string formatted = levelToString(level) + message;

        std::ostream& out = (level == LogLevel::Error) ? std::cerr : std::cout;
        out << formatted << std::endl;

        if (fileLoggingEnabled && logFile.is_open()) {
            logFile << formatted << std::endl;
        }
    }
}

// Inside the Logger class (public section)
void Logger::silente_log(const std::string& message) {
    std::lock_guard<std::mutex> lock(mutex);
    if (fileLoggingEnabled && logFile.is_open()) {
        logFile << message << std::endl;  // INFO level by default
    }
}

void Logger::log_assert(LogLevel level, bool assertion, const std::string& message) {
    if (assertion ) {
        return;
    }

    switch (level) {
        case LogLevel::Error: error(message);
        case LogLevel::Warning: warn(message);
        case LogLevel::Info: info(message);
        case LogLevel::Debug: debug(message);
        case LogLevel::Trace: trace(message);
        default: error("Unknown LogLevel " + levelToString(level));
    }
}

void Logger::setLogFile(const std::string& filename) {
    std::lock_guard<std::mutex> lock(mutex);

    logFile.open(filename, std::ios::out | std::ios::trunc);
    if (!logFile.is_open()) {
        std::cerr << "Logger Error: Failed to open log file: " << filename << std::endl;
        return;
    }

    fileLoggingEnabled = true;
}

void Logger::debug(const std::string& msg) { log(LogLevel::Debug, msg); }
void Logger::info(const std::string& msg)  { log(LogLevel::Info, msg); }
void Logger::warn(const std::string& msg)  { log(LogLevel::Warning, msg); }
void Logger::error(const std::string& msg) { log(LogLevel::Error, msg); }
void Logger::trace(const std::string& msg) { log(LogLevel::Trace, msg); }
void Logger::silente(const std::string& msg) { silente_log(msg); }

// Do the same thing with stringstreams
void Logger::log(LogLevel level, const std::stringstream& message) { log(level, message.str()); }
void Logger::debug(const std::stringstream& msg) { debug(msg.str()); }
void Logger::info(const std::stringstream& msg)  { info(msg.str()); }
void Logger::warn(const std::stringstream& msg)  { warn(msg.str()); }
void Logger::error(const std::stringstream& msg) { error(msg.str()); }
void Logger::trace(const std::stringstream& msg) { trace(msg.str()); }
void Logger::silente(const std::stringstream& msg) { silente(msg.str()); }

} // end namespace
