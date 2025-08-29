#include "log.hpp"

// EXAMPLE :
// stoat::LOG_INFO("Program started");
// stoat::LOG_DEBUG("Loaded " + std::to_string(node_count) + " nodes");
// stoat::LOG_WARN("Using fallback parameter");
// stoat::LOG_ERROR("Cannot open file");
// stoat::LOG_FATAL("File not found: " + filename);
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
        std::ostream& out = (level == LogLevel::Error) ? std::cerr : std::cout;
        out << levelToString(level) << message << std::endl;
    }
}

void Logger::log_assert(LogLevel level, bool assertion, const std::string& message) {
    if (assertion ) {
        return;
    }
    switch (level) {
        // TODO: There isn't a log level that will stop the run so use fatal if the level is error
        case LogLevel::Error: fatal(message);
        case LogLevel::Warning: warn(message);
        case LogLevel::Info: info(message);
        case LogLevel::Debug: debug(message);
        case LogLevel::Trace: trace(message);
        default: error( "unknown LogLevel " + levelToString(level));
    }
}

void Logger::debug(const std::string& msg) { log(LogLevel::Debug, msg); }
void Logger::info(const std::string& msg)  { log(LogLevel::Info, msg); }
void Logger::warn(const std::string& msg)  { log(LogLevel::Warning, msg); }
void Logger::error(const std::string& msg) { log(LogLevel::Error, msg); }
void Logger::trace(const std::string& msg) { log(LogLevel::Trace, msg); }


void Logger::fatal(const std::string& msg) {
    log(LogLevel::Error, msg);
    throw std::runtime_error("Fatal error. Exiting.");
    // log(LogLevel::Error, "Fatal error. Exiting.");
    // std::exit(EXIT_FAILURE);
}

// Do the same thing with stringstreams
void Logger::log(LogLevel level, const std::stringstream& message) { log(level, message.str()); }
void Logger::debug(const std::stringstream& msg) { debug(msg.str()); }
void Logger::info(const std::stringstream& msg)  { info(msg.str()); }
void Logger::warn(const std::stringstream& msg)  { warn(msg.str()); }
void Logger::error(const std::stringstream& msg) { error(msg.str()); }
void Logger::trace(const std::stringstream& msg) { trace(msg.str()); }
void Logger::fatal(const std::stringstream& msg) { fatal(msg.str()); }

} // end namespace
