#include "log.hpp"

// EXAMPLE :
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
        std::ostream& out = (level == LogLevel::Error) ? std::cerr : std::cout;
        out << levelToString(level) << message << std::endl;
    }
}

void Logger::debug(const std::string& msg) { log(LogLevel::Debug, msg); }
void Logger::info(const std::string& msg)  { log(LogLevel::Info, msg); }
void Logger::warn(const std::string& msg)  { log(LogLevel::Warning, msg); }
void Logger::error(const std::string& msg) { log(LogLevel::Error, msg); }
void Logger::trace(const std::string& msg) { log(LogLevel::Trace, msg); }

} // end namespace