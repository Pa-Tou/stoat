#ifndef BANNER_HPP
#define BANNER_HPP

#include <iostream>
#include "log.hpp"

namespace stoat {
    void print_ascii_banner(const std::string &version);
    void print_banner(const std::string &version);
}

#endif //end banner hpp
