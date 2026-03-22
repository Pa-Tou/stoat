#include <iostream>
#include <string>

#include "banner.hpp"

namespace stoat {

void print_full_banner() {
    std::cerr << R"(--------------------------------------------------
                         (\ /)
 ___|‾|_ ___   __ _|‾|_  (. .)
/ __| __/ _ \ / _\ | __|  >o<\
\__ \ || (_) | (_| | |_    \  `--------.,=-.,_
|___/\__\___/ \__ _|\__|    `,,________,,)  `---.
                             ""        ""
--------------------------------------------------
)";
}

void print_small_banner(const std::string& version) {
    std::cerr << "\n--------------------------------------------------\n"
                 "Stoat version " << version << "\n"
                 "--------------------------------------------------\n";
}

} // namespace stoat
