#include <iostream>
#include <string>
#include <algorithm>

#include "banner.hpp"

namespace stoat {

void print_ascii_banner(const std::string& version) {
        std::cout << R"(-------------------------------------------------
                         (\ /)
 ___|‾|_ ___   __ _|‾|_  (. .)
/ __| __/ _ \ / _\ | __|  >o<\
\__ \ || (_) | (_| | |_    \  `--------.,=-.,
|___/\__\___/ \__ _|\__|    `,,________,,)  `--.
                            ,,        ,,
-------------------------------------------------
)";
        stoat::LOG_INFO("Stoat " + version);
}

// Modern minimal banner :
// Stoat v0.0.4
void print_banner(const std::string& version) {
        stoat::LOG_INFO("Stoat " + version);
}

} // namespace stoat
