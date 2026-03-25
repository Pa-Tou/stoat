#include <iostream>
#include <string>
#include <algorithm>

#include "banner.hpp"

namespace stoat {

void print_ascii_banner(const std::string& version) {
        std::cerr << R"(-------------------------------------------------
                         (\ /)
 ___|‾|_ ___   __ _|‾|_  (. .)
/ __| __/ _ \ / _\ | __|  >o<\
\__ \ || (_) | (_| | |_    \  `--------.,=-.,
|___/\__\___/ \__ _|\__|    `,,________,,)  `--.
                            ,,        ,,
Stoat )" << version << R"( - Pangenome GWAS toolkit
Run 'stoat --help' for full documentation.
-------------------------------------------------
)";
}

// Modern minimal banner :
// Stoat v0.0.4
void print_banner(const std::string& version) {
    std::cerr << "\nStoat v" << version << "\n";
}

} // namespace stoat
