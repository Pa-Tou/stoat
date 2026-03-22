#include <iostream>
#include <string>
#include <algorithm>

#include "banner.hpp"

namespace stoat {

void print_ascii_banner(const std::string& version) {
        std::cerr << 
    R"(-------------------------------------------------
                             (\ /)
     ___|‾|_ ___   __ _|‾|_  (. .)
    / __| __/ _ \ / _\ | __|  >o<\
    \__ \ || (_) | (_| | |_    \  `--------.,=-.,_
    |___/\__\___/ \__ _|\__|    `,,________,,)  `---.
                                ,,        ,,
    Stoat )" << version << R"( - Pangenome GWAS toolkit
    Run 'stoat --help' for full documentation.
    -------------------------------------------------
    )";
}

// Boxed banner:

// ┌──────────────────────────────────────┐
// │ Stoat VERSION                        │
// │ Pangenome GWAS toolkit               │
// │                                      │
// │ Type --help for usage information    │
// └──────────────────────────────────────┘
void print_banner_v1(const std::string& version) {
    const std::string tool_name = "Stoat " + version;
    const std::string tagline = "Pangenome GWAS toolkit";
    const std::string hint = "Type --help for usage information";

    // Determine box width based on longest line
    size_t width = std::max({tool_name.size(), tagline.size(), hint.size()}) + 4; // padding

    // Top border
    std::cerr << "\n┌" << std::string(width, '─') << "┐\n";

    // Content lines
    std::cerr << "│ " << tool_name
              << std::string(width - tool_name.size() - 2, ' ')
              << "│\n";
    std::cerr << "│ " << tagline
              << std::string(width - tagline.size() - 2, ' ')
              << "│\n";
    std::cerr << "│ " << std::string(width - 2, ' ') << "│\n"; // empty line
    std::cerr << "│ " << hint
              << std::string(width - hint.size() - 2, ' ')
              << "│\n";

    // Bottom border
    std::cerr << "└" << std::string(width, '─') << "┘\n";
}

// Modern minimal banner :

// Stoat 0.0.4
// Pangenome GWAS toolkit
//
// Usage: stoat [options] <input>
//
// Run 'stoat --help' for full documentation.
void print_banner_v2(const std::string& version) {
    std::cerr << "\nStoat " << version << "\n";
    std::cerr << "Pangenome GWAS toolkit\n\n";
    std::cerr << "Usage: stoat [options] <input>\n\n";
    std::cerr << "Run 'stoat --help' for full documentation.\n";
}

} // namespace stoat
