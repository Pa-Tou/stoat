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

// Boxed banner:

// ┌──────────────────────────────────────┐
// | Stoat VERSION                        |
// | Pangenome GWAS toolkit               |
// |                                      |
// | Type --help for usage information    |
// └──────────────────────────────────────┘
void print_banner_v1(const std::string& version) {
    const std::string tool_name = "Stoat " + version;
    const std::string tagline = "Pangenome GWAS toolkit";
    const std::string hint = "Type --help for usage information";

    // Determine the content width (longest string)
    size_t content_width = std::max({tool_name.size(), tagline.size(), hint.size()});
    size_t width = content_width + 2; // 1 space padding on each side

    // Top border
    std::cerr << "┌" << std::string(width, '-') << "┐\n";

    // Helper to print a line with proper padding
    auto print_line = [width](const std::string& text) {
        size_t padding = width - text.size();
        std::cerr << "| " << text << std::string(padding - 1, ' ') << "|\n";
    };

    print_line(tool_name);
    print_line(tagline);
    print_line("");       // empty line
    print_line(hint);

    // Bottom border
    std::cerr << "└" << std::string(width, '-') << "┘\n";
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
