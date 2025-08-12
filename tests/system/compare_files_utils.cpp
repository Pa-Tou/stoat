#include "compare_files_utils.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <string>

namespace fs = std::filesystem;

// === UTILITY FUNCTIONS ===

void clean_output_dir(const std::string& output_dir) {
    if (fs::exists(output_dir))
        fs::remove_all(output_dir);
    fs::create_directory(output_dir);
}

bool files_equal(const std::string& file1, const std::string& file2) {
    std::unordered_map<std::string, std::string> map1, map2;

    auto load_file = [](const std::string& path, std::unordered_map<std::string, std::string>& map) {
        std::ifstream infile(path);
        std::string line;
        int snarl_column = -1;
        bool header_processed = false;

        while (std::getline(infile, line)) {
            if (line.empty()) continue;

            std::istringstream ss(line);
            std::string token;
            std::vector<std::string> columns;

            while (std::getline(ss, token, '\t')) {
                columns.push_back(token);
            }

            // First non-empty line is the header
            if (!header_processed) {
                for (size_t i = 0; i < columns.size(); ++i) {
                    if (columns[i] == "SNARL") {
                        snarl_column = static_cast<int>(i);
                        break;
                    }
                }

                if (snarl_column == -1) {
                    std::cerr << "SNARL column not found in header of file: " << path << std::endl;
                    exit(1);
                }

                header_processed = true;
                continue; // skip header line
            }

            if (snarl_column >= columns.size()) {
                std::cerr << "Invalid line (too few columns) in file: " << path << "\nLine: " << line << std::endl;
                exit(1);
            }

            std::string snarl_id = columns[snarl_column];
            map[snarl_id] = line;
        }
    };

    load_file(file1, map1);
    load_file(file2, map2);

    // Compare file1 against file2
    for (const auto& [snarl, line1] : map1) {
        auto it = map2.find(snarl);
        if (it == map2.end()) {
            std::cerr << "Missing SNARL in file2: " << snarl << std::endl;
            return false;
        } else if (line1 != it->second) {
            std::cerr << "Mismatch for SNARL " << snarl << ":\n"
                      << "File1: " << line1 << "\n"
                      << "File2: " << it->second << "\n";
            return false;
        }
    }

    // Compare file2 against file1
    for (const auto& [snarl, _] : map2) {
        if (map1.find(snarl) == map1.end()) {
            std::cerr << "Missing SNARL in file1: " << snarl << std::endl;
            return false;
        }
    }

    return true;
}

bool compare_output_dirs(const std::string& output_dir, const std::string& expected_dir) {
    for (const auto& file : fs::directory_iterator(expected_dir)) {
        auto expected_file = file.path();
        auto output_file = fs::path(output_dir) / expected_file.filename();

        if (!fs::exists(output_file)) {
            std::cerr << "Missing output file: " << output_file << "\n";
            return false;
        }
        if (!files_equal(expected_file, output_file)) {
            std::cerr << "Mismatch in file: " << expected_file.filename() << "\n";
            return false;
        }
    }
    return true;
}
