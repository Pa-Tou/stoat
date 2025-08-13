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
using namespace std;

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

bool fasta_equal(const std::string& file, const std::vector<std::tuple<size_t, std::string, std::string>>& fasta_records) {
    std::unordered_map<std::string, std::pair<size_t, std::string>> header_to_sequence;
    size_t record_count = 0;
    for (auto& record : fasta_records) {
        record_count = std::max(record_count, std::get<0>(record));
        header_to_sequence.emplace(std::string(std::get<1>(record)), std::make_pair(std::get<0>(record), std::get<2>(record)));
    }

    //Check that each equivalence class is represented in the file 
    std::vector<bool> has_record(record_count, false);
    size_t actual_record_count = 0;

    std::ifstream infile(file);
    std::string line;

    while (std::getline(infile, line)) {
        std::string header = line;
        std::string seq = "";
        while (infile.peek() != '>' && infile.peek() != EOF) {
            std::getline(infile, line);
            seq += line;
        }
        if (!header_to_sequence.count(header)) {
            std::cerr << "FASTA output contains unknown header: " << header << std::endl;
            return false;
        }
        const std::pair<size_t, std::string>& truth = header_to_sequence.at(header);
        if (seq != truth.second) {
            std::cerr << "FASTA output with header" << std::endl << header << std::endl;
            std::cerr << "contains different sequence" << std::endl;
            std::cerr << "FASTA: " << seq << std::endl;
            std::cerr << "Truth: " << truth.second << std::endl;

            return false;
        }

        if (has_record[truth.first]) {
            std::cerr << "FASTA output has duplicate header " << header << std::endl;
        }

        has_record[truth.first] = true;
    }
    infile.close();
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
