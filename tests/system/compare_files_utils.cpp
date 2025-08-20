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
#include <regex>

namespace fs = std::filesystem;
using namespace std;

// === UTILITY FUNCTIONS ===

void clean_output_dir(const std::string& output_dir) {
    if (fs::exists(output_dir))
        fs::remove_all(output_dir);
    fs::create_directory(output_dir);
}

std::tuple<int, int> parse_header_eqtl(
    const std::string& header_line, 
    const std::string& file_name) {

    std::istringstream ss(header_line);
    std::string token;
    size_t index = 0;
    int snarl_column_index = -1;
    int gene_column_index = -1;

    while (std::getline(ss, token, '\t')) {

        if (token == "SNARL") {
            snarl_column_index = index;
        } else if (token == "GENE") {
            gene_column_index = index;
        }

        if (snarl_column_index != -1 &&
            gene_column_index != -1) {
            break;
        }

        ++index;
    }

    if (snarl_column_index == -1 || gene_column_index == -1) {
        std::cerr << "SNARL and/or GENE column not found in header of file: " << file_name << std::endl;
        std::exit(1);
    }

    return {snarl_column_index, gene_column_index};
}

void process_tsv_line_eqtl(const std::string& line,
    std::unordered_map<std::string, std::string>& map,
    const int& snarl_column, const int& gene_column,
    const std::string& file_name) {
    
    std::istringstream ss(line);
    std::string token;
    std::vector<std::string> columns;

    while (std::getline(ss, token, '\t')) {
        columns.push_back(token);
    }

    if (snarl_column >= columns.size() || gene_column >= columns.size()) {
        std::cerr << "Invalid line (too few columns) in file: " << file_name
                  << "\nLine: " << line << std::endl;
        std::exit(1);
    }

    const std::string& snarl_id = columns[snarl_column];
    const std::string& gene_id  = columns[gene_column];

    const std::string composite_key = snarl_id + "_" + gene_id;
    map[composite_key] = line;
}

int parse_header(const std::string& header_line, 
    const std::string& file_name) {

    std::istringstream ss(header_line);
    std::string token;
    int snarl_column_index = 0;

    while (std::getline(ss, token, '\t')) {
        if (token == "SNARL") {
            return snarl_column_index;
        }
        ++snarl_column_index;
    }

    std::cerr << "SNARL column not found in header of file: " << file_name << std::endl;
    std::exit(1);
}

void process_tsv_line(const std::string& line,
    std::unordered_map<std::string, std::string>& map,
    const int& snarl_column,
    const std::string& file_name) {
    
    std::istringstream ss(line);
    std::string token;
    std::vector<std::string> columns;

    while (std::getline(ss, token, '\t')) {
        columns.push_back(token);
    }

    if (snarl_column >= columns.size()) {
        std::cerr << "Invalid line (too few columns) in file: " << file_name
                  << "\nLine: " << line << std::endl;
        std::exit(1);
    }

    const std::string& snarl_id = columns[snarl_column];
    map[snarl_id] = line;
}

void load_tsv_file_eqtl(const std::string& path,
    std::unordered_map<std::string, std::string>& map) {

    std::ifstream infile(path);
    if (!infile.is_open()) {
        std::cerr << "Failed to open file: " << path << std::endl;
        std::exit(1);
    }

    std::string line;

    // --- Parse header ---
    if (!std::getline(infile, line) || line.empty()) {
        std::cerr << "Empty file or missing header in file: " << path << std::endl;
        std::exit(1);
    }

    auto [snarl_column, gene_column] = parse_header_eqtl(line, path);

    // --- Process data lines ---
    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        process_tsv_line_eqtl(line, map, snarl_column, gene_column, path);
    }
}

void load_tsv_file(const std::string& path,
    std::unordered_map<std::string, std::string>& map) {

    std::ifstream infile(path);
    if (!infile.is_open()) {
        std::cerr << "Failed to open file: " << path << std::endl;
        std::exit(1);
    }

    std::string line;
    int snarl_column = -1;

    // --- Parse header ---
    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        snarl_column = parse_header(line, path);
        break;
    }

    // --- Process data lines ---
    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        process_tsv_line(line, map, snarl_column, path);
    }
}

void load_tsv_file(const std::vector<std::string>& lines, 
    std::unordered_map<std::string, std::string>& map) {

    int snarl_column = -1;
    bool found_header = false;

    for (const std::string& line : lines) {
        if (line.empty()) continue;

        if (!found_header) {
            snarl_column = parse_header(line, "supposed truth vector");
            found_header = true;
            continue; // skip header line
        }

        process_tsv_line(line, map, snarl_column, "supposed truth vector");
    }
}

bool files_equal_eqtl(const std::string& file1, const std::string& file2) {
    std::unordered_map<std::string, std::string> map1, map2;

    load_tsv_file_eqtl(file1, map1);
    load_tsv_file_eqtl(file2, map2);

    return compare_file(map1, map2);
}

bool files_equal(const std::string& file1, const std::string& file2) {
    std::unordered_map<std::string, std::string> map1, map2;

    load_tsv_file(file1, map1);
    load_tsv_file(file2, map2);

    return compare_file(map1, map2);
}

bool files_equal(const std::string& file, const std::vector<std::string>& lines) {
    std::unordered_map<std::string, std::string> map1, map2;

    load_tsv_file(file, map1);
    load_tsv_file(lines, map2);

    return compare_file(map1, map2);
}

bool compare_file(const std::unordered_map<std::string, std::string>& map1, 
    const std::unordered_map<std::string, std::string>& map2) {

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

bool is_valid_fasta(const std::string& file) {

    std::ifstream infile(file);
    std::string line;
    while (std::getline(infile, line)) {
        if (line[0] == '>') {
            continue;
        }
        // Check that this is a valid sequence line
        if (!std::regex_match(line, std::regex("[ACGTN]*"))) {
            cerr << "Invalid FASTA line: " << line << endl;
            return false;
        }
        if (infile.peek() != '>' && infile.peek() != EOF) {
            if (line.size() > 80) {
                cerr << "FASTA sequence line longer than 80 characters" << endl;
                return false;
            }
        }
    }
    return true;
}

bool fasta_equal(const std::string& file, 
    const std::vector<std::tuple<size_t, std::string, std::string>>& fasta_records) {

    std::unordered_map<std::string, std::pair<size_t, std::string>> header_to_sequence;
    std::unordered_map<size_t, std::vector<std::string>> set_to_header;
    size_t record_count = 0;
    for (auto& record : fasta_records) {
        size_t set_id = std::get<0>(record);
        record_count = std::max(record_count, set_id);
        header_to_sequence.emplace(std::string(std::get<1>(record)), std::make_pair(set_id, std::get<2>(record)));

        if (!set_to_header.count(set_id)) {
            set_to_header.emplace(set_id, std::vector<std::string>());
        }
        set_to_header.at(set_id).emplace_back(std::get<1>(record));
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
        actual_record_count ++;
    }
    infile.close();
    for (size_t i = 1 ; i < record_count+1 ; i++) {
        if (!has_record[i]) {
            cerr << i << endl;
            cerr << "Fasta output should have one of the following records" << endl;
            for (const auto& header : set_to_header.at(i)) {
                cerr << "\t" << header << endl;
            }
            return false;
        }
    }

    if (actual_record_count != record_count) {
        cerr << "Fasta output has " << actual_record_count << " records, should have " << record_count << endl;
        return false;
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

        const std::string filename = expected_file.filename().string();

        if (filename.find("eqtl") != std::string::npos) {
            if (!files_equal_eqtl(expected_file, output_file)) {
                std::cerr << "Mismatch in eQTL file: " << filename << "\n";
                return false;
            }
        } else {
            if (!files_equal(expected_file, output_file)) {
                std::cerr << "Mismatch in file: " << filename << "\n";
                return false;
            }
        }
    }
    return true;
}
