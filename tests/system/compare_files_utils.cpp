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
#include "../../src/snarl_data_collection.hpp"

namespace fs = std::filesystem;
using namespace std;

// === UTILITY FUNCTIONS ===

void clean_output_dir(const std::string& output_dir) {
    if (fs::exists(output_dir))
        fs::remove_all(output_dir);
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

        if (token == "START_NODE") {
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
        std::cerr << "START_NODE and/or GENE column not found in header of file: " << file_name << std::endl;
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

    const std::string& snarl_start_id = columns[snarl_column];
    const std::string& snarl_end_id = columns[snarl_column+1];
    const std::string& gene_id  = columns[gene_column];

    const std::string composite_key = snarl_start_id + snarl_end_id + "_" + gene_id;
    map[composite_key] = line;
}

int parse_header(const std::string& header_line, 
    const std::string& file_name) {

    std::istringstream ss(header_line);
    std::string token;
    int snarl_column_index = 0;

    while (std::getline(ss, token, '\t')) {
        if (token == "START_NODE") {
            return snarl_column_index;
        }
        ++snarl_column_index;
    }

    std::cerr << "START_NODE column not found in header of file: " << file_name << std::endl;
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

    const std::string& snarl_start_id = columns[snarl_column];
    const std::string& snarl_end_id = columns[snarl_column+1];
    map[snarl_start_id + snarl_end_id] = line;
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
            std::cerr << "Missing START_NODE/END_NODE in file2: " << snarl << std::endl;
            return false;
        } else if (line1 != it->second) {
            std::cerr << "Mismatch for START_NODE/END_NODE " << snarl << ":\n"
                      << "File1: " << line1 << "\n"
                      << "File2: " << it->second << "\n";
            return false;
        }
    }

    // Compare file2 against file1
    for (const auto& [snarl, _] : map2) {
        if (map1.find(snarl) == map1.end()) {
            std::cerr << "Missing START_NODE/END_NODE in file1: " << snarl << std::endl;
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

    std::ifstream reinfile(file);
    for (size_t i = 1 ; i < record_count+1 ; i++) {
        if (!has_record[i]) {
            cerr << i << endl;
            cerr << "Fasta output should have one of the following records" << endl;
            for (const auto& header : set_to_header.at(i)) {
                cerr << "\t" << header << endl;
            }
            cerr << "Fasta output:" << endl;
            while (std::getline(reinfile, line)) {
                cerr << line << endl;
            }
            return false;
        }
    }
    reinfile.close();

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
        } else if (filename.find("assoc.pvalues") != std::string::npos) {
            if (!is_equivalent_assoc_file(expected_file, output_file)) {
                std::cerr << "Mismatch in assoc.pvalues file: " << expected_file << "and " << output_file << std::endl;
                return false;
            }
        } else if (filename.find("snarl_genotypes") == std::string::npos || filename.find("snarl_info") == std::string::npos) {
            // Ignore the snarl file because it gets done separately
            if (!is_equivalent_snarl_collection_file(expected_file, output_file)) {
                std::cerr << "Mismatch in file: " << filename << "\n";
                return false;
            }
        }
    }
    return true;
}

assoc_vals_t load_assoc_line(stoat::phenotype_type_t phenotype_type, const std::string& line) {
    assoc_vals_t vals;
    std::stringstream linestream(line);
    std::string part;

    //Chr
    std::getline(linestream, part, '\t');
    vals.chr = part;

    //start offset
    std::getline(linestream, part, '\t');
    vals.start_offset = part == "." ? std::numeric_limits<size_t>::max() : std::stoull(part);

    //end offset
    std::getline(linestream, part, '\t');
    vals.end_offset =part == "." ? std::numeric_limits<size_t>::max() :  std::stoull(part);
    
    //Snarl start node traversal
    std::getline(linestream, part, '\t');
    vals.start_node = part;

    //Snarl end node traversal
    std::getline(linestream, part, '\t');
    vals.end_node = part;

    //Allele lengths
    std::getline(linestream, part, '\t');
    if (part != ".") {
        std::stringstream lengthstream(part);
        std::string length;
        while (std::getline(lengthstream, length, ',')){
            vals.allele_lengths.emplace_back(std::move(length));
        }
    }

    // For eqtl, this is the gene name next
    if (phenotype_type == stoat::EQTL) {
        std::getline(linestream, part, '\t');
        vals.gene_name = part;
    }

    //p value
    std::getline(linestream, part, '\t');
    vals.p_value = part == "." ? std::numeric_limits<double>::max() : std::stod(part);

    if (phenotype_type == stoat::BINARY) {
        // Binary has the chi2 p value then allele count per phenotype
        std::getline(linestream, part, '\t');
        vals.p_value_chi2 = part == "." ? std::numeric_limits<double>::max() : std::stod(part);

        std::getline(linestream, part, '\t');
        if (part != ".") {
            std::stringstream countstream(part);
            std::string counts;
            while (std::getline(countstream, counts, ',')) {
                vals.allele_counts_per_pheno.emplace_back();
                std::stringstream per_phenostream(counts);
                std::string perpheno;
                while (std::getline(per_phenostream, perpheno, ':')) {
                    vals.allele_counts_per_pheno.back().emplace_back(part == "." ? std::numeric_limits<size_t>::max() : std::stoull(perpheno));
                }
            }
        }
    } else {
        // all others have allele count
        std::getline(linestream, part, '\t');
        if (part != ".") {
            std::stringstream countstream(part);
            std::string counts;
            while (std::getline(countstream, counts, ',')) {
                vals.allele_counts.emplace_back(part == "." ? std::numeric_limits<size_t>::max() : std::stoull(counts));
            }
        }
    }

    // depth
    std::getline(linestream, part, '\t');
    vals.depth = part == "." ? std::numeric_limits<size_t>::max() : std::stoull(part);

    assert(!std::getline(linestream, part, '\t'));

    return vals;
}

bool is_equivalent_assoc(stoat::phenotype_type_t phenotype_type, assoc_vals_t& vals1, assoc_vals_t& vals2) {
    if (! ((vals1.start_node == vals2.start_node && vals1.end_node == vals2.end_node) ||
           (vals1.start_node == vals2.end_node && vals1.end_node == vals2.start_node))) {
        return false;
    }

    if (vals1.chr != vals2.chr || vals1.start_offset != vals2.start_offset || vals1.end_offset != vals2.end_offset) {
        std::cerr << "For snarl " << vals1.start_node << vals1.end_node << ": mismatch in coordinates" << std::endl;
        std::cerr << "\t" << vals1.chr << ":" << vals1.start_offset << "-" << vals1.end_offset << std::endl;
        std::cerr << "\t" << vals2.chr << ":" << vals2.start_offset << "-" << vals2.end_offset << std::endl;
        return false;
    }
    if (vals1.p_value != vals2.p_value) {
        std::cerr << "For snarl " << vals1.start_node << vals1.end_node << ": mismatch in p-values" << std::endl;
        std::cerr << "For snarl " << vals2.start_node << vals2.end_node << ": mismatch in p-values" << std::endl;
        std::cerr << "\t" << vals1.p_value << std::endl;
        std::cerr << "\t" << vals2.p_value << std::endl;
        return false;
    }

    if (phenotype_type == stoat::BINARY) {
        //check chi2 p value
        if (vals1.p_value_chi2 != vals2.p_value_chi2) {
            std::cerr << "For snarl " << vals1.start_node << vals1.end_node << ": mismatch in p-values" << std::endl;
            std::cerr << "\t" << vals1.p_value_chi2 << std::endl;
            std::cerr << "\t" << vals2.p_value_chi2 << std::endl;
            return false;
        }
        if  (vals1.allele_counts_per_pheno.size() != vals2.allele_counts_per_pheno.size()) {
            std::cerr << "For snarl " << vals1.start_node << vals1.end_node << ": mismatch in number of alleles per pheno" << std::endl;
        }
    } else {
        if (phenotype_type == stoat::EQTL && vals1.gene_name != vals2.gene_name) {
            std::cerr << "For snarl " << vals1.start_node << vals1.end_node << ": mismatch in gene name" << std::endl;
            std::cerr << "\t" << vals1.gene_name << std::endl;
            std::cerr << "\t" << vals2.gene_name << std::endl;
            return false;
            
        }
        if  (vals1.allele_counts.size() != vals2.allele_counts.size()) {
            std::cerr << "For snarl " << vals1.start_node << vals1.end_node << ": mismatch in number of alleles counts" << std::endl;
        }
    }
    if  (vals1.allele_lengths.size() != vals2.allele_lengths.size()) {
        std::cerr << "For snarl " << vals1.start_node << vals1.end_node << ": mismatch in number of allele lengths" << std::endl;
    }

    // Now check allele_count/allele_counts_per_pheno and allele_lengths

    // For each allele in vals1, what is the matching allele in vals2?
    std::vector<size_t> allele_of_1_in_2 (vals1.allele_lengths.size(), std::numeric_limits<size_t>::max());
    std::vector<size_t> allele_of_2_in_1 (vals1.allele_lengths.size(), std::numeric_limits<size_t>::max());
    if (vals1.allele_counts.size() > 0 || vals1.allele_counts_per_pheno.size() > 0) {
        for (size_t allele_num1 = 0 ; allele_num1 < vals1.allele_lengths.size() ; allele_num1++) {
            // For each allele in 1, try to find an allele in 2 that matches it and hasn't already been assigned to another allele in 1
            for (size_t allele_num2 = 0 ; allele_num2 < vals2.allele_lengths.size() ; allele_num2++) {
                if (allele_of_2_in_1[allele_num2] == std::numeric_limits<size_t>::max() &&
                    vals1.allele_lengths[allele_num1] == vals2.allele_lengths[allele_num2]) {
                    // If the alleles have the same length, check the counts
                    bool match = false;
                    if (phenotype_type == stoat::BINARY) {
                        // If this is a binary phenotype, then there is a vector of counts per allele
                        std::sort(vals1.allele_counts_per_pheno[allele_num1].begin(), vals1.allele_counts_per_pheno[allele_num1].end());
                        std::sort(vals2.allele_counts_per_pheno[allele_num2].begin(), vals2.allele_counts_per_pheno[allele_num2].end());
                        if (vals1.allele_counts_per_pheno[allele_num1] == vals2.allele_counts_per_pheno[allele_num2]) {
                            allele_of_1_in_2[allele_num1] = allele_num2;
                            allele_of_2_in_1[allele_num2] = allele_num1;
                            break;
                        }
                    } else {
                        if (vals1.allele_counts[allele_num1] == vals2.allele_counts[allele_num2]) {
                            allele_of_1_in_2[allele_num1] = allele_num2;
                            allele_of_2_in_1[allele_num2] = allele_num1;
                            break;
                        }
                    }
                }
            }
        }

        // Each allele from vals1 should match exactly one allele in vals2 
        for (size_t assignment1 : allele_of_1_in_2) {
            if (assignment1 == std::numeric_limits<size_t>::max()) {
                return false;
            }
        }
        for (size_t assignment2 : allele_of_2_in_1) {
            if (assignment2 == std::numeric_limits<size_t>::max()) {
                return false;
            }
        }
    }
    return true;
}
bool is_equivalent_assoc_file(const std::string& file1, const std::string& file2){
    // Read each file and save contents as a map from start node to assoc_vals_t
    // Since the bounds might be flipped, save both start and end node to map

    std::unordered_map<std::string, assoc_vals_t> vals1;
    std::unordered_map<std::string, assoc_vals_t> vals2; 


    std::ifstream infile1(file1);
    std::string line;
    std::getline(infile1, line);
    // Get the phenotype type based on the contents of the header
    stoat::phenotype_type_t phenotype_type = stoat::QUANTITATIVE;
    if (line.find("P_CHI2") != std::string::npos) {
        phenotype_type = stoat::BINARY;
    } else if (line.find("GENE") != std::string::npos) {
        phenotype_type = stoat::EQTL;
    }
    while (std::getline(infile1, line)) {
        assoc_vals_t vals = load_assoc_line(phenotype_type, line);
        assert(vals1.count(vals.start_node+vals.gene_name) == 0);
        assert(vals1.count(vals.end_node+vals.gene_name) == 0);
        vals1.insert({vals.start_node+vals.gene_name, vals});
        vals1.insert({vals.end_node+vals.gene_name, std::move(vals)});
    }
    infile1.close();

    std::ifstream infile2(file2);
    std::getline(infile2, line);

    // verify the phenotype type based on the contents of the header
    if (line.find("P_CHI2") != std::string::npos) {
        if (phenotype_type != stoat::BINARY) {
            std::cerr << "Files come from two different types of phenotype" << std::endl;
            return false;
        }
    } else if (line.find("GENE") != std::string::npos) {
        if (phenotype_type != stoat::EQTL) {
            std::cerr << "Files come from two different types of phenotype" << std::endl;
            return false;
        }
    } else {
        // This doesn't matter since the other two are the only ones with special fields
        if (phenotype_type != stoat::QUANTITATIVE) {
            std::cerr << "Files come from two different types of phenotype" << std::endl;
            return false;
        }
    }
    while (std::getline(infile2, line)) {
        assoc_vals_t vals = load_assoc_line(phenotype_type, line);
        assert(vals2.count(vals.start_node+vals.gene_name) == 0);
        assert(vals2.count(vals.end_node+vals.gene_name) == 0);
        vals2.insert({vals.start_node+vals.gene_name, vals});
        vals2.insert({vals.end_node+vals.gene_name, std::move(vals)});
    }
    infile2.close();

    // Now check that for each snarl in 1, there is an equivalent snarl in 2
    for (auto& snarl1 : vals1) {
        if (vals2.count(snarl1.first) == 0) {
            std::cerr << "Snarl " << snarl1.first << ": " << snarl1.second.start_node << snarl1.second.end_node << " missing in second file" << std::endl;
            return false;
        }
        if (!is_equivalent_assoc(phenotype_type, snarl1.second, vals2[snarl1.first])) {
            return false;
        }
    }
    // And for each snarl in 2, there is an equivalent snarl in 1
    for (auto& snarl2 : vals2) {
        if (vals1.count(snarl2.first) == 0) {
            std::cerr << "Snarl " << snarl2.first << ": " << snarl2.second.start_node << snarl2.second.end_node << " missing in first file" << std::endl;
            return false;
        }
        if (!is_equivalent_assoc(phenotype_type, snarl2.second, vals1[snarl2.first])) {
            return false;
        }
    }

    return true;
}

bool is_equivalent_snarl_collection_file(const std::string& file1, const std::string& file2){

    std::shared_ptr<SnarlCoordinates> snarl_coords1 (new SnarlCoordinates);
    SnarlDataCollection snarls1(snarl_coords1, 0);
    StdReader reader1(file1);
    snarls1.load_snarl_data_collection(reader1);
    reader1.close();
    
    std::shared_ptr<SnarlCoordinates> snarl_coords2 (new SnarlCoordinates);
    SnarlDataCollection snarls2(snarl_coords2, 0);
    StdReader reader2(file2);
    snarls2.load_snarl_data_collection(reader2);
    reader2.close();
    
    return SnarlDataCollection::is_equivalent(snarls1, snarls2);

}
