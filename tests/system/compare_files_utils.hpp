#include <iostream>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <regex>
#include "../../src/types_and_structs.hpp"

namespace fs = std::filesystem;

// === UTILITY FUNCTIONS ===

// rm -rf this directory
void clean_output_dir(const std::string& output_dir);

// Check that the two directories have the same files and that the files match.
// This only knows how to compare snarl_info, snarl_genotypes, and assoc.pvalues files using functions is_equivalent_snarl_collection_file() and is_equivalent_assoc_file
bool compare_output_dirs(const std::string& output_dir, const std::string& expected_dir);

// Check if a fasta file is equivalent to a set of fasta records.
// Because there may be multiple options for which path is represented in the fasta (eg two paths that take the same walk), this must allow different headers in an equivalence class (walk through a snarl).
// fasta_records is a tuple of <equivalence class, header, sequence>
// There should be as many lines are there are equivalence classes
// This doesn't check the values of equivalence classes so they are assumed to start at 1 and increase by 1
bool fasta_equal(const std::string& file, const std::vector<std::tuple<size_t, std::string, std::string>>& fasta_records);

// Check if a fasta is valid: lines are only headers or sequences and sequence lines are less than 80 characters
bool is_valid_fasta(const std::string& file);

/// Load and compare two SnarlDataCollections (snarl_info or snarl_genotypes)
/// Returns true if they are equivalent and false otherwise.
/// Values must be the same but the order of alleles, orientation of paths, etc may change
bool is_equivalent_snarl_collection_file(const std::string& file1, const std::string& file2);

///////////////////////////////// Load and compare an assoc file (output of stoat test)

/// Load and compare two assoc.pvalues files, returns true if they are equivalent and false otherwise
/// All values must be the same but the order of alleles may change
bool is_equivalent_assoc_file(const std::string& file1, const std::string& file2);

struct assoc_vals_t {
    std::string start_node;
    std::string end_node;
    std::string chr;
    size_t start_offset;
    size_t end_offset;
    std::vector<std::string> allele_lengths;//Kept as a string to compare min/max
    double p_value;
    double p_value_chi2;//For binary, p_value is fisher and this is chi2
    std::vector<size_t> allele_counts;
    std::vector<std::vector<size_t>> allele_counts_per_pheno;
    size_t depth;
    std::string gene_name;
};

assoc_vals_t load_assoc_line(stoat::phenotype_type_t phenotype_type, const std::string& line);

// Are the assoc vals equivalent? True for equivalent, false for mismatch
bool is_equivalent_assoc(stoat::phenotype_type_t phenotype_type, assoc_vals_t& vals1, assoc_vals_t& vals2);

