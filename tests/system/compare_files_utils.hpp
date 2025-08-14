#include <iostream>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <regex>

namespace fs = std::filesystem;

// === UTILITY FUNCTIONS ===

void clean_output_dir(const std::string& output_dir);

// Check if two TSV files contain the same lines. 
// This assumes that the TSVs have a column for a unique snarl identifier, and checks that the snarls are the same and have the same lines
bool files_equal(const std::string& file1, const std::string& file2);

// Check if two tsv files contain the same lines, but this time one is in the form of a vector of lines instead of a file
bool files_equal(const std::string& file1, const std::vector<std::string>& lines);


// Helper function for other files_equal() functions. Each of the others transforms the input into maps and calls this 
bool files_equal(const std::unordered_map<std::string, std::string>& map1, const std::unordered_map<std::string, std::string>& map2);


// Check if a fasta file is equivalent to a set of fasta records.
// Because there may be multiple options for which path is represented in the fasta (eg two paths that take the same walk), this must allow different headers in an equivalence class (walk through a snarl).
// fasta_records is a tuple of <equivalence class, header, sequence>
// There should be as many lines are there are equivalence classes
// This doesn't check the values of equivalence classes so they are assumed to start at 1 and increase by 1
bool fasta_equal(const std::string& file, const std::vector<std::tuple<size_t, std::string, std::string>>& fasta_records);

// Check if a fasta is valid: lines are only headers or sequences and sequence lines are less than 80 characters
bool is_valid_fasta(const std::string& file);

bool compare_output_dirs(const std::string& output_dir, const std::string& expected_dir);

// Helper function to load a tsv from a file into a map from snarl to line
void load_tsv_file (const std::string& path, std::unordered_map<std::string, std::string>& map);

// Helper function to load a tsv from a vector of lines into a map from snarl to line
void load_tsv_file (const std::vector<std::string>& lines, std::unordered_map<std::string, std::string>& map);

// Helper function for the two load_tsv_file() functions. Called per line
void process_tsv_line(const std::string& line, std::unordered_map<std::string, std::string>& map, int& snarl_column, bool& header_processed, const std::string& file_name);
