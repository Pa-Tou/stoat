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

// Check if a fasta file is equivalent to a set of fasta records.
// Because there may be multiple options for which path is represented in the fasta (eg two paths that take the same walk), this must allow different headers in an equivalence class (walk through a snarl).
// fasta_records is a tuple of <equivalence class, header, sequence>
// There should be as many lines are there are equivalence classes
// This doesn't check the values of equivalence classes so they are assumed to start at 1 and increase by 1
bool fasta_equal(const std::string& file, const std::vector<std::tuple<size_t, std::string, std::string>>& fasta_records);

bool compare_output_dirs(const std::string& output_dir, const std::string& expected_dir);
