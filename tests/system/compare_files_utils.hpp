#include <iostream>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

namespace fs = std::filesystem;

// === UTILITY FUNCTIONS ===

void clean_output_dir(const std::string& output_dir);

bool files_equal(const std::string& file1, const std::string& file2);

bool compare_output_dirs(const std::string& output_dir, const std::string& expected_dir);
