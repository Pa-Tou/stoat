#include <iostream>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <regex>

struct binary_table_values_t{
    std::string chr;
    size_t start_pos;
    size_t end_pos;
    std::pair<size_t, size_t> snarl_ids;
    std::vector<std::string> path_lengths; // min/max lengths of paths per allele
    std::string p_fishers;
    std::string p_chi2;
    std::vector<std::string> group_paths; // count of each phenotype group per allele
    size_t depth;

    // Everything must be equal, and the path lengths and group paths must correspond to each other and have the same values, but 
    // they may be in a different order
    bool operator==(const binary_table_values_t& other) const;
};

//load info from a line from the binary tsv file 
//CHR START_POS END_POS SNARL PATH_LENGTHS P_FISHER P_CHI2 GROUP_PATHS DEPTH
binary_table_values_t load_binary_snarl_line(std::string);
