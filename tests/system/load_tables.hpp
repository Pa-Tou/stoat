#include <iostream>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <regex>


////////////////////////////////////////// snarl_genotype.tsv
struct snarl_genotype_values_t {
    std::string start_node;
    std::string end_node;
    std::string reference;
    size_t start_offset;
    size_t end_offset;
    size_t depth;
    std::vector<std::string> allele_lengths;
    std::vector<std::string> walks;
    std::vector<std::string> sequences;
    std::unordered_map<std::string, size_t> sample_to_allele;

    bool operator==(const snarl_genotype_values_t& other) const;
};

inline std::ostream& operator<<(std::ostream& out, const snarl_genotype_values_t& vals) {
    std::ostringstream rest;
    rest << "\tallele_lengths: ";
    for (std::string len : vals.allele_lengths) {
        rest << len << ","; 
    }
    rest << std::endl << "\twalks: ";
    for (std::string w : vals.walks) {
        rest << w << ",";
    }
    rest << std::endl << "\tsequences: ";
    for (std::string w : vals.sequences) {
        rest << w << ",";
    }
    rest << std::endl << "\tallele assignments: " << std::endl;
    for (const auto& x : vals.sample_to_allele) {
        rest << "\t\t" << x.first << ": " << x.second << std::endl;
    }
    return out << "Snarl between nodes " << vals.start_node << " and " << vals.end_node << std::endl
               << "\ton reference " << vals.reference << ":" << vals.start_offset << "-" << vals.end_offset << std::endl
               << "\tat depth " << vals.depth << std::endl
               << rest.str() << std::endl;
}
// Given the file name, load the snarl genotypes into a vector
// Throw an error if the file seems bad
std::vector<snarl_genotype_values_t> load_genotype_file(const std::string& infile);


/////////////////////////////////////////// Binary table from stoat test
struct binary_table_values_t{
    std::string chr;
    size_t start_pos;
    size_t end_pos;
    std::pair<std::string, std::string> snarl_ids;
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


