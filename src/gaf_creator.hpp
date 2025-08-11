#ifndef GAF_CREATOR_HPP
#define GAF_CREATOR_HPP

#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <regex>
#include <fstream>

#include <bdsg/hash_graph.hpp>
#include <bdsg/packed_graph.hpp>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>

#include "snarl_data_t.hpp"
#include "utils.hpp"
#include "snarl_analyzer.hpp"
#include "matrix.hpp"
#include "binary_table.hpp"

using namespace std;

namespace stoat_vcf {

std::pair<double, double> calcul_proportion_signi(size_t number_ind_group0, size_t number_ind_group1, double p_value);
std::string addSuffixToFilename(const std::string& filename, const std::string& suffix);
std::string add_suffix_to_filename(const std::string& filename, const std::string& suffix);
void write_gaf_lines(const std::string& sequence_name, const std::string& path, int length, double prop, ofstream& outfile);

std::vector<size_t> decompose_snarl(const std::string& snarl);
int calcul_path_length(bdsg::PackedGraph& pg, const std::string& snarl);
void write_gaf_lines(const std::string& sequence_name, const std::string& path, int length, double prop, ofstream& outfile);
void gaf_creation(const std::string& input_file, std::unordered_map<std::string, std::vector<stoat::Snarl_data_t>>& snarl_chr,
    bdsg::PackedGraph& pg, const std::string& output_file);

} //end stoat namespace

#endif // GAF_CREATOR_HPP