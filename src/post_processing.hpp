#ifndef POST_PROCESSING_HPP
#define POST_PROCESSING_HPP

#include <vector>
#include <tuple>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdio>
#include "utils.hpp"

using namespace std;

namespace stoat{

// Given a vector of <p-value, 1.0, line index from the input file>, fill in the vector with the adjusted p-value
// and sort the vector by adjusted p-value
void adjust_pvalues_with_BH(std::vector<std::tuple<double, double, size_t>>& data);

// Read a tsv from input_file, collect the p-values from the correct column (depending on phenotype_type), 
// and write the same file plus a BH-adjusted p-value to outupt_file_significant.
void add_BH_adjusted_column(
    const std::string& input_file,
    const std::string& output_dir,
    const std::string& output_file_significant,
    const stoat::phenotype_type_t& phenotype_type);

// The same, but specify the column number (0-indexed) of the p-value and the adjusted p-value
void add_BH_adjusted_column(
    const std::string& input_file,
    const std::string& output_dir,
    const std::string& output_file_significant,
    size_t p_col, size_t adjusted_col_index);

} // namespace stoat_vcf

#endif // ADJUSTED_PVALUE_HPP
