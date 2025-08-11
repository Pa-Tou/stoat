#ifndef quantitative_table_HPP
#define quantitative_table_HPP

#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <unordered_set>
#include <tuple>
#include <iomanip>

#include "arg_parser.hpp"
#include "matrix.hpp"
#include "snarl_analyzer.hpp"
#include "utils.hpp"
#include "stats_test.hpp"

using namespace std;

namespace stoat_vcf {

// Return a tuple of genotypes, index_used, allele_paths
std::tuple<std::vector<std::vector<double>>, std::unordered_set<size_t>, std::vector<size_t>>
process_table_quantitative(
    const size_t& number_samples,
    const std::vector<stoat::Path_traversal_t>& column_headers,
    const stoat_vcf::EdgeBySampleMatrix& matrix);

// Given the number of samples (length_sample), the paths through the snarl (column_headers), the binary or quantitative phenotype of each sample (phenotype)
// and a matrix of edges in each sample
// Return a tuple of 
// - genotypes_filtered: a matrix where each row is a sample, each column is an allele (from column_headers), counts divided by the sum of each row
// - phenotype_filtered: the phenotypes for each genotype 
// - allele_paths: the number of samples that take each path through the snarl (per column) 
template <typename T>
std::tuple<std::vector<std::vector<double>>, std::vector<T>, std::vector<size_t>>
create_quantitative_table(
    const size_t& number_samples,
    const std::vector<stoat::Path_traversal_t>& column_headers,
    const std::vector<T>& phenotype,
    const stoat_vcf::EdgeBySampleMatrix& matrix);

// Given the number of samples (length_sample), the paths through the snarl (column_headers), and a matrix of edges in each sample,
// Return a tuple of 
// - genotypes_filtered: a matrix where each row is a sample, each column is an allele (from column_headers), counts divided by the sum of each row
// - index_used: row (samples) indices that were filled in
// - allele_paths: the number of samples that take each path through the snarl (per column) 
std::tuple<std::vector<std::vector<double>>, std::unordered_set<size_t>, std::vector<size_t>>
create_eqtl_table(
    const size_t& number_samples,
    const std::vector<stoat::Path_traversal_t>& column_headers,
    const stoat_vcf::EdgeBySampleMatrix& matrix);

} // namespace stoat

#endif
