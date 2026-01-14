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
#include <set>
#include <tuple>
#include <iomanip>

// moved from (old) binary_table.hpp
#include <unordered_map>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include "arg_parser.hpp"
#include "matrix.hpp"
#include "snarl_analyzer.hpp"
#include "utils.hpp"
#include "stats_test.hpp"


namespace stoat_vcf {

// Return a tuple of genotypes, index_used, allele_paths
std::tuple<std::vector<std::vector<double>>, std::set<size_t>, std::vector<size_t>>
process_table_quantitative(
    const size_t& number_samples,
    const std::vector<stoat::PathTraversal>& column_headers,
    const stoat_vcf::EdgeBySampleMatrix& matrix);

// Given the number of samples (length_sample), the paths through the snarl (column_headers), the binary or quantitative phenotype of each sample (phenotype)
// and a matrix of edges in each sample
// Return a tuple of 
// - genotypes_filtered: a matrix where each row is a sample, each column is an allele (from column_headers), counts divided by the sum of each row
// - phenotype_filtered: the phenotypes for each genotype 
// - allele_paths: the number of samples that take each path through the snarl (per column) 
template <typename T>
std::tuple<std::vector<std::vector<double>>, std::vector<T>, std::vector<std::string>, std::vector<size_t>>
create_quantitative_table(
    const size_t& number_samples,
    const std::vector<stoat::PathTraversal>& column_headers,
    const std::vector<T>& phenotype,
    const stoat_vcf::EdgeBySampleMatrix& matrix);

// Given the number of samples (length_sample), the paths through the snarl (column_headers), and a matrix of edges in each sample,
// Return a tuple of 
// - genotypes_filtered: a matrix where each row is a sample, each column is an allele (from column_headers), counts divided by the sum of each row
// - index_used: row (samples) indices that were filled in
// - allele_paths: the number of samples that take each path through the snarl (per column) 
std::tuple<std::vector<std::vector<double>>, std::set<size_t>, std::vector<std::string>, std::vector<size_t>>
create_eqtl_table(
    const size_t& number_samples,
    const std::vector<stoat::PathTraversal>& column_headers,
    const stoat_vcf::EdgeBySampleMatrix& matrix);

// ------------------------ Binary table ------------------------

// Write a std::string of: g0[0]:g1[1],g0[1]:g1[1],g0[2]:g1[2]...
std::string format_group_paths(const std::vector<size_t>& g0, const std::vector<size_t>& g1);

/// Given two (empty) vectors of genotypes representing two groups (with length number_paths), fill them in with counts of the number of times each path is seen  
/// g0 and g1 can be used in binary_stat_test()
/// Return the number of individuals represented in the table
size_t create_binary_table(
    std::vector<size_t>& g0, std::vector<size_t>& g1,
    const std::vector<bool>& binary_phenotype, 
    const std::vector<stoat::PathTraversal>& list_path_snarl, 
    const stoat_vcf::EdgeBySampleMatrix& matrix);

    
} // namespace stoat


#endif
