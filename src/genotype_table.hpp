#ifndef genotype_table_HPP
#define genotype_table_HPP

#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <iomanip>
#include <sstream>

#include "matrix.hpp"
#include "snarl_analyzer.hpp"
#include "utils.hpp"

using namespace std;

namespace stoat_vcf {

    // Given two vectors of genotypes representing two groups (with length number_paths), fill them in with counts of the number of times each path is seen
    // g0 and g1 can be used in binary_stat_test()
    std::vector<std::vector<char>> create_genotype_table(
        const size_t &number_samples,
        const std::vector<stoat::Path_traversal_t> &column_headers,
        const stoat_vcf::EdgeBySampleMatrix &matrix);

} // namespace stoat

#endif
