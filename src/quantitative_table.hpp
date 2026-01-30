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

// Write a std::string of: g0[0]:g1[1],g0[1]:g1[1],g0[2]:g1[2]...
// JEAN copied to feature table, remove this one when the migration is over
std::string format_group_paths(const std::vector<size_t>& g0, const std::vector<size_t>& g1);
    
} // namespace stoat


#endif
