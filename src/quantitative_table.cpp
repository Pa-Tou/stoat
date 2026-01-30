#include "quantitative_table.hpp"

namespace stoat_vcf {

// JEAN remove once the feature table migration is over
std::string format_group_paths(const std::vector<size_t>& g0, const std::vector<size_t>& g1) {

    std::string result;
    size_t numb_col = g0.size();
    for (size_t index_col = 0; index_col < numb_col; ++index_col) {
        result += std::to_string(g0[index_col]) + ":" + std::to_string(g1[index_col]);
        if (index_col < numb_col - 1) {
            result += ","; // Separate row pairs with ','
        }
    }
    return result;
}

} // namespace stoat
