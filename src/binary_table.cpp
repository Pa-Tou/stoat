#include "binary_table.hpp"

// ------------------------ Binary table & stats ------------------------
namespace stoat_vcf {

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

std::pair<size_t, size_t> create_binary_table(
    std::vector<size_t>& g0, std::vector<size_t>& g1,
    const std::vector<bool>& binary_phenotype, 
    const std::vector<stoat::Path_traversal_t>& list_path_snarl, 
    const size_t& number_paths,
    const size_t& number_samples,
    const stoat_vcf::EdgeBySampleMatrix& matrix) {

    size_t total_sum = 0;
    std::vector<bool> sample_included(number_samples, false);
    for (size_t idx_g = 0; idx_g < number_paths; ++idx_g) {
        const stoat::Path_traversal_t& path_snarl = list_path_snarl[idx_g];
        std::vector<stoat::Edge_t> list_edge_path = stoat_vcf::decompose_path_to_edges(path_snarl);
        std::vector<size_t> idx_srr_save = stoat_vcf::identify_path(list_edge_path, matrix, number_samples * 2);

        for (size_t idx : idx_srr_save) {
            bool group = binary_phenotype[idx / 2];
            sample_included[idx / 2] = true;
            if (group) {
                g1[idx_g] += 1;
            } else {
                g0[idx_g] += 1;
            }
            total_sum++;
        }
    }

    // Count the number of individuals included
    size_t individuals_included = 0;
    for (bool included : sample_included) {
        if (included) individuals_included++;
    }

    return {total_sum, individuals_included};
}

} // namespace stoat
