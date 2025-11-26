#include "genotype_table.hpp"

// ------------------------ VCF table & stats ------------------------
namespace stoat_vcf {

std::vector<std::vector<char>> create_genotype_table(
        const size_t &number_samples,
        const std::vector<stoat::PathTraversal> &column_headers,
        const stoat_vcf::EdgeBySampleMatrix &matrix) {

    // Each sample can have up to 2 alleles → initialize with two empty strings
    std::vector<std::vector<char>> genotypes(number_samples, std::vector<char>(2, '.'));

    // Loop over all columns / paths
    for (size_t col_idx = 0; col_idx < column_headers.size(); ++col_idx) {
        // Identify samples that use this path
        const auto& path_snarl = column_headers[col_idx];
        std::vector<size_t> idx_sample_save = matrix.get_samples_on_path(path_snarl);

        // Fill genotype matrix
        for (size_t idx : idx_sample_save) {
            size_t sample_idx = idx / 2;
            size_t allele_idx = idx % 2;  // 0 for first allele, 1 for second allele
            genotypes[sample_idx][allele_idx] = static_cast<char>('0' + (col_idx % 10));
        }
    }
    return genotypes;
}


} // namespace stoat
