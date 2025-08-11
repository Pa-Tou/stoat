#include "quantitative_table.hpp"

using namespace std;

namespace stoat_vcf {

// Explicit template instantiations
template std::tuple<std::vector<std::vector<double>>, std::vector<double>, std::vector<size_t>>
    create_quantitative_table<double>(
        const size_t&,
        const std::vector<stoat::Path_traversal_t>&,
        const std::vector<double>&,
        const stoat_vcf::EdgeBySampleMatrix&);

template std::tuple<std::vector<std::vector<double>>, std::vector<bool>, std::vector<size_t>>
    create_quantitative_table<bool>(
        const size_t&,
        const std::vector<stoat::Path_traversal_t>&,
        const std::vector<bool>&,
        const stoat_vcf::EdgeBySampleMatrix&);

std::tuple<std::vector<std::vector<double>>,  
std::unordered_set<size_t>, std::vector<size_t>> process_table_quantitative(
        const size_t& number_samples,
        const std::vector<stoat::Path_traversal_t>& column_headers,
        const stoat_vcf::EdgeBySampleMatrix& matrix) {

    size_t length_column = column_headers.size();

    std::vector<size_t> allele_paths(length_column, 0);

    std::vector<std::vector<double>> genotypes(number_samples);
    for (auto& row : genotypes)
        row.reserve(length_column);

    std::vector<size_t> kept_columns; // Indices of valid columns
    std::unordered_set<size_t> index_used;

    // Loop over all columns
    for (size_t col_idx = 0; col_idx < length_column; ++col_idx) {
        const stoat::Path_traversal_t& path_snarl = column_headers[col_idx];
        std::vector<stoat::Edge_t> list_edge_path = stoat_vcf::decompose_path_to_edges(path_snarl);

        //Get the indices of all samples that take this path
        std::vector<size_t> idx_srr_save = identify_path(list_edge_path, matrix, number_samples * 2);

        if (idx_srr_save.empty())
            continue; // Skip if column is empty

        kept_columns.push_back(col_idx); // Valid column

        // Ensure rows have space for new column
        for (size_t i = 0; i < number_samples; ++i) {
            if (genotypes[i].size() < kept_columns.size())
                genotypes[i].resize(kept_columns.size(), 0.0);
        }

        size_t numb_all = idx_srr_save.size();
        allele_paths[col_idx] = numb_all;

        // Fill genotype matrix
        for (size_t idx : idx_srr_save) {
            size_t srr_idx = idx / 2;
            genotypes[srr_idx][kept_columns.size() - 1] += 1.0;
            index_used.insert(srr_idx);
        }
    }

    return {genotypes, index_used, allele_paths};   
}

// Function template definition
template<typename T>
std::tuple<std::vector<std::vector<double>>, std::vector<T>,  std::vector<size_t>> create_quantitative_table(
    const size_t& number_samples,
    const std::vector<stoat::Path_traversal_t>& column_headers,
    const std::vector<T>& phenotype,
    const stoat_vcf::EdgeBySampleMatrix& matrix) {

    const auto& [genotypes, index_used, allele_paths] = 
    process_table_quantitative(number_samples, column_headers, matrix);

    std::vector<std::vector<double>> genotypes_filtered;
    genotypes_filtered.reserve(index_used.size());

    std::vector<T> phenotype_filtered;
    phenotype_filtered.reserve(index_used.size());

    for (size_t i : index_used) {
        const auto& row = genotypes[i];
        double row_sum = std::accumulate(row.begin(), row.end(), 0.0);

        std::vector<double> normalized_row;
        size_t max_col = row.size();
        normalized_row.reserve(max_col);

        for (size_t j = 0; j < max_col; ++j) {
            normalized_row.push_back(row[j] > 0.0 ? row[j] / row_sum : 0.0);
        }

        genotypes_filtered.push_back(std::move(normalized_row));
        phenotype_filtered.push_back(phenotype[i]);
    }

    return {genotypes_filtered, phenotype_filtered, allele_paths};
}

std::tuple<std::vector<std::vector<double>>, std::unordered_set<size_t>, std::vector<size_t>> create_eqtl_table(
    const size_t& number_samples,
    const std::vector<stoat::Path_traversal_t>& column_headers,
    const stoat_vcf::EdgeBySampleMatrix& matrix) {

    const auto& [genotypes, index_used, allele_paths] = 
    process_table_quantitative(number_samples, column_headers, matrix);

    std::vector<std::vector<double>> genotypes_filtered;
    genotypes_filtered.reserve(index_used.size());

    for (size_t i : index_used) {
        const auto& row = genotypes[i];
        double row_sum = std::accumulate(row.begin(), row.end(), 0.0);

        std::vector<double> normalized_row;
        size_t max_col = row.size();
        normalized_row.reserve(max_col);

        for (size_t j = 0; j < max_col; ++j) {
            normalized_row.push_back(row[j] > 0.0 ? row[j] / row_sum : 0.0);
        }

        genotypes_filtered.push_back(std::move(normalized_row));
    }

    return {genotypes_filtered, index_used, allele_paths};
}

} // namespace stoat
