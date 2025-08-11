#include "snarl_analyzer.hpp"
#include "matrix.hpp"
#include "binary_table.hpp"
#include "quantitative_table.hpp"
#include "utils.hpp"
#include "arg_parser.hpp"
#include "writer.hpp"
#include "omp.h"

namespace stoat_vcf {

SnarlAnalyzer::SnarlAnalyzer(
    const std::unordered_map<std::string, std::vector<stoat::Snarl_data_t>>& chr_to_snarl_data, 
    EdgeBySampleMatrix& edge_matrix,
    const std::vector<std::string>& list_samples, 
    const std::vector<std::vector<double>>& covariate, 
    const double& maf_threshold, 
    const double& table_threshold,
    const size_t& min_individuals,
    const size_t& min_haplotypes,
    const std::string& regression_dir) :

        chr_to_snarl_data(chr_to_snarl_data), 
        edge_matrix(edge_matrix),
        list_samples(list_samples), 
        covariate(covariate), 
        maf_threshold(maf_threshold), 
        table_threshold(table_threshold),
        min_individuals(min_individuals),
        min_haplotypes(min_haplotypes),
        regression_dir(regression_dir)
        {};

BinarySnarlAnalyzer::BinarySnarlAnalyzer(
    const std::unordered_map<std::string, std::vector<stoat::Snarl_data_t>>& chr_to_snarl_data,
    EdgeBySampleMatrix& edge_matrix,
    const std::vector<std::string>& list_samples, 
    const double& maf_threshold,
    const double& table_threshold,
    const std::vector<bool>& binary_phenotype,
    const size_t& min_individuals,
    const size_t& min_haplotypes,
    const std::string& regression_dir) :

        SnarlAnalyzer(chr_to_snarl_data, edge_matrix, list_samples, {}, maf_threshold, table_threshold, min_individuals, min_haplotypes, regression_dir), 
        binary_phenotype(binary_phenotype), fk() {};

BinaryCovarSnarlAnalyzer::BinaryCovarSnarlAnalyzer(
    const std::unordered_map<std::string, std::vector<stoat::Snarl_data_t>>& chr_to_snarl_data,
    EdgeBySampleMatrix& edge_matrix,
    const std::vector<std::string>& list_samples, 
    const std::vector<std::vector<double>>& covariate, 
    const double& maf_threshold, 
    const double& table_threshold,
    const std::vector<bool>& binary_phenotype,
    const size_t& min_individuals,
    const size_t& min_haplotypes,
    const std::string& regression_dir) :

        SnarlAnalyzer(chr_to_snarl_data, edge_matrix, list_samples, covariate, maf_threshold, table_threshold, min_individuals, min_haplotypes, regression_dir), 
        binary_phenotype(binary_phenotype), lr() {};

QuantitativeSnarlAnalyzer::QuantitativeSnarlAnalyzer(
    const std::unordered_map<std::string, std::vector<stoat::Snarl_data_t>>& chr_to_snarl_data, 
    EdgeBySampleMatrix& edge_matrix,
    const std::vector<std::string>& list_samples, 
    const std::vector<std::vector<double>>& covariate,
    const double& maf_threshold, 
    const double& table_threshold,
    const std::vector<double>& quantitative_phenotype,
    const size_t& min_individuals,
    const size_t& min_haplotypes,
    const std::string& regression_dir) :

        SnarlAnalyzer(chr_to_snarl_data, edge_matrix, list_samples, covariate, maf_threshold, table_threshold, min_individuals, min_haplotypes, regression_dir), 
        quantitative_phenotype(quantitative_phenotype), lr() {};

EQTLSnarlAnalyzer::EQTLSnarlAnalyzer(
    const std::unordered_map<std::string, std::vector<stoat::Snarl_data_t>>& chr_to_snarl_data, 
    EdgeBySampleMatrix& edge_matrix,
    const std::vector<std::string>& list_samples, 
    const std::vector<std::vector<double>>& covariate, 
    const double& maf_threshold, 
    const double& table_threshold,
    const std::unordered_map<std::string, std::vector<stoat_vcf::Qtl_data>>& eqtl_map,
    const size_t& windows_gene_threshold,
    const size_t& min_individuals,
    const size_t& min_haplotypes,
    const std::string& regression_dir) :

        SnarlAnalyzer(chr_to_snarl_data, edge_matrix, list_samples, covariate, maf_threshold, table_threshold, min_individuals, min_haplotypes, regression_dir), 
        eqtl_map(eqtl_map), windows_gene_threshold(windows_gene_threshold), lr() {};

void BinarySnarlAnalyzer::write_header(std::ofstream& outf) {
    stoat::write_binary_header(outf);
}

void BinaryCovarSnarlAnalyzer::write_header(std::ofstream& outf) {
    stoat::write_binary_covar_header(outf);
}

void QuantitativeSnarlAnalyzer::write_header(std::ofstream& outf) {
    stoat::write_quantitative_header(outf);
}

void EQTLSnarlAnalyzer::write_header(std::ofstream& outf) {
    stoat::write_eqtl_header(outf);
}

void SnarlAnalyzer::process_snarls_by_chromosome_chunk(
    htsFile* &ptr_vcf,
    bcf_hdr_t* &hdr,
    bcf1_t* &rec, 
    const std::string& output_filename) {

    std::ofstream outf(output_filename, std::ios::binary);

    // Write the header
    write_header(outf);

    // Go through the vcf and get chunks by chromosome. 
    while (bcf_read(ptr_vcf, hdr, rec) >= 0) {

        std::string chr = bcf_hdr_id2name(hdr, rec->rid);

        // Skip chromosomes not in chr_to_snarl_data
        while (chr_to_snarl_data.find(chr) == chr_to_snarl_data.end()) {
            stoat::LOG_WARN("Chromosome " + chr + " not found in snarl paths file. Skipping.");

            bool found_new_chr = false;
            while (bcf_read(ptr_vcf, hdr, rec) >= 0) {
                std::string chr_next = bcf_hdr_id2name(hdr, rec->rid);
                if (chr_next != chr) {
                    chr = chr_next;  // Update to the new chromosome
                    found_new_chr = true;
                    break;
                }
            }

            if (!found_new_chr) {
                return;  // exit if no more records are available
            }
        }

        auto start_1 = std::chrono::high_resolution_clock::now();

        stoat::LOG_INFO("Analysing chr : " + chr);
        size_t size_chr = chr_to_snarl_data.at(chr).size();

        // Make genotype matrix by chromosome    
        auto [ptr_vcf_new, hdr_new, rec_new] = make_edge_matrix(ptr_vcf, hdr, rec, chr, size_chr);
        ptr_vcf = ptr_vcf_new;
        hdr = hdr_new;
        rec = rec_new;

        auto end_1 = std::chrono::high_resolution_clock::now();
        stoat::LOG_INFO("Matrix time analysis in chr " + chr + " : " + std::to_string(std::chrono::duration<double>(end_1 - start_1).count()) + " s");

        const auto& snarls = chr_to_snarl_data.at(chr);
        
        auto start_2 = std::chrono::high_resolution_clock::now();

        #pragma omp parallel for schedule(static)

        // Make the snarl test analysis
        // Iterate over each snarl
        for (const stoat::Snarl_data_t& snarl_data_s : snarls) {
            analyze_and_write_snarl(snarl_data_s, chr, outf);
        }

        auto end_2 = std::chrono::high_resolution_clock::now();

        stoat::LOG_INFO("Snarl time analysis in chr " + chr + " : " + std::to_string(std::chrono::duration<double>(end_2 - start_2).count()) + " s");
        stoat::LOG_INFO("Total time for chr " + chr + " : " + std::to_string(std::chrono::duration<double>(end_2 - start_1).count()) + " s");

    }

    // Cleanup
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(ptr_vcf);
}

std::tuple<htsFile*, bcf_hdr_t*, bcf1_t*> SnarlAnalyzer::make_edge_matrix(htsFile *ptr_vcf, bcf_hdr_t *hdr, bcf1_t *rec, std::string &chr, size_t &num_paths_chr) {

    edge_matrix.reset(list_samples, num_paths_chr*4, list_samples.size() * 2);

    // loop over the VCF file for each line and stop where chr is different
    do {
        bcf_unpack(rec, BCF_UN_STR);

        // Check the INFO field for LV (Level Variant) and skip if LV != 0
        int32_t *lv = nullptr;
        int n_lv = 0;

        // Extract LV field from INFO skip if variant is lv != 0 to avoid duplication paths/snarl variant analysis
        if (bcf_get_info_int32(hdr, rec, "LV", &lv, &n_lv) > 0) {
            if (lv[0] != 0) {
                free(lv);
                continue;
            }
        }
        free(lv);

        // Extract genotypes (GT)
        int ngt = 0;
        int32_t *gt = nullptr;
        ngt = bcf_get_genotypes(hdr, rec, &gt, &ngt);

        // Extract AT field from INFO
        char *at_str = nullptr;
        int nat = 0;
        nat = bcf_get_info_string(hdr, rec, "AT", &at_str, &nat);
        std::vector<std::string> path_list;

        std::string at_value(at_str);  // Convert C-string to C++ std::string
        free(at_str);  // Free HTSlib-allocated memory

        // Split by comma
        std::stringstream ss(at_value);
        std::string item;
        while (std::getline(ss, item, ',')) {
            path_list.push_back(item);
        }

        // Decompose snarl paths [vector std::string] into [vector vectorstoat::Edge_t]
        // paths : >123>213<234,>123<234,>123<234<345
        // list_paths_edge : [[Edge_t(123, 213),stoat::Edge_t(213, 234)], [...]]
        const std::vector<std::vector<stoat::Edge_t>> list_paths_edge = decompose_path_list_str(path_list);

        for (int i = 0; i < rec->n_sample; ++i) {
            int idex_path_allele_1 = bcf_gt_allele(gt[i * 2]);
            int idex_path_allele_2 = bcf_gt_allele(gt[i * 2 + 1]);
            size_t col_idx = i * 2;

            if (idex_path_allele_1 != -1) { // Handle missing genotypes
                for (const auto &edge_path_1 : list_paths_edge[idex_path_allele_1]) {
                    edge_matrix.push_matrix(edge_path_1, col_idx);
                }
            }

            if (idex_path_allele_2 != -1) { // Handle missing genotypes
                for (const auto &edge_path_2 : list_paths_edge[idex_path_allele_2]) {
                    edge_matrix.push_matrix(edge_path_2, col_idx + 1);
                }
            }
        }
        free(gt);

    } while ((bcf_read(ptr_vcf, hdr, rec) >= 0) && (chr == bcf_hdr_id2name(hdr, rec->rid)));

    edge_matrix.shrink();
    return std::make_tuple(ptr_vcf, hdr, rec);
}

// Decompose path stoat::Path_traversal_t to vectorstoat::Edge_t
std::vector<stoat::Edge_t> decompose_path_to_edges(const stoat::Path_traversal_t& list_paths) {
    std::vector<stoat::Edge_t> edges;
    const std::vector<stoat::Node_traversal_t>& list_nodes = list_paths.get_paths();
    size_t length_s = list_nodes.size();
    edges.reserve(length_s - 1); // Reserve memory

    for (size_t i = 0; i < length_s - 1; ++i) {
        edges.emplace_back(list_nodes[i], list_nodes[i + 1]);
    }

    return edges;
}

// Decompose path std::string to vectorstoat::Edge_t
std::vector<stoat::Edge_t> decompose_path_str_to_edge(const std::string& s) {
    std::vector<stoat::Edge_t> edges;
    std::vector<stoat::Node_traversal_t> nodes;

    size_t i = 0;
    while (i < s.size()) {
        if (s[i] == '>' || s[i] == '<') {
            bool is_rev = (s[i] == '<');
            ++i;

            size_t node_id = 0;
            while (i < s.size() && isdigit(s[i])) {
                node_id = node_id * 10 + (s[i] - '0');
                ++i;
            }
            nodes.emplace_back(node_id, is_rev);
        } else {
            ++i; // Skip invalid characters
        }
    }

    for (size_t j = 0; j + 1 < nodes.size(); ++j) {
        edges.emplace_back(nodes[j], nodes[j + 1]);
    }

    return edges;
}

// Decompose a list of paths str into a vector ofstoat::Edge_t
const std::vector<std::vector<stoat::Edge_t>> decompose_path_list_str(const std::vector<std::string>& list_paths) {
    std::vector<std::vector<stoat::Edge_t>> paths_snarl;
    for (const auto& path : list_paths) {
        paths_snarl.push_back(decompose_path_str_to_edge(path));
    }
    return paths_snarl;
}

// Function to identify the path in the edge matrix
std::vector<size_t> identify_path(
    const std::vector<stoat::Edge_t>& list_edge_path,
    const stoat_vcf::EdgeBySampleMatrix& edge_matrix,
    const size_t num_cols) {

    std::vector<size_t> rows_to_check;
    rows_to_check.reserve(list_edge_path.size());

    // Map snarl names to row indices
    for (const stoat::Edge_t& edge : list_edge_path) {
        const auto& [node_id_1, node_id_2] = edge.print_pair_edge(); // Convertstoat::Edge_t to std::pair<size_t, size_t>
        
        // Skip if snarl contains '*' (here * == 0) aka complex path
        if (node_id_1 == 0 || node_id_2 == 0) {
            continue;
        }
        size_t row_index = edge_matrix.find_edge(edge);
        if (row_index != std::numeric_limits<size_t>::max()) {
            rows_to_check.push_back(row_index);
        } else {
            return {}; // If any snarl isn't found, abort early
        }
    }

    std::vector<size_t> idx_srr_save;
    idx_srr_save.reserve(num_cols);

    // Loop columns first (better cache locality if matrix is column-major or similar)
    for (size_t col = 0; col < num_cols; ++col) {
        bool all_ones = true;
        for (size_t row : rows_to_check) {
            if (!edge_matrix(row, col)) {
                all_ones = false;
                break;
            }
        }
        if (all_ones) {
            idx_srr_save.push_back(static_cast<int>(col));
        }
    }
    return idx_srr_save;
}

void BinarySnarlAnalyzer::analyze_and_write_snarl(
    const stoat::Snarl_data_t& snarl_data_s, const std::string& chr, std::ofstream& outf) {

    std::ostringstream oss;

    for (size_t i = 0; i < snarl_data_s.type_variants.size(); ++i) {
        if (i != 0) oss << ",";
        oss << snarl_data_s.type_variants[i];
    }

    std::string type_var_str = oss.str();

    size_t paths_number = snarl_data_s.snarl_paths.size();
    std::vector<size_t> g0(paths_number, 0);
    std::vector<size_t> g1(paths_number, 0);

    size_t total_sum = stoat_vcf::create_binary_table(g0, g1, binary_phenotype, snarl_data_s.snarl_paths, paths_number, list_samples.size(), edge_matrix);
    remove_empty_columns_binary_table(g0, g1);
    bool filtration = filtration_binary_table(g0, g1, total_sum, min_individuals, min_haplotypes, maf_threshold);

    // Binary analysis single test
    if (!filtration) { // good table
        auto group_paths = format_group_paths(g0, g1);

        auto [fastfisher_p_value, chi2_p_value] = fk.fisher_khi2(g0, g1);

        # pragma omp critical (outf) 
        {
            stoat::write_binary(outf, chr, snarl_data_s, type_var_str, fastfisher_p_value, chi2_p_value, "",  group_paths);
        }
    }
}

void BinaryCovarSnarlAnalyzer::analyze_and_write_snarl( 
    const stoat::Snarl_data_t& snarl_data_s, const std::string& chr, std::ofstream& outf) {

    std::ostringstream oss;

    for (size_t i = 0; i < snarl_data_s.type_variants.size(); ++i) {
        if (i != 0) oss << ",";
        oss << snarl_data_s.type_variants[i];
    }

    std::string type_var_str = oss.str();

    auto [df, phenotype_filtered, allele_paths] = create_quantitative_table(list_samples.size(), snarl_data_s.snarl_paths, binary_phenotype, edge_matrix);
    remove_empty_columns_quantitative_table(df);
    
    bool filtration = filtration_quantitative_table(df, min_individuals, min_haplotypes, maf_threshold);
    remove_last_columns_quantitative_table(df);

    if (!filtration) { // filtred snarl
        // logistic regression with covariates if not empty
        const auto& [p_value, beta, se, r2] = lr.logistic_regression(df, phenotype_filtered, covariate);

        // Plot regression table
        if (table_threshold != -1 && stoat::isPValueSignificant(table_threshold, p_value)) {
            std::string variant_file_name = regression_dir + "/" + stoat::pairToString(snarl_data_s.snarl_ids) + ".tsv";
            stoat::writeSignificantTableToTSV(df,stoat::stringToVector<std::string>(stoat::vectorPathToString(snarl_data_s.snarl_paths)), edge_matrix.sampleNames, variant_file_name);
        }
        # pragma omp critical (outf) 
        {
            stoat::write_binary_covar(outf, chr, snarl_data_s, type_var_str, p_value, "", r2, beta, se, allele_paths);
        }
    }
}

// Quantitative Table Generation
void QuantitativeSnarlAnalyzer::analyze_and_write_snarl(
    const stoat::Snarl_data_t& snarl_data_s, const std::string& chr, std::ofstream& outf) {

    auto [df, phenotype_filtered, allele_paths] = create_quantitative_table(list_samples.size(), snarl_data_s.snarl_paths, quantitative_phenotype, edge_matrix);
    remove_empty_columns_quantitative_table(df);
    
    bool filtration = filtration_quantitative_table(df, min_individuals, min_haplotypes, maf_threshold);
    remove_last_columns_quantitative_table(df);

    // make a std::string separated by ',' from a vector of std::string
    std::ostringstream oss;
    for (size_t i = 0; i < snarl_data_s.type_variants.size(); ++i) {
        if (i != 0) oss << ","; // Add comma before all elements except the first
        oss << snarl_data_s.type_variants[i];
    }

    std::string type_var_str = oss.str();
    std::stringstream data;
    
    if (!filtration) { // filtred snarl
        auto [p_value, beta, se, r2] = lr.linear_regression(df, phenotype_filtered, covariate);
        
        if (table_threshold != -1 && stoat::isPValueSignificant(table_threshold, p_value)) {
            std::string variant_file_name = regression_dir + "/" + stoat::pairToString(snarl_data_s.snarl_ids) + ".tsv";
            stoat::writeSignificantTableToTSV(df,stoat::stringToVector<std::string>(stoat::vectorPathToString(snarl_data_s.snarl_paths)), edge_matrix.sampleNames, variant_file_name);
        }
        
        #pragma omp critical (outf)
        {
            stoat::write_quantitative(outf, chr, snarl_data_s, type_var_str, p_value, "", r2, beta, se, allele_paths);
        }
    }
}

// Identify genes index that will be tested for this snarl by matching position
// eqtl : <gene_name, gene_expression, start_pos, end_pos>
std::vector<size_t> found_gene_snarl(
    const std::vector<Qtl_data>& gene_position, 
    const size_t& start_pos, 
    const size_t& end_pos,
    const size_t& windows_gene_threshold) {

    std::vector<size_t> gene_index;
    size_t start_pos_threshold = (start_pos > windows_gene_threshold) ? start_pos - windows_gene_threshold : 0;
    size_t end_pos_threshold = end_pos + windows_gene_threshold;

    for (size_t i = 0; i < gene_position.size(); ++i) {
        size_t gene_start = gene_position[i].start_pos;
        size_t gene_end = gene_position[i].end_pos;

        // Check if the gene overlaps with the snarl region
        if (!(gene_end < start_pos_threshold || gene_start > end_pos_threshold)) {
            gene_index.push_back(i);
        }
    }
    return gene_index;
}

void EQTLSnarlAnalyzer::analyze_and_write_snarl(
    const stoat::Snarl_data_t& snarl_data_s, const std::string& chr, std::ofstream& outf) {

    std::vector<size_t> list_gene_index = found_gene_snarl(eqtl_map.at(chr), snarl_data_s.start_positions, snarl_data_s.end_positions, windows_gene_threshold);
    auto [df, index_filtered, allele_paths] = stoat_vcf::create_eqtl_table(list_samples.size(), snarl_data_s.snarl_paths, edge_matrix);
    remove_empty_columns_quantitative_table(df);
    
    bool filtration = filtration_quantitative_table(df, min_individuals, min_haplotypes, maf_threshold);
    remove_last_columns_quantitative_table(df);

    for (size_t i = 0; i < list_gene_index.size(); ++i) {
        size_t gene_idx = list_gene_index[i];
        std::string gene_name = eqtl_map.at(chr)[gene_idx].geneName;
        std::vector<double> gene_expression = eqtl_map.at(chr)[gene_idx].sampleExpresion;
        stoat::retain_indices(gene_expression, index_filtered);

        // make a std::string separated by ',' from a vector of std::string
        std::ostringstream oss;
        for (size_t i = 0; i < snarl_data_s.type_variants.size(); ++i) {
            if (i != 0) oss << ","; // Add comma before all elements except the first
            oss << snarl_data_s.type_variants[i];
        }

        std::string type_var_str = oss.str();
        std::stringstream data;

        if (!filtration) { // filtred snarl
            auto [p_value, beta, se, r2] = lr.linear_regression(df, gene_expression, covariate);

            if (table_threshold != -1 && stoat::isPValueSignificant(table_threshold, p_value)) {
                std::string variant_file_name = regression_dir + "/" + stoat::pairToString(snarl_data_s.snarl_ids) + ".tsv";
                stoat::writeSignificantTableToTSV(df,stoat::stringToVector<std::string>(stoat::vectorPathToString(snarl_data_s.snarl_paths)), edge_matrix.sampleNames, variant_file_name);
            }

            #pragma omp critical (outf)
            {
                stoat::write_eqtl(outf, chr, snarl_data_s, type_var_str, gene_name, p_value, "", r2, beta, se, allele_paths);
            }
        }
    }
}

// Return true when snarl must be filtered and false if not
bool filtration_quantitative_table(
    const std::vector<std::vector<double>>& df,
    const size_t& min_individuals,
    const size_t& min_haplotypes, 
    const double& maf_threshold) {
    
    // number of path < 2 OR not enougth individuals
    if (df.empty() || df[0].size() < 2 || df.size() < min_individuals) {
        return true; // Not enough data → filter out
    }

    size_t numPaths = df[0].size();
    std::vector<double> table(numPaths, 0.0);
    double totalSum = 0.0;

    // Compute column sums and total sum
    for (const auto& row : df) {
        for (size_t i = 0; i < numPaths; ++i) {
            table[i] += row[i];
            totalSum += row[i];
        }
    }

    if (totalSum < min_haplotypes) { // not good because 0.5 only in row can append
        return true;
    }

    int count_above_threshold = 0;
    for (size_t i = 0; i < numPaths; ++i) {
        double freq = table[i] / totalSum;
        double maf = std::min(freq, 1.0 - freq);
        if (maf > maf_threshold) {
            ++count_above_threshold;
        }
    }

    return count_above_threshold < 2;
}

void remove_empty_columns_quantitative_table(
    std::vector<std::vector<double>>& df) {

    if (df.empty()) return;

    size_t num_rows = df.size();
    size_t num_cols = df[0].size();

    // Identify non-empty columns
    std::vector<bool> keep_column(num_cols, false);

    for (size_t col = 0; col < num_cols; ++col) {
        for (size_t row = 0; row < num_rows; ++row) {
            double val = df[row][col];
            if (val != 0.0 && !std::isnan(val)) {
                keep_column[col] = true;
                break;
            }
        }
    }

    // Create filtered df
    std::vector<std::vector<double>> df_filtered;
    df_filtered.reserve(num_rows);

    for (size_t row = 0; row < num_rows; ++row) {
        std::vector<double> new_row;
        for (size_t col = 0; col < num_cols; ++col) {
            if (keep_column[col]) {
                new_row.push_back(df[row][col]);
            }
        }
        df_filtered.push_back(std::move(new_row));
    }

    // Replace original df with filtered one
    df = std::move(df_filtered);
}

void remove_last_columns_quantitative_table(std::vector<std::vector<double>>& df) {
    if (df.empty() || df[0].empty()) return;

    for (auto& row : df) {
        if (!row.empty()) {
            row.pop_back(); // Remove last column from each row
        }
    }
}

void remove_empty_columns_binary_table(
    std::vector<size_t>& g0, 
    std::vector<size_t>& g1) {

    std::vector<size_t> g0_filtered;
    std::vector<size_t> g1_filtered;

    for (size_t i = 0; i < g0.size(); ++i) {
        if (g0[i] + g1[i] != 0) {
            g0_filtered.push_back(g0[i]);
            g1_filtered.push_back(g1[i]);
        }
    }

    g0 = std::move(g0_filtered);
    g1 = std::move(g1_filtered);
}

// true : filtration on; false : no filtration
bool filtration_binary_table(
    std::vector<size_t>& g0, 
    std::vector<size_t>& g1,
    const size_t& totalSum, 
    const size_t& min_individuals,
    const size_t& min_haplotypes,
    const double& maf_threshold) {

    // Not enougth individuals OR not enougth haplotypes OR number of paths < 2
    if (totalSum/2 < min_individuals || totalSum < min_haplotypes || g0.size() < 2) { // not good don't take account on 1 allele only in a nested snarl
        return true; // Empty or invalid input → filter
    }

    int count_above_threshold = 0;

    for (size_t i = 0; i < g0.size(); ++i) {
        size_t columnSum = g0[i] + g1[i];

        double freq1 = static_cast<double>(g1[i]) / columnSum;
        double maf = std::min(freq1, 1.0 - freq1);

        if (maf > maf_threshold) {
            ++count_above_threshold;
        }
    }

    return count_above_threshold < 2; // Keep if at least two MAFs path > MAF threshold
}

} // end namespace stoat
