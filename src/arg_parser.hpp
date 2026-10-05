#ifndef ARG_PARSER_HPP
#define ARG_PARSER_HPP

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <htslib/vcf.h>
#include <htslib/hts.h>

#include "log.hpp"
#include "types_and_structs.hpp"
#include "feature_tables.hpp"


namespace stoat_vcf {

// Predicts which statistical methods to used based on the values in
// the phenotype file (when a model is not specified by the user)
//   binary without covariate -> chi2
//   binary with covariate -> logistic
//   other -> linear
std::string methods_stats_prediction(const std::string& file_path, const bool& covariate);
 
// Parse a binary phenotype file and return a BinaryPhenotypeTable
// if the sample-to-index map provided is be empty, it will be filled with the samples in the file.
// otherwise the Table will use that sample-to-index map
std::unique_ptr<stoat::BinaryPhenotypeTable> parse_binary_pheno_table(const std::string& file_path, std::unordered_map<std::string, size_t>& sample_to_index);

// Parse a quantitative phenotype file and return a QuantitativePhenotypeTable
// if the sample-to-index map provided is be empty, it will be filled with the samples in the file.
// otherwise the Table will use that sample-to-index map
std::unique_ptr<stoat::QuantitativePhenotypeTable> parse_quantitative_pheno_table(const std::string& file_path, std::unordered_map<std::string, size_t>& sample_to_index);

// Parse a pair of gene expression and position files covariate file and return a GeneExpressionTable
// if the sample-to-index map provided is empty, it will be filled with the samples in the file.
// otherwise the Table will use that sample-to-index map
// the gene_to_index map is expected to be empty and will be filled
std::unique_ptr<stoat::GeneExpressionTable> parse_gene_expression_table(const std::string& gene_expression_path, const std::string& gene_position_path, std::unordered_map<std::string, size_t>& sample_to_index, std::unordered_map<std::string, size_t>& gene_to_index);

// Parse a covariate file and return a CovariateTable
// if the sample-to-index map provided is empty, it will be filled with the samples in the file.
// otherwise the Table will use that sample-to-index map
// specify the names of the covariables to use with the covar_to_index map
// the covar_to_index map can't be empty (although it could be useful to allow this when we want all covariables)
std::unique_ptr<stoat::CovariateTable> parse_covariate_table(const std::string& file_path, std::unordered_map<std::string, size_t>& sample_to_index, std::unordered_map<std::string, size_t>& covar_to_index);

std::tuple<htsFile*, bcf_hdr_t*, bcf1_t*> parse_vcf(const std::string& vcf_path);

std::tuple<std::vector<std::string>, htsFile*, bcf_hdr_t*, bcf1_t*> parseHeader(const std::string& vcf_path);

std::unordered_set<std::string> parse_chromosome_reference(const std::string& file_path);

template <typename T>
void check_match_samples(const std::unordered_map<std::string, T>& map, const std::vector<std::string>& keys);

void check_file(const std::string& file_path);

void check_methods(const std::string& method_name);

} //end stoat namespace

#endif
