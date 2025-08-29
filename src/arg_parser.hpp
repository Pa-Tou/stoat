#ifndef ARG_PARSER_HPP
#define ARG_PARSER_HPP

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <sstream>
#include <vector>
#include <tuple>
#include <map>
#include <algorithm>
#include <limits>
#include <filesystem>
#include <stdexcept>
#include <regex>

#include <Eigen/Dense>
#include <htslib/vcf.h>
#include <htslib/hts.h>

#include "log.hpp"
#include "snarl_data_t.hpp"

using namespace std;

namespace stoat_vcf {

struct KinshipMatrix {
    std::vector<std::string> ids;
    std::vector<std::vector<double>> matrix;

    // Default constructor
    KinshipMatrix() = default;

    // Parameterized constructor
    KinshipMatrix(const std::vector<std::string>& ids,
                  const std::vector<std::vector<double>>& matrix)
        : ids(ids), matrix(matrix) {}

    void parseKinshipMatrix(const std::string& filename);

    const bool empty() const {
        return ids.empty() || matrix.empty();
    } 
};

struct Qtl_data {
    std::string geneName;
    std::vector<double> sampleExpresion;
    size_t start_pos;
    size_t end_pos;

    Qtl_data(const std::string& geneName_,
        const std::vector<double>& sampleExpresion_,
        const size_t& start_pos_,
        const size_t& end_pos_) : 
        geneName(geneName_),
        start_pos(start_pos_),
        end_pos(end_pos_)
        {sampleExpresion = std::move(sampleExpresion_);}
};

std::unordered_map<std::string, std::vector<double>> parse_qtl_file(
    const std::string& filename, const std::vector<std::string>& list_samples);

std::unordered_map<std::string, std::tuple<std::string, size_t, size_t>> parse_gene_positions(
    const std::string& filename);

std::unordered_map<std::string, std::vector<stoat_vcf::Qtl_data>> parse_qtl_gene_file(
    const std::string& eqtl_path, 
    const std::string& gene_position_path, 
    const std::vector<std::string>& list_samples);

std::vector<std::vector<double>> parse_covariates(
    const std::string& filename, 
    const std::vector<std::string>& covar_names,
    const std::vector<std::string>& list_samples);

// Parse a binary phenotype file, formatted FID, IID, phenotype
// If list_samples is given, then it is const and the samples in the phenotype file will be checked against it.
// If list_samples is empty, then fill it in with the samples in the phenotype file 
std::vector<bool> parse_binary_pheno(
    const std::string& file_path,
    std::vector<std::string>& list_samples);

// Parses the phenotype file and returns a map with IID as keys and PHENO as float values.
std::vector<double> parse_quantitative_pheno(
    const std::string& file_path, 
    const std::vector<std::string>& list_samples);

std::tuple<htsFile*, bcf_hdr_t*, bcf1_t*> parse_vcf(const std::string& vcf_path);

std::tuple<std::vector<std::string>, htsFile*, bcf_hdr_t*, bcf1_t*> parseHeader(const std::string& vcf_path);

std::unordered_set<std::string> parse_chromosome_reference(const std::string& file_path);

template <typename T>
void check_match_samples(const std::unordered_map<std::string, T>& map, const std::vector<std::string>& keys);

void check_file(const std::string& file_path);

} //end stoat namespace

#endif
