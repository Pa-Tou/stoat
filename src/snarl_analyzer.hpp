#ifndef snarl_analyzer_HPP
#define snarl_analyzer_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <unordered_set>
#include <sstream>
#include <cstdlib>
#include <utility>
#include <iostream>
#include <thread>
#include <mutex>
#include <future>
#include <chrono>
#include <htslib/vcf.h>
#include <htslib/hts.h>

#include "arg_parser.hpp"
#include "stats_test.hpp"
#include "matrix.hpp"
#include "snarl_data_t.hpp"
#include "binary_table.hpp"
#include "quantitative_table.hpp"
#include "utils.hpp"
#include "log.hpp"

using namespace std;

namespace stoat_vcf {

class SnarlAnalyzer {
public:
    SnarlAnalyzer(
        const std::unordered_map<std::string, std::vector<stoat::Snarl_data_t>>& chr_to_snarl_data,
        EdgeBySampleMatrix& edge_matrix,
        const std::vector<std::string>& list_samples, 
        const std::vector<std::vector<double>>& covariate,
        const double& maf_threshold,
        const double& table_threshold,
        const size_t& min_individuals,
        const size_t& min_haplotypes,
        const std::string& regression_dir);

    ~SnarlAnalyzer()=default;

    /// Go through the vcf by chromosome, parse it to get a matrix of genotypes (either binary, quantitative, or eqtl, depending on the phenotype type),
    /// then write the output (also depending on the phenotype type).
    /// This calls write_header() to write the appropriate output header and analyze_and_write_snarl() for each snarl
    void process_snarls_by_chromosome_chunk(
        htsFile* &ptr_vcf, 
        bcf_hdr_t* &hdr, 
        bcf1_t* &rec,
        const std::string& output_filename);

    /// Update the EdgeBySampleMatrix representing the genotypes in a vcf and the pointers to the vcf but advanced to the end of the chromosome?
    std::tuple<htsFile*, bcf_hdr_t*, bcf1_t*> make_edge_matrix(
        htsFile *ptr_vcf, 
        bcf_hdr_t *hdr, 
        bcf1_t *rec, 
        std::string &chr, 
        size_t &num_paths_ch);

    /// For the given snarl, analyze the snarl and write it to outf
    virtual bool analyze_and_write_snarl(const stoat::Snarl_data_t& snarl_data, const std::string& chr, std::ofstream& outf) = 0;

    /// Write the header of the output tsv file
    /// This should ideally call a write_header() function from writer.hpp to keep things consistent
    virtual void write_header(std::ofstream& outf) = 0;

//////////////// Private data members
protected:
    
    // Map chromosome name to a vector of snarl_data_t
    const std::unordered_map<std::string, std::vector<stoat::Snarl_data_t>>& chr_to_snarl_data;

    // A list of sample names
    const std::vector<std::string>& list_samples;

    // Covariate matrix
    const std::vector<std::vector<double>>& covariate;

    // Matrix of edges in each sample/haplotype
    // This generally is a per-chromosome or per-chunk matrix, so it must be updated for each new chunk being analyzed 
    EdgeBySampleMatrix& edge_matrix;
    const double& maf_threshold; 
    const double& table_threshold;
    const size_t& min_individuals;
    const size_t& min_haplotypes;
    const std::string& regression_dir;
    std::ofstream* outf;
};

class BinarySnarlAnalyzer : public SnarlAnalyzer {

public:
    
    BinarySnarlAnalyzer(
        const std::unordered_map<std::string, std::vector<stoat::Snarl_data_t>>& chr_to_snarl_data,
        EdgeBySampleMatrix& edge_matrix,
        const std::vector<std::string>& list_samples, 
        const double& maf_threshold,
        const double& table_threshold,
        const std::vector<bool>& binary_phenotype,
        const size_t& min_individuals,
        const size_t& min_haplotypes,
        const std::string& regression_dir);

    bool analyze_and_write_snarl(const stoat::Snarl_data_t& snarl_data, const std::string& chr, std::ofstream& outf);

    void write_header(std::ofstream &outf);

/////////////////// Private data members
protected:

    const std::vector<bool>& binary_phenotype;
    FisherKhi2 fk;
};

class BinaryCovarSnarlAnalyzer : public SnarlAnalyzer {

public:
    
    BinaryCovarSnarlAnalyzer(
        const std::unordered_map<std::string, std::vector<stoat::Snarl_data_t>>& chr_to_snarl_data,
        EdgeBySampleMatrix& edge_matrix,
        const std::vector<std::string>& list_samples, 
        const std::vector<std::vector<double>>& covariate, 
        const double& maf_threshold, 
        const double& table_threshold,
        const std::vector<bool>& binary_phenotype,
        const size_t& min_individuals,
        const size_t& min_haplotypes,
        const std::string& regression_dir);

    bool analyze_and_write_snarl(const stoat::Snarl_data_t& snarl_data, const std::string& chr, std::ofstream& outf);

    void write_header(std::ofstream &outf);

/////////////////// Private data members
protected:

    const std::vector<bool>& binary_phenotype;
    LogisticRegression lr;
};

class QuantitativeSnarlAnalyzer : public SnarlAnalyzer {

public:
    
    QuantitativeSnarlAnalyzer(
        const std::unordered_map<std::string, std::vector<stoat::Snarl_data_t>>& chr_to_snarl_data, 
        EdgeBySampleMatrix& edge_matrix,
        const std::vector<std::string>& list_samples, 
        const std::vector<std::vector<double>>& covariate, 
        const double& maf_threshold, 
        const double& table_threshold,
        const std::vector<double>& quantitative_phenotype,
        const size_t& min_individuals,
        const size_t& min_haplotypes,
        const std::string& regression_dir);

    bool analyze_and_write_snarl(const stoat::Snarl_data_t& snarl_data, const std::string& chr, std::ofstream& outf) ;

    void write_header(std::ofstream &outf);

/////////////////// Private data members
protected:

    const std::vector<double>& quantitative_phenotype;
    LinearRegression lr;
};

class EQTLSnarlAnalyzer : public SnarlAnalyzer {

public:
    
    EQTLSnarlAnalyzer(
        const std::unordered_map<std::string, std::vector<stoat::Snarl_data_t>>& chr_to_snarl_data, 
        EdgeBySampleMatrix& edge_matrix,
        const std::vector<std::string>& list_samples, 
        const std::vector<std::vector<double>>& covariate, 
        const double& maf_threshold, 
        const double& table_threshold,
        const std::unordered_map<std::string, std::vector<Qtl_data>>& eqtl_map,
        const size_t& windows_gene_threshold,
        const size_t& min_individuals,
        const size_t& min_haplotypes,
        const std::string& regression_dir);

    bool analyze_and_write_snarl(const stoat::Snarl_data_t& snarl_data, const std::string& chr, std::ofstream& outf);

    void write_header(std::ofstream &outf);

/////////////////// Private data members
protected:

    // TODO idk what these are 
    // Maps something to something else?
    // Matis ans : eqtl_map is an {chr name : std::vector<Qtl_data>}
    // is organise like that in the first place to optimize edge_matrix / eqtl linking
    // but now we can just use std::vector<Qtl_data> because we already know the chr
    // that we gonna use
    const std::unordered_map<std::string, std::vector<Qtl_data>>& eqtl_map;
    const size_t& windows_gene_threshold;
    LinearRegression lr;
};

void combine_identical_columns_quantitative_table(
    std::vector<std::vector<double>>& df);

void remove_empty_columns_binary_table(
    std::vector<size_t>& g0, 
    std::vector<size_t>& g1);

void remove_empty_columns_quantitative_table(
    std::vector<std::vector<double>>& df);

void remove_last_columns_quantitative_table(
    std::vector<std::vector<double>>& df);

/// Return true if snarl must be filtered
bool filtration_quantitative_table(
    const std::vector<std::vector<double>>& df,
    const size_t& min_individuals,
    const size_t& min_haplotypes,
    const double& maf);

bool filtration_binary_table(
    std::vector<size_t>& g0, 
    std::vector<size_t>& g1,
    const size_t& totalSum, 
    const size_t& min_individuals,
    const size_t& min_haplotypes,
    const double& maf);

std::vector<size_t> found_gene_snarl(
    const std::vector<Qtl_data>& gene_position, 
    const size_t& start_pos, 
    const size_t& end_pos,
    const size_t& windows_gene_threshold);

// Decompose path stoat::Path_traversal_t to vectorstoat::Edge_t
std::vector<stoat::Edge_t> decompose_path_to_edges(const stoat::Path_traversal_t& s);

// Decompose a list of paths std::string into a vector ofstoat::Edge_t
const std::vector<std::vector<stoat::Edge_t>> decompose_path_list_str(const std::vector<std::string>& list_paths);

// Decompose path std::string to vectorstoat::Edge_t
std::vector<stoat::Edge_t> decompose_path_str_to_edge(const std::string& s);

/// Given a path through the snarl, a matrix of edges for each sample/haplotype, and the number of columns (samples/haplotypes),
/// return the columns for which all edges (rows) in the path are set, i.e. the haplotypes with the given path.
std::vector<size_t> identify_path(
    const std::vector<stoat::Edge_t>& list_edge_path,
    const EdgeBySampleMatrix& matrix,
    const size_t num_cols);

} //end stoat namespace

#endif
