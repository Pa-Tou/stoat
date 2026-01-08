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
#include "types_and_structs.hpp"
#include "quantitative_table.hpp"
#include "snarl_data_collection.hpp"
#include "utils.hpp"
#include "log.hpp"


namespace stoat_vcf {

class SnarlAnalyzer {
public:
    SnarlAnalyzer(
        const stoat::SnarlDataCollection& snarl_collection,
        const std::unordered_set<std::string>& ref_chrs,
        const std::vector<std::string>& list_samples, 
        const std::vector<std::vector<double>>& covariate,
        const double maf_threshold,
        const size_t min_individuals);

    ~SnarlAnalyzer()=default;

    /// Go through the vcf by chromosome, parse it to get a matrix of genotypes (either binary, quantitative, or eqtl, depending on the phenotype type),
    /// then test the association and write the output (also depending on the phenotype type).
    /// This calls write_header() to write the appropriate output header and test_and_write_snarl() for each snarl
    void genotype_test_snarls_by_chr_from_vcf(
        htsFile* &ptr_vcf, 
        bcf_hdr_t* &hdr, 
        bcf1_t* &rec,
        const std::string output_dir);

    /// Update the EdgeBySampleMatrix representing the genotypes in a vcf and the pointers to the vcf but advanced to the end of the chromosome?
    std::tuple<htsFile*, bcf_hdr_t*, bcf1_t*> make_edge_matrix(
        htsFile *ptr_vcf, 
        bcf_hdr_t *hdr, 
        bcf1_t *rec, 
        std::string &chr);

    // get the type of phenotype, in case we do specific things outside of this class (although we should try to avoid it)
    stoat::phenotype_type_t get_phenotype_type() const;
    
    /// For the given snarl, get the genotypes and test the snarl, then write results to outf
    virtual bool test_and_write_snarl(const stoat::snarl_info_t& snarl_data, const std::string chr, std::ofstream& outf) = 0;

    /// Write the header of the output tsv file
    /// This should ideally call a write_header() function from writer.hpp to keep things consistent
    virtual void write_header(std::ofstream& outf);
    
//////////////// Private data members
protected:
    
    // a collection of all snarls
    const stoat::SnarlDataCollection& snarl_collection;

    // The reference path names
    const std::unordered_set<std::string>& ref_chrs;

    // A list of sample names
    const std::vector<std::string>& list_samples;

    // Covariate matrix
    const std::vector<std::vector<double>>& covariate;

    // save the type of phenotype for that SnarlAnalyzer (e.g. BINARY, QUANTITATIVE)
    stoat::phenotype_type_t phenotype_type;

    // Matrix of edges in each sample/haplotype
    // This generally is a per-chromosome or per-chunk matrix, so it must be updated for each new chunk being analyzed 
    EdgeBySampleMatrix edge_matrix;

    // threshold used to filter snarls 
    const double maf_threshold; 
    const size_t min_individuals;
};

class BinarySnarlAnalyzer : public SnarlAnalyzer {

public:
    
    BinarySnarlAnalyzer(
        const stoat::SnarlDataCollection& snarl_collection,
        const std::unordered_set<std::string>& ref_chrs,
        const std::vector<std::string>& list_samples, 
        const double maf_threshold,
        const std::vector<bool>& binary_phenotype,
        const size_t min_individuals);

    bool test_and_write_snarl(const stoat::snarl_info_t& snarl_data, const std::string chr, std::ofstream& outf);

/////////////////// Private data members
protected:

    const std::vector<bool>& binary_phenotype;
    stoat::FisherChi2 fchi;
};

class BinaryCovarSnarlAnalyzer : public SnarlAnalyzer {

public:
    
    BinaryCovarSnarlAnalyzer(
        const stoat::SnarlDataCollection& snarl_collection,
        const std::unordered_set<std::string>& ref_chrs,
        const std::vector<std::string>& list_samples, 
        const std::vector<std::vector<double>>& covariate, 
        const double maf_threshold, 
        const std::vector<bool>& binary_phenotype,
        const size_t min_individuals);

    bool test_and_write_snarl(const stoat::snarl_info_t& snarl_data, const std::string chr, std::ofstream& outf);

/////////////////// Private data members
protected:

    const std::vector<bool>& binary_phenotype;
    stoat::LogisticRegression lr;
};

class QuantitativeSnarlAnalyzer : public SnarlAnalyzer {

public:
    
    QuantitativeSnarlAnalyzer(
        const stoat::SnarlDataCollection& snarl_collection, 
        const std::unordered_set<std::string>& ref_chrs,
        const std::vector<std::string>& list_samples, 
        const std::vector<std::vector<double>>& covariate, 
        const double maf_threshold, 
        const std::vector<double>& quantitative_phenotype,
        const size_t min_individuals);

    bool test_and_write_snarl(const stoat::snarl_info_t& snarl_data, const std::string chr, std::ofstream& outf) ;

/////////////////// Private data members
protected:

    const std::vector<double>& quantitative_phenotype;
    stoat::LinearRegression lr;
};

class EQTLSnarlAnalyzer : public SnarlAnalyzer {

public:
    
    EQTLSnarlAnalyzer(
        const stoat::SnarlDataCollection& snarl_collection, 
        const std::unordered_set<std::string>& ref_chrs,
        const std::vector<std::string>& list_samples, 
        const std::vector<std::vector<double>>& covariate, 
        const double maf_threshold, 
        const std::unordered_map<std::string, std::vector<Qtl_data>>& eqtl_map,
        const size_t windows_gene_threshold,
        const size_t min_individuals);

    bool test_and_write_snarl(const stoat::snarl_info_t& snarl_data, const std::string chr, std::ofstream& outf);

/////////////////// Private data members
protected:

    // TODO idk what these are 
    // Maps something to something else?
    // Matis ans: eqtl_map is an {chr name : std::vector<Qtl_data>}
    // is organise like that in the first place to optimize edge_matrix / eqtl linking
    // but now we can just use std::vector<Qtl_data> because we already know the chr
    // that we gonna use
    const std::unordered_map<std::string, std::vector<Qtl_data>>& eqtl_map;
    const size_t windows_gene_threshold;
    stoat::LinearRegression lr;
};

// find which genes are close enough to the specified position
std::vector<size_t> get_genes_around_pos(
    const std::vector<Qtl_data>& gene_position, 
    const size_t start_pos, 
    const size_t end_pos,
    const size_t max_distance);

// Decompose path std::string to vectorstoat::edge_t
std::vector<stoat::edge_t> decompose_path_str_to_edge(const std::string s);

} //end stoat namespace

#endif
