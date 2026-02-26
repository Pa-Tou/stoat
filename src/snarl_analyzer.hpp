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
        const stoat::CovariateTable& covariate,
        const double maf_threshold,
        const size_t min_individuals);

    SnarlAnalyzer(
        const stoat::SnarlDataCollection& snarl_collection,
        const std::unordered_set<std::string>& ref_chrs,
        const double maf_threshold,
        const size_t min_individuals);

    ~SnarlAnalyzer()=default;

    /// Go through the SnarlCollection (which contains genotypes) by chromosome
    /// then test the association with phenotype and write the output (also depending on the phenotype type).
    /// This calls write_header() to write the appropriate output header and test_and_write_snarl() for each snarl
    void genotype_test_snarls_by_chr(const std::string output_dir);

    // Go throught the snarls in a file and test the association with the phenotype.
    // Avoids loading the entire snarl collection with all the genotypes at once.
    void test_snarls_from_file(const std::string snarl_genotype_path, const std::string output_dir);
    
    // get the type of phenotype, in case we do specific things outside of this class (although we should try to avoid it)
    stoat::phenotype_type_t get_phenotype_type() const;
    
    /// For the given snarl, get the genotypes and test the snarl, then write results to outf
    virtual bool test_and_write_snarl(stoat::snarl_info_t& snarl_data, std::ofstream& outf) = 0;

    /// Write the header of the output tsv file
    /// This should ideally call a write_header() function from writer.hpp to keep things consistent
    virtual void write_header(std::ofstream& outf);
    
//////////////// Private data members
protected:
    
    // a collection of all snarls
    const stoat::SnarlDataCollection& snarl_collection;

    // The reference path names
    const std::unordered_set<std::string>& ref_chrs;

    // Covariate matrix
    const stoat::CovariateTable& covariate;

    // save the type of phenotype for that SnarlAnalyzer (e.g. BINARY, QUANTITATIVE)
    stoat::phenotype_type_t phenotype_type;

    // threshold used to filter snarls 
    const double maf_threshold; 
    const size_t min_individuals;
};

class BinarySnarlAnalyzer : public SnarlAnalyzer {

public:
    
    BinarySnarlAnalyzer(
        const stoat::SnarlDataCollection& snarl_collection,
        const std::unordered_set<std::string>& ref_chrs,
        const double maf_threshold,
        const stoat::BinaryPhenotypeTable& phenotype,
        const size_t min_individuals);

    bool test_and_write_snarl(stoat::snarl_info_t& snarl_data, std::ofstream& outf);

protected:

    const stoat::BinaryPhenotypeTable& phenotype;
    stoat::FisherChi2 fchi;
};

class BinaryCovarSnarlAnalyzer : public SnarlAnalyzer {

public:
    
    BinaryCovarSnarlAnalyzer(
        const stoat::SnarlDataCollection& snarl_collection,
        const std::unordered_set<std::string>& ref_chrs,
        const stoat::CovariateTable& covariate,
        const double maf_threshold, 
        const stoat::BinaryPhenotypeTable& phenotype,
        const size_t min_individuals);

    bool test_and_write_snarl(stoat::snarl_info_t& snarl_data, std::ofstream& outf);

/////////////////// Private data members
protected:

    const stoat::BinaryPhenotypeTable& phenotype;
    stoat::LogisticRegression lr;
};

class QuantitativeSnarlAnalyzer : public SnarlAnalyzer {

public:
    
    QuantitativeSnarlAnalyzer(
        const stoat::SnarlDataCollection& snarl_collection, 
        const std::unordered_set<std::string>& ref_chrs,
        const stoat::CovariateTable& covariate,
        const double maf_threshold, 
        const stoat::QuantitativePhenotypeTable& phenotype,
        const size_t min_individuals);

    bool test_and_write_snarl(stoat::snarl_info_t& snarl_data, std::ofstream& outf) ;

/////////////////// Private data members
protected:

    const stoat::QuantitativePhenotypeTable& phenotype;
    stoat::LinearRegression lr;
};

class EQTLSnarlAnalyzer : public SnarlAnalyzer {

public:
    
    EQTLSnarlAnalyzer(const stoat::SnarlDataCollection& snarl_collection, 
                      const std::unordered_set<std::string>& ref_chrs,
                      const stoat::CovariateTable& covariate,
                      const double maf_threshold, 
                      const stoat::GeneExpressionTable& gene_expression,
                      const size_t max_gene_dist,
                      const size_t min_individuals);
    
    bool test_and_write_snarl(stoat::snarl_info_t& snarl_data, std::ofstream& outf);
    
protected:

    // the Table with gene expression and positions
    const stoat::GeneExpressionTable& gene_expression;

    // the maximum distance allowed between a snarl and a gene
    const size_t max_gene_dist;

    // an object from the linear regression helper class
    stoat::LinearRegression lr;
};

// Decompose path std::string to vectorstoat::edge_t
std::vector<stoat::edge_t> decompose_path_str_to_edge(const std::string s);

} //end stoat namespace

#endif
