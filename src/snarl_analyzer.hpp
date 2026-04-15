#ifndef snarl_analyzer_HPP
#define snarl_analyzer_HPP

#include <string>
#include <vector>
#include <unordered_set>

#include "snarl_data_collection.hpp"
#include "writer.hpp"
#include "feature_tables.hpp"

namespace stoat_vcf {

class SnarlAnalyzer {
public:
    SnarlAnalyzer(
        const stoat::SnarlDataCollection& snarl_collection,
        const stoat::CovariateTable& covariate,
        const double maf_threshold,
        const size_t min_individuals);

    SnarlAnalyzer(
        const stoat::SnarlDataCollection& snarl_collection,
        const double maf_threshold,
        const size_t min_individuals);

    ~SnarlAnalyzer()=default;

    // Go throught the snarls in a file and test the association with the phenotype.
    // Avoids loading the entire snarl collection with all the genotypes at once.
    void test_snarls_from_file(stoat::Reader& gt_reader, stoat::Writer& out_writer);
    
    // get the type of phenotype, in case we do specific things outside of this class (although we should try to avoid it)
    stoat::phenotype_type_t get_phenotype_type() const;
    
    /// For the given snarl, test and return formatted output string (empty if filtered)
    virtual std::string test_and_format_snarl(stoat::snarl_info_t& snarl_data, stoat::Writer& formatter) = 0;
    
//////////////// Private data members
protected:
    
    // a collection of all snarls
    const stoat::SnarlDataCollection& snarl_collection;

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
        const double maf_threshold,
        const stoat::BinaryPhenotypeTable& phenotype,
        const size_t min_individuals);

    std::string test_and_format_snarl(stoat::snarl_info_t& snarl_data, stoat::Writer& formatter);

protected:

    const stoat::BinaryPhenotypeTable& phenotype;
    stoat::FisherChi2 fchi;
};

class ExactBinarySnarlAnalyzer : public SnarlAnalyzer {

public:
    
    ExactBinarySnarlAnalyzer(
        const stoat::SnarlDataCollection& snarl_collection,
        const double maf_threshold,
        const stoat::BinaryPhenotypeTable& phenotype,
        const size_t min_individuals);

    std::string test_and_format_snarl(stoat::snarl_info_t& snarl_data, stoat::Writer& formatter);

protected:
    std::pair<std::set<std::string>, std::set<std::string>> sample_sets;
    
};

class BinaryCovarSnarlAnalyzer : public SnarlAnalyzer {

public:
    
    BinaryCovarSnarlAnalyzer(
        const stoat::SnarlDataCollection& snarl_collection,
        const stoat::CovariateTable& covariate,
        const double maf_threshold,
        const stoat::BinaryPhenotypeTable& phenotype,
        const size_t min_individuals);

    std::string test_and_format_snarl(stoat::snarl_info_t& snarl_data, stoat::Writer& formatter);

/////////////////// Private data members
protected:

    const stoat::BinaryPhenotypeTable& phenotype;
    stoat::LogisticRegression lr;
};

class QuantitativeSnarlAnalyzer : public SnarlAnalyzer {

public:
    
    QuantitativeSnarlAnalyzer(
        const stoat::SnarlDataCollection& snarl_collection,
        const stoat::CovariateTable& covariate,
        const double maf_threshold,
        const stoat::QuantitativePhenotypeTable& phenotype,
        const size_t min_individuals);

    std::string test_and_format_snarl(stoat::snarl_info_t& snarl_data, stoat::Writer& formatter);

/////////////////// Private data members
protected:

    const stoat::QuantitativePhenotypeTable& phenotype;
    stoat::LinearRegression lr;
};

class EQTLSnarlAnalyzer : public SnarlAnalyzer {

public:
    
    EQTLSnarlAnalyzer(const stoat::SnarlDataCollection& snarl_collection,
                      const stoat::CovariateTable& covariate,
                      const double maf_threshold,
                      const stoat::GeneExpressionTable& gene_expression,
                      const size_t max_gene_dist,
                      const size_t min_individuals);

    std::string test_and_format_snarl(stoat::snarl_info_t& snarl_data, stoat::Writer& formatter);
    
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
