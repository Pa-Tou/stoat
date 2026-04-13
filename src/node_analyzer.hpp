#ifndef NODE_ANALYZER_HPP
#define NODE_ANALYZER_HPP

#include <string>
#include <vector>
#include <unordered_set>

#include "node_data_collection.hpp"
#include "writer.hpp"
#include "feature_tables.hpp"

namespace stoat_vcf {

class NodeAnalyzer {
public:
    NodeAnalyzer(
        const stoat::NodeDataCollection& node_collection,
        const stoat::CovariateTable& covariate,
        const double maf_threshold,
        const size_t min_individuals);

    NodeAnalyzer(
        const stoat::NodeDataCollection& node_collection,
        const double maf_threshold,
        const size_t min_individuals);

    ~NodeAnalyzer()=default;

    // Go throught the nodes in a file and test the association with the phenotype.
    // Avoids loading the entire node collection with all the genotypes at once.
    void test_nodes_from_file(stoat::Reader& gt_reader, stoat::Writer& out_writer);
    
    // get the type of phenotype, in case we do specific things outside of this class (although we should try to avoid it)
    stoat::phenotype_type_t get_phenotype_type() const;
    
    /// For the given node, get the genotypes and test the node, then write results to outf
    virtual bool test_and_write_node(stoat::node_info_t& node_data, stoat::Writer& out_writer) = 0;
    
//////////////// Private data members
protected:
    
    // a collection of all nodes
    const stoat::NodeDataCollection& node_collection;

    // Covariate matrix
    const stoat::CovariateTable& covariate;

    // save the type of phenotype for that NodeAnalyzer (e.g. BINARY, QUANTITATIVE)
    stoat::phenotype_type_t phenotype_type;

    // threshold used to filter nodes 
    const double maf_threshold; 
    const size_t min_individuals;
};

class BinaryNodeAnalyzer : public NodeAnalyzer {

public:
    
    BinaryNodeAnalyzer(
        const stoat::NodeDataCollection& node_collection,
        const double maf_threshold,
        const stoat::BinaryPhenotypeTable& phenotype,
        const size_t min_individuals);

    bool test_and_write_node(stoat::node_info_t& node_data, stoat::Writer& out_writer);

protected:

    const stoat::BinaryPhenotypeTable& phenotype;
    stoat::FisherChi2 fchi;
};

class ExactBinaryNodeAnalyzer : public NodeAnalyzer {

public:
    
    ExactBinaryNodeAnalyzer(
        const stoat::NodeDataCollection& node_collection,
        const double maf_threshold,
        const stoat::BinaryPhenotypeTable& phenotype,
        const size_t min_individuals);

    bool test_and_write_node(stoat::node_info_t& node_data, stoat::Writer& out_writer);

protected:
    std::pair<std::set<std::string>, std::set<std::string>> sample_sets;
    
};

class BinaryCovarNodeAnalyzer : public NodeAnalyzer {

public:
    
    BinaryCovarNodeAnalyzer(
        const stoat::NodeDataCollection& node_collection,
        const stoat::CovariateTable& covariate,
        const double maf_threshold, 
        const stoat::BinaryPhenotypeTable& phenotype,
        const size_t min_individuals);

    bool test_and_write_node(stoat::node_info_t& node_data, stoat::Writer& out_writer);

/////////////////// Private data members
protected:

    const stoat::BinaryPhenotypeTable& phenotype;
    stoat::LogisticRegression lr;
};

class QuantitativeNodeAnalyzer : public NodeAnalyzer {

public:
    
    QuantitativeNodeAnalyzer(
        const stoat::NodeDataCollection& node_collection, 
        const stoat::CovariateTable& covariate,
        const double maf_threshold, 
        const stoat::QuantitativePhenotypeTable& phenotype,
        const size_t min_individuals);

    bool test_and_write_node(stoat::node_info_t& node_data, stoat::Writer& out_writer) ;

/////////////////// Private data members
protected:

    const stoat::QuantitativePhenotypeTable& phenotype;
    stoat::LinearRegression lr;
};

class EQTLNodeAnalyzer : public NodeAnalyzer {

public:
    
    EQTLNodeAnalyzer(const stoat::NodeDataCollection& node_collection, 
                      const stoat::CovariateTable& covariate,
                      const double maf_threshold, 
                      const stoat::GeneExpressionTable& gene_expression,
                      const size_t max_gene_dist,
                      const size_t min_individuals);
    
    bool test_and_write_node(stoat::node_info_t& node_data, stoat::Writer& out_writer);
    
protected:

    // the Table with gene expression and positions
    const stoat::GeneExpressionTable& gene_expression;

    // the maximum distance allowed between a node and a gene
    const size_t max_gene_dist;

    // an object from the linear regression helper class
    stoat::LinearRegression lr;
};

// Decompose path std::string to vectorstoat::edge_t
std::vector<stoat::edge_t> decompose_path_str_to_edge(const std::string s);

} //end stoat namespace

#endif
