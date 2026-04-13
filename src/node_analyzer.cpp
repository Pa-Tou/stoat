#include "node_analyzer.hpp"
#include "omp.h"
#include "log.hpp"

using namespace stoat;
namespace stoat_vcf {

    NodeAnalyzer::NodeAnalyzer(
        const NodeDataCollection &node_collection,
        const stoat::CovariateTable& covariate,
        const double maf_threshold,
        const size_t min_individuals) :
        node_collection(node_collection),
        covariate(covariate),
        maf_threshold(maf_threshold),
        min_individuals(min_individuals),
        phenotype_type(stoat::BINARY) {};

    NodeAnalyzer::NodeAnalyzer(
        const NodeDataCollection &node_collection,
        const double maf_threshold,
        const size_t min_individuals) :
        node_collection(node_collection),
        covariate(covariate),
        maf_threshold(maf_threshold),
        min_individuals(min_individuals),
        phenotype_type(stoat::BINARY) {};

    BinaryNodeAnalyzer::BinaryNodeAnalyzer(
        const NodeDataCollection &node_collection,
        const double maf_threshold,
        const stoat::BinaryPhenotypeTable& phenotype,
        const size_t min_individuals) :
        NodeAnalyzer(node_collection, maf_threshold, min_individuals),
        phenotype(phenotype), fchi() {
        phenotype_type = stoat::BINARY;
    };

    ExactBinaryNodeAnalyzer::ExactBinaryNodeAnalyzer(
        const NodeDataCollection &node_collection,
        const double maf_threshold,
        const stoat::BinaryPhenotypeTable& phenotype,
        const size_t min_individuals) :
        NodeAnalyzer(node_collection, maf_threshold, min_individuals) {
        phenotype_type = stoat::BINARY;

        // fill the sample sets
        for (std::string sample: phenotype.get_sample_names()) {
            if (phenotype.get_value_for_sample(sample)) {
                sample_sets.first.emplace(sample);
            } else {
                sample_sets.second.emplace(sample);
            }
        }
    };

    BinaryCovarNodeAnalyzer::BinaryCovarNodeAnalyzer(
        const NodeDataCollection &node_collection,
        const stoat::CovariateTable& covariate,
        const double maf_threshold,
        const stoat::BinaryPhenotypeTable& phenotype,
        const size_t min_individuals) :

        NodeAnalyzer(node_collection, covariate, maf_threshold, min_individuals),
        phenotype(phenotype), lr() {
        phenotype_type = stoat::BINARY_COVAR;
    };

    QuantitativeNodeAnalyzer::QuantitativeNodeAnalyzer(
        const NodeDataCollection &node_collection,
        const stoat::CovariateTable& covariate,
        const double maf_threshold,
        const stoat::QuantitativePhenotypeTable& phenotype,
        const size_t min_individuals) :

        NodeAnalyzer(node_collection, covariate, maf_threshold, min_individuals),
        phenotype(phenotype), lr() {
        phenotype_type = stoat::QUANTITATIVE;
    };

    EQTLNodeAnalyzer::EQTLNodeAnalyzer(
        const NodeDataCollection &node_collection,
        const stoat::CovariateTable& covariate,
        const double maf_threshold,
        const stoat::GeneExpressionTable& gene_expression,
        const size_t max_gene_dist,
        const size_t min_individuals) :
        NodeAnalyzer(node_collection, covariate, maf_threshold, min_individuals),
        gene_expression(gene_expression), max_gene_dist(max_gene_dist), lr() {
        phenotype_type = stoat::EQTL;
    };

    stoat::phenotype_type_t NodeAnalyzer::get_phenotype_type() const {
        return (phenotype_type);
    }

void NodeAnalyzer::test_nodes_from_file(stoat::Reader& gt_reader, stoat::Writer& out_writer) {

    // prepare node collection that will stream the nodes and open connection to the file
    stoat::NodeDataCollection node_collection_stream(node_collection.get_sample_to_index_copy());

    // Write the header of the output file
    out_writer.write_node_stoat_output_header(phenotype_type);
 
    // read each node and test it
    // JEAN parallelize here?
    size_t number_node_filtered = 0;
    node_collection_stream.for_each_node_in_file(gt_reader, [&](node_info_t& node_info) {
        bool filtered = test_and_write_node(node_info, out_writer);
        number_node_filtered += (filtered ? 1 : 0);
    });

    stoat::LOG_INFO("Total number of node filtered : " + std::to_string(number_node_filtered));
}

bool BinaryNodeAnalyzer::test_and_write_node(stoat::node_info_t &node_data, stoat::Writer& out_writer) {

    // link to the phenotype
    node_data.genotypes.link_to_binary_phenotype(phenotype);
    
    // remove non-variable allele, e.g. absent in both groups
    node_data.genotypes.remove_noncovered_samples();
    node_data.genotypes.remove_constant_predictors();

    // prepare an output objet and init to NA
    test_result_t test_res;
    test_res.pv = std::nan("");
    test_res.second_pv = std::nan("");

    // should we test this node?
    if (node_data.genotypes.passes_filters(maf_threshold, min_individuals)){

        // fill up the contingency table (one vector per group)
        std::vector<size_t> g0;
        std::vector<size_t> g1;
        node_data.genotypes.fill_contingency_table(g0, g1);

        // performs the test
        auto fc_res = fchi.fisher_chi2(g0, g1);
        test_res.pv = fc_res.first;
        test_res.second_pv = fc_res.second;
        test_res.group_paths = stoat::format_group_paths(g0, g1);
    } else {
        // JEAN what should we do here? returning NA for now
        stoat::LOG_DEBUG("Filtered: didn't pass the filters");
    }
    
    // test this node using Fisher exact test and chi-squared test
    if (std::isnan(test_res.pv) && std::isnan(test_res.second_pv)) {
        stoat::LOG_DEBUG("filtered: " + node_data.node.to_string());
        return true;
    }
    
    out_writer.write_node_binary(node_data, test_res);
    return false;
}

bool ExactBinaryNodeAnalyzer::test_and_write_node(stoat::node_info_t &node_data, stoat::Writer& out_writer) {
    // This test checks if all members of one of the phenotype groups has the same allele that no other sample has.

    // prepare an output objet and init to NA
    test_result_t test_res;
    test_res.pv = std::nan("");
    test_res.second_pv = std::nan("");
    test_res.group_paths = "NA";

    // From the genotype matrix, make sets of sample that have the same genotype and compare to the sets of phenotype groups.
    std::unordered_map<std::string, std::set<std::string>> genotype_to_sample_set;
    for (const sample_hap_t& sample_hap : node_data.all_sample_haplotypes) {
        std::string genotype_str = node_data.genotypes.get_genotype_as_string(sample_hap.sample);

        if (genotype_to_sample_set.count(genotype_str) == 0) {
            genotype_to_sample_set.emplace(genotype_str, std::set<std::string>());
        }

        genotype_to_sample_set.at(genotype_str).emplace(sample_hap.sample);
    }

    // If we only want to know if there is one allele that matches exactly one of the phenotype groups
    bool write_output = false;
    for (const auto& genotype_samples : genotype_to_sample_set) {
        const std::set<std::string>& partition = genotype_samples.second;

        // If one partition exactly matches one group we want, then all other partitions combined (including things not in the node) will match
        // the other.
        // TODO: This could be better but I don't think it's worth working on it yet
        if (partition == sample_sets.first || partition == sample_sets.second) {
            // There was an exact match, we'll want to write that one
            write_output = true;
            break;
        }
    }

    // skip if we don't want to write that one
    if (!write_output) {
        return true;
    }
    
    out_writer.write_node_binary(node_data, test_res);
    return false;
}

bool BinaryCovarNodeAnalyzer::test_and_write_node(stoat::node_info_t &node_data, stoat::Writer& out_writer) {

    // link the phenotype
    node_data.genotypes.link_to_binary_phenotype(phenotype);
    node_data.genotypes.link_to_covariates(covariate);

    // remove non-variable predictors, e.g. alleles absent in all samples
    node_data.genotypes.remove_noncovered_samples();
    node_data.genotypes.remove_constant_predictors();    

    // prepare an output objet and init to NA
    test_result_t test_res;
    test_res.pv = std::nan("");
    
    // should we test this node?
    if (node_data.genotypes.passes_filters(maf_threshold, min_individuals)){

        // add covariate with the number of alleles (if necessary) to correct for the parent node effect (or normalize)
        node_data.genotypes.add_total_allele_count_covariable();
    
        // before performing the regression, try to reduce potential colinearity
        node_data.genotypes.remove_duplicated_predictors();    
        node_data.genotypes.remove_one_allele();
    
        // prepare the matrices, fit the logistic model and test effect of alleles
        Eigen::MatrixXd X = node_data.genotypes.make_matrixXd_features();
        Eigen::VectorXd Y = node_data.genotypes.make_vectorxd_phenotype();
        test_res.pv = lr.logistic_regression(X, Y, node_data.genotypes.get_n_active_alleles());
    } else {
        stoat::LOG_DEBUG("Filtered: didn't pass the filters");
    }
    
    if (std::isnan(test_res.pv)) {
        stoat::LOG_DEBUG("filtered: " + node_data.node.to_string());
        return true;
    }
 
    out_writer.write_node_binary_covar(node_data, test_res);
    return false;
}

// Quantitative Table Generation
bool QuantitativeNodeAnalyzer::test_and_write_node(stoat::node_info_t &node_data, stoat::Writer& out_writer) {

    // link the phenotype
    node_data.genotypes.link_to_quantitative_phenotype(phenotype);
    node_data.genotypes.link_to_covariates(covariate);

    // remove non-variable allele, e.g. absent in both groups
    node_data.genotypes.remove_noncovered_samples();    
    node_data.genotypes.remove_constant_predictors();

    // prepare an output objet and init to NA
    test_result_t test_res;
    test_res.pv = std::nan("");

    // should we test this node?
    if (node_data.genotypes.passes_filters(maf_threshold, min_individuals)){

        // add covariate with the number of alleles (if necessary) to correct for the parent node effect (or normalize)
        node_data.genotypes.add_total_allele_count_covariable();

        // before performing the regression, try to reduce potential colinearity
        node_data.genotypes.remove_duplicated_predictors();
        node_data.genotypes.remove_one_allele();

        Eigen::MatrixXd X = node_data.genotypes.make_matrixXd_features();
        Eigen::VectorXd Y = node_data.genotypes.make_vectorxd_phenotype();

        test_res.pv = lr.linear_regression(X, Y, node_data.genotypes.get_n_active_alleles());
    } else {
        // JEAN what should we do here? returning NA for now
        stoat::LOG_DEBUG("Filtered: didn't pass the filters");
    }

    if (std::isnan(test_res.pv)) {
        stoat::LOG_DEBUG("filtered: " + node_data.node.to_string());
        return true;
    }
 
    out_writer.write_node_quantitative(node_data, test_res);
    return false;
}

bool EQTLNodeAnalyzer::test_and_write_node(stoat::node_info_t &node_data, stoat::Writer& out_writer) {

    // get genes near node
    std::vector<std::string> genes_near = gene_expression.get_genes_around_pos(node_data.ref_path, node_data.position, max_gene_dist);
    bool filtered = true;

    // test against the expression of each nearby gene
    for (std::string gene_name: genes_near) {

        // make a QuantitativePhenotypeTable for this gene
        QuantitativePhenotypeTable gene_phenotype(gene_expression.get_sample_to_index());
        for (std::string sample_name: gene_expression.get_sample_names()) {
            gene_phenotype.set_value_for_sample(sample_name, gene_expression.get_value_for_sample_and_feature(sample_name, gene_name));
        }

        // test the gene
        // reinitialize the genotype object (remove masks etc set when testing other genes)
        node_data.genotypes.clear();

        // link the phenotype
        node_data.genotypes.link_to_quantitative_phenotype(gene_phenotype);
        node_data.genotypes.link_to_covariates(covariate);

        // remove non-variable allele, e.g. absent in both groups
        node_data.genotypes.remove_noncovered_samples();    
        node_data.genotypes.remove_constant_predictors();

        // prepare an output objet and init to NA
        test_result_t test_res;
        test_res.pv = std::nan("");

        // should we test this node?
        if (node_data.genotypes.passes_filters(maf_threshold, min_individuals)){

            // add covariate with the number of alleles (if necessary) to correct for the parent node effect (or normalize)
            node_data.genotypes.add_total_allele_count_covariable();

            // before performing the regression, try to reduce potential colinearity
            node_data.genotypes.remove_duplicated_predictors();
            node_data.genotypes.remove_one_allele();
    
            Eigen::MatrixXd X = node_data.genotypes.make_matrixXd_features();
            Eigen::VectorXd Y = node_data.genotypes.make_vectorxd_phenotype();

            test_res.pv = lr.linear_regression(X, Y, node_data.genotypes.get_n_active_alleles());
        } else {
            // JEAN what should we do here? returning NA for now
            stoat::LOG_DEBUG("Filtered: didn't pass the filters");
        }

        if (std::isnan(test_res.pv)) {
            stoat::LOG_DEBUG("filtered: gene " + gene_name + ", " + node_data.node.to_string());
            continue;
        }
  
        out_writer.write_node_eqtl(node_data, gene_name, test_res);

        // at least this test was not filtered
        filtered = false;
    }

    // return if the node was filtered in all considered genes
    return filtered;
}

} // end namespace stoat
