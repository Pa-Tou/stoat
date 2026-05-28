#include "snarl_analyzer.hpp"
#include "omp.h"
#include "log.hpp"

using namespace stoat;
namespace stoat_vcf {

    SnarlAnalyzer::SnarlAnalyzer(
        const SnarlDataCollection &snarl_collection,
        const stoat::CovariateTable& covariate,
        const double maf_threshold,
        const size_t min_individuals) :
        snarl_collection(snarl_collection),
        covariate(covariate),
        maf_threshold(maf_threshold),
        min_individuals(min_individuals),
        phenotype_type(stoat::BINARY) {};

    SnarlAnalyzer::SnarlAnalyzer(
        const SnarlDataCollection &snarl_collection,
        const double maf_threshold,
        const size_t min_individuals) :
        snarl_collection(snarl_collection),
        covariate(covariate),
        maf_threshold(maf_threshold),
        min_individuals(min_individuals),
        phenotype_type(stoat::BINARY) {};

    BinarySnarlAnalyzer::BinarySnarlAnalyzer(
        const SnarlDataCollection &snarl_collection,
        const double maf_threshold,
        const stoat::BinaryPhenotypeTable& phenotype,
        const size_t min_individuals) :
        SnarlAnalyzer(snarl_collection, maf_threshold, min_individuals),
        phenotype(phenotype), fchi() {
        phenotype_type = stoat::BINARY;
    };

    ExactBinarySnarlAnalyzer::ExactBinarySnarlAnalyzer(
        const SnarlDataCollection &snarl_collection,
        const double maf_threshold,
        const stoat::BinaryPhenotypeTable& phenotype,
        const size_t min_individuals) :
        SnarlAnalyzer(snarl_collection, maf_threshold, min_individuals) {
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

    BinaryCovarSnarlAnalyzer::BinaryCovarSnarlAnalyzer(
        const SnarlDataCollection &snarl_collection,
        const stoat::CovariateTable& covariate,
        const double maf_threshold,
        const stoat::BinaryPhenotypeTable& phenotype,
        const size_t min_individuals) :

        SnarlAnalyzer(snarl_collection, covariate, maf_threshold, min_individuals),
        phenotype(phenotype), lr() {
        phenotype_type = stoat::BINARY_COVAR;
    };

    QuantitativeSnarlAnalyzer::QuantitativeSnarlAnalyzer(
        const SnarlDataCollection &snarl_collection,
        const stoat::CovariateTable& covariate,
        const double maf_threshold,
        const stoat::QuantitativePhenotypeTable& phenotype,
        const size_t min_individuals) :

        SnarlAnalyzer(snarl_collection, covariate, maf_threshold, min_individuals),
        phenotype(phenotype), lr() {
        phenotype_type = stoat::QUANTITATIVE;
    };

    EQTLSnarlAnalyzer::EQTLSnarlAnalyzer(
        const SnarlDataCollection &snarl_collection,
        const stoat::CovariateTable& covariate,
        const double maf_threshold,
        const stoat::GeneExpressionTable& gene_expression,
        const size_t max_gene_dist,
        const size_t min_individuals) :
        SnarlAnalyzer(snarl_collection, covariate, maf_threshold, min_individuals),
        gene_expression(gene_expression), max_gene_dist(max_gene_dist), lr() {
        phenotype_type = stoat::EQTL;
    };

    stoat::phenotype_type_t SnarlAnalyzer::get_phenotype_type() const {
        return (phenotype_type);
    }

void SnarlAnalyzer::test_snarls_from_file(stoat::Reader& gt_reader, stoat::Writer& out_writer) {

    // prepare snarl collection that will stream the snarls and open connection to the file
    stoat::SnarlDataCollection snarl_collection_stream(0, 0, 0);

    // Write the header of the output file
    out_writer.write_stoat_output_header(phenotype_type);
 
    // read each snarl and test it
    // JEAN parallelize here?
    size_t number_snarl_filtered = 0;
    snarl_collection_stream.for_each_snarl_in_file(gt_reader, [&](snarl_info_t& snarl_info) {
        bool filtered = test_and_write_snarl(snarl_info, out_writer);
        number_snarl_filtered += (filtered ? 1 : 0);
    });

    stoat::LOG_INFO("Total number of snarl filtered : " + std::to_string(number_snarl_filtered));
}

// Decompose path std::string to vectorstoat::edge_t
std::vector<stoat::edge_t> decompose_path_str_to_edge(const std::string s) {
    std::vector<stoat::edge_t> edges;
    stoat::PathTraversal nodes;

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
            nodes.add_node_traversal_t({node_id, is_rev});
        } else {
            // JEAN should we throw an error here? What are invalid characters?
            ++i; // Skip invalid characters
        }
    }

    // we try flipping the path here to avoid most inconsistencies with vg call's ATs
    // inconsistencies are still possibly because this is potentially a very long path
    // traversing the top-level snarl only while the ones used exploring the snarl tree
    // and preparing the snarl paths are the "simplified"/net versions
    nodes.check_path_flip();

    for (size_t j = 0; j + 1 < nodes.get_path().size(); ++j) {
        stoat::edge_t edge(nodes.get_path()[j], nodes.get_path()[j + 1]);
        edges.emplace_back(edge);
    }
    return edges;
}

bool BinarySnarlAnalyzer::test_and_write_snarl(stoat::snarl_info_t &snarl_data, stoat::Writer& out_writer) {

    // link to the phenotype
    snarl_data.genotypes.link_to_binary_phenotype(phenotype);
    
    // remove non-variable allele, e.g. absent in both groups
    snarl_data.genotypes.remove_noncovered_samples();
    snarl_data.genotypes.remove_constant_predictors();

    // Which alleles are actually used? Used for making sure that we only write the alleles we used
    std::vector<bool> active_alleles = snarl_data.genotypes.get_active_alleles();

    // prepare an output objet and init to NA
    test_result_t test_res;
    test_res.pv = std::nan("");
    test_res.second_pv = std::nan("");

    // should we test this snarl?
    if (snarl_data.genotypes.passes_filters(maf_threshold, min_individuals)){
        // fill up the contingency table (one vector per group)
        std::vector<size_t> g0;
        std::vector<size_t> g1;
        snarl_data.genotypes.fill_contingency_table(g0, g1);

        // performs the test
        auto fc_res = fchi.fisher_chi2(g0, g1);
        test_res.pv = fc_res.first;
        test_res.second_pv = fc_res.second;
        test_res.group_paths = stoat::format_group_paths(g0, g1);
    } else {
        // JEAN what should we do here? returning NA for now
        stoat::LOG_DEBUG("Filtered: didn't pass the filters");
    }
    
    // test this snarl using Fisher exact test and chi-squared test
    if (std::isnan(test_res.pv) && std::isnan(test_res.second_pv)) {
        stoat::LOG_DEBUG("filtered: " + snarl_data.start_node.to_string() + snarl_data.end_node.to_string());
        return true;
    }
    
    out_writer.write_binary(snarl_data, test_res, active_alleles);
    return false;
}

bool ExactBinarySnarlAnalyzer::test_and_write_snarl(stoat::snarl_info_t &snarl_data, stoat::Writer& out_writer) {
    // This test checks if all members of one of the phenotype groups has the same allele that no other sample has.

    // prepare an output objet and init to NA
    test_result_t test_res;
    test_res.pv = std::nan("");
    test_res.second_pv = std::nan("");
    test_res.group_paths = "NA";

    // From the genotype matrix, make sets of sample that have the same genotype and compare to the sets of phenotype groups.
    std::unordered_map<std::string, std::set<std::string>> genotype_to_sample_set;
    for (const sample_hap_t& sample_hap : snarl_data.all_sample_haplotypes) {
        std::string genotype_str = snarl_data.genotypes.get_genotype_as_string(sample_hap.sample);

        if (genotype_to_sample_set.count(genotype_str) == 0) {
            genotype_to_sample_set.emplace(genotype_str, std::set<std::string>());
        }

        genotype_to_sample_set.at(genotype_str).emplace(sample_hap.sample);
    }

    // If we only want to know if there is one allele that matches exactly one of the phenotype groups
    bool write_output = false;
    for (const auto& genotype_samples : genotype_to_sample_set) {
        const std::set<std::string>& partition = genotype_samples.second;

        // If one partition exactly matches one group we want, then all other partitions combined (including things not in the snarl) will match
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
    
    out_writer.write_binary(snarl_data, test_res, snarl_data.genotypes.get_active_alleles());
    return false;
}
    
bool BinaryCovarSnarlAnalyzer::test_and_write_snarl(stoat::snarl_info_t &snarl_data, stoat::Writer& out_writer) {

    // link the phenotype
    snarl_data.genotypes.link_to_binary_phenotype(phenotype);
    snarl_data.genotypes.link_to_covariates(covariate);

    // remove non-variable predictors, e.g. alleles absent in all samples
    snarl_data.genotypes.remove_noncovered_samples();
    snarl_data.genotypes.remove_constant_predictors();    

    // prepare an output objet and init to NA
    test_result_t test_res;
    test_res.pv = std::nan("");

    // Which alleles are actually used? Used for making sure that we only write the alleles we used
    std::vector<bool> active_alleles;
    
    // should we test this snarl?
    if (snarl_data.genotypes.passes_filters(maf_threshold, min_individuals)){
        // add the allele path info to include in the output
        test_res.allele_paths = snarl_data.genotypes.allele_paths_as_str();
    
        // add covariate with the number of alleles (if necessary) to correct for the parent snarl effect (or normalize)
        snarl_data.genotypes.add_total_allele_count_covariable();
    
        // before performing the regression, try to reduce potential colinearity
        snarl_data.genotypes.remove_duplicated_predictors();    
        active_alleles = snarl_data.genotypes.get_active_alleles();
        snarl_data.genotypes.remove_one_allele();
    
        // prepare the matrices, fit the logistic model and test effect of alleles
        Eigen::MatrixXd X = snarl_data.genotypes.make_matrixXd_features();
        if (X.size() == 0) {
            stoat::LOG_DEBUG("Filtered: Genotype matrix is empty");
            return true;
        }

        Eigen::VectorXd Y = snarl_data.genotypes.make_vectorxd_phenotype();
        if (Y.size() == 0) {
            stoat::LOG_DEBUG("Filtered: Phenotype matrix is empty");
            return true;
        }
        test_res.pv = lr.logistic_regression(X, Y, snarl_data.genotypes.get_n_active_alleles());
    } else {
        stoat::LOG_DEBUG("Filtered: didn't pass the filters");
    }
    
    if (std::isnan(test_res.pv)) {
        stoat::LOG_DEBUG("filtered: " + snarl_data.start_node.to_string() + snarl_data.end_node.to_string());
        return true;
    }
 
    out_writer.write_binary_covar(snarl_data, test_res, active_alleles);
    return false;
}

// Quantitative Table Generation
bool QuantitativeSnarlAnalyzer::test_and_write_snarl(stoat::snarl_info_t &snarl_data, stoat::Writer& out_writer) {

    // link the phenotype
    snarl_data.genotypes.link_to_quantitative_phenotype(phenotype);
    snarl_data.genotypes.link_to_covariates(covariate);

    // remove non-variable allele, e.g. absent in both groups
    snarl_data.genotypes.remove_noncovered_samples();    
    snarl_data.genotypes.remove_constant_predictors();

    // prepare an output objet and init to NA
    test_result_t test_res;
    test_res.pv = std::nan("");

    // Which alleles are actually used? Used for making sure that we only write the alleles we used
    std::vector<bool> active_alleles;

    // should we test this snarl?
    if (snarl_data.genotypes.passes_filters(maf_threshold, min_individuals)){
        // add the allele path info to include in the output
        test_res.allele_paths = snarl_data.genotypes.allele_paths_as_str();

        // add covariate with the number of alleles (if necessary) to correct for the parent snarl effect (or normalize)
        snarl_data.genotypes.add_total_allele_count_covariable();

        // before performing the regression, try to reduce potential colinearity
        snarl_data.genotypes.remove_duplicated_predictors();
        active_alleles = snarl_data.genotypes.get_active_alleles();
        snarl_data.genotypes.remove_one_allele();

        Eigen::MatrixXd X = snarl_data.genotypes.make_matrixXd_features();
        if (X.size() == 0) {
            stoat::LOG_DEBUG("Filtered: Genotype matrix is empty");
            return true;
        }
        Eigen::VectorXd Y = snarl_data.genotypes.make_vectorxd_phenotype();
        if (Y.size() == 0) {
            stoat::LOG_DEBUG("Filtered: Phenotype matrix is empty");
            return true;
        }

        test_res.pv = lr.linear_regression(X, Y, snarl_data.genotypes.get_n_active_alleles());
    } else {
        // JEAN what should we do here? returning NA for now
        stoat::LOG_DEBUG("Filtered: didn't pass the filters");
    }

    if (std::isnan(test_res.pv)) {
        stoat::LOG_DEBUG("filtered: " + snarl_data.start_node.to_string() + snarl_data.end_node.to_string());
        return true;
    }
 
    out_writer.write_quantitative(snarl_data, test_res, active_alleles);
    return false;
}

bool EQTLSnarlAnalyzer::test_and_write_snarl(stoat::snarl_info_t &snarl_data, stoat::Writer& out_writer) {

    // get genes near snarl
    std::vector<std::string> genes_near = gene_expression.get_genes_around_pos(snarl_data.ref_path, snarl_data.start_position, snarl_data.end_position, max_gene_dist);
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
        snarl_data.genotypes.clear();

        // link the phenotype
        snarl_data.genotypes.link_to_quantitative_phenotype(gene_phenotype);
        snarl_data.genotypes.link_to_covariates(covariate);

        // remove non-variable allele, e.g. absent in both groups
        snarl_data.genotypes.remove_noncovered_samples();    
        snarl_data.genotypes.remove_constant_predictors();

        // prepare an output objet and init to NA
        test_result_t test_res;
        test_res.pv = std::nan("");

    // Which alleles are actually used? Used for making sure that we only write the alleles we used
        std::vector<bool> active_alleles;

        // should we test this snarl?
        if (snarl_data.genotypes.passes_filters(maf_threshold, min_individuals)){

            // add the allele path info to include in the output
            test_res.allele_paths = snarl_data.genotypes.allele_paths_as_str();

            // add covariate with the number of alleles (if necessary) to correct for the parent snarl effect (or normalize)
            snarl_data.genotypes.add_total_allele_count_covariable();

            // before performing the regression, try to reduce potential colinearity
            snarl_data.genotypes.remove_duplicated_predictors();
            active_alleles = snarl_data.genotypes.get_active_alleles();
            snarl_data.genotypes.remove_one_allele();
    
            Eigen::MatrixXd X = snarl_data.genotypes.make_matrixXd_features();
            if (X.size() == 0) {
                stoat::LOG_DEBUG("Filtered: Genotype matrix is empty");
                continue;
            }
            Eigen::VectorXd Y = snarl_data.genotypes.make_vectorxd_phenotype();
            if (Y.size() == 0) {
                stoat::LOG_DEBUG("Filtered: Phenotype matrix is empty");
                continue;
            }

            test_res.pv = lr.linear_regression(X, Y, snarl_data.genotypes.get_n_active_alleles());
        } else {
            // JEAN what should we do here? returning NA for now
            stoat::LOG_DEBUG("Filtered: didn't pass the filters");
        }

        if (std::isnan(test_res.pv)) {
            stoat::LOG_DEBUG("filtered: gene " + gene_name + ", " + snarl_data.start_node.to_string() + snarl_data.end_node.to_string());
            continue;
        }
  
        out_writer.write_eqtl(snarl_data, gene_name, test_res, active_alleles);

        // at least this test was not filtered
        filtered = false;
    }

    // return if the snarl was filtered in all considered genes
    return filtered;
}

} // end namespace stoat
