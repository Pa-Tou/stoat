#include "snarl_analyzer.hpp"
#include "matrix.hpp"
#include "quantitative_table.hpp"
#include "utils.hpp"
#include "arg_parser.hpp"
#include "writer.hpp"
#include "omp.h"

using namespace stoat;
namespace stoat_vcf {

    SnarlAnalyzer::SnarlAnalyzer(
        const SnarlDataCollection &snarl_collection,
        const std::unordered_set<std::string>& ref_chrs,
        const stoat::CovariateTable& covariate,
        const double maf_threshold,
        const size_t min_individuals) :
        snarl_collection(snarl_collection),
        ref_chrs(ref_chrs),
        covariate(covariate),
        maf_threshold(maf_threshold),
        min_individuals(min_individuals),
        phenotype_type(stoat::BINARY) {};

    SnarlAnalyzer::SnarlAnalyzer(
        const SnarlDataCollection &snarl_collection,
        const std::unordered_set<std::string>& ref_chrs,
        const double maf_threshold,
        const size_t min_individuals) :
        snarl_collection(snarl_collection),
        ref_chrs(ref_chrs),
        covariate(covariate),
        maf_threshold(maf_threshold),
        min_individuals(min_individuals),
        phenotype_type(stoat::BINARY) {};

    BinarySnarlAnalyzer::BinarySnarlAnalyzer(
        const SnarlDataCollection &snarl_collection,
        const std::unordered_set<std::string>& ref_chrs,
        const double maf_threshold,
        const stoat::BinaryPhenotypeTable& phenotype,
        const size_t min_individuals) :
        SnarlAnalyzer(snarl_collection, ref_chrs, maf_threshold, min_individuals),
        phenotype(phenotype), fchi() {
        phenotype_type = stoat::BINARY;
    };

    BinaryCovarSnarlAnalyzer::BinaryCovarSnarlAnalyzer(
        const SnarlDataCollection &snarl_collection,
        const std::unordered_set<std::string>& ref_chrs,
        const stoat::CovariateTable& covariate,
        const double maf_threshold,
        const stoat::BinaryPhenotypeTable& phenotype,
        const size_t min_individuals) :

        SnarlAnalyzer(snarl_collection, ref_chrs, covariate, maf_threshold, min_individuals),
        phenotype(phenotype), lr() {
        phenotype_type = stoat::BINARY_COVAR;
    };

    QuantitativeSnarlAnalyzer::QuantitativeSnarlAnalyzer(
        const SnarlDataCollection &snarl_collection,
        const std::unordered_set<std::string>& ref_chrs,
        const stoat::CovariateTable& covariate,
        const double maf_threshold,
        const stoat::QuantitativePhenotypeTable& phenotype,
        const size_t min_individuals) :

        SnarlAnalyzer(snarl_collection, ref_chrs, covariate, maf_threshold, min_individuals),
        phenotype(phenotype), lr() {
        phenotype_type = stoat::QUANTITATIVE;
    };

    EQTLSnarlAnalyzer::EQTLSnarlAnalyzer(
        const SnarlDataCollection &snarl_collection,
        const std::unordered_set<std::string>& ref_chrs,
        const stoat::CovariateTable& covariate,
        const double maf_threshold,
        const stoat::GeneExpressionTable& gene_expression,
        const size_t max_gene_dist,
        const size_t min_individuals) :
        SnarlAnalyzer(snarl_collection, ref_chrs, covariate, maf_threshold, min_individuals),
        gene_expression(gene_expression), max_gene_dist(max_gene_dist), lr() {
        phenotype_type = stoat::EQTL;
    };

    void SnarlAnalyzer::write_header(std::ofstream &outf)
    {
        stoat::write_stoat_output_header(outf, phenotype_type);
    }

    stoat::phenotype_type_t SnarlAnalyzer::get_phenotype_type() const {
        return (phenotype_type);
    }

void SnarlAnalyzer::genotype_test_snarls_by_chr(const std::string output_dir) {
    
    // JEAN if we want to keep track of what was run, we might as well include a header in the file with the full info (all parameters, input files, etc)
    std::string output_filename = output_dir + "/stoat.assoc.pvalues.tsv";
    std::ofstream outf(output_filename, std::ios::binary);
    
    // Write the header of the output file
    write_header(outf);

    // for the log
    size_t total_number_snarl_filtered = 0;

    for (std::string chrom : ref_chrs) {
        // start analyzing this chromosome
        stoat::LOG_INFO("Analysing chr : " + chrom);
        auto timer_start_chr = std::chrono::high_resolution_clock::now();
        size_t chr_number_snarl_filtered = 0;

        // JEAN parallelize here?    
        snarl_collection.for_each_snarl([&](snarl_info_t& snarl_info) {
            if (snarl_info.ref_path == chrom) {
                bool filtered = test_and_write_snarl(snarl_info, outf);
                chr_number_snarl_filtered += (filtered ? 1 : 0);
            }
        });

        total_number_snarl_filtered += chr_number_snarl_filtered;
        auto timer_end_chr = std::chrono::high_resolution_clock::now();

        stoat::LOG_INFO("Number of snarl filtered in chr " + chrom + " : " + std::to_string(chr_number_snarl_filtered));
        stoat::LOG_INFO("Total time for chr " + chrom + " : " + std::to_string(std::chrono::duration<double>(timer_end_chr - timer_start_chr).count()) + " s");
    }
    
    stoat::LOG_INFO("Total number of snarl filtered : " + std::to_string(total_number_snarl_filtered));
}

void SnarlAnalyzer::test_snarls_from_file(const std::string snarl_genotype_path, const std::string output_dir) {

    // prepare snarl collection that will stream the snarls and open connection to the file
    stoat::SnarlDataCollection snarl_collection_stream(0, 0, 0);
    std::ifstream snarl_stream;
    snarl_stream.open(snarl_genotype_path);

    std::string output_filename = output_dir + "/stoat.assoc.pvalues.tsv";
    std::ofstream outf(output_filename, std::ios::binary);
    
    // Write the header of the output file
    write_header(outf);

    // count snarls filterd for the log
    size_t total_number_snarl_filtered = 0;

    // read each snarl and test it
    // JEAN parallelize here?
    size_t number_snarl_filtered = 0;
    snarl_collection_stream.for_each_snarl_in_file(snarl_stream, [&](snarl_info_t& snarl_info) {
        bool filtered = test_and_write_snarl(snarl_info, outf);
        number_snarl_filtered += (filtered ? 1 : 0);
    });

    snarl_stream.close();
    outf.close();
    
    stoat::LOG_INFO("Total number of snarl filtered : " + std::to_string(number_snarl_filtered));
}

// Decompose path std::string to vectorstoat::edge_t
std::vector<stoat::edge_t> decompose_path_str_to_edge(const std::string s) {
    std::vector<stoat::edge_t> edges;
    stoat::PathTraversal nodes;

    size_t i = 0;
    while (i < s.size())
        {
            if (s[i] == '>' || s[i] == '<')
                {
                    bool is_rev = (s[i] == '<');
                    ++i;

                    size_t node_id = 0;
                    while (i < s.size() && isdigit(s[i]))
                        {
                            node_id = node_id * 10 + (s[i] - '0');
                            ++i;
                        }
                    nodes.add_node_traversal_t({node_id, is_rev});
                }
            else
                {
                    // JEAN should we throw an error here? What are invalid characters?
                    ++i; // Skip invalid characters
                }
        }

    // we try flipping the path here to avoid most inconsistencies with vg call's ATs
    // inconsistencies are still possiblt because this is potentially a very long path
    // traversing the top-level snarl only while the ones used exploring the snarl tree
    // and preparing the snarl paths are the "simplified"/net versions
    nodes.check_path_flip();

    for (size_t j = 0; j + 1 < nodes.get_path().size(); ++j)
        {
            stoat::edge_t edge(nodes.get_path()[j], nodes.get_path()[j + 1]);
            edges.emplace_back(edge);
        }
    return edges;
}

bool BinarySnarlAnalyzer::test_and_write_snarl(stoat::snarl_info_t &snarl_data, std::ofstream &outf) {
    // test this snarl using Fisher exact test and chi-squared test
    test_result_t test_res = fchi.fisher_chi2(phenotype, snarl_data.genotypes, maf_threshold, min_individuals);
    if (test_res.pv == "NA" && test_res.second_pv == "NA") {
        stoat::LOG_DEBUG("filtered: " + snarl_data.start_node.to_string() + snarl_data.end_node.to_string());
        return true;
    }
    
#pragma omp critical(outf)
    {
        stoat::write_binary(outf, snarl_data, test_res);
    }
    return false;
}
    
bool BinaryCovarSnarlAnalyzer::test_and_write_snarl(stoat::snarl_info_t &snarl_data, std::ofstream &outf) {
    test_result_t test_res = lr.logistic_regression(phenotype, snarl_data.genotypes, covariate, maf_threshold, min_individuals);
    if (test_res.pv == "NA") {
        stoat::LOG_DEBUG("filtered: " + snarl_data.start_node.to_string() + snarl_data.end_node.to_string());
        return true;
    }
 
#pragma omp critical(outf)
    {
        stoat::write_binary_covar(outf, snarl_data, test_res);
    }
    return false;
}

// Quantitative Table Generation
bool QuantitativeSnarlAnalyzer::test_and_write_snarl(stoat::snarl_info_t &snarl_data, std::ofstream &outf) {
    test_result_t test_res = lr.linear_regression(phenotype, snarl_data.genotypes, covariate, maf_threshold, min_individuals);
    if (test_res.pv == "NA") {
        stoat::LOG_DEBUG("filtered: " + snarl_data.start_node.to_string() + snarl_data.end_node.to_string());
        return true;
    }
 
#pragma omp critical(outf)
    {
        stoat::write_quantitative(outf, snarl_data, test_res);
    }
    return false;
}

bool EQTLSnarlAnalyzer::test_and_write_snarl(stoat::snarl_info_t &snarl_data, std::ofstream &outf) {
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
        test_result_t test_res = lr.linear_regression(gene_phenotype, snarl_data.genotypes, covariate, maf_threshold, min_individuals);
        if (test_res.pv == "NA") {
            stoat::LOG_DEBUG("filtered: gene " + gene_name + ", " + snarl_data.start_node.to_string() + snarl_data.end_node.to_string());
            continue;
        }
 
#pragma omp critical(outf)
        {
            stoat::write_eqtl(outf, snarl_data, gene_name, test_res);
        }
        // at least this test was not filtered
        filtered = false;
    }

    // return if the snarl was filtered in all considered genes
    return filtered;
}

} // end namespace stoat
