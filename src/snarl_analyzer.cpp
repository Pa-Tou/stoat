#include "snarl_analyzer.hpp"
#include "matrix.hpp"
#include "quantitative_table.hpp"
#include "genotype_table.hpp"
#include "utils.hpp"
#include "arg_parser.hpp"
#include "writer.hpp"
#include "omp.h"

namespace stoat_vcf {

    SnarlAnalyzer::SnarlAnalyzer(
        const std::unordered_map<std::string, std::vector<stoat::Snarl_data_t>> &chr_to_snarl_data,
        const std::vector<std::string> &list_samples,
        const std::vector<std::vector<double>> &covariate,
        const double maf_threshold,
        const size_t min_individuals) :
        chr_to_snarl_data(chr_to_snarl_data),
        list_samples(list_samples),
        edge_matrix(list_samples, 0),
        covariate(covariate),
        maf_threshold(maf_threshold),
        min_individuals(min_individuals),
        phenotype_type(stoat::BINARY) {};

    BinarySnarlAnalyzer::BinarySnarlAnalyzer(
        const std::unordered_map<std::string, std::vector<stoat::Snarl_data_t>> &chr_to_snarl_data,
        const std::vector<std::string> &list_samples,
        const double maf_threshold,
        const std::vector<bool> &binary_phenotype,
        const size_t min_individuals) :
        SnarlAnalyzer(chr_to_snarl_data, list_samples, {}, maf_threshold, min_individuals),
        binary_phenotype(binary_phenotype), fchi() {
        phenotype_type = stoat::BINARY;
    };

    BinaryCovarSnarlAnalyzer::BinaryCovarSnarlAnalyzer(
        const std::unordered_map<std::string, std::vector<stoat::Snarl_data_t>> &chr_to_snarl_data,
        const std::vector<std::string> &list_samples,
        const std::vector<std::vector<double>> &covariate,
        const double maf_threshold,
        const std::vector<bool> &binary_phenotype,
        const size_t min_individuals) :

        SnarlAnalyzer(chr_to_snarl_data, list_samples, covariate, maf_threshold, min_individuals),
        binary_phenotype(binary_phenotype), lr() {
        phenotype_type = stoat::BINARY_COVAR;
    };

    QuantitativeSnarlAnalyzer::QuantitativeSnarlAnalyzer(
        const std::unordered_map<std::string, std::vector<stoat::Snarl_data_t>> &chr_to_snarl_data,
        const std::vector<std::string> &list_samples,
        const std::vector<std::vector<double>> &covariate,
        const double maf_threshold,
        const std::vector<double> &quantitative_phenotype,
        const size_t min_individuals) :

        SnarlAnalyzer(chr_to_snarl_data, list_samples, covariate, maf_threshold, min_individuals),
        quantitative_phenotype(quantitative_phenotype), lr() {
        phenotype_type = stoat::QUANTITATIVE;
    };

    EQTLSnarlAnalyzer::EQTLSnarlAnalyzer(
        const std::unordered_map<std::string, std::vector<stoat::Snarl_data_t>> &chr_to_snarl_data,
        const std::vector<std::string> &list_samples,
        const std::vector<std::vector<double>> &covariate,
        const double maf_threshold,
        const std::unordered_map<std::string, std::vector<stoat_vcf::Qtl_data>> &eqtl_map,
        const size_t windows_gene_threshold,
        const size_t min_individuals) :
        SnarlAnalyzer(chr_to_snarl_data, list_samples, covariate, maf_threshold, min_individuals),
        eqtl_map(eqtl_map), windows_gene_threshold(windows_gene_threshold), lr() {
        phenotype_type = stoat::EQTL;
    };

    void SnarlAnalyzer::write_header(std::ofstream &outf)
    {
        stoat::write_stoat_output_header(outf, phenotype_type);
    }

    stoat::phenotype_type_t SnarlAnalyzer::get_phenotype_type() const {
        return (phenotype_type);
    }

    void SnarlAnalyzer::genotype_test_snarls_by_chr_from_vcf(
        htsFile *&ptr_vcf,
        bcf_hdr_t *&hdr,
        bcf1_t *&rec,
        const std::string output_dir)
    {

        // JEAN if we want to keep track of what was run, we might as well include a header in the file with the full info (all parameters, input files, etc)
        std::string output_filename = output_dir + "/stoat.assoc.pvalues.tsv";
        std::ofstream outf(output_filename, std::ios::binary);

        // Write the header of the output file
        write_header(outf);

        // for the log
        size_t total_number_snarl_filtered = 0;

        // Go through the vcf and get chunks by chromosome.
        while (bcf_read(ptr_vcf, hdr, rec) >= 0) {
            
            std::string chr = bcf_hdr_id2name(hdr, rec->rid);

            // Skip chromosomes not in chr_to_snarl_data
            while (chr_to_snarl_data.find(chr) == chr_to_snarl_data.end()) {
                stoat::LOG_WARN("Chromosome " + chr + " not found in snarl paths file. Skipping.");

                bool found_new_chr = false;
                while (bcf_read(ptr_vcf, hdr, rec) >= 0)
                {
                    std::string chr_next = bcf_hdr_id2name(hdr, rec->rid);
                    if (chr_next != chr)
                    {
                        chr = chr_next; // Update to the new chromosome
                        found_new_chr = true;
                        break;
                    }
                }

                if (!found_new_chr)
                {
                    return; // exit if no more records are available
                }
            }

            // start analyzing this chromosome chr
            stoat::LOG_INFO("Analysing chr : " + chr);
            auto timer_start_chr = std::chrono::high_resolution_clock::now();

            // prepare the edge matrix for this chromosome by reading the VCF
            auto [ptr_vcf_new, hdr_new, rec_new] = make_edge_matrix(ptr_vcf, hdr, rec, chr);
            ptr_vcf = ptr_vcf_new;
            hdr = hdr_new;
            rec = rec_new;

            auto timer_end_matrix = std::chrono::high_resolution_clock::now();
            stoat::LOG_INFO("Edge matrix construction for chr " + chr + " : " + std::to_string(std::chrono::duration<double>(timer_end_matrix - timer_start_chr).count()) + " s");

            const auto &snarls = chr_to_snarl_data.at(chr);
            size_t chr_number_snarl_filtered = 0;

            #pragma omp parallel for schedule(static)

            // Make the snarl test analysis
            // Iterate over each snarl
            for (const stoat::Snarl_data_t &snarl_data_s : snarls) {
                bool filtered = test_and_write_snarl(snarl_data_s, chr, outf);
                chr_number_snarl_filtered += (filtered ? 1 : 0);
            }

            total_number_snarl_filtered += chr_number_snarl_filtered;
            auto timer_end_chr = std::chrono::high_resolution_clock::now();

            stoat::LOG_INFO("Number of snarl filtered in chr " + chr + " : " + std::to_string(chr_number_snarl_filtered));
            stoat::LOG_INFO("Snarl analysis in chr " + chr + " : " + std::to_string(std::chrono::duration<double>(timer_end_chr - timer_end_matrix).count()) + " s");
            stoat::LOG_INFO("Total time for chr " + chr + " : " + std::to_string(std::chrono::duration<double>(timer_end_chr - timer_start_chr).count()) + " s");
        }

        stoat::LOG_INFO("Total number of snarl filtered : " + std::to_string(total_number_snarl_filtered));

        // Cleanup
        bcf_destroy(rec);
        bcf_hdr_destroy(hdr);
        bcf_close(ptr_vcf);
    }

    std::tuple<htsFile *, bcf_hdr_t *, bcf1_t *> SnarlAnalyzer::make_edge_matrix(htsFile *ptr_vcf, bcf_hdr_t *hdr, bcf1_t *rec, std::string &chr)
    {

        // init the edge matrix, allocating about 4 edges per snarl to start
        size_t n_snarls = chr_to_snarl_data.at(chr).size();
        edge_matrix.clear_edges(n_snarls * 4);

        // loop over the VCF file for each line and stop where chr is different
        do
        {
            bcf_unpack(rec, BCF_UN_STR);

            // check the INFO field for LV (Level in the snarl tree)
            // skip if variant is lv != 0 to avoid duplication paths/snarl variant analysis
            int32_t *lv = nullptr;
            int n_lv = 0;
            if (bcf_get_info_int32(hdr, rec, "LV", &lv, &n_lv) > 0)
            {
                if (lv[0] != 0)
                {
                    free(lv);
                    continue;
                }
            }
            free(lv);

            // extract genotypes GT
            int ngt = 0;
            int32_t *gt = nullptr;
            ngt = bcf_get_genotypes(hdr, rec, &gt, &ngt);

            if (ngt <= 0 || gt == nullptr)
            {
                throw std::invalid_argument("GT field is missing in VCF at position " + std::to_string(rec->pos + 1));
            }

            // extract AT field from INFO
            char *at = nullptr;
            int nat = 0;
            nat = bcf_get_info_string(hdr, rec, "AT", &at, &nat);
            if (nat <= 0 || !at)
            {
                // AT field is mandatory, throw an error
                throw std::invalid_argument("AT field is missing in VCF at position " + std::to_string(rec->pos + 1));
            }
            std::string at_str(at); // convert to C++ std::string
            free(at);

            // split by comma and save as a vector of edge lists [vector vector stoat::Edge_t]
            // from: ">123>213<234", ">123<234", ">123<234<345"
            // to: [[Edge_t(123, 213),stoat::Edge_t(213, 234)], [...]]
            std::vector<std::vector<stoat::Edge_t>> list_paths_edge;
            std::stringstream at_ss(at_str);
            std::string item;
            while (std::getline(at_ss, item, ','))
            {
                std::vector<stoat::Edge_t> edge_vec = decompose_path_str_to_edge(item);
                list_paths_edge.push_back(edge_vec);
            }

            for (int i = 0; i < rec->n_sample; ++i)
            {
                int ploidy = 2;
                size_t col_idx = i * 2;
                for (int al_idx = 0; al_idx < ploidy; ++al_idx){
                    // allele al_idx of that sample
                    // JEAN here we are assuming diploid genotypes. check how to make sure we're really/always getting the genotype for sample i with bcf_gt_allele
                    int idx_path_allele = bcf_gt_allele(gt[col_idx + al_idx]);
                    if (idx_path_allele != -1)
                        { // Handle missing genotypes
                            for (const auto &edge_path : list_paths_edge[idx_path_allele])
                                {
                                    edge_matrix.add_sample_edge(edge_path, col_idx + al_idx);
                                }
                        }
                }
            }
            free(gt);

        } while ((bcf_read(ptr_vcf, hdr, rec) >= 0) && (chr == bcf_hdr_id2name(hdr, rec->rid)));

        edge_matrix.shrink();
        return std::make_tuple(ptr_vcf, hdr, rec);
    }

    // Decompose path std::string to vectorstoat::Edge_t
    std::vector<stoat::Edge_t> decompose_path_str_to_edge(const std::string s)
    {
        std::vector<stoat::Edge_t> edges;
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
            stoat::Edge_t edge(nodes.get_path()[j], nodes.get_path()[j + 1]);
            edges.emplace_back(edge);
        }

        return edges;
    }

    bool BinarySnarlAnalyzer::test_and_write_snarl(const stoat::Snarl_data_t &snarl_data_s, 
        const std::string chr, std::ofstream &outf) {

            std::string al_lens = vectorPathToString(snarl_data_s.paths, true);
            size_t paths_number = snarl_data_s.paths.size();
            std::vector<size_t> g0(paths_number, 0);
            std::vector<size_t> g1(paths_number, 0);

            // JEAN for later when we use SnarlDataCollecions and Tables
            // extract genotypes for this snarl
            // GenotypeTable snarl_gt = SNARLCOLLECTION.FUNCTION(snarl_data_s.???)
            // auto [chi2_p_value, fastfisher_p_value] = fchi.fisher_chi2(phenotype_tab, snarl_gt, maf_threshold, min_individuals);
            // if (chi2_p_value == "NA" && fastfisher_p_value == "NA") {
            //     stoat::LOG_DEBUG("filtered: " + stoat::pairToString(snarl_data_s.ids));
            //     return true;
            // }
            // JEAN still need to put group_paths somewhere (in fisher_chi2?)
            
            size_t individuals_included = stoat_vcf::create_binary_table(g0, g1, binary_phenotype, snarl_data_s.paths, edge_matrix);
            stoat::remove_empty_columns_binary_table(g0, g1);
            bool to_filter = stoat::filter_binary_table(g0, g1, individuals_included, min_individuals, maf_threshold);

            if (to_filter)
            {
                stoat::LOG_DEBUG("filtered: " + stoat::pairToString(snarl_data_s.ids));
                return to_filter;
            }

            // Binary analysis single test
            auto group_paths = format_group_paths(g0, g1);
            auto [chi2_p_value, fastfisher_p_value] = fchi.fisher_chi2(g0, g1);

            #pragma omp critical(outf)
            {
                stoat::write_binary(outf, chr, snarl_data_s, al_lens, fastfisher_p_value, chi2_p_value, group_paths);
            }

            return to_filter;
    }

    bool BinaryCovarSnarlAnalyzer::test_and_write_snarl(
        const stoat::Snarl_data_t &snarl_data_s, const std::string chr, std::ofstream &outf)
    {
        std::string al_lens = vectorPathToString(snarl_data_s.paths, true);
        

        auto [X, Y, samples_name, allele_paths] = create_quantitative_table(list_samples.size(), snarl_data_s.paths, binary_phenotype, edge_matrix);
        stoat::remove_empty_columns_quantitative_table(X);

        bool to_filter = stoat::filter_quantitative_table(X, min_individuals, maf_threshold);

        if (to_filter)
        {
            stoat::LOG_DEBUG("filtered: " + stoat::pairToString(snarl_data_s.ids) + " -> filter_quantitative_table");
            return to_filter;
        }

        stoat::combine_identical_columns_quantitative_table(X);
        stoat::remove_last_columns_quantitative_table(X);

        to_filter = stoat::check_last_columns_quantitative_table(X);

        if (to_filter)
        {
            stoat::LOG_DEBUG("filtered: " + stoat::pairToString(snarl_data_s.ids) + " -> check_last_columns_quantitative_table");
            return to_filter;
        }
        
        // logistic regression with covariates if not empty
        const auto &[p_value, beta, se] = lr.logistic_regression(X, Y, covariate);

        #pragma omp critical(outf)
        {
            stoat::write_binary_covar(outf, chr, snarl_data_s, al_lens, p_value, allele_paths);
        }
        return to_filter;
    }

    // Quantitative Table Generation
    bool QuantitativeSnarlAnalyzer::test_and_write_snarl(
        const stoat::Snarl_data_t &snarl_data_s,
        const std::string chr,
        std::ofstream &outf)
    {

        bool to_filter = false;
        auto [X, Y, samples_name, allele_paths] = create_quantitative_table(list_samples.size(), snarl_data_s.paths, quantitative_phenotype, edge_matrix);
        stoat::remove_empty_columns_quantitative_table(X);

        to_filter = stoat::filter_quantitative_table(X, min_individuals, maf_threshold);

        if (to_filter)
        {
            stoat::LOG_DEBUG("filtered: " + stoat::pairToString(snarl_data_s.ids) + " -> filter_quantitative_table");
            return to_filter;
        }

        stoat::combine_identical_columns_quantitative_table(X);
        stoat::remove_last_columns_quantitative_table(X);
        to_filter = stoat::check_last_columns_quantitative_table(X);

        if (to_filter)
        {
            stoat::LOG_DEBUG("filtered: " + stoat::pairToString(snarl_data_s.ids) + " -> check_last_columns_quantitative_table");
            return to_filter;
        }

        std::string al_lens = vectorPathToString(snarl_data_s.paths, true);

        auto [p_value, r2] = lr.linear_regression(X, Y, covariate);

        #pragma omp critical(outf)
        {
            stoat::write_quantitative(outf, chr, snarl_data_s, al_lens, p_value, r2, allele_paths);
        }
        return to_filter;
    }

    // Identify genes index that will be tested for this snarl by matching position
    // Qtl_data : <gene_name, gene_expression, start_pos, end_pos>
    std::vector<size_t> get_genes_around_pos(
        const std::vector<Qtl_data> &gene_position,
        const size_t start_pos,
        const size_t end_pos,
        const size_t max_distance)
    {
        // we look for genes in the specified range +- the maximum distance
        size_t start_pos_threshold = (start_pos > max_distance) ? start_pos - max_distance : 0;
        size_t end_pos_threshold = end_pos + max_distance;

        // loop over all genes and save the ones that are close enough
        std::vector<size_t> gene_index;
        for (size_t i = 0; i < gene_position.size(); ++i)
        {
            size_t gene_start = gene_position[i].start_pos;
            size_t gene_end = gene_position[i].end_pos;

            // Check if the gene overlaps with the snarl region
            if (!(gene_end < start_pos_threshold || gene_start > end_pos_threshold))
            {
                gene_index.push_back(i);
            }
        }
        return gene_index;
    }

    bool EQTLSnarlAnalyzer::test_and_write_snarl(
        const stoat::Snarl_data_t &snarl_data_s, const std::string chr, std::ofstream &outf)
    {

        bool to_filter = false;
        std::vector<size_t> list_gene_index = get_genes_around_pos(eqtl_map.at(chr), snarl_data_s.start_position, snarl_data_s.end_position, windows_gene_threshold);
        auto [X, index_filtered, samples_name, allele_paths] = stoat_vcf::create_eqtl_table(list_samples.size(), snarl_data_s.paths, edge_matrix);

        stoat::remove_empty_columns_quantitative_table(X);
        to_filter = stoat::filter_quantitative_table(X, min_individuals, maf_threshold);

        if (to_filter)
        {
            stoat::LOG_DEBUG("filtered: " + stoat::pairToString(snarl_data_s.ids) + " -> filter_quantitative_table");
            return to_filter;
        }

        stoat::combine_identical_columns_quantitative_table(X);
        stoat::remove_last_columns_quantitative_table(X);
        to_filter = stoat::check_last_columns_quantitative_table(X);

        if (to_filter)
        {
            stoat::LOG_DEBUG("filtered: " + stoat::pairToString(snarl_data_s.ids) + " -> check_last_columns_quantitative_table");
            return to_filter;
        }

        for (size_t i = 0; i < list_gene_index.size(); ++i)
        {
            size_t gene_idx = list_gene_index[i];
            std::string gene_name = eqtl_map.at(chr)[gene_idx].geneName;
            std::vector<double> gene_expression = eqtl_map.at(chr)[gene_idx].sampleExpresion;
            // keep only the desired samples (in index_filtered)
            // used to be "retain_indices"
            // JEAN will be replaced by the new tables I think
            size_t write_idx = 0;
            for (size_t read_idx = 0; read_idx < gene_expression.size(); ++read_idx) {
                if (index_filtered.count(read_idx)) {
                    gene_expression[write_idx++] = gene_expression[read_idx];
                }
            }
            gene_expression.resize(write_idx);

            
            std::string al_lens = vectorPathToString(snarl_data_s.paths, true);

            auto [p_value, r2] = lr.linear_regression(X, gene_expression, covariate);

            #pragma omp critical(outf)
            {
                stoat::write_eqtl(outf, chr, snarl_data_s, al_lens, gene_name, p_value, r2, allele_paths);
            }
        }
        return to_filter;
    }

} // end namespace stoat
