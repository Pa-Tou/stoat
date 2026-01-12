#include "vcf.hpp"

#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <cstdlib>
#include <getopt.h>
#include <omp.h>

#include "../log.hpp"
#include "../snarl_analyzer.hpp"
#include "../arg_parser.hpp"
#include "../matrix.hpp"
#include "../gaf_creator.hpp"
#include "../post_processing.hpp"
#include "../io/register_io.hpp"
#include "../feature_tables.hpp"

//#define USE_CALLGRIND

#ifdef USE_CALLGRIND
    #include <valgrind/callgrind.h>
#endif


namespace stoat_command {

void print_help_vcf() {
    std::cerr << "Usage: stoat vcf [options]\n\n"
              << "  -g, --graph FILE                Path to the graph file (only Packed Graph works for now)\n"
              << "  -d, --dist FILE                 Path to the distance index file\n"
              << "  -v, --vcf FILE                  Path to the VCF file\n"
              << "  -s, --snarl FILE                Path to the snarl file\n"
              << "  -r, --reference-chrs FILE       Path to the chromosome reference file, one path name per line\n"
              << "  -b, --binary FILE               Path to the binary phenotype group file\n"
              << "  -q, --quantitative FILE         Path to the quantitative phenotype file\n"
              << "  -e, --gene-expression FILE      Path to the gene expression file (for eQTL analysis)\n"
              << "  -c, --covariate FILE            Path to the covariate file\n"
              << "  -C, --covar-name NAME           Covariate column name(s) used for GWAS\n"
              << "  -G, --gaf                       Generate a GAF file with the GWAS results\n"
              << "  -I, --min-individuals INT       Minimum number of individuals per snarl [0]\n"
              << "  -i, --children INT              Max number of children per snarl in decomposition [50]\n"
              << "  -y, --cycle INT                 Max number of authorized cycles in snarl decomposition [1]\n"
              << "  -l, --path-length INT           Max number of nodes in paths during snarl decomposition [50]\n"
              << "  -P, --gene-position FILE        Path to the gene position file\n"
              << "  -w, --max-gene-distance INT     Include snarls up to this distance from the gene when looking for eQTLs [1000000]\n"
              << "  -M, --maf FLOAT                 Minimum allele frequency threshold [0.05]\n"
              << "  -t, --threads INT               Number of threads to use [1]\n"
              << "  -V, --verbose INT               Verbosity level (0=error, 1=warn, 2=info, 3=debug, 4=trace) [2]\n"
              << "  -o, --output FILE               Output directory name\n"
              << "  -h, --help                      Print this help message\n";
}

int main_stoat_vcf(int argc, char* argv[]) {
    
    // Declare variables to hold argument values
    // JEAN maybe avoid having so many different input files. Couldn't we guess the type of phenotype/analysis or just add a "mode" argument?
    std::string vcf_path, snarl_path, graph_path, dist_path, 
        chromosome_path, binary_path, quantitative_path, 
        eqtl_path, covariate_path, gene_position_path;

    size_t phenotype = 0;
    size_t cycle_threshold = 1;
    size_t children_threshold = 50;
    size_t min_individuals = 0;
    // JEAN this threshold is a bit redundant with children_threshold and cycle_threshold but I guess could be useful if we want to set it lower than (children_threshold * (cycle_threshold+1))
    size_t path_length_threshold = 50;
    size_t max_gene_dist = 1000000;

    double maf_threshold = 0.05;
    std::string output_dir = "output";
    
    bool gaf = false;
    bool only_snarl_parsing = false;

    std::vector<std::string> covar_names;

    // Parse arguments
    int c;

    static struct option long_options[] = {
        {"vcf", required_argument, 0, 'v'},
        {"snarl", required_argument, 0, 's'},
        {"graph", required_argument, 0, 'g'},
        {"dist", required_argument, 0, 'd'},
        {"reference-chrs", required_argument, 0, 'r'},
        {"binary", required_argument, 0, 'b'},
        {"quantitative", required_argument, 0, 'q'},
        {"gene-expression", required_argument, 0, 'e'},
        {"covariate", required_argument, 0, 'c'},
        {"covar-name", required_argument, 0, 'C'},
        {"gaf", no_argument, 0, 'G'},
        {"min-individuals", required_argument, 0, 'I'},
        {"children", required_argument, 0, 'i'},
        {"cycle", required_argument, 0, 'y'},
        {"path-length", required_argument, 0, 'l'},
        {"gene-position", required_argument, 0, 'P'},
        {"max-gene-distance", required_argument, 0, 'w'},
        {"maf", required_argument, 0, 'M'},
        {"thread", required_argument, 0, 't'},
        {"verbose", required_argument, 0, 'V'},
        {"output", required_argument, 0, 'o'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    while ((c = getopt_long(argc, argv, "v:s:g:d:r:b:q:e:c:C:i:y:l:P:w:M:t:V:I:o:Gh", long_options, nullptr)) != -1) {
        switch (c) {
            case 'v': vcf_path = optarg; stoat_vcf::check_file(vcf_path); break;
            case 's': snarl_path = optarg; stoat_vcf::check_file(snarl_path); break;
            case 'g': graph_path = optarg; stoat_vcf::check_file(graph_path); break;
            case 'd': dist_path = optarg; stoat_vcf::check_file(dist_path); break;
            case 'r': chromosome_path = optarg; stoat_vcf::check_file(chromosome_path); break;
            case 'b': binary_path = optarg; phenotype++; stoat_vcf::check_file(binary_path); break;
            case 'q': quantitative_path = optarg; phenotype++; stoat_vcf::check_file(quantitative_path); break;
            case 'e': eqtl_path = optarg; phenotype++; stoat_vcf::check_file(eqtl_path); break;
            case 'c': covariate_path = optarg; stoat_vcf::check_file(covariate_path); break;
            case 'C': {
                std::stringstream ss(optarg);
                std::string token;
                while (std::getline(ss, token, ',')) covar_names.push_back(token);
                break;
            }
            case 'G': gaf = true; break;
            case 'I':
                min_individuals = std::stoi(optarg);
                if (min_individuals < 2) {
                    throw std::runtime_error("Error: [stoat vcf] min_individuals threshold must be > 1");
                }
                break;
            case 'i':
                children_threshold = std::stoi(optarg);
                if (children_threshold < 2) {
                    throw std::runtime_error("Error: [stoat vcf] Children threshold must be > 1");
                }
                break;
            case 'y':
                cycle_threshold = std::stoi(optarg);
                if (cycle_threshold < 1) {
                    throw std::runtime_error("Error: [stoat vcf] Cycle threshold must be > 0");
                }
                break;
            case 'l':
                path_length_threshold = std::stoi(optarg);
                if (path_length_threshold < 2) {
                    throw std::runtime_error("Error: [stoat vcf] Path length threshold must be > 1");
                }
                break;
            case 'P': gene_position_path = optarg; stoat_vcf::check_file(gene_position_path); break;
            case 'w':
                max_gene_dist = std::stoi(optarg);
                if (max_gene_dist < 1) {
                    throw std::runtime_error("Error: [stoat vcf] Maximum gene distance (-w, --max-gene-distance) must be > 0");
                }
                break;
            case 'M':
                maf_threshold = std::stod(optarg);
                if (maf_threshold < 0 || maf_threshold > 1) {
                    throw std::runtime_error("Error: [stoat vcf] MAF must be in [0,1]");
                }
                break;
            case 't':
                if (std::stoi(optarg) < 1) {
                    throw std::runtime_error("Error: [stoat vcf] Number of threads must be > 0");
                }
                omp_set_num_threads(std::stoi(optarg));
                break;
            case 'V': 
                {
                int level = std::stoi(optarg);
                if (level < 0 || level > 4) {
                    throw std::runtime_error("Error: [stoat vcf] Invalid verbosity level. Use 0=Error, 1=Warn, 2=Info, 3=Debug, 4=Trace");
                }
                stoat::LogLevel logLevel = static_cast<stoat::LogLevel>(level);
                stoat::Logger::instance().setLevel(logLevel);                
                break;
                }
            case 'o': output_dir = optarg; break;
            case 'h': 
                print_help_vcf(); 
                return EXIT_SUCCESS; 
            default:
                stoat::LOG_ERROR("[stoat vcf] Unknown argument");
                print_help_vcf();
                return EXIT_FAILURE;
        }
    }

    if (argc == 2) {
        print_help_vcf();
        return EXIT_FAILURE;
    }
    
    if (!covariate_path.empty() && covar_names.empty()) {
        stoat::LOG_ERROR("[stoat vcf] If --covariate path is provided you must add the column name(s), using --covar-name");
        print_help_vcf();
        return EXIT_FAILURE;
    }

    if ((!eqtl_path.empty() && gene_position_path.empty()) || (eqtl_path.empty() && !gene_position_path.empty())) {
        stoat::LOG_ERROR("[stoat vcf] eqtl phenotype file and gene position file must be provided together");
        print_help_vcf();
        return EXIT_FAILURE;
    }

    std::filesystem::create_directory(output_dir);
    std::unordered_set<std::string> ref_chrs = (!chromosome_path.empty()) ? stoat_vcf::parse_chromosome_reference(chromosome_path) : std::unordered_set<std::string>{};
    stoat::Logger::instance().setLogFile(output_dir + "/stoat.log");

    // add command launch in log file
    std::stringstream ss;
    ss << "stoat ";
    for (int i = 0; i < argc; ++i) ss << argv[i] << " ";
    stoat::LOG_SILENTE(ss.str());

    // Enforce valid argument combinations
    if ((!snarl_path.empty() || (!graph_path.empty() && !dist_path.empty())) && !vcf_path.empty() && phenotype == 1) {
        //stoat::LOG_TRACE("Case Gwas");
        // Case 2: graph_path + dist_path + vcf_path + phenotype
        // Case 3: snarl_path + vcf_path + phenotype
    } else if (!graph_path.empty() && !dist_path.empty() && vcf_path.empty() && snarl_path.empty() && phenotype == 0) {
        //stoat::LOG_TRACE("Case Snarl path decomposition");
        // Case 1: Only graph_path + dist_path
        only_snarl_parsing = true;
    } else {
        stoat::LOG_ERROR("[stoat vcf] " +
            std::string("Invalid argument combination provided.\n") +
            "There are only 3 ways to launch stoat vcf:\n" +
            "Case 1 (snarl path decomposition): graph_path + dist_path\n" +
            "Case 2 (GWAS + snarl path decomposition): graph_path + dist_path + vcf_path + phenotype (+ optional file)\n" +
            "Case 3 (GWAS only): snarl_path + vcf_path + phenotype (+ optional file)"
        );
        print_help_vcf();
        return EXIT_FAILURE;
    }

    if ((gaf == true && binary_path.empty()) || (gaf == true && graph_path.empty())) {
        stoat::LOG_ERROR("[stoat vcf] GAF file can be generated only with binary phenotype AND providing a graph (-g)");
        print_help_vcf();
        return EXIT_FAILURE;
    }

    auto start_total_timer = std::chrono::high_resolution_clock::now();

    // start reading the VCF to get the sample list
    std::vector<std::string> list_samples;
    htsFile* ptr_vcf;
    bcf_hdr_t* hdr;
    bcf1_t* rec;

    if (!only_snarl_parsing) {
        stoat::LOG_TRACE("Parsing header VCF file");
        std::tie(list_samples, ptr_vcf, hdr, rec) = stoat_vcf::parseHeader(vcf_path); 
    }

    //////////////////// Load the phenotypes and covariate matrix from files

    // JEAN these phenotypes/covariate could be 1-2 objects, well-defined, and keeping the sample names/order somewhere
    std::vector<bool> binary_phenotype;
    std::vector<double> quantitative_phenotype;
    unique_ptr<stoat::BinaryPhenotypeTable> binary_phenotype_table;

    // dict chr:string : vector{geneName:string, sample_expression:vector<double>, start_pos:size_t, end_pos:size_t}
    std::unordered_map<std::string, std::vector<stoat_vcf::Qtl_data>> eqtl_phenotype;
    std::vector<std::vector<double>> covariate;

    if (!covariate_path.empty()) {
        stoat::LOG_TRACE("Parsing covariate file");
        covariate = stoat_vcf::parse_covariates(covariate_path, covar_names, list_samples);
    }

    if (!binary_path.empty()) {
        stoat::LOG_TRACE("Parsing binary phenotype file");
        binary_phenotype = stoat_vcf::parse_binary_pheno(binary_path, list_samples);
    } else if (!quantitative_path.empty()) {
        stoat::LOG_TRACE("Parsing quantitative phenotype file");
        quantitative_phenotype = stoat_vcf::parse_quantitative_pheno(quantitative_path, list_samples);

    } else if (!eqtl_path.empty() && !gene_position_path.empty()) {
        stoat::LOG_TRACE("Parsing eqtl phenotype file");
        eqtl_phenotype = stoat_vcf::parse_qtl_gene_file(eqtl_path, gene_position_path, list_samples);
    }

    // Make an empty SnarlDataCollection, to be filled in or loaded
    // TODO: Double check that these thresholds are doing the right thing
    stoat::SnarlDataCollection snarl_collection(0, children_threshold, path_length_threshold);

    unique_ptr<bdsg::SnarlDistanceIndex> distance_index;
    unique_ptr<handlegraph::PathHandleGraph> graph;

    bdsg::PathPositionOverlayHelper overlay_helper;
    bdsg::PathPositionHandleGraph* path_position_graph;

    // handlegraph::net_handle_t root;

// Start tracking with callgrind
#ifdef USE_CALLGRIND
    CALLGRIND_START_INSTRUMENTATION;
#endif


    if (!snarl_path.empty()){ // If we have already saved the paths in snarls, load them
        stoat::LOG_TRACE("Reading snarl path file");
        ifstream snarls_in;
        snarls_in.open(snarl_path);
        snarl_collection.load_snarl_data_collection(snarls_in);
        snarls_in.close();

    } else { // Otherwise, find them from the graph and snarl tree
        stoat::LOG_INFO("Starting snarl decomposition... ");
        auto start_dec_timer = std::chrono::high_resolution_clock::now();

        // Load the snarl tree and graph
        std::tie(distance_index, graph) = stoat::load_graph_tree(graph_path, dist_path);
        path_position_graph =  overlay_helper.apply(graph.get());

        // Check if chr present in chr file is present in the graph
        for (const auto& chr : ref_chrs) {
            stoat::LOG_TRACE("Sequence name not found in -r/--chr file: " + chr);
            if (!graph->has_path(chr)) {
                throw std::runtime_error("Reference chromosome: " + chr + " not present in graph");
            }
        }

        // The snarl collection requires sample_haplotypes instead of samples so make copy sample_hap_t's with empty haplotypes
        std::vector<stoat::sample_hap_t> sample_haplotypes;
        for (string sample : list_samples) {
            sample_haplotypes.emplace_back(sample, "");
        }

        snarl_collection.fill_in_snarl_info(*path_position_graph, *distance_index, sample_haplotypes,
            true, //find_alleles_first, doesn't matter in this case
            true, // walks_requested
            [&] (const net_handle_t& snarl, const snarl_info_t& snarl_data,std::vector<PathTraversal>& walks) { // function to fill in walks
                SnarlDataCollection::get_all_walks_through_snarl(*path_position_graph, *distance_index, snarl, snarl_data, walks, cycle_threshold); //TODO: Use Matis's version and write the skipped snarls somewhere
            },
            false, //alleles_requested
            [&] (const net_handle_t& snarl, const snarl_info_t& snarl_data, const std::vector<stoat::sample_hap_t>& all_sample_haplotypes) { //function to find alleles
                return std::vector<size_t>();
            }, 
            false, // sequence_requested 
            ref_chrs, // reference 
            false //check distances
            ); 

        // Always write the snarls
        ofstream snarls_out;
        snarls_out.open(output_dir + "/snarl_info.tsv");
        snarl_collection.write_snarl_data_collection(snarls_out);
        snarls_out.close();

        auto end_dec_timer = std::chrono::high_resolution_clock::now();
        stoat::LOG_INFO("Snarl decomposition took " + std::to_string(std::chrono::duration<double>(end_dec_timer - start_dec_timer).count()) + " s");

        if (only_snarl_parsing) {
            return EXIT_SUCCESS;
        }

        // Clean up unique_ptr except graph
        distance_index.reset();
    }

    // If there were no references given, fill them in with the references from the snarl collection
    if (ref_chrs.empty()) {
        for (const std::string& ref : snarl_collection.get_reference_names()) {
            ref_chrs.insert(ref);
        }
    }

    //////////////////////////////////////// Go through the vcf, do the analysis, and write the output
    stoat::LOG_INFO("Starting GWAS analysis...");
    auto start_gwas_timer = std::chrono::high_resolution_clock::now();
    std::shared_ptr<stoat_vcf::SnarlAnalyzer> snarl_analyzer;

    // Decide which type of SnarlAnalyzer we want
    if (!binary_path.empty()) {
        // binary
        if (!covariate.empty()){
            // Binary covariate
            snarl_analyzer.reset(new stoat_vcf::BinaryCovarSnarlAnalyzer(snarl_collection, ref_chrs, list_samples, covariate, maf_threshold, binary_phenotype, min_individuals));
        } else {
            // Binary without covariate
            snarl_analyzer.reset(new stoat_vcf::BinarySnarlAnalyzer(snarl_collection, ref_chrs, list_samples, maf_threshold,
                                                                    binary_phenotype, min_individuals));
        }
    } else if (!quantitative_path.empty()) {
        // Quantitative
        snarl_analyzer.reset(new stoat_vcf::QuantitativeSnarlAnalyzer(snarl_collection, ref_chrs, list_samples, covariate, maf_threshold, 
                                                                      quantitative_phenotype, min_individuals));
    } else if (!eqtl_path.empty()) {
        // EQTL
        snarl_analyzer.reset(new stoat_vcf::EQTLSnarlAnalyzer(snarl_collection, ref_chrs, list_samples, covariate, maf_threshold, 
                                                              eqtl_phenotype, max_gene_dist, min_individuals));
    }

    // read the VCF by chromosomome, genotype each snarl and perform the association test
    snarl_analyzer->genotype_test_snarls_by_chr_from_vcf(ptr_vcf, hdr, rec, output_dir);

    // eventually, make a GAF to visualize the tested snarls and their association signal
    if (snarl_analyzer->get_phenotype_type() == stoat::BINARY && gaf) {
        stoat::LOG_TRACE("Create GAF");
        std::string output_gaf = output_dir + "/stoat.assoc.gaf";
        stoat_vcf::gaf_creation(output_dir + "/stoat.assoc.pvalues.tsv", snarl_collection, *graph, output_gaf);
    }

    auto end_total_timer = std::chrono::high_resolution_clock::now();
    stoat::LOG_INFO("GWAS analysis took " + std::to_string(std::chrono::duration<double>(end_total_timer - start_gwas_timer).count()) + " s");
    stoat::LOG_INFO("Total time: " + std::to_string(std::chrono::duration<double>(end_total_timer - start_total_timer).count()) + " s");
    return EXIT_SUCCESS;
}

} // end stoat
