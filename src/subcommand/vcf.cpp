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
#include "../path_partitioner.hpp"

// #define USE_CALLGRIND

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
              << "  -i, --children INT              Max number of children per snarl in decomposition [50]\n"
              << "  -y, --cycle INT                 Max number of authorized cycles in snarl decomposition [1]\n"
              << "  -l, --path-length INT           Max number of nodes in paths during snarl decomposition [50]\n"
              << "  -R, --resolve-vcf               Resolve conflicting calls in the VCF that may arise in nested snarls. This may be slow\n"
              << "  -t, --threads INT               Number of threads to use [1]\n"
              << "  -V, --verbose INT               Verbosity level (0=error, 1=warn, 2=info, 3=debug, 4=trace) [2]\n"
              << "  -o, --output FILE               Output directory name\n"
              << "  -h, --help                      Print this help message\n";
}

int main_stoat_vcf(int argc, char* argv[]) {
    
    // Declare variables to hold argument values
    std::string vcf_path, snarl_path, graph_path, dist_path, chromosome_path;

    size_t cycle_threshold = 1;
    size_t children_threshold = 50;
    size_t min_individuals = 0;
    // JEAN this threshold is a bit redundant with children_threshold and cycle_threshold but I guess could be useful if we want to set it lower than (children_threshold * (cycle_threshold+1))
    size_t path_length_threshold = 50;
    std::string output_dir = "output";
    bool only_prepare_snarls = false;
    bool resolve_vcf = false;

    // Parse arguments
    int c;

    static struct option long_options[] = {
        {"vcf", required_argument, 0, 'v'},
        {"snarl", required_argument, 0, 's'},
        {"graph", required_argument, 0, 'g'},
        {"dist", required_argument, 0, 'd'},
        {"reference-chrs", required_argument, 0, 'r'},
        {"children", required_argument, 0, 'i'},
        {"cycle", required_argument, 0, 'y'},
        {"path-length", required_argument, 0, 'l'},
        {"resolve-vcf", no_argument, 0, 'R'},
        {"thread", required_argument, 0, 't'},
        {"verbose", required_argument, 0, 'V'},
        {"output", required_argument, 0, 'o'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    while ((c = getopt_long(argc, argv, "v:s:g:d:r:i:y:l:Rt:V:o:h", long_options, nullptr)) != -1) {
        switch (c) {
            case 'v': vcf_path = optarg; stoat_vcf::check_file(vcf_path); break;
            case 's': snarl_path = optarg; stoat_vcf::check_file(snarl_path); break;
            case 'g': graph_path = optarg; stoat_vcf::check_file(graph_path); break;
            case 'd': dist_path = optarg; stoat_vcf::check_file(dist_path); break;
            case 'r': chromosome_path = optarg; stoat_vcf::check_file(chromosome_path); break;
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
            case 'R':
                resolve_vcf=true;
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
    
    std::filesystem::create_directory(output_dir);
    std::unordered_set<std::string> ref_chrs = (!chromosome_path.empty()) ? stoat_vcf::parse_chromosome_reference(chromosome_path) : std::unordered_set<std::string>{};
    stoat::Logger::instance().setLogFile(output_dir + "/stoat.vcf.log");

    // add command launch in log file
    std::stringstream ss;
    ss << "stoat ";
    for (int i = 0; i < argc; ++i) ss << argv[i] << " ";
    stoat::LOG_SILENTE(ss.str());

    // Enforce valid argument combinations
    if (!graph_path.empty() && !dist_path.empty() && vcf_path.empty() && snarl_path.empty()) {
        //stoat::LOG_TRACE("Case Snarl path decomposition");
        // Case 1: Only graph_path + dist_path
        only_prepare_snarls = true;
    } else if ( vcf_path.empty() || (snarl_path.empty() && (graph_path.empty() || dist_path.empty())) ) {
        stoat::LOG_ERROR("[stoat vcf] " +
            std::string("Invalid argument combination provided.\n") +
            "There are only 3 ways to launch stoat vcf:\n" +
            "Case 1 (snarl path decomposition): graph_path + dist_path\n" +
            "Case 2 (snarl path decomposition and genotyping): graph_path + dist_path + vcf_path\n" +
            "Case 3 (snarl genotyping): snarl_path + vcf_path"
        );
        print_help_vcf();
        return EXIT_FAILURE;
    }

    auto start_total_timer = std::chrono::high_resolution_clock::now();

    // Make an empty SnarlDataCollection, to be filled in or loaded
    // TODO: Double check that these thresholds are doing the right thing
    stoat::SnarlDataCollection snarl_collection(0, children_threshold, path_length_threshold);
    
    // Start tracking with callgrind
#ifdef USE_CALLGRIND
    CALLGRIND_START_INSTRUMENTATION;
#endif

    // we might need the graph at the end also, so this will point to it and ensure it's kept until the end
    std::unique_ptr<handlegraph::PathHandleGraph> graph;
        
    //////////////////////////////////////// Enumerate/load the snarls in the pangenome
    if (!snarl_path.empty()){ // if we have already saved the snarls info, load it
        stoat::LOG_INFO("Loading snarls from " + snarl_path);
        snarl_collection.load_snarl_data_collection(snarl_path);

    } else { // otherwise, find them from the pangenome graph and snarl tree
        stoat::LOG_INFO("Starting snarl decomposition... ");
        auto start_dec_timer = std::chrono::high_resolution_clock::now();

        // Load the snarl tree and graph
        // Tell the IO library about libvg types.
        if (!stoat::io::register_libvg_io()) {
            throw std::runtime_error("error[stoat vgio]: Could not register libvg types with libvgio");
        }

        // Load the graph and make it a PathPositionHandleGraph
        graph = std::move(vg::io::VPKG::load_one<handlegraph::PathHandleGraph>(graph_path));
        bdsg::PathPositionOverlayHelper overlay_helper;
        bdsg::PathPositionHandleGraph* path_position_graph;
        path_position_graph =  overlay_helper.apply(graph.get());

        // Load the distance index
        std::unique_ptr<bdsg::SnarlDistanceIndex> distance_index = std::make_unique<bdsg::SnarlDistanceIndex>();
        distance_index->deserialize(dist_path);

        // Check if chromosomes specified in the --chr file are present in the graph
        for (const auto& chr : ref_chrs) {
            stoat::LOG_TRACE("Sequence name not found in -r/--chr file: " + chr);
            if (!graph->has_path(chr)) {
                throw std::runtime_error("Reference chromosome: " + chr + " not present in graph");
            }
        }
        // JEAN maybe here add to find the reference paths from the graph if ref_chrs was empty

        // The snarl collection requires a sample_haplotypes but we don't want to work on the pangenome's haplotypes here,
        // because we'll work on the samples from the VCF later, so we use an empty set of haplotypes
        std::vector<stoat::sample_hap_t> sample_haplotypes;

        // equivalent to what was done before in stoat vcf: enumerate all walks through a snarl
        snarl_collection.fill_in_snarl_info(*path_position_graph, *distance_index, sample_haplotypes,
            true, //find_alleles_first, doesn't matter in this case
            true, // walks_requested
            [&] (const net_handle_t& snarl, const snarl_info_t& snarl_data, std::vector<PathTraversal>& walks) { // function to fill in walks
                SnarlDataCollection::get_all_walks_through_snarl(*path_position_graph, *distance_index, snarl, snarl_data, walks, cycle_threshold); //TODO: Use Matis's version and write the skipped snarls somewhere
            },
            false, //alleles_requested
            [&] (const net_handle_t& snarl, const snarl_info_t& snarl_data, const std::vector<stoat::sample_hap_t>& all_sample_haplotypes) { //function to find alleles
                return std::vector<size_t>();
            }, 
            false, // sequence_requested 
            ref_chrs, // reference 
            false, //check distances
            "", // Filename to write the intermediate files
            true // Keep the snarls in the collection?
            ); 

        // Always write the snarls
        ofstream snarls_out;
        snarls_out.open(output_dir + "/snarl_info.tsv");
        snarl_collection.write_snarl_data_collection(snarls_out);
        snarls_out.close();

        auto end_dec_timer = std::chrono::high_resolution_clock::now();
        stoat::LOG_INFO("Snarl decomposition took " + std::to_string(std::chrono::duration<double>(end_dec_timer - start_dec_timer).count()) + " s");

        if (only_prepare_snarls) {
            return EXIT_SUCCESS;
        }
    }

    // If there were no references given, fill them in with the references from the snarl collection
    if (ref_chrs.empty()) {
        for (const std::string& ref : snarl_collection.get_reference_names()) {
            ref_chrs.insert(ref);
        }
    }

    //////////////////////////////////////// Go through the vcf, genotype all the snarls ( and save the intermediate file)
    if (!only_prepare_snarls) {
        stoat::LOG_INFO("Retrieving genotypes for all snarls...");
        auto start_gt_timer = std::chrono::high_resolution_clock::now();

        // start reading the VCF to get the sample list
        stoat_vcf::VCFParser vcf_parser(resolve_vcf);
        std::vector<std::string> list_samples = vcf_parser.initialize_parser(vcf_path);

        snarl_collection.genotype_snarls_by_chr_from_vcf(list_samples, vcf_parser);

        // We are done reading through the vcf file so close it
        vcf_parser.close_vcf();

        auto end_gt_timer = std::chrono::high_resolution_clock::now();
        stoat::LOG_INFO("Retrieving snarl genotypes took " + std::to_string(std::chrono::duration<double>(end_gt_timer - start_gt_timer).count()) + " s");
              
        // write the genotypes
        ofstream snarls_out;
        // JEAN would be nice to write in a bgzip file directly, otherwise it might get very big
        stoat::LOG_INFO("Writing genotypes in " + output_dir + "/snarl_genotypes.tsv");
        snarls_out.open(output_dir + "/snarl_genotypes.tsv");
        // JEAN would reduce memory to write the collection while genotyping the snarls, one chr at a time, appending to the output file (or in separate chr files).
        snarl_collection.write_snarl_data_collection(snarls_out);
        snarls_out.close();
    }
    
    auto end_total_timer = std::chrono::high_resolution_clock::now();
    stoat::LOG_INFO("Total time: " + std::to_string(std::chrono::duration<double>(end_total_timer - start_total_timer).count()) + " s");
    return EXIT_SUCCESS;
}

} // end stoat
