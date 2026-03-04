#include "graph.hpp"
#include <iostream>
#include <string>
#include <getopt.h>
#include <omp.h>
#include <filesystem>

#include <bdsg/snarl_distance_index.hpp>
#include <bdsg/overlays/overlay_helper.hpp>
#include <handlegraph/path_handle_graph.hpp>
#include <vg/io/vpkg.hpp>

#include "../snarl_data_collection.hpp"
#include "../path_partitioner.hpp"
#include "../io/register_io.hpp"
#include "../arg_parser.hpp"

//#define USE_CALLGRIND

#ifdef USE_CALLGRIND
    #include <valgrind/callgrind.h>
#endif

namespace stoat_command {

void print_help_graph() {
    std::cerr << "usage: stoat graph -g [graph] -d [distance index] [options]" << std::endl
        << "Retrieves snarl genotypes based on the haplotype paths present in the graph"<< std::endl
        << std::endl
        << "input:" << std::endl
        << "  -g, --graph FILE                   Use this graph (only hash graph works for now) (required)" << std::endl
        << "  -d, --distance-index FILE          Use this distance index (required if -s is not given)" << std::endl
        << std::endl
        << "output:" << std::endl
        << "  -o, --output DIR                   Output directory name [stoat_output]" << std::endl
        << "  -L, --allele-lengths               Find the lengths of alleles (they will be NA without this flag). This makes stoat slow." << std::endl
        << "  -u, --no-bgzip                     Don't compress the output file with bgzip\n"
        << std::endl
        << "options:" << std::endl
        << "  -t, --threads N                    Number of threads to use" << std::endl
        << "  -V, --verbose INT                  Verbosity level (0=error, 1=warn, 2=info, 3=debug, 4=trace)" << std::endl
        << "  -l, --allele-size-limit INT        Don't report variants smaller than this [0]" << std::endl
        << "  -r, --reference-chrs FILE          Path to the chromosome reference file, one path name per line. These paths must be REFERENCE- or GENERIC-sense paths (check with vg paths -M)." << std::endl
        << "                                     If not given, use any reference-sense paths in the graph as the references" << std::endl
        << "  -h, --help                         Print this help message" << std::endl;
}

int main_stoat_graph(int argc, char *argv[]) {

    if (argc <= 1) {
        print_help_graph();
        return EXIT_FAILURE;
    }

    std::string graph_name;
    std::string distance_name;
    size_t allele_size_limit = 0;
    std::string reference_path;
    std::vector<std::string> samples;
    std::string output_dir="stoat_output";
    bool bgzip_output = true;

    bool find_allele_lengths = false; 

    int c = 0;
    optind = 1;
    while (true) {
        static struct option long_options[] =
            {
                {"graph", required_argument, 0, 'g'},
                {"distance-index", required_argument, 0, 'd'},
                {"allele-size-limit", required_argument, 0, 'l'},
                {"threads", required_argument, 0, 't'},
                {"reference-chrs", required_argument, 0, 'r'},
                {"output", required_argument, 0, 'o'},
                {"allele-length", no_argument, 0, 'L'},
                {"verbose", required_argument, 0, 'V'},
                {"no-bgz", no_argument, 0, 'u'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "g:d:l:t:r:V:o:Luh",
                        long_options, &option_index); 
        if (c == -1) {
            break;
        }
        switch (c) {
            case 'g':
                graph_name = optarg;
                break;
            case 'd':
                distance_name = optarg;
                break;
            case 'l':
                if (std::stoi(optarg) < 0) {
                    stoat::LOG_ERROR("[stoat graph] Number of allele size limit must be >= 0");
                    return EXIT_FAILURE;
                }
                allele_size_limit = std::stoi(optarg);
                break;
            case 't':
                if (std::stoi(optarg) < 1) {
                    stoat::LOG_ERROR("[stoat graph] Number of threads must be > 0");
                    return EXIT_FAILURE;
                }
                omp_set_num_threads(std::stoi(optarg));
                break;
            case 'V':
                {
                int level = std::stoi(optarg);
                if (level < 0 || level > 4) {
                    stoat::LOG_ERROR("[stoat graph] Invalid verbosity level. Use 0=Error, 1=Warn, 2=Info, 3=Debug, 4=Trace");
                    return EXIT_FAILURE;
                }
                stoat::LogLevel logLevel = static_cast<stoat::LogLevel>(level);
                stoat::Logger::instance().setLevel(logLevel);                
                break;
                }
            case 'r':
                reference_path = optarg;
                break;
            case 'o':
                output_dir = optarg;
                break;
            case 'L':
                find_allele_lengths = true;
                break;
            case 'u':
                bgzip_output = false;
                break;
            case 'h':
                print_help_graph();
                return EXIT_SUCCESS;
            default:
                stoat::LOG_ERROR("[stoat graph] Unknown argument");
                print_help_graph();
                return EXIT_FAILURE;
        }
    }

    std::string snarls_filename = output_dir + "/snarl_genotypes.tsv";
    if (bgzip_output) {
        snarls_filename += ".gz";
    }

    //////////////////////////////////////////////////// Check the inputs and outputs and logs

    // Check that the inputs are ok
    if (graph_name.empty()) {
        stoat::LOG_ERROR("[stoat graph] stoat graph requires a graph file");
        return EXIT_FAILURE; 
    }

    if (distance_name.empty()) {
        stoat::LOG_ERROR("[stoat graph] stoat graph requires a distance index file");
        return EXIT_FAILURE; 
    }

    // Make the output directory
    std::filesystem::create_directory(output_dir);
    stoat::Logger::instance().setLogFile(output_dir + "/stoat.graph.log");

    // add command launch in log file
    std::stringstream ss;
    ss << "stoat ";
    for (int i = 0; i < argc; ++i) ss << argv[i] << " ";
    stoat::LOG_SILENTE(ss.str());

    // Tell the IO library about libvg types.
    if (!stoat::io::register_libvg_io()) {
        stoat::LOG_ERROR("[stoat vgio] Could not register libvg types with libvgio");
        return EXIT_FAILURE;
    }

    ///////////////////////////////////////////////// Load the graph and stuff

    auto start_1 = std::chrono::high_resolution_clock::now();
    stoat::LOG_INFO("Loading graph and preparing indexes...");

    // Load the graph and make it a PathPositionHandleGraph
    std::unique_ptr<handlegraph::PathHandleGraph> handle_graph = vg::io::VPKG::load_one<handlegraph::PathHandleGraph>(graph_name);

    /// For the PathPositionHandleGraph, haplotypes are not indexed automatically so we need to give additional path names
    /// that we want to be included in the index.
    /// We used to select specific samples/haplotypes but now include all paths in this "genotype retrieval" subcommand.

    // Get a list of paths to include in the path position overlay
    std::unordered_set<std::string> paths_set;

    // A set of the samples+haplotypes in the graph
    std::vector<stoat::sample_hap_t> all_sample_haplotypes;

    // go through all paths in the pangenome and save them
    handle_graph->for_each_path_matching(nullptr, nullptr, nullptr, [&] (handlegraph::path_handle_t path) {
        std::string sample_name = stoat::get_sample_name_from_path(*handle_graph, path);
        // Get the sample haplotypes that we want
        // TODO: For now, if we are saving the snarls to be used again, force all samples to be included. This may change when this is a separate subcommand
        // JEAN right, maybe we do want to be able to say which ones to include (new input file or prefix selection?)
        paths_set.emplace(handle_graph->get_path_name(path));
        all_sample_haplotypes.emplace_back(stoat::sample_hap_t(*handle_graph, path));
        return true;
    });

    bdsg::PathPositionOverlayHelper overlay_helper;
    bdsg::PathPositionHandleGraph* path_position_graph = overlay_helper.apply(handle_graph.get(), paths_set);

    // Load the distance index
    bdsg::SnarlDistanceIndex distance_index;
    if (!distance_name.empty()) {
        // Load the distance index
        distance_index.deserialize(distance_name);
    }

    // Get the reference sample names
    std::unordered_set<std::string> reference_samples = (!reference_path.empty()) ? stoat_vcf::parse_chromosome_reference(reference_path) : std::unordered_set<std::string>{};

    //////////////////////////////// Make the snarls file and load it if possible
    // If it is being built, it will count towards the time of analysis

    // TODO: Get these from the command line, infinite for now (which may also be fine)
    size_t snarl_child_limit = std::numeric_limits<size_t>::max();
    size_t walk_steps_limit = std::numeric_limits<size_t>::max();
    SnarlDataCollection snarl_collection(allele_size_limit, snarl_child_limit, walk_steps_limit);
    
    ////////////////////////////////////////////////// Start doing work

    auto end_1 = std::chrono::high_resolution_clock::now();
    stoat::LOG_INFO("Graph loading time : " + std::to_string(std::chrono::duration<double>(end_1 - start_1).count()) + " s");
    stoat::LOG_INFO("Start finding alleles in snarls...");
    auto start_2 = std::chrono::high_resolution_clock::now();

#ifdef USE_CALLGRIND
    CALLGRIND_START_INSTRUMENTATION;
#endif

    // prepare a Writer for the collection
    std::shared_ptr<stoat::Writer> snarl_writer;
    if ((snarls_filename.compare(snarls_filename.length()-3, 3, ".gz") == 0) ||
        (snarls_filename.compare(snarls_filename.length()-4, 4, ".bgz") == 0)) {
        snarl_writer.reset(new BgzWriter(snarls_filename));
    } else {
        snarl_writer.reset(new StdWriter(snarls_filename));
    }

    // Find and write the information (inc. paths and genotypes) for each snarl in the index
    snarl_collection.fill_in_snarl_info(*path_position_graph, distance_index, all_sample_haplotypes, 
                                        true, // Find the sets of samples in each allele (walk through the snarl) before finding the walks themselves
                                        find_allele_lengths, // find walks (used for sequences or for lengths)
                                        [&] (const net_handle_t& snarl, const snarl_info_t& snarl_data, //Function to find the walks
                                             std::vector<PathTraversal>& walks) {
                                            // If we actually need the walks to get the sequences or the lengths
                                            return SnarlDataCollection::get_walks_from_alleles(*path_position_graph, distance_index, snarl, snarl_data, walks);
                                        },
                                        true, // find the alleles
                                        // Function to find the alleles 
                                        [&] (const net_handle_t& snarl, const snarl_info_t& snarl_data,
                                             const std::vector<stoat::sample_hap_t>& sample_haplotypes) {
                                            return stoat_graph::partition_embedded_paths_in_snarl(*path_position_graph, distance_index, snarl, sample_haplotypes);
                                        },
                                        false, // find the sequences, only for fasta format
                                        reference_samples,
                                        distance_index.has_distances(),
                                        *snarl_writer, // Filename to write the snarls to
                                        false // Keep the snarls in the collection?
                                        );

    // we're done writing the snarl info, close Writer
    snarl_writer->close();
    
    auto end_2 = std::chrono::high_resolution_clock::now();
    stoat::LOG_INFO("Snarl parsing time : " + std::to_string(std::chrono::duration<double>(end_2 - start_2).count()) + " s");
    stoat::LOG_INFO("Total time : " + std::to_string(std::chrono::duration<double>(end_2 - start_1).count()) + " s");
    return EXIT_SUCCESS;
}
} //end namespace
