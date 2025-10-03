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

#include "../graph_path_association_finder.hpp"
#include "../io/register_io.hpp"
#include "../post_processing.hpp"
#include "../arg_parser.hpp"

//#define USE_CALLGRIND

#ifdef USE_CALLGRIND
    #include <valgrind/callgrind.h>
#endif

using namespace std;
namespace stoat_command {

void print_help_graph() {
    std::cerr << "usage: stoat graph -g [graph] -d [distance index] -b [phenotype file] [options]" << endl
        << "Find associated variants based on the haplotype paths in the graph"<< endl
        << "Computing the snarls from the distance index may be slow, so they can be saved with -s." << endl
        << "Requires either -b to compute the associations, or -s to save the snarls in the graph. Or both to do both" << endl 
        << endl
        << "input:" << endl
        << "  -g, --graph FILE                   Use this graph (only hash graph works for now) (required)" << endl
        << "  -d, --distance-index FILE          Use this distance index (required if -s is not given)" << endl
        << "  -b, --binary-pheno FILE            A tsv of \"FID IID PHENO\" for family id, sample name, and phenotype (1 or 2), one per line (required if -s is not given)" << endl
        << "                                     If this is not give, then -s is required to save the snarls." << endl
        << "  -s, --snarls FILE                  If this is file is empty, then save the snarl paths in the graph to this file. (required if -b is not given) " << endl
        << "                                     If this is file is not empty, then load the snarl paths from this file instead of recomputing them. " << endl
        << endl
        << "output:" << endl
        << "  -o, --output DIR                   Output directory name [output]" << endl
        << "  -O, --output-format NAME           The format of the output (tsv / fasta) [tsv]" << endl
        << "                                     Output will be written to DIR/binary_table_graph.tsv or DIR/binary_output.fasta" << endl
        << "options:" << endl
        << "  -t, --threads N                    Number of threads to use" << endl
        << "  -T, --test NAME                    Which test will be used to determine association (exact / chi2) [exact]" << endl
        //<< "  -p, --p-value-threshold FLOAT    What is the threshold p-value to be considered significant? [0.05]" << endl
        //<< "                                   When used with multiple testing, discard any p-value above this threshold without doing multiple testing" << endl
        << "  -V, --verbose INT                  Verbosity level (0=error, 1=warn, 2=info, 3=debug, 4=trace)" << endl
        //<< "  -m, --method NAME                What method is used to find associations? (paths) [paths]" << endl
        << "  -l, --allele-size-limit INT        Don't report variants smaller than this [0]" << endl
        << "  -r, --reference-sample NAME        If there is no reference in the graph, use this sample as the reference" << endl
        << "  -h, --help                         Print this help message" << endl;
}

int main_stoat_graph(int argc, char *argv[], stoat::LogLevel &verbosity) {

    if (argc <= 1) {
        print_help_graph();
        return EXIT_FAILURE;
    }

    std::string graph_name;
    std::string distance_name;
    size_t allele_size_limit = 0;
    //double p_value = 0.05;
    std::string method_name = "paths";
    std::string test_method = "exact";
    std::string reference_sample;
    std::string samples_filename;
    std::string snarls_filename;
    std::vector<std::string> samples;
    std::string output_format= "tsv";
    std::string output_dir="output";

    int c = 0;
    optind = 1;
    while (true) {
        static struct option long_options[] =
            {
                {"graph", required_argument, 0, 'g'},
                {"distance-index", required_argument, 0, 'd'},
                {"allele-size-limit", required_argument, 0, 'l'},
                {"threads", required_argument, 0, 't'},
                {"test", required_argument, 0, 'T'},
                //{"p-value", required_argument, 0, 'p'},
                //{"method", required_argument, 0, 'm'},
                {"reference-sample", required_argument, 0, 'r'},
                {"binary-pheno", required_argument, 0, 'b'},
                {"snarls", required_argument, 0, 's'},
                {"output", required_argument, 0, 'o'},
                {"output-format", required_argument, 0, 'O'},
                {"verbose", required_argument, 0, 'V'},
                {"skip-bh-correction", no_argument, 0, 'B'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "g:d:l:t:T:r:b:s:V:o:O:h",
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
            case 'T':
                test_method = optarg;
                break;
            //case 'p':
            //    p_value = std::stof(optarg);
            //    break;
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
            //case 'm':
            //    method_name = optarg;
            //    break;
            case 'r':
                reference_sample = optarg;
                break;
            case 'b':
                samples_filename = optarg;
                break;
            case 's':
                snarls_filename = optarg;
                break;
            case 'o':
                output_dir = optarg;
                break;
            case 'O':
                output_format = optarg;
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

    // Check that the inputs are ok
    if (graph_name.empty()) {
        stoat::LOG_ERROR("[stoat graph] stoat graph requires a graph file");
        return EXIT_FAILURE; 
    }

    if (distance_name.empty() && snarls_filename.empty()) {
        stoat::LOG_ERROR("[stoat graph] stoat graph requires a distance index file or saved snarls");
        return EXIT_FAILURE; 
    }

    if (output_format != "tsv" && output_format != "fasta") {
        stoat::LOG_ERROR("[stoat graph] invalid output format " + output_format);
        return EXIT_FAILURE; 
    }

    // Make the output directory
    std::filesystem::create_directory(output_dir);
    stoat::Logger::instance().setLogFile(output_dir + "/stoat_graph.log");

    // Load the samples from a file
    if (samples_filename.empty() && snarls_filename.empty()) {
        stoat::LOG_ERROR("[stoat graph]: stoat graph requires samples of interest to compute associations or an empty snarls file to pre-compute the snarls");
        return EXIT_FAILURE; 
    }


    // If we want to save the snarls, then double check that the file is empty
    // If we want to use the snarls, then double check that the file is not empty.
    // Either way, tell the user what we're doing
    bool save_snarls = false;
    bool load_snarls = false;
    if (!snarls_filename.empty()) {
        if (std::filesystem::exists(snarls_filename)) {
            if (samples_filename.empty()) {
                // If this file already exists but we wanted to write it (no -b), then tell the user and error
                stoat::LOG_ERROR("[stoat graph] --binary-pheno was not given but --snarls is not empty.\n\tstoat will not overwrite the snarls file. Please delete it if you want to re-write it.");
            } else {
                // Otherwise, we are going to use the precomputed snarls
                stoat::LOG_INFO("[stoat graph] using precomputed snarl paths from " + snarls_filename);
                load_snarls = true;
            }
        } else {
            //conversely, if the file doesn't exist and we don't have a distance index then also error
            if (distance_name.empty()) {
                stoat::LOG_ERROR("[stoat graph] --snarls file is empty. --distance-index is required to compute the snarls");
            } else {
                // Otherwise, we are going to write snarls
                stoat::LOG_INFO("[stoat graph] write computed snarl paths to " + snarls_filename);
                save_snarls = true;
                // TODO: Make sure that the directory exists
            }
        }
    }

    
    std::vector<bool> phenotypes;
    if (!samples_filename.empty()) {
        phenotypes = stoat_vcf::parse_binary_pheno(samples_filename, samples);
    }
    std::pair<std::set<std::string>, std::set<std::string>> sample_sets;

    for (size_t i = 0 ; i < samples.size() ; i++ ) {
        if (phenotypes[i]) {
            sample_sets.first.emplace(samples[i]);
        } else {
            sample_sets.second.emplace(samples[i]);
        }
    }

    //Output trace info
    stoat::LOG_TRACE("Truth sample set 1: ");
    for (auto& sample : sample_sets.first) {
        stoat::LOG_TRACE("\t" + sample);
    }
    stoat::LOG_TRACE("Truth sample set 2: ");
    for (auto& sample : sample_sets.second) {
        stoat::LOG_TRACE("\t" + sample);
    }


    // Tell the IO library about libvg types.
    if (!stoat::io::register_libvg_io()) {
        stoat::LOG_ERROR("[stoat vgio] Could not register libvg types with libvgio");
        return EXIT_FAILURE;
    }

    auto start_1 = std::chrono::high_resolution_clock::now();
    stoat::LOG_INFO("Start Sample haplotype analysis...");

    // Load the graph and make it a PathPositionHandleGraph
    unique_ptr<handlegraph::PathHandleGraph> path_graph = vg::io::VPKG::load_one<handlegraph::PathHandleGraph>(graph_name);


    /// For the PathPositionHandleGraph, haplotypes are not indexed automatically so we need to give additional path names
    /// that we want to be included in the index.
    /// We also want a list of all samples and haplotypes that we include in the analysis
    /// Get both at the same time by going through all paths in the graph and checking if the sample is one of our samples of interest

    // Get a list of paths to include in the path position overlay
    std::unordered_set<std::string> paths_set;

    // A set of the samples+haplotypes in the graph that match the ones from the phenotype file
    std::set<stoat::sample_hap_t> all_sample_haplotypes;

    path_graph->for_each_path_matching(nullptr, nullptr, nullptr, [&] (handlegraph::path_handle_t path) {
        std::string sample_name = stoat::get_sample_name_from_path(*path_graph, path);
        if (samples_filename.empty() || sample_sets.first.count(sample_name) == 1 || sample_sets.second.count(sample_name) == 1) {
            paths_set.emplace(path_graph->get_path_name(path));
            all_sample_haplotypes.emplace(stoat::get_sample_and_haplotype(*path_graph, path));
        }
        return true;
    });

    if (stoat::Logger::instance().at_level(stoat::LogLevel::Error)) {
        stoat::Logger::instance().log_assert(stoat::LogLevel::Error, 
                                             all_sample_haplotypes.size() >= sample_sets.first.size() + sample_sets.second.size(), 
                                             "there are more samples given than haplotypes in the graph");
    }

    bdsg::PathPositionOverlayHelper overlay_helper;
    bdsg::PathPositionHandleGraph* graph = overlay_helper.apply(path_graph.get(), paths_set);


    bdsg::SnarlDistanceIndex* distance_index_ptr = nullptr;
    bdsg::SnarlDistanceIndex distance_index;
    if (!distance_name.empty()) {
        // Load the distance index
        cerr << " LOAD DISTANCE INDEX" << endl;
        distance_index.deserialize(distance_name);
        distance_index_ptr = &distance_index;
    }

    auto end_1 = std::chrono::high_resolution_clock::now();
    stoat::LOG_INFO("Sample haplotype time : " + std::to_string(std::chrono::duration<double>(end_1 - start_1).count()) + " s");
    stoat::LOG_INFO("Start GWAS analysis...");
    auto start_2 = std::chrono::high_resolution_clock::now();

    string filename = output_dir + "/";
    if (output_format == "tsv") {
        filename += "binary_table_graph.tsv";
    } else if (output_format == "fasta") {
        filename += "binary_output.fasta";
    }

    // Get the out streams
    std::ofstream out_stream;
    out_stream.open(filename);


    // Make the partitioner
    std::shared_ptr<stoat_graph::SnarlTraverserAndPartitioner> partitioner;
    if (method_name == "paths") {
        partitioner.reset(new stoat_graph::SnarlTraverserAndPathPartitioner(all_sample_haplotypes, distance_index_ptr, reference_sample, save_snarls));
        if (load_snarls) {
            partitioner->deserialize(snarls_filename, *graph);
        }
    } else {
        stoat::LOG_ERROR("[stoat graph] unknown method " + method_name);
        return EXIT_FAILURE; 
    }

#ifdef USE_CALLGRIND
    CALLGRIND_START_INSTRUMENTATION;
#endif

    if (!samples_filename.empty()) {
        // If we actually want to do the analysis
        // TODO: Double check that the distance index isn't actually used
        stoat_graph::AssociationFinder af (*graph, 
                                       partitioner,
                                       sample_sets,
                                       reference_sample,
                                       test_method,
                                       output_format,
                                       allele_size_limit,
                                       out_stream);
        af.test_snarls();
    } else if (save_snarls) {
        // If we just want to precompute the snarls
        // The first function says which snarls to include - all of them
        // The second function says what to do with the partitions - nothing because it saves them automatically 
        partitioner->for_each_snarl_partition(*graph, 
            [](const bdsg::SnarlDistanceIndex& dist_index, const net_handle_t& net) {return true;},
            [&](const auto& snarl_info) { return; });
    }

    if (save_snarls) {
        // Save the snarls
        partitioner->serialize(snarls_filename);
    }

    //Close streams
    out_stream.close();

    auto end_2 = std::chrono::high_resolution_clock::now();
    stoat::LOG_INFO("GWAS time analysis : " + std::to_string(std::chrono::duration<double>(end_2 - start_2).count()) + " s");
    stoat::LOG_INFO("Total time : " + std::to_string(std::chrono::duration<double>(end_2 - start_1).count()) + " s");
    return EXIT_SUCCESS;
}
} //end namespace
