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
#include "../writer.hpp"
#include "../binary_table.hpp"
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
        << "Find associated variants based on the haplotype paths present in the graph"<< endl
        << "Computing the snarls from the distance index may be slow, so they can be saved or loaded with -s." << endl
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
        << "  -T, --test NAME                    Which test will be used to determine association (exact / chi2) [chi2]" << endl
        //<< "  -p, --p-value-threshold FLOAT    What is the threshold p-value to be considered significant? [0.05]" << endl
        //<< "                                   When used with multiple testing, discard any p-value above this threshold without doing multiple testing" << endl
        << "  -V, --verbose INT                  Verbosity level (0=error, 1=warn, 2=info, 3=debug, 4=trace)" << endl
        //<< "  -m, --method NAME                  What method is used to find associations? (paths) [paths]" << endl
        << "  -M, --maf FLOAT                    Only consider a snarl if the allele frequencies of at least two alleles are greater than FLOAT [0.05]" << endl
        << "  -I, --min-individuals INT          If there are fewer than INT individuals/samples in a snarl, then ignore the snarl [1]\n"
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
    std::string test_method = "chi2";
    std::string reference_sample;
    std::string samples_filename;
    std::string snarls_filename;
    std::vector<std::string> samples;
    std::string output_format= "tsv";
    std::string output_dir="output";

    double maf_threshold = 0.05;
    size_t min_individuals = 1;

    int c = 0;
    optind = 1;
    while (true) {
        static struct option long_options[] =
            {
                {"graph", required_argument, 0, 'g'},
                {"distance-index", required_argument, 0, 'd'},
                {"allele-size-limit", required_argument, 0, 'l'},
                {"maf", required_argument, 0, 'M'},
                {"min-individuals", required_argument, 0, 'I'},
                {"threads", required_argument, 0, 't'},
                {"test", required_argument, 0, 'T'},
                // {"p-value", required_argument, 0, 'p'},
                // {"method", required_argument, 0, 'm'},
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
        c = getopt_long(argc, argv, "g:d:l:t:T:r:b:s:V:I:M:o:O:h",
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
            case 'M':
                maf_threshold = std::stod(optarg);
                if (maf_threshold < 0 || maf_threshold > 1) {
                    throw std::runtime_error("Error: [stoat graph] MAF must be in [0,1]");
                }
                break;
             case 'I':
                 min_individuals = std::stoi(optarg);
                 if (min_individuals < 2) {
                     throw std::runtime_error("Error: [stoat graph] min_individuals threshold must be > 1");
                 }
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

    //////////////////////////////////////////////////// Check the inputs and outputs and logs

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

    // add command launch in log file
    std::stringstream ss;
    ss << "stoat ";
    for (int i = 0; i < argc; ++i) ss << argv[i] << " ";
    stoat::LOG_SILENTE(ss.str());

    // Load the samples from a file
    if (samples_filename.empty() && snarls_filename.empty()) {
        stoat::LOG_ERROR("[stoat graph]: stoat graph requires samples of interest to compute associations or an empty snarls file to pre-compute the snarls");
        return EXIT_FAILURE; 
    }

    // Decide if we want to save or load the files
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
                return EXIT_FAILURE;
            } else {
                // Otherwise, we are going to use the precomputed snarls
                stoat::LOG_INFO("[stoat graph] using precomputed snarl paths from " + snarls_filename);
                load_snarls = true;
            }
        } else {
            //conversely, if the file doesn't exist and we don't have a distance index then also error
            if (distance_name.empty()) {
                stoat::LOG_ERROR("[stoat graph] --snarls file is empty. --distance-index is required to compute the snarls");
                return EXIT_FAILURE;
            } else {
                // Otherwise, we are going to write snarls
                stoat::LOG_INFO("[stoat graph] write computed snarl paths to " + snarls_filename);
                save_snarls = true;
                // TODO: Make sure that the directory exists
            }
        }
    }

    
    // Load the phenotypes
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

    if (method_name != "paths") {
        stoat::LOG_ERROR("[stoat graph] unknown method " + method_name);
        return EXIT_FAILURE; 
    }

    ///////////////////////////////////////////////// Load the graph and stuff

    auto start_1 = std::chrono::high_resolution_clock::now();
    stoat::LOG_INFO("Loading graph and preparing indexes...");

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


    // Load the distance index
    bdsg::SnarlDistanceIndex* distance_index_ptr = nullptr;
    bdsg::SnarlDistanceIndex distance_index;
    if (!distance_name.empty()) {
        // Load the distance index
        distance_index.deserialize(distance_name);
        distance_index_ptr = &distance_index;
    }

    string filename = output_dir + "/";
    if (output_format == "tsv") {
        filename += "binary_table_graph.tsv";
    } else if (output_format == "fasta") {
        filename += "binary_output.fasta";
    }

    // Get the out streams
    std::ofstream out_stream;
    out_stream.open(filename);

    //////////////////////////////// Make the snarls file and load it if possible
    // If it is being built, it will count towards the time of analysis

    // TODO: Get these from the command line, infinite for now (which may also be fine)
    size_t snarl_child_limit = std::numeric_limits<size_t>::max();
    size_t walk_steps_limit = std::numeric_limits<size_t>::max();
    SnarlDataCollection snarl_collection(allele_size_limit, snarl_child_limit, walk_steps_limit);
    if (load_snarls) {
        std::ifstream in_snarls;
        in_snarls.open(snarls_filename);
        snarl_collection.load_snarl_data_collection(in_snarls);
        in_snarls.close();
    }

    ////////////////////////////////////////////////// Start doing work


    auto end_1 = std::chrono::high_resolution_clock::now();
    stoat::LOG_INFO("Graph loading time : " + std::to_string(std::chrono::duration<double>(end_1 - start_1).count()) + " s");
    stoat::LOG_INFO("Start GWAS analysis...");
    auto start_2 = std::chrono::high_resolution_clock::now();

#ifdef USE_CALLGRIND
    CALLGRIND_START_INSTRUMENTATION;
#endif

    ////////////////// If we didn't load the snarls, then we need to calculate them here
    if (!load_snarls) {
        snarl_collection.fill_in_snarl_info(*graph, distance_index, 
                                            true, // Find the sets of samples in each allele (walk through the snarl) before finding the walks themselves
                                            true, // find walks
                                            SnarlDataCollection::get_walks_from_sample_sets, // Function to find the walks 
                                            true, // find the sample sets
                                            // Function to find the sample sets 
                                            [&] (const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                                                 const net_handle_t& snarl, const snarl_info_t& snarl_data,
                                                 std::vector<std::set<sample_hap_t>>& sample_sets_by_allele) {
                                                stoat_graph::partition_embedded_paths_in_snarl(graph, distance_index, snarl, all_sample_haplotypes, sample_sets_by_allele);
                                            },
                                            true, // find the sequences
                                            reference_sample,
                                            distance_index.has_distances());
        if (save_snarls) {
            ofstream out_snarls;
            out_snarls.open(snarls_filename);
            snarl_collection.write_snarl_data_collection(out_snarls);
            out_snarls.close();
        }
    }

    ////////////////////////////////// Now do the stastistics and write the output

    // Make a tester
    stoat::FisherKhi2 fisher_chi2_tester;

    if (!samples_filename.empty()) {
        // If we actually want to do the analysis

        ///////////// Write the header, if necessary
        if (output_format == "tsv") {
            stoat::write_binary_header(out_stream);
        }

        snarl_collection.for_each_snarl([&](const snarl_info_t& snarl_info){
            // For each snarl, get the genotype/phenotype matrix, do the statistics, and write the output


            // Declare a bunch of strings that are needed for the output
            string group_paths = "NA";
            string fastfisher_p_value = "NA";
            string chi2_p_value = "NA"; 

            // Should we write the output? Don't always if the exact test fails
            bool write_output = false;

            ////////////// Do statistics
            if (test_method == "exact") {
                // If we only want to know if there is one allele that matches exactly one of the phenotype groups
                for (const std::set<sample_hap_t>& sample_hap_partition : snarl_info.sample_sets_by_allele) {
                
                    // Make a set of just the sample names
                    // TODO: This isn't super efficient
                    std::set<std::string> partition;
                    for (const sample_hap_t& sample : sample_hap_partition) {
                        partition.emplace(sample.sample);
                    }
                
                    // If one partition exactly matches one group we want, then all other partitions combined (including things not in the snarl) will match
                    // the other.
                    // TODO: This could be better but I don't think it's worth working on it yet
                    if (partition == sample_sets.first || partition == sample_sets.second) {
                        // There was an exact match
                        write_output = true;
                    }
                }

            } else if (test_method == "chi2") {


                // Fill in the phenotype/genotype count vectors. Each item in these vectors is an allele (/walk through the snarl/partition of samples)
                // One vector for each phenotype
                std::vector<size_t> sample_count_by_allele1(snarl_info.sample_sets_by_allele.size(), 0);
                std::vector<size_t> sample_count_by_allele2(snarl_info.sample_sets_by_allele.size(), 0);
                
                // How many individuals/samples are included?
                std::unordered_set<std::string> seen_samples;

                // Put the sample sets into genotype vectors (one entry in the vector is one allele/path/genotype)
                for (size_t i = 0 ; i < snarl_info.sample_sets_by_allele.size() ; i++) {
                    // For each allele
                    for (const sample_hap_t& sample : snarl_info.sample_sets_by_allele[i]) {
                        // Count each sample with this allele as being in phenotype group 1 or 2
                        if (sample_sets.first.count(sample.sample) == 1) {
                            sample_count_by_allele1[i]++;
                        } else if (sample_sets.second.count(sample.sample) == 1) {
                            sample_count_by_allele2[i]++;
                        }
                        seen_samples.insert(sample.sample);
                    }
                }
                if (!stoat::filtration_binary_table(sample_count_by_allele1, sample_count_by_allele2, seen_samples.size(), min_individuals, maf_threshold)) {
                    // TODO: This could do what pangwas was doing to keep track of only good p-values instead of writing everything

                    //Get a bunch of strings that get used for the output
                    // TODO: This function should probably be part of the output function
                
                    //Get a string representing the number of samples of each phenotype with each allele
                    group_paths = stoat_vcf::format_group_paths(sample_count_by_allele1, sample_count_by_allele2);
                
                    // Run the statistical test
                    std::tie(chi2_p_value, fastfisher_p_value) = fisher_chi2_tester.fisher_khi2(sample_count_by_allele1, sample_count_by_allele2);

                    write_output = true;
                
                }
            }
            /////////////////////////////// Write the output
            if (write_output) {
                if (output_format == "tsv") {
                    # pragma omp critical (out_associated) 
                    {
                        // Leave adjusted p-value blank, to be filled in later
                        stoat::write_binary(out_stream, snarl_info, fastfisher_p_value, chi2_p_value, group_paths);
                    }
                } else if (output_format == "fasta") {
                
                    # pragma omp critical (out_associated) 
                    {
                        stoat::write_fasta(out_stream, *graph, distance_index, snarl_info);
                    }
                }
            }
            return;

        });
    }

    //Close streams
    out_stream.close();

    auto end_2 = std::chrono::high_resolution_clock::now();
    stoat::LOG_INFO("GWAS time analysis : " + std::to_string(std::chrono::duration<double>(end_2 - start_2).count()) + " s");
    stoat::LOG_INFO("Total time : " + std::to_string(std::chrono::duration<double>(end_2 - start_1).count()) + " s");
    return EXIT_SUCCESS;
}
} //end namespace
