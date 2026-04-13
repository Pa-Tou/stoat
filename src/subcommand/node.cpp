#include "node.hpp"

#include <string>
#include <unordered_set>
#include <chrono>
#include <getopt.h>
#include <omp.h>
#include <iostream>
#include <filesystem>

#include <bdsg/overlays/overlay_helper.hpp>
#include <vg/io/vpkg.hpp>

#include "../banner.hpp"
#include "../log.hpp"
#include "../arg_parser.hpp"
#include "../io/register_io.hpp"
#include "../node_data_collection.hpp"
#include "../writer.hpp"
#include "../vcf_parser.hpp"

// #define USE_CALLGRIND

#ifdef USE_CALLGRIND
    #include <valgrind/callgrind.h>
#endif

namespace stoat_command {

// STOAT_VERSION is defined in the CMAKELIST file
void print_help_node() {
    stoat::print_banner(std::string(STOAT_VERSION));
    std::cerr << "Usage: stoat node [options]\n\n"
              << "  -v, --node FILE                  Path to the node file\n"
              << "  -f, --resolve-node               Resolve conflicting calls in the node that may arise in nested nodes. This may be slow\n"
              << "  -t, --threads INT               Number of threads to use [1]\n"
              << "  -V, --verbose INT               Verbosity level (0=error, 1=warn, 2=info, 3=debug, 4=trace) [2]\n"
              << "  -o, --output FILE               Output directory name [stoat_output]\n"
              << "  -u, --no-bgzip                  Don't compress the output file with bgzip\n"
              << "  -a, --ascii                     Print the STOAT ascii art banner\n"
              << "  -h, --help                      Print this help message\n";
}

int main_stoat_node(int argc, char* argv[]) {

    // Declare variables to hold argument values
    std::string vcf_path;

    std::string output_dir = "stoat_output";
    bool resolve_vcf = false;
    bool bgzip_output = true;
    bool ascii = false;

    // Parse arguments
    int c;

    static struct option long_options[] = {
        {"vcf", required_argument, 0, 'v'},
        {"resolve-vcf", no_argument, 0, 'f'},
        {"thread", required_argument, 0, 't'},
        {"verbose", required_argument, 0, 'V'},
        {"output", required_argument, 0, 'o'},
        {"no-bgzip", no_argument, 0, 'u'},
        {"ascii", no_argument, 0, 'a'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    while ((c = getopt_long(argc, argv, "v:ft:V:o:uah", long_options, nullptr)) != -1) {
        switch (c) {
            case 'v': vcf_path = optarg; stoat_vcf::check_file(vcf_path); break;
            case 'a': ascii = true; break;
            case 'f':
                resolve_vcf=true;
                break;
            case 't':
                if (std::stoi(optarg) < 1) {
                    throw std::runtime_error("Error: [stoat node] Number of threads must be > 0");
                }
                omp_set_num_threads(std::stoi(optarg));
                break;
            case 'V': 
                {
                int level = std::stoi(optarg);
                if (level < 0 || level > 4) {
                    throw std::runtime_error("Error: [stoat node] Invalid verbosity level. Use 0=Error, 1=Warn, 2=Info, 3=Debug, 4=Trace");
                }
                stoat::LogLevel logLevel = static_cast<stoat::LogLevel>(level);
                stoat::Logger::instance().setLevel(logLevel);                
                break;
                }
            case 'o': output_dir = optarg; break;
            case 'u':
                bgzip_output = false;
                break;
            case 'h': 
                print_help_node(); 
                return EXIT_SUCCESS; 
            default:
                stoat::LOG_ERROR("[stoat node] Unknown argument");
                print_help_node();
                return EXIT_FAILURE;
        }
    }

    if (argc == 2) {
        print_help_node();
        return EXIT_FAILURE;
    }

    // Enforce valid argument combinations
    // Either we just want to prepare the nodes (from pangenome index files) or genotype those nodes from a node (or both)
    if (!vcf_path.empty()) {
        stoat::LOG_ERROR("[stoat node] " +
            std::string("Invalid argument : vcf_path must be specify"));
        print_help_node();
        return EXIT_FAILURE;
    }

    // create output directory and prepare the log file
    std::filesystem::create_directory(output_dir);
    stoat::Logger::instance().setLogFile(output_dir + "/stoat.node.log");

    // add command launch in log file + print banner
    if (ascii) {
        print_ascii_banner(std::string(STOAT_VERSION));
    } else {
        print_banner(std::string(STOAT_VERSION));
    }

    std::stringstream ss;
    ss << "stoat ";
    for (int i = 0; i < argc; ++i) ss << argv[i] << " ";
    stoat::LOG_SILENTE(ss.str());

    // start the overall timer
    auto start_total_timer = std::chrono::high_resolution_clock::now();

    // Start tracking with callgrind
    #ifdef USE_CALLGRIND
        CALLGRIND_START_INSTRUMENTATION;
    #endif

    stoat::LOG_INFO("Retrieving genotypes for all nodes...");
    auto start_gt_timer = std::chrono::high_resolution_clock::now();

    // start reading the node to get the sample list
    stoat_vcf::vcfParser vcf_parser(resolve_vcf);
    std::vector<std::string> list_samples = vcf_parser.initialize_parser(vcf_path);

    // retrieve genotypes one chromosome at a time
    stoat::NodeDataCollection node_collection();
    node_collection.genotype_node_by_chr_from_vcf(list_samples, vcf_parser);

    // We are done reading through the node file so close it
    vcf_parser.close_node();

    auto end_gt_timer = std::chrono::high_resolution_clock::now();
    stoat::LOG_INFO("Retrieving node genotypes took " + std::to_string(std::chrono::duration<double>(end_gt_timer - start_gt_timer).count()) + " s");

    // write the genotypes
    std::string genotype_path = output_dir + "/node_genotypes.tsv";
    if (bgzip_output) {
        genotype_path += ".gz";
    }

    std::shared_ptr<stoat::Writer> gt_writer;
    if ((genotype_path.compare(genotype_path.length()-3, 3, ".gz") == 0) ||
        (genotype_path.compare(genotype_path.length()-4, 4, ".bgz") == 0)) {
        gt_writer.reset(new BgzWriter(genotype_path));
    } else {
        gt_writer.reset(new StdWriter(genotype_path));
    }

    stoat::LOG_INFO("Writing genotypes in " + genotype_path);

    // JEAN would reduce memory to write the collection while genotyping the nodes, one chr at a time, appending to the output file (or in separate chr files).
    auto start_writegt_timer = std::chrono::high_resolution_clock::now();
    node_collection.write_node_data_collection(*gt_writer);
    gt_writer->close();

    auto end_writegt_timer = std::chrono::high_resolution_clock::now();
    stoat::LOG_INFO("Writing genotypes took " + std::to_string(std::chrono::duration<double>(end_writegt_timer - start_writegt_timer).count()) + " s");
    
    auto end_total_timer = std::chrono::high_resolution_clock::now();
    stoat::LOG_INFO("Total time: " + std::to_string(std::chrono::duration<double>(end_total_timer - start_total_timer).count()) + " s");
    return EXIT_SUCCESS;
}

} // end stoat_command