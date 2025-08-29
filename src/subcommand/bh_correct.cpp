#include "bh_correct.hpp"

#include <iostream>
#include <string>
#include <getopt.h>
#include <omp.h>
#include <filesystem>

#include "../post_processing.hpp"

using namespace std;
namespace stoat_command {

void print_help_bh_correct() {
    std::cerr << "usage: stoat BHcorrect [options] " << endl
              << endl
              << "options:" << endl
              << "  -t, --tsv FILE                  The TSV file to be processed" << endl
              << "  -p, --p-index N                 The column of the p-value in the tsv (1-indexed)" << endl
              << "  -a, --adjusted-p-index N        The column number of the adjusted p-value in the tsv (1-indexed)" << endl
              << "  -v, --top-variant-file FILE     Write the most significant variants to this file" << endl
              << "  -V, --verbose INT               Verbosity level (0=error, 1=warn, 2=info, 3=debug, 4=trace)" << endl
              << "  -o, --output-directory DIR      Put the output in this directory" << endl;
}

int main_stoat_bh_correct(int argc, char *argv[], stoat::LogLevel &verbosity) {

    if (argc <= 1) {
        print_help_bh_correct();
        return EXIT_FAILURE;
    }

    std::string tsv_name;
    size_t p_index = std::numeric_limits<size_t>::max();
    size_t adjusted_p_index = std::numeric_limits<size_t>::max();
    std::string top_variant;
    std::string output_dir;

    int c = 0;
    optind = 1;
    while (true) {
        static struct option long_options[] =
            {
                {"tsv", required_argument, 0, 't'},
                {"p-index", required_argument, 0, 'p'},
                {"adjusted-p-index", required_argument, 0, 'a'},
                {"top-variant-file", required_argument, 0, 'v'},
                {"verbose", required_argument, 0, 'V'},
                {"output-directory", required_argument, 0, 'o'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "t:p:a:v:V:o:h",
                        long_options, &option_index); 
        if (c == -1) {
            break;
        }
        switch (c) {
            case 't':
                tsv_name = optarg;
                break;
            case 'p':
                p_index = std::stoi(optarg);
                break;
            case 'a':
                adjusted_p_index = std::stoi(optarg);
                break;
            case 'v':
                top_variant = optarg;
                break;
            case 'V': 
                {
                int level = std::stoi(optarg);
                if (level < 0 || level > 4) {
                    throw std::runtime_error("Invalid verbosity level. Use 0=Error, 1=Warn, 2=Info, 3=Debug, 4=Trace");
                }
                stoat::LogLevel logLevel = static_cast<stoat::LogLevel>(level);
                stoat::Logger::instance().setLevel(logLevel);                
                break;
                }
            case 'o':
                output_dir = optarg;
                break;
            case 'h':
                print_help_bh_correct();
                return EXIT_FAILURE;
            default:
                stoat::LOG_ERROR("[stoat BHcorrect] Unknown argument");
                print_help_bh_correct();
                return EXIT_FAILURE;
        }
    }

    if (tsv_name.empty()) {
        stoat::LOG_ERROR("[stoat BHcorrect] stoat BHcorrect requires an input tsv");
        print_help_bh_correct();
        return EXIT_FAILURE;
    }
    if (top_variant.empty()) {
        stoat::LOG_ERROR("[stoat BHcorrect] stoat BHcorrect requires a top variant file");
        print_help_bh_correct();
        return EXIT_FAILURE;
    }
    if (output_dir.empty()) {
        stoat::LOG_ERROR("[stoat BHcorrect] stoat BHcorrect requires an output directory");
        print_help_bh_correct();
        return EXIT_FAILURE;
    }
    if (p_index == std::numeric_limits<size_t>::max()) {
        stoat::LOG_ERROR("[stoat BHcorrect] stoat BHcorrect requires a p-value column");
        print_help_bh_correct();
        return EXIT_FAILURE;
    }
    if (adjusted_p_index == std::numeric_limits<size_t>::max()) {
        stoat::LOG_ERROR("[stoat BHcorrect] stoat BHcorrect requires an adjusted p-value column");
        print_help_bh_correct();
        return EXIT_FAILURE;
    }
    // Add the BH adjusted column
    // Indices are 1-indexed by the subcommand, 0-indexed by the actual function
    stoat::add_BH_adjusted_column(tsv_name, output_dir, top_variant, p_index-1, adjusted_p_index-1);

    return EXIT_SUCCESS;
}
} //end namespace
