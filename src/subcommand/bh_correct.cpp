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
         << "  -o, --output-directory DIR      Put the output in this directory" << endl;
}

int main_stoat_bh_correct(int argc, char *argv[]) {

    if (argc <= 1) {
        print_help_bh_correct();
        return 1;
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
                {"output-directory", required_argument, 0, 'o'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "t:p:a:v:o:h",
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
            case 'o':
                output_dir = optarg;
                break;
            case 'h':
                print_help_bh_correct();
                exit(1);
                break;
            default:
                abort();
        }
    }

    if (tsv_name.empty()) {
        cerr << "error [stoat BHcorrect]: stoat BHcorrect requires an input tsv" << endl;
        return 1;
    }
    if (top_variant.empty()) {
        cerr << "error [stoat BHcorrect]: stoat BHcorrect requires a top variant file" << endl;
        return 1;
    }
    if (output_dir.empty()) {
        cerr << "error [stoat BHcorrect]: stoat BHcorrect requires an output directory" << endl;
        return 1;
    }
    if (p_index == std::numeric_limits<size_t>::max()) {
        cerr << "error [stoat BHcorrect]: stoat BHcorrect requires a p-value column" << endl;
        return 1;
    }
    if (adjusted_p_index == std::numeric_limits<size_t>::max()) {
        cerr << "error [stoat BHcorrect]: stoat BHcorrect requires an adjusted p-value column" << endl;
        return 1;
    }
    // Add the BH adjusted column
    // Indices are 1-indexed by the subcommand, 0-indexed by the actual function
    stoat::add_BH_adjusted_column(tsv_name, output_dir, top_variant, p_index-1, adjusted_p_index-1);

    return 0;
}
} //end namespace
