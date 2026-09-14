#include "bh_correct.hpp"

#include <iostream>
#include <string>
#include <getopt.h>
#include <omp.h>
#include <filesystem>

#include "../banner.hpp"
#include "../post_processing.hpp"

using namespace std;
namespace stoat_command {

void print_help_bh_correct() {
    stoat::print_banner(std::string(STOAT_VERSION));
    std::cerr << "Usage: stoat BHcorrect [options] " << endl << endl
              << "options:" << endl
              << "  -t, --tsv FILE                  The TSV file to be processed" << endl
              << "  -p, --p-index N                 The column of the p-value in the tsv (1-indexed)" << endl
              << "                                  If a header is present then this will be the column with the appropriate label by default" << endl
              << "  -V, --verbose INT               Verbosity level (0=error, 1=warn, 2=info, 3=debug, 4=trace)" << endl;
}

int main_stoat_bh_correct(int argc, char *argv[]) {

    if (argc <= 1) {
        print_help_bh_correct();
        return EXIT_FAILURE;
    }

    std::string tsv_name;
    size_t p_index = std::numeric_limits<size_t>::max();

    int c = 0;
    optind = 1;
    while (true) {
        static struct option long_options[] =
            {
                {"tsv", required_argument, 0, 't'},
                {"p-index", required_argument, 0, 'p'},
                {"verbose", required_argument, 0, 'V'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "t:p:V:h",
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

    // print banner
    stoat::print_banner(std::string(STOAT_VERSION));

    std::string temp_file_name = tsv_name + ".tmp";

    std::shared_ptr<stoat::Reader> reader;
    std::shared_ptr<stoat::Reader> reader_copy;
    std::shared_ptr<stoat::Writer> writer;

    if ((tsv_name.compare(tsv_name.length()-3, 3, ".gz") == 0) ||
        (tsv_name.compare(tsv_name.length()-4, 4, ".bgz") == 0)) {
        reader.reset(new stoat::BgzReader(tsv_name));
        reader_copy.reset(new stoat::BgzReader(tsv_name));
        writer.reset(new stoat::BgzWriter(temp_file_name));
    } else {
        reader.reset(new stoat::StdReader(tsv_name));
        reader_copy.reset(new stoat::StdReader(tsv_name));
        writer.reset(new stoat::StdWriter(temp_file_name));
    }


    // Add the BH adjusted column
    // Indices are 1-indexed by the subcommand, 0-indexed by the actual function
    stoat::add_BH_adjusted_column(reader, reader_copy, writer, p_index == std::numeric_limits<size_t>::max() ? p_index : p_index-1);

    reader->close();
    reader_copy->close();
    writer->close();

    // Move the temporary file
    std::remove(tsv_name.c_str());
    std::rename(temp_file_name.c_str(), tsv_name.c_str());


    return EXIT_SUCCESS;
}
} //end namespace
