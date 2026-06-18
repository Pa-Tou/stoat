#include "change_reference.hpp"

#include <iostream>
#include <string>
#include <getopt.h>
#include <omp.h>
#include <filesystem>

#include <bdsg/snarl_distance_index.hpp>
#include <handlegraph/path_handle_graph.hpp>
#include <bdsg/overlays/overlay_helper.hpp>
#include <vg/io/vpkg.hpp>

#include "../banner.hpp"
#include "../io/register_io.hpp"
#include "../post_processing.hpp"

using namespace std;
namespace stoat_command {

void print_help_change_reference() {
    stoat::print_banner(std::string(STOAT_VERSION));
    std::cerr << "Usage: stoat change-ref -T [stoat.assoc.pvalues.tsv] -g [graph] -d [distance-index] -r [reference name] > [renamed_tsv]" << endl << endl
              << "options:" << endl
              << "  -T, --tsv FILE                  The TSV file to be processed, the output file of stoat" << endl
              << "  -o, --out-tsv FILE              The output file to be written" << endl
              << "  -g, --graph FILE                The graph used to find coordinates" << endl
              << "  -d, --distance-index FILE       The distance index" << endl
              << "  -r, --reference-names FILE      Rewrite the reference coordinates to be relative to the paths listed in FILE (one per line)" << endl
              << "                                  Paths in the graph can be found with vg paths" << endl
              << "  -t, --threads N                 Number of threads to use" << endl;
}

int main_stoat_change_reference(int argc, char *argv[]) {

    if (argc <= 1) {
        print_help_change_reference();
        return EXIT_FAILURE;
    }

    std::string in_tsv_name;
    std::string out_tsv_name;
    std::string graph_name;
    std::string dist_name;
    std::string reference_file;

    omp_set_num_threads(1);

    int c = 0;
    optind = 1;
    while (true) {
        static struct option long_options[] =
            {
                {"tsv", required_argument, 0, 'T'},
                {"out-tsv", required_argument, 0, 'o'},
                {"graph", required_argument, 0, 'g'},
                {"dist_name", required_argument, 0, 'd'},
                {"reference-prefix", required_argument, 0, 'r'},
                {"threads", required_argument, 0, 't'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "T:o:g:d:r:t:h",
                        long_options, &option_index); 
        if (c == -1) {
            break;
        }
        switch (c) {
            case 'T':
                in_tsv_name = optarg;
                break;
            case 'o':
                out_tsv_name = optarg;
                break;
            case 'g':
                graph_name = optarg;
                break;
            case 'd':
                dist_name = optarg;
                break;
            case 'r':
                reference_file = optarg;
                break;
            case 't':
                if (std::stoi(optarg) < 1) {
                    stoat::LOG_ERROR("[stoat graph] Number of threads must be > 0");
                    return EXIT_FAILURE;
                }
                omp_set_num_threads(std::stoi(optarg));
                break;

            case 'h':
                print_help_change_reference();
                return EXIT_FAILURE;
            default:
                stoat::LOG_ERROR("[stoat change-ref] Unknown argument");
                print_help_change_reference();
                return EXIT_FAILURE;
        }
    }

    if (in_tsv_name.empty()) {
        stoat::LOG_ERROR("[stoat change-ref] stoat change-ref requires an input tsv");
        print_help_change_reference();
        return EXIT_FAILURE;
    }

    if (out_tsv_name.empty()) {
        stoat::LOG_ERROR("[stoat change-ref] stoat change-ref requires an output tsv");
        print_help_change_reference();
        return EXIT_FAILURE;
    }
    if (graph_name.empty()) {
        stoat::LOG_ERROR("[stoat change-ref] stoat change-ref requires a graph");
        print_help_change_reference();
        return EXIT_FAILURE;
    }
    if (dist_name.empty()) {
        stoat::LOG_ERROR("[stoat change-ref] stoat change-ref requires a distance index");
        print_help_change_reference();
        return EXIT_FAILURE;
    }
    if (reference_file.empty()) {
        stoat::LOG_ERROR("[stoat change-ref] stoat change-ref requires a reference prefix");
        print_help_change_reference();
        return EXIT_FAILURE;
    }

    // Tell the IO library about libvg types.
    if (!stoat::io::register_libvg_io()) {
        stoat::LOG_ERROR("[stoat vgio] Could not register libvg types with libvgio");
        return EXIT_FAILURE;
    }


    // Load the graph and make it a PathPositionHandleGraph
    std::unique_ptr<handlegraph::PathHandleGraph> graph = vg::io::VPKG::load_one<handlegraph::PathHandleGraph>(graph_name);
    bdsg::PathPositionOverlayHelper overlay_helper;
    bdsg::PathPositionHandleGraph* path_position_graph =  overlay_helper.apply(graph.get());

    // print banner
    stoat::print_banner(std::string(STOAT_VERSION));

    // Load the distance index
    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize(dist_name);

    // Get the reference paths
    std::unordered_set<std::string> reference_names;
    ifstream ref_stream;
    ref_stream.open(reference_file);
    std::string ref_name;
    while (std::getline(ref_stream, ref_name)) {
        reference_names.emplace(ref_name);
    }
    ref_stream.close();

    std::shared_ptr<stoat::Reader> reader;
    std::shared_ptr<stoat::Writer> writer;

    if ((in_tsv_name.compare(in_tsv_name.length()-3, 3, ".gz") == 0) ||
        (in_tsv_name.compare(in_tsv_name.length()-4, 4, ".bgz") == 0)) {
        reader.reset(new stoat::BgzReader(in_tsv_name));
        writer.reset(new stoat::BgzWriter(out_tsv_name));
    } else {
        reader.reset(new stoat::StdReader(in_tsv_name));
        writer.reset(new stoat::StdWriter(out_tsv_name));
    }

    // Write the new file to stdout
    stoat::change_reference(*path_position_graph, distance_index, reader, writer, reference_names);

    reader->close();
    writer->close();

    return EXIT_SUCCESS;
}
} //end namespace
