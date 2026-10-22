#include "get_tree.hpp"

#include <iostream>
#include <string>
#include <getopt.h>
#include <omp.h>
#include <filesystem>

#include <bdsg/overlays/overlay_helper.hpp>
#include <vg/io/vpkg.hpp>

#include "../banner.hpp"
#include "../io/register_io.hpp"
#include "../post_processing.hpp"

using namespace std;
namespace stoat_command {

void print_help_get_tree() {
    stoat::print_banner(std::string(STOAT_VERSION));
    std::cerr << "Usage: stoat get-tree -g [graph] -d [distance-index] -p [reference range] > [tsv]" << endl << endl
              << "options:" << endl
              << "  -g, --graph FILE                The graph" << endl
              << "  -d, --distance-index FILE       The distance index" << endl
              << "  -p, --reference-range FILE      The reference range, written as path_name:start-end" << endl;
}

int main_stoat_get_tree(int argc, char *argv[]) {

    if (argc <= 1) {
        print_help_get_tree();
        return EXIT_FAILURE;
    }

    std::string graph_name;
    std::string dist_name;
    std::string reference_range;

    omp_set_num_threads(1);

    int c = 0;
    optind = 1;
    while (true) {
        static struct option long_options[] =
            {
                {"graph", required_argument, 0, 'g'},
                {"dist_name", required_argument, 0, 'd'},
                {"reference-range", required_argument, 0, 'p'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "g:d:p:h",
                        long_options, &option_index); 
        if (c == -1) {
            break;
        }
        switch (c) {
            case 'g':
                graph_name = optarg;
                break;
            case 'd':
                dist_name = optarg;
                break;
            case 'p':
                reference_range = optarg;
                break;

            case 'h':
                print_help_get_tree();
                return EXIT_FAILURE;
            default:
                stoat::LOG_ERROR("[stoat get-tree] Unknown argument");
                print_help_get_tree();
                return EXIT_FAILURE;
        }
    }

    if (graph_name.empty()) {
        stoat::LOG_ERROR("[stoat get-tree] stoat get-tree requires a graph");
        print_help_get_tree();
        return EXIT_FAILURE;
    }
    if (dist_name.empty()) {
        stoat::LOG_ERROR("[stoat get-tree] stoat get-tree requires a distance index");
        print_help_get_tree();
        return EXIT_FAILURE;
    }
    if (reference_range.empty()) {
        stoat::LOG_ERROR("[stoat get-tree] stoat get-tree requires a reference range");
        print_help_get_tree();
        return EXIT_FAILURE;
    }

    // Tell the IO library about libvg types.
    if (!stoat::io::register_libvg_io()) {
        stoat::LOG_ERROR("[stoat vgio] Could not register libvg types with libvgio");
        return EXIT_FAILURE;
    }


    // Load the graph and make it a PathPositionHandleGraph

    std::unique_ptr<handlegraph::PathHandleGraph> graph = vg::io::VPKG::load_one<handlegraph::PathHandleGraph>(graph_name);

    // print banner
    stoat::print_banner(std::string(STOAT_VERSION));

    // Load the distance index
    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize(dist_name);

    // Get the reference path

    size_t colon_i = reference_range.rfind(":"); 
    std::string ref_path = reference_range.substr(0, colon_i).c_str();
    size_t dash_i = reference_range.find("-", colon_i); 
    size_t start = std::atoi(reference_range.substr(colon_i+1, dash_i-colon_i).c_str());
    size_t end = std::atoi(reference_range.substr(dash_i+1, std::string::npos).c_str());


    std::unordered_set<std::string> reference_names;
    reference_names.emplace(ref_path);

    bdsg::ReferencePathOverlayHelper overlay_helper;
    bdsg::PathPositionHandleGraph* path_position_graph =  overlay_helper.apply(graph.get(), reference_names);

    // Write the tsv to stdout
    print_tree_tsv(*path_position_graph, distance_index, ref_path, start, end);

    return EXIT_SUCCESS;
}


// The function to print the snarl, parent tsv. I'm being lazy and writing it here because it won't get used anywhere else
void print_tree_tsv(const bdsg::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index, std::string& ref_path, size_t start, size_t end) {
    // Start by getting the bounding nodes/snarls in the lowest common ancestor chain for the reference path
    handlegraph::path_handle_t path_handle = graph.get_path_handle(ref_path);
    if ( start > end || start >= graph.get_path_length(path_handle) || end >= graph.get_path_length(path_handle)) {
        std::cerr << "Range " << start << "-" << end << " is invalid for path " << ref_path << std::endl;
        throw std::runtime_error("error: invalid coordinate range");
    }
    handlegraph::net_handle_t start_node = distance_index.get_net(graph.get_handle_of_step(graph.get_step_at_position(path_handle, start)), &graph);
    handlegraph::net_handle_t end_node = distance_index.get_net(graph.get_handle_of_step(graph.get_step_at_position(path_handle, end)), &graph);

    std::unordered_set<handlegraph::net_handle_t> ancestors_of_start;
    handlegraph::net_handle_t ancestor = distance_index.start_end_traversal_of(distance_index.get_parent(start_node));
    while (!distance_index.is_root(ancestor)) {
        ancestors_of_start.emplace(ancestor);
        handlegraph::net_handle_t next = distance_index.start_end_traversal_of(distance_index.get_parent(ancestor));
        ancestor = next;
    }
    handlegraph::net_handle_t end_ancestor = distance_index.start_end_traversal_of(end_node);
    ancestor = distance_index.start_end_traversal_of(distance_index.get_parent(end_node));
    while (!distance_index.is_root(ancestor)) {
        if (ancestors_of_start.count(ancestor)) {
            break;
        } else {
            end_ancestor = ancestor;
            ancestor = distance_index.start_end_traversal_of(distance_index.get_parent(end_ancestor));
        }
    }

    // ancestor is now the lowest common ancestor
    // end_ancestor is now the ancestor of start_node that is a child of the lowest common ancestor
    // Find start_ancestor
    handlegraph::net_handle_t start_ancestor = distance_index.start_end_traversal_of(distance_index.get_parent(start_node));
    while (start_ancestor != ancestor) {
        start_node = start_ancestor;
        start_ancestor = distance_index.start_end_traversal_of(distance_index.get_parent(start_node));
    }
    start_ancestor = distance_index.start_end_traversal_of(start_node);
    if (distance_index.is_root(ancestor)) {
        throw std::runtime_error("error: the offsets are on different connected components or a weird top-level snarl");
        //TODO: Could try to deal with top-level snarls
    }

    // Now we want to walk across the lowest common ancestor chain and get all snarls
    // If the lowest common ancestor is a snarl, then just get that snarl  
    if (distance_index.is_snarl(ancestor)) {
        start_ancestor = ancestor;
        end_ancestor = ancestor;
    }

    // Get the start and end in the right orientation
    // Both start and end ancestors are currently in their local start-to-end orientation, which may not be pointing forward in the chain
    bool traverse_forward = distance_index.is_ordered_in_chain(start_ancestor, end_ancestor);
    if (traverse_forward == distance_index.is_reversed_in_parent(start_ancestor)) {
        handlegraph::net_handle_t flipped = distance_index.flip(start_ancestor);
        start_ancestor = flipped;
    }
    if (traverse_forward == distance_index.is_reversed_in_parent(end_ancestor)) {
        handlegraph::net_handle_t flipped = distance_index.flip(end_ancestor);
        end_ancestor = flipped;
    }
    handlegraph::net_handle_t stop_ancestor = distance_index.get_root();
    distance_index.follow_net_edges(end_ancestor, &graph, false, [&](net_handle_t next_child) {
        stop_ancestor = next_child;
        return true;
    });

    // Now walk through the chain from start to end ancestors and emit all snarls contained in the chain
     while (start_ancestor != stop_ancestor) {

        if (distance_index.is_snarl(start_ancestor)) {

            // Put all snarls we need to print in to_explore
            // We take this snarl as the parent and print all the child snarls, then put the children in to_explore
            std::vector<std::pair<handlegraph::net_handle_t, std::string>> to_explore;

            // Get the bounds of the snarl facing in for printing later
            handlegraph::handle_t current_start_in = distance_index.get_handle(distance_index.get_bound(start_ancestor, false, true), &graph);
            handlegraph::handle_t current_end_in = distance_index.get_handle(distance_index.get_bound(start_ancestor, true, true), &graph);
            
            stoat::node_traversal_t current_start_trav(graph.get_id(current_start_in),
                                                            graph.get_is_reverse(current_start_in));
            stoat::node_traversal_t current_end_trav(graph.get_id(current_end_in),
                                                          graph.get_is_reverse(current_end_in));
            to_explore.emplace_back(start_ancestor, current_start_trav.to_string() + current_end_trav.to_string());

            // Write this snarl with no ancestor
            std::cout << current_start_trav.to_string() << current_end_trav.to_string() << "\t." << std::endl;

            // TODO: Check the order/orientation of snarl children to put them in a roughly linear order

            while (!to_explore.empty()) {
                handlegraph::net_handle_t current_snarl = to_explore.back().first;
                std::string current_snarl_str = to_explore.back().second;

                to_explore.pop_back();
                
                distance_index.for_each_child(current_snarl, [&] (net_handle_t child_chain) {
                    distance_index.for_each_child(child_chain, [&] (net_handle_t child_snarl) {
                        if (distance_index.is_snarl(child_snarl)) {
                            // For each grandchild snarl of the current snarl

                            // Get the bounds as strings
                            handlegraph::handle_t child_start_in = distance_index.get_handle(distance_index.get_bound(child_snarl, false, true), &graph);
                            handlegraph::handle_t child_end_in = distance_index.get_handle(distance_index.get_bound(child_snarl, true, true), &graph);
                            
                            stoat::node_traversal_t child_start_trav(graph.get_id(child_start_in),
                                                                            graph.get_is_reverse(child_start_in));
                            stoat::node_traversal_t child_end_trav(graph.get_id(child_end_in),
                                                                          graph.get_is_reverse(child_end_in));
                            std::string child_str = child_start_trav.to_string() + child_end_trav.to_string();

                            // Write this grandchild snarl
                            std::cout << child_str << "\t" << current_snarl_str << std::endl;

                            // Remember to print its children
                            to_explore.emplace_back(child_snarl, child_str);
                        }
                        return true; });
                    return true;
                });
            }
        }


        handlegraph::net_handle_t next = distance_index.get_root();
        distance_index.follow_net_edges(start_ancestor, &graph, false, [&](net_handle_t next_child) {
            next = next_child;
            return true;
        });
        start_ancestor = next;
    }

}
} //end namespace
