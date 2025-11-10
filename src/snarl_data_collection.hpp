#ifndef STOAT_SNARL_DATA_COLLECTION_HPP_INCLUDED
#define STOAT_SNARL_DATA_COLLECTION_HPP_INCLUDED

#include <iostream>
#include <handlegraph/path_position_handle_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include "utils.hpp"
#include "snarl_data_t.hpp"

using namespace std;
using namespace stoat;

namespace stoat {

// This holds snarl information from the distance index and graph: id of the snarl, reference coordinates, walks, and (optionally) partitions of embedded paths
struct snarl_data_t {
    public:
        // Constructor definition
        snarl_data_t() {};

        // Constructor from a snarl, calculate everything from the graph and distance index 
        snarl_data_t(bdsg::net_handle_t snarl, const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index, bool find_partitions=false);

        // Constructor from elements
        snarl_data_t(std::vector<std::string> type_variants, std::vector<Path_traversal_t> snarl_paths, stoat::Node_traversal_t start_node,
                     stoat::Node_traversal_t end_node, std::string ref_path, size_t start_position, size_t end_position, size_t depth,
                     std::vector<std::set<sample_hap_t>> partitions) :
                     type_variants(type_variants), snarl_paths(snarl_paths), start_node(start_node), end_node(end_node), 
                     ref_path(ref_path), start_position(start_position), end_position(end_position), depth(depth), partitions(partitions)
                     {};

        // Type variant for each path: either 0 for deletion, length of the path, or min/max length of the path
        // TODO: Change name to path_lengths or something
        std::vector<std::string> type_variants;
        // All possible walks through the snarl
        std::vector<Path_traversal_t> snarl_paths;
        // Start and end nodes, both pointing into the snarl
        stoat::Node_traversal_t start_node;
        stoat::Node_traversal_t end_node;

        std::string ref_path; //TODO: I think this could get pretty big, might want to save it as an index into a list of reference paths
        // Start and end offset along the reference path
        size_t start_position;
        size_t end_position;
        size_t depth;

        std::vector<std::set<sample_hap_t>> partitions;
};
}

#endif
