#ifndef snarl_data_t
#define snarl_data_t

#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <iostream>
#include <array>
#include <chrono>
#include <cassert>
#include <regex>
#include <cstddef>
#include <stdexcept>
#include <utility>

#include <bdsg/hash_graph.hpp>
#include <bdsg/packed_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include <bdsg/overlays/packed_path_position_overlay.hpp>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>
#include <bdsg/overlays/overlay_helper.hpp>
#include <vg/io/vpkg.hpp>

#include "log.hpp"
#include "utils.hpp"

#include <filesystem>
#include "io/register_io.hpp"

using namespace std;

using handlegraph::step_handle_t;
using handlegraph::handle_t;
using handlegraph::net_handle_t;

namespace stoat {

// Define a Node_traversal_t structure to represent a node with orientation
struct Node_traversal_t { // 64 bits per node 
    private:
        size_t node_id : 63; // 63 bits for node ID
        bool is_reverse : 1; // 1 bit for orientation (true for reverse, false for forward)

    public:
        Node_traversal_t(const size_t &id, const bool &rev);
        
        // Getters
        size_t get_node_id() const;
        bool get_is_reverse() const;

        // Convert to std::string representation
        std::string to_string() const;

        bool operator==(const Node_traversal_t& other) const;
};

// Define a Edge_t structure to represent an edge between two Node_traversal_t nodes
struct Edge_t { // 128 bits per edge 
    private:
        std::pair<Node_traversal_t, Node_traversal_t> edge;

    public:
        Edge_t(const Node_traversal_t &node_traversal_1, const Node_traversal_t &node_traversal_2);
        
        // Converter
        std::pair<size_t, size_t> print_pair_edge() const;
        std::string print_string_edge() const;

        // Accessor to edge, useful for hashing and comparison
        const std::pair<Node_traversal_t, Node_traversal_t>& get_edge() const;

        // Comparison operator
        bool operator==(const Edge_t &other) const;
};

// Define a Path_traversal_t structure to represent a path through the graph
struct Path_traversal_t {
    private:
        std::vector<Node_traversal_t> paths; // Nodes in the path

    public:
        // add a node traversal to the path
        Path_traversal_t() = default;
        void add_node_traversal_t(const Node_traversal_t &paths);

        // Getters
        const std::vector<Node_traversal_t>& get_paths() const;
        
        // convert to std::string representation
        std::string to_string() const;
};

// Define a Snarl_data_t structure to hold snarl information
struct Snarl_data_t {
    public:
        // Constructor definition
        Snarl_data_t(bdsg::net_handle_t snarl_, const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index);
        Snarl_data_t(net_handle_t snarl_,
                    std::pair<size_t, size_t> snarl_ids_,
                    std::vector<Path_traversal_t> snarl_paths_,
                    const size_t start_positions_, const size_t end_positions_,
                    std::vector<std::string> type_variants_,
                    size_t depth);  // Assuming path_nodes correspond to type_variants

        std::vector<std::string> type_variants;
        std::vector<Path_traversal_t> snarl_paths;
        net_handle_t snarl; // handlegraph::subrange_t Snarl_data_t::snarl_id
        std::pair<size_t, size_t> snarl_ids;
        size_t start_positions;
        size_t end_positions;
        size_t depth;
};

// Converter
std::string pairToString(const std::pair<size_t, size_t>& name);
std::pair<size_t, size_t> stringToPair(const std::string& str);
std::string vectorPathToString(const std::vector<Path_traversal_t>& vec_paths);
std::vector<Path_traversal_t> stringToVectorPath(std::string& str);

// A class representing a path as a vector of strings representing nodes
class Path {
private:
    std::vector<size_t> nodes;
    std::vector<bool> orients;

public:
    // Constructor
    Path();

    // Add a node with known orientation
    void addNode(const size_t& node, bool orient);

    // Add a node handle and extract information using the std::string representation
    bool addNodeHandle(const handlegraph::net_handle_t& node_h, const bdsg::SnarlDistanceIndex& stree);

    // Get the std::string representation of the path
    Path_traversal_t print() const;

    // Flip the path orientation
    void flip();

    // Get the size of the path
    size_t size() const;

    // Count the number of reversed nodes
    size_t nreversed() const;
};

// Parses the snarl path file and returns a map with snarl as keys and paths as a list of strings.
std::unordered_map<std::string, std::vector<Snarl_data_t>> parse_snarl_path(const std::string& path_file);

void write_snarl_data_output(std::ostream& outstream);
void write_snarl_data_fail(std::ostream& outstream);

// Load the distance index and graph and return unique_ptrs to them
std::tuple<std::unique_ptr<bdsg::SnarlDistanceIndex>, 
           std::unique_ptr<bdsg::PackedGraph>, 
           handlegraph::net_handle_t, 
           std::unique_ptr<handlegraph::PathHandleGraph>,
           std::unique_ptr<bdsg::PackedPositionOverlay>>
parse_graph_tree(const std::string& pg_file, const std::string& dist_file);

// Function to calculate the type of variant
// Given a vector of <size node 2, min length of the snarl, max length of the snarl, path length, sum_path, is_complex)
// TODO : change sum_path to definition using the length of the path including in the boundary nodes
// Matis ans : i don t know how to do it
std::vector<std::string> calcul_pos_type_variant(const std::vector<std::tuple<size_t, size_t, size_t>>& list_length_paths);

// Function to follow edges
void follow_edges(bdsg::SnarlDistanceIndex& stree,
    std::vector<std::vector<handlegraph::net_handle_t>>& finished_paths,
    const std::vector<handlegraph::net_handle_t>& path,
    std::vector<std::vector<handlegraph::net_handle_t>>& paths,
    bdsg::PackedGraph& pg,
    const bool& cycle);

// Function to save snarls
std::vector<std::tuple<handlegraph::net_handle_t, std::string, size_t, size_t, bool>> save_snarls(
                            bdsg::SnarlDistanceIndex& stree, 
                            handlegraph::net_handle_t& root,
                            bdsg::PackedGraph& pg, 
                            unordered_set<std::string>& ref_paths,
                            bdsg::PackedPositionOverlay& ppo);

// Function to fill pretty paths
tuple<std::vector<Path_traversal_t>, std::vector<std::string>> fill_pretty_paths(
                            bdsg::SnarlDistanceIndex& stree, 
                            bdsg::PackedGraph& pg, 
                            std::vector<std::vector<handlegraph::net_handle_t>>& finished_paths);

// Function to loop over snarls and write output to output_file
// Output is a tsv of <chromosome, start pos, end pos, snarl, paths, variant type, reference>
// Returns a map from chromosome name to a vector of <snarl name, paths, start position, end position, variant type>
std::unordered_map<std::string, std::vector<Snarl_data_t>> loop_over_snarls_write(
                            bdsg::SnarlDistanceIndex& stree, 
                            std::vector<std::tuple<handlegraph::net_handle_t, std::string, size_t, size_t, bool>>& snarls, 
                            bdsg::PackedGraph& pg, 
                            const std::string& output_file, 
                            const std::string& output_snarl_not_analyse, 
                            const size_t& children_treshold,
                            const size_t& path_length_threshold,
                            const size_t& cycle_threshold);

} // end namespace stoat

// Hash functions for Node_traversal_t
namespace std {
    template <>
    struct hash<stoat::Node_traversal_t> {
        size_t operator()(const stoat::Node_traversal_t& node) const {
            // Simple way: Shift node_id and pack is_reverse into the lower bit
            return (node.get_node_id() << 1) | static_cast<size_t>(node.get_is_reverse());
        }
    };

    // Hash function for Edge_t
    template <>
    struct hash<stoat::Edge_t> {
        size_t operator()(const stoat::Edge_t& edge) const {
            const auto& pair = edge.get_edge();
            size_t h1 = hash<stoat::Node_traversal_t>()(pair.first);
            size_t h2 = hash<stoat::Node_traversal_t>()(pair.second);

            // Combines two hash values (h1 and h2) into a single hash using bitwise operations.
            // 0x9e3779b9 is a large prime constant (from the golden ratio) used to improve distribution.
            // (h1 << 6) and (h1 >> 2) add additional mixing by shifting bits left and right.
            // This reduces hash collisions by ensuring small changes in input produce different hashes.
            return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
        }
    };

} // end namespace std

#endif
