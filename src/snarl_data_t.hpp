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
#include <bdsg/snarl_distance_index.hpp>
#include <bdsg/overlays/path_position_overlays.hpp>
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
        Node_traversal_t(const std::string &str);
        
        // Setter
        void set_is_reverse(const bool &rev) { is_reverse = rev; };
        void set_node_id(const size_t &id) { node_id = id; };

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
    // JEAN why is that a pair? can't we just have two nodes in that struct?
        std::pair<Node_traversal_t, Node_traversal_t> edge;

    public:
        Edge_t(const Node_traversal_t &node_traversal_1, const Node_traversal_t &node_traversal_2);
        
        // Converter
        std::pair<size_t, size_t> get_node_pair() const;
        std::string to_string() const;

        // Accessor to edge, useful for hashing and comparison
        const std::pair<Node_traversal_t, Node_traversal_t>& get_edge() const;

    // return a flipped version of this edge
    Edge_t get_flipped() const;

        // Comparison operator
        bool operator==(const Edge_t &other) const;
};

// Define a PathTraversal structure to represent a path through the graph
class PathTraversal {
private:
    std::vector<Node_traversal_t> path; // Nodes in the path
    size_t min_allele_len;
    size_t max_allele_len;
    
public:
    // add a node traversal to the path
    PathTraversal() : min_allele_len(0), max_allele_len(0) {}
    void add_node(const size_t& node, bool is_rev);
    void add_node_handle(const handlegraph::net_handle_t& node_h, const bdsg::SnarlDistanceIndex& distance_index);
    void add_node_traversal_t(const Node_traversal_t &node_trav);
    
    // Check and flip the Path if necessary to ensure consistent orientation
    void check_path_flip();
    void path_flip();

    void add_min_allele_len(size_t len);
    void add_max_allele_len(size_t len);
    void set_allele_length_from_string(std::string al_len_str);

    // TODO : change sum_path to definition using the length of the path including in the boundary nodes
    // Matis ans : i don t know how to do it
    std::string get_allele_length() const;
        
    // Getters
    const std::vector<Node_traversal_t>& get_path() const;
    size_t get_max_allele_length() const { return max_allele_len; }
    size_t get_min_allele_length() const { return min_allele_len; }
    
    // convert to std::string representation
    std::string to_string() const;
    size_t size() const;
};

// TODO: Make std::pair<size_t, size_t> a struct?
//struct snarl_id {
//    size_t start_id;
//    size_t end_id;
//}
    
//TODO: This is eventually going to be replaced with snarl_info_t from snarl_data_collection.hpp
// Define a Snarl_data_t structure to hold snarl information
struct Snarl_data_t {
    public:
        // Constructor definition
        Snarl_data_t() {};
        Snarl_data_t(bdsg::net_handle_t net_, const bdsg::SnarlDistanceIndex& distance_index);
        Snarl_data_t(net_handle_t net_,
                    std::pair<size_t, size_t> ids_,
                    std::vector<PathTraversal> paths_,
                    const size_t start_position_, const size_t end_position_,
                    size_t depth);
    // constructor when the input information is in string form (when reading a file for example)
        Snarl_data_t(net_handle_t net_,
                     std::string snarl_ids_,
                     std::string paths_,
                     std::string allele_lengths_,
                     const size_t start_position_, const size_t end_position_,
                    size_t depth);

        std::vector<PathTraversal> paths;
        net_handle_t net; // handlegraph::subrange_t Snarl_data_t::snarl_id
        //TODO: Make this a node_traversal_t
        std::pair<size_t, size_t> ids;
        size_t start_position;
        size_t end_position;
        size_t depth;
        bool is_on_ref;
};

// Convert a pair of size_t, for example defining a snarl ID to a string of them separated by an underscore
std::string pairToString(const std::pair<size_t, size_t>& name);

// convert a vector of path traversals to a string, either with node or allele length information
std::string vectorPathToString(const std::vector<PathTraversal>& vec_paths, bool allele_lengths = false);

// Get a vector of path traversals from their string representation
std::vector<stoat::PathTraversal> string_to_path_traversals(const std::string& path_string, const std::string& path_lengths_string);

// Parses the snarl path file and returns a map with snarl as keys and paths as a list of strings.
std::unordered_map<std::string, std::vector<Snarl_data_t>> read_snarl_path(const std::string& path_file);

// Load the distance index and graph and return unique_ptrs to them
std::tuple<
    unique_ptr<bdsg::SnarlDistanceIndex>,
    unique_ptr<handlegraph::PathHandleGraph>>
load_graph_tree(const std::string& graph_file, const std::string& dist_file);

// Function to list all the snarls and their position on the reference paths (if possible)
// group the snarls by chromosome/reference path in the output map (chr -> snarls)
std::unordered_map<std::string, std::vector<Snarl_data_t>> list_all_snarls_with_pos(
                            const bdsg::SnarlDistanceIndex& distance_index, 
                            handlegraph::PathHandleGraph& graph, 
                            unordered_set<std::string>& ref_paths);

// convert paths from the simple vector of net handles to the PathTraversal object
std::vector<PathTraversal> convert_path_traversals(
                            const bdsg::SnarlDistanceIndex& distance_index, 
                            const handlegraph::PathHandleGraph& graph, 
                            std::vector<std::vector<handlegraph::net_handle_t>>& finished_paths);

// Function to loop over snarls and write output to output_file
// Output is a tsv of <chromosome, start pos, end pos, snarl, paths, variant type, reference>
// Returns a map from chromosome name to a vector of <snarl name, paths, start position, end position, variant type>
std::unordered_map<std::string, std::vector<Snarl_data_t>> write_snarls_with_paths(
                            const bdsg::SnarlDistanceIndex& distance_index, 
                            std::unordered_map<std::string, std::vector<Snarl_data_t>>& chr_to_snarls, 
                            handlegraph::PathHandleGraph& graph, 
                            const std::string& output_dir,
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
