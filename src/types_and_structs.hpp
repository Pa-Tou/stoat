#ifndef types_and_structs
#define types_and_structs

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

#include <filesystem>
#include "io/register_io.hpp"

using namespace std;

using handlegraph::step_handle_t;
using handlegraph::handle_t;
using handlegraph::net_handle_t;

namespace stoat {

// Define a node_traversal_t structure to represent a node with orientation
struct node_traversal_t { // 64 bits per node 
    private:
        size_t node_id : 63; // 63 bits for node ID
        bool is_reverse : 1; // 1 bit for orientation (true for reverse, false for forward)

    public:
        node_traversal_t(const size_t &id, const bool &rev);
        node_traversal_t(const std::string &str);
        
        // Setter
        void set_is_reverse(const bool &rev) { is_reverse = rev; };
        void set_node_id(const size_t &id) { node_id = id; };

        // Getters
        size_t get_node_id() const;
        bool get_is_reverse() const;

        // Convert to std::string representation
        std::string to_string() const;

        bool operator==(const node_traversal_t& other) const;
};

// Define a edge_t structure to represent an edge between two node_traversal_t nodes
struct edge_t { // 128 bits per edge 
    private:
    // JEAN why is that a pair? can't we just have two nodes in that struct?
        std::pair<node_traversal_t, node_traversal_t> edge;

    public:
        edge_t(const node_traversal_t &node_traversal_1, const node_traversal_t &node_traversal_2);
        
        // Converter
        std::pair<size_t, size_t> get_node_pair() const;
        std::string to_string() const;

        // Accessor to edge, useful for hashing and comparison
        const std::pair<node_traversal_t, node_traversal_t>& get_edge() const;

    // return a flipped version of this edge
    edge_t get_flipped() const;

        // Comparison operator
        bool operator==(const edge_t &other) const;
};

// Define a PathTraversal structure to represent a path through the graph
class PathTraversal {
private:
    std::vector<node_traversal_t> path; // Nodes in the path
    size_t min_allele_len;
    size_t max_allele_len;
    
public:
    // add a node traversal to the path
    PathTraversal() : min_allele_len(0), max_allele_len(0) {}
    void add_node(const size_t& node, bool is_rev);
    void add_node_handle(const handlegraph::net_handle_t& node_h, const bdsg::SnarlDistanceIndex& distance_index);
    void add_node_traversal_t(const node_traversal_t &node_trav);
    
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
    const std::vector<node_traversal_t>& get_path() const;
    size_t get_max_allele_length() const { return max_allele_len; }
    size_t get_min_allele_length() const { return min_allele_len; }
    
    // convert to std::string representation
    std::string to_string() const;
    size_t size() const;
};

// Convert a pair of size_t, for example defining a snarl ID to a string of them separated by an underscore
std::string pairToString(const std::pair<size_t, size_t>& name);

// convert a vector of path traversals to a string, either with node or allele length information
std::string vectorPathToString(const std::vector<PathTraversal>& vec_paths, bool allele_lengths = false);

// Get a vector of path traversals from their string representation
std::vector<stoat::PathTraversal> string_to_path_traversals(const std::string& path_string, const std::string& path_lengths_string);

// Load the distance index and graph and return unique_ptrs to them
std::tuple<
    unique_ptr<bdsg::SnarlDistanceIndex>,
    unique_ptr<handlegraph::PathHandleGraph>>
load_graph_tree(const std::string& graph_file, const std::string& dist_file);

// convert paths from the simple vector of net handles to the PathTraversal object
std::vector<PathTraversal> convert_path_traversals(
                            const bdsg::SnarlDistanceIndex& distance_index, 
                            const handlegraph::PathHandleGraph& graph, 
                            std::vector<std::vector<handlegraph::net_handle_t>>& finished_paths);

// This stores a sample name and haplotype identifier
struct sample_hap_t {
    std::string sample;
    std::string haplotype;

    sample_hap_t() {};
    sample_hap_t(const handlegraph::PathHandleGraph& graph, const handlegraph::path_handle_t& path);

    sample_hap_t(std::string path_name){
        // Get the sample name up to #
        std::stringstream stream(path_name);
        if (std::getline(stream, sample, '#')) {
            std::getline(stream, haplotype);
        } else {
            haplotype = "";
        }

    };
    std::string to_string() const {
        return sample + "#" + haplotype;
    }
    sample_hap_t(std::string samp, std::string hap) :
        sample(std::move(samp)), haplotype(std::move(hap)) {};

    const inline bool operator==(const sample_hap_t& other) const {
        return (sample==other.sample && haplotype==other.haplotype);
    }
    const inline bool operator<(const sample_hap_t& other) const {
        if (sample == other.sample) {
            return haplotype < other.haplotype;
        } else {
            return sample < other.sample;
        }
    }
};

inline std::ostream& operator<<(std::ostream& out, const sample_hap_t& sample) {
    return out << sample.sample << "#" << sample.haplotype;
}

// A struct for holding a range along the path
struct path_range_t {
    handlegraph::step_handle_t start;
    handlegraph::step_handle_t end;
};


} // end namespace stoat

// Hash functions for node_traversal_t
namespace std {
    template <>
    struct hash<stoat::node_traversal_t> {
        size_t operator()(const stoat::node_traversal_t& node) const {
            // Simple way: Shift node_id and pack is_reverse into the lower bit
            return (node.get_node_id() << 1) | static_cast<size_t>(node.get_is_reverse());
        }
    };

    // Hash function for edge_t
    template <>
    struct hash<stoat::edge_t> {
        size_t operator()(const stoat::edge_t& edge) const {
            const auto& pair = edge.get_edge();
            size_t h1 = hash<stoat::node_traversal_t>()(pair.first);
            size_t h2 = hash<stoat::node_traversal_t>()(pair.second);

            // Combines two hash values (h1 and h2) into a single hash using bitwise operations.
            // 0x9e3779b9 is a large prime constant (from the golden ratio) used to improve distribution.
            // (h1 << 6) and (h1 >> 2) add additional mixing by shifting bits left and right.
            // This reduces hash collisions by ensuring small changes in input produce different hashes.
            return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
        }
    };

    // Define hash for sample_hap_t
    template <>
    struct hash<stoat::sample_hap_t> {
        size_t operator()(const stoat::sample_hap_t& sample_hap) const {
            return std::hash<std::string>()(sample_hap.sample + "#" + sample_hap.haplotype);
        }
    };


} // end namespace std

#endif
