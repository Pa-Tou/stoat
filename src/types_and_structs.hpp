#ifndef types_and_structs
#define types_and_structs

#include <vector>
#include <string>
#include <iostream>

#include <bdsg/snarl_distance_index.hpp>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>

#include "log.hpp"
#include "feature_tables.hpp"

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


        // Get the same node but flipped
        node_traversal_t get_flipped() { return node_traversal_t(node_id, !is_reverse); }

        bool operator==(const node_traversal_t& other) const;
        bool operator!=(const node_traversal_t& other) const;
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
    std::string get_allele_length(size_t& count_path_lengths_warn) const;
        
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
// This will return an empty string if vec_paths is empty
// If is_allele_included is not empty, then only included those paths for which the corresponding entry in is_allele_included is true.
std::string vectorPathToString(const std::vector<PathTraversal>& vec_paths, bool allele_lengths = false, const std::vector<bool>& is_allele_included = {});

// Get a vector of path traversals from their string representation (from vectorPathToString())
std::vector<stoat::PathTraversal> string_to_path_traversals(const std::string& path_string, const std::string& path_lengths_string);

// Get a vector of node_traversal_t's representing a single path given as a string
std::vector<stoat::node_traversal_t> string_to_path_node_traversal(const std::string& path_string);

// Get a string representing a path of node_traversal_t's
std::string path_node_traversal_to_string(const std::vector<stoat::node_traversal_t>& path);

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

// Used to store the allele assignment for a vector of samples. 
// This struct only makes sense when it is paired with a vector of samples or sample_hap_t's (as a member of a snarl_info_t). 
// Each entry in alleles is the identifier of an allele
// The allele identifiers are also used as indices into other per-allele vectors, so they should start from 0 and go up to the number of alleles-1.
// allele_count is the number of non-inf alleles, or 1 + max value in alleles_by_sample. This is not enforced though
struct allele_by_sample_t{
    size_t allele_count;
    std::vector<size_t> alleles;

    allele_by_sample_t(size_t count, std::vector<size_t> alleles) :
        allele_count(count),
        alleles(std::move(alleles)) {}
    allele_by_sample_t() :
        allele_count(0), alleles(0) {}
};

// This holds snarl information from the distance index and graph: id and reference coordinates of the snarl, genotypes, reference coordinates, walks, 
// the sequences of the walks and sets samples/haplotypes taking each walk
// The latter four may be empty if they were not filled in by the SnarlDataCollection
struct snarl_info_t {
    public:

        // Constructor from elements
        snarl_info_t(stoat::node_traversal_t start_node, 
                     stoat::node_traversal_t end_node, 
                     std::string ref_path, 
                     size_t start_position, 
                     size_t end_position, 
                     size_t depth,
                     GenotypeTable& genotypes, 
                     const std::vector<stoat::sample_hap_t>& all_sample_haplotypes,
                     const allele_by_sample_t& alleles_by_sample,
                     const std::vector<PathTraversal>& walks_by_allele, 
                     const std::vector<std::string>& sequences_by_allele) :

                     start_node(start_node), end_node(end_node), 
                     ref_path(ref_path), start_position(start_position), end_position(end_position), depth(depth),
                     genotypes(genotypes), all_sample_haplotypes(all_sample_haplotypes), alleles_by_sample(alleles_by_sample),
                     walks_by_allele(walks_by_allele), sequences_by_allele(sequences_by_allele)
                     {};

        // Start and end nodes, both pointing into the snarl
        stoat::node_traversal_t start_node;
        stoat::node_traversal_t end_node;

        // The reference chromosome/path
        std::string ref_path; 

        // Start and end offset along the reference path
        size_t start_position;
        size_t end_position;

        // The depth of the snarl in the snarl tree
        size_t depth;

        // A genotype table of the counts of each allele for each sample (see feature_table.hpp).
        // Alleles have the same numbering as walks_by_allele and sequences_by_allele
        GenotypeTable& genotypes;

        // This stores all the sample/haplotypes 
        const std::vector<stoat::sample_hap_t>& all_sample_haplotypes;

        // For each haplotype in all_sample_haplotypes, an assignment to an allele number
        const allele_by_sample_t& alleles_by_sample;

        // For each allele, the walk through the snarl
        // The walk includes the boundary nodes of the snarl
        // Nested chains are included in the walk as the boundary node of the chain going into the chain, an empty node 0 going forward, 
        //    and the other bound of the chain going out.
        // When a walk leaves the snarl and comes back, it will include the boundary node of the snarl leaving it, an empty node 0 going backward, and the 
        //    boundary of the snarl going back in
        const std::vector<PathTraversal>& walks_by_allele;

        // For each allele, what is its sequence?
        // The sequences don't include sequences of the boundary nodes
        const std::vector<std::string>& sequences_by_allele;

};

struct node_info_t {
    public:

        // Constructor from elements
        node_info_t(stoat::node_traversal_t node, 
                     std::string ref_path, 
                     size_t position, 
                     size_t depth,
                     GenotypeTable& genotypes, 
                     const std::vector<stoat::sample_hap_t>& all_sample_haplotypes,
                     const allele_by_sample_t& alleles_by_sample,
                     const std::string& sequences_by_allele) :

                     node(node), position(position), ref_path(ref_path), 
                     depth(depth), genotypes(genotypes), all_sample_haplotypes(all_sample_haplotypes), 
                     alleles_by_sample(alleles_by_sample), sequences_by_allele(sequences_by_allele) {};

        // Start and end nodes, both pointing into the node
        stoat::node_traversal_t node;

        // The reference chromosome/path
        std::string ref_path; 

        // Start and end offset along the reference path
        size_t position;

        // The depth of the node in the node tree
        size_t depth;

        // A genotype table of the counts of each allele for each sample (see feature_table.hpp).
        // Alleles have the same numbering as walks_by_allele and sequences_by_allele
        GenotypeTable& genotypes;

        // This stores all the sample/haplotypes 
        const std::vector<stoat::sample_hap_t>& all_sample_haplotypes;

        // For each haplotype in all_sample_haplotypes, an assignment to an allele number
        const allele_by_sample_t& alleles_by_sample;

        const std::string& sequences_by_allele

};

enum phenotype_type_t { BINARY = 1, BINARY_COVAR, QUANTITATIVE, EQTL };

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
