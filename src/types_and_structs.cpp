#include "types_and_structs.hpp"

//#define DEBUG_STRUCTS

namespace stoat {

// node_traversal_t
node_traversal_t::node_traversal_t(const size_t &id, const bool &rev)
        : node_id(id), is_reverse(rev) {}

// node_traversal_t from a string
node_traversal_t::node_traversal_t(const std::string& str) {
    is_reverse = str.at(0) == '<';
    node_id = std::stoi(std::string(str.begin()+1, str.end()));
}

// Convert node_traversal_t to node + path representation [string]
std::string node_traversal_t::to_string() const {
    return (is_reverse ? "<" : ">") + std::to_string(node_id);
}

// Getters for node_traversal_t
size_t node_traversal_t::get_node_id() const { return node_id; }
bool node_traversal_t::get_is_reverse() const { return is_reverse; }

bool node_traversal_t::operator==(const node_traversal_t& other) const {
    return node_id == other.node_id && is_reverse == other.is_reverse;
}

// edge_t
edge_t::edge_t(const node_traversal_t &node_traversal_1, 
               const node_traversal_t &node_traversal_2) :
    edge(std::make_pair(node_traversal_1, node_traversal_2)) {}

// return a flipped version of this edge
edge_t edge_t::get_flipped() const {
    node_traversal_t first_node = node_traversal_t(edge.first.get_node_id(),
                                                   !edge.first.get_is_reverse());
    node_traversal_t second_node = node_traversal_t(edge.second.get_node_id(),
                                                   !edge.second.get_is_reverse());
    return (edge_t(second_node, first_node));
}
    
// Convert edge_t to std::pair<size_t, size_t>
std::pair<size_t, size_t> edge_t::get_node_pair() const {
    return std::make_pair(edge.first.get_node_id(), edge.second.get_node_id());
}

// Convert edge_t to std::string
std::string edge_t::to_string() const {
    return edge.first.to_string() + edge.second.to_string();
}

// Accessor to edge, useful for hashing and comparison
const std::pair<node_traversal_t, node_traversal_t>& edge_t::get_edge() const {
    return edge;
}

bool edge_t::operator==(const edge_t &other) const {
    return edge == other.edge;
}

// Add a node with known orientation
void PathTraversal::add_node(const size_t& node, bool is_rev) {
    node_traversal_t node_traversal(node, is_rev); 
    // Add node to path
    this->path.push_back(node_traversal);
}

    
// Add a node handle and extract information using the std::string representation
void PathTraversal::add_node_handle(const handlegraph::net_handle_t& node_h, const bdsg::SnarlDistanceIndex& distance_index) {
    // find the node orientation
    bool is_rev = distance_index.ends_at(node_h) != bdsg::SnarlDistanceIndex::END;
    node_traversal_t node_traversal(distance_index.node_id(node_h), is_rev); 
    // Add node to path
    this->path.push_back(node_traversal);
}
    
// add a node traversal to the path
void PathTraversal::add_node_traversal_t(const node_traversal_t &node_trav) {
    this->path.push_back(node_trav);
}

void PathTraversal::add_min_allele_len(size_t len){
    min_allele_len += len;
}

void PathTraversal::add_max_allele_len(size_t len){
    max_allele_len += len;
}

void PathTraversal::set_allele_length_from_string(std::string al_len_str){
    // parse the string, could be either one size or MIN/MAX (see get_allele_length below)
    std::istringstream al_len_stream(al_len_str);
    std::string al_len;
    std::vector<size_t> min_max_al_lens;
    while (std::getline(al_len_stream, al_len, '/')) {
        min_max_al_lens.push_back(std::stoull(al_len));
    }
    // set the allele length range
    min_allele_len = min_max_al_lens[0];
    if (min_max_al_lens.size() == 1) {
        max_allele_len = min_max_al_lens[0];
    } else {
        max_allele_len = min_max_al_lens[1];
    }
}

    
// TODO : change sum_path to definition using the length of the path including in the boundary nodes
// Matis ans : i don t know how to do it
std::string PathTraversal::get_allele_length () const {
    if (path.size() >= 3) {
        // If there is at least one node other than the boundaries
        if (min_allele_len != max_allele_len) {
            // a "complex" variant with no fixed size because of nested variants
            // return a range of possible lengths
            return std::to_string(min_allele_len) + "/" + std::to_string(max_allele_len);
        } else {
            // one length when simple variant or when nested variants are SNPs or balanced MNPs
            return std::to_string(min_allele_len);
        }        
    } else if (path.size() == 2) {
        // if only the boundary nodes, it's a deletion
        return "0";
    } else {
        // This should probably never happen right ?
        stoat::LOG_WARN("path_lengths is empty");
        return "NA";
    }
}
    
// Check and flip the path if necessary to ensure consistent orientation
void PathTraversal::check_path_flip() {
    // Check if the path is already in the good orientation (aka min ID >> max ID)
    // JEAN this is quite bad, needs better condition to flip or not, and also deal with potential >0 added previously
    // JEAN maybe slightly better if we make sure the first bound is first on the reference (if on the reference path)
    
    if (path[0].get_node_id() > path.back().get_node_id()) {
        // flip the path
        path_flip();
    }
}

// Flip the PathTraversal
void PathTraversal::path_flip() {
    std::reverse(path.begin(), path.end());

    for (size_t i = 0; i < path.size(); ++i) {
        // JEAN maybe here never flip >0?
        path[i].set_is_reverse(!path[i].get_is_reverse());    
    }
}

// convert PathTraversal to path representation
std::string PathTraversal::to_string() const {
    std::string result;
    for (const auto& node : path) {
        result += node.to_string();
    }
    return result;
}

const std::vector<node_traversal_t>& PathTraversal::get_path() const { 
    return path;
};

// Get the size of the path
size_t PathTraversal::size() const {
    return path.size();
}

    // JEAN add to Snarl_data_t
std::string pairToString(const std::pair<size_t, size_t>& name) {
    std::ostringstream oss;
    oss << name.first << "_" << name.second;
    return oss.str();
}
    
std::string vectorPathToString(const std::vector<stoat::PathTraversal>& vec_paths, bool allele_lengths, const std::vector<bool>& is_allele_included) {
    std::ostringstream oss;
    for (size_t i = 0; i < vec_paths.size(); ++i) {
        if (is_allele_included.size() == 0 || is_allele_included[i]) {
            if (i > 0) oss << ",";
            if (allele_lengths) {
                oss << vec_paths[i].get_allele_length();
            } else {
                oss << vec_paths[i].to_string();
            }
        }
    }
    return oss.str();
}

std::vector<stoat::PathTraversal> string_to_path_traversals(const std::string& path_string, const std::string& path_lengths_string) {

    std::vector<stoat::PathTraversal> paths;
    if (path_string.size() != 0 && path_string != ".") {
        // extract the paths from the input string (comma separated)
        std::istringstream paths_iss(path_string);
        std::string path_str;
        while (std::getline(paths_iss, path_str, ',')) {
            stoat::PathTraversal path;
            size_t i = 0;
            const size_t len = path_str.size();
            while (i < len) {
                // Assume direction is always present and correct
                bool is_reverse = (path_str[i] == '<');
                ++i; // Move past '<' or '>'
                // Parse node_id
                size_t node_id = 0;
                while (i < len && path_str[i] >= '0' && path_str[i] <= '9') {
                    node_id = node_id * 10 + (path_str[i++] - '0');
                }
                path.add_node_traversal_t(stoat::node_traversal_t(node_id, is_reverse));
            }
            paths.emplace_back(std::move(path));
        }
    }
    if (path_lengths_string.size() != 0 && path_lengths_string != ".") {

        // extract the allele length information (also comma separated)
        std::istringstream al_lens_stream(path_lengths_string);
        std::string al_lens;
        int path_idx = 0;
        while (std::getline(al_lens_stream, al_lens, ',')) {
            if (path_idx > paths.size()-1) {
                paths.emplace_back();
            }
            paths[path_idx].set_allele_length_from_string(al_lens);
            path_idx++;
        }    
    }

    return paths;
}

std::tuple<
    std::unique_ptr<bdsg::SnarlDistanceIndex>,
    std::unique_ptr<handlegraph::PathHandleGraph>>
    load_graph_tree(
        const std::string& graph_file, 
        const std::string& dist_file) {

    // Tell the IO library about libvg types.
    if (!stoat::io::register_libvg_io()) {
        throw std::runtime_error("error[stoat vgio]: Could not register libvg types with libvgio");
    }

    // Load the graph and make it a PathPositionHandleGraph
    std::unique_ptr<handlegraph::PathHandleGraph> graph = vg::io::VPKG::load_one<handlegraph::PathHandleGraph>(graph_file);

    // Load the distance index
    std::unique_ptr<bdsg::SnarlDistanceIndex> distance_index = std::make_unique<bdsg::SnarlDistanceIndex>();
    distance_index->deserialize(dist_file);

    return std::make_tuple(
        std::move(distance_index),
        std::move(graph)
    );
}

std::vector<stoat::PathTraversal> convert_path_traversals(
    const bdsg::SnarlDistanceIndex& distance_index, 
    const handlegraph::PathHandleGraph& graph, 
    std::vector<std::vector<handlegraph::net_handle_t>>& finished_paths) {

    // list of paths
    std::vector<stoat::PathTraversal> path_travs;

    for (const auto& path : finished_paths) {
        PathTraversal path_trav;
        // saves the size of each node/element in the path
        std::vector<size_t> size_node;
        size_node.resize(path.size(), 0);

        for (int i=0; i<path.size(); i++) {
            handlegraph::net_handle_t net = path[i];

            if (distance_index.is_sentinel(net)) {
                net = distance_index.get_node_from_sentinel(net);
            }

            // Node case
            if (distance_index.is_node(net)) {
                path_trav.add_node_handle(net, distance_index);
                handlegraph::nid_t node_start_id = distance_index.node_id(net);
                handlegraph::handle_t node_handle = graph.get_handle(node_start_id);
                size_node[i] = graph.get_length(node_handle);
            }

            // Trivial chain case
            else if (distance_index.is_trivial_chain(net)) {
                path_trav.add_node_handle(net, distance_index);
                // it's a trivial chain so there is just one node inside. get it's node ID and size
                handlegraph::nid_t triv_node_id = distance_index.node_id(distance_index.get_bound(net, false, true));
                handlegraph::handle_t triv_node_h = graph.get_handle(triv_node_id);
                size_node[i] = graph.get_length(triv_node_h);
            }

            // Chain case (can be nested snarl or just chain nodes)
            else if (distance_index.is_chain(net)) {
                handlegraph::net_handle_t nodl, nodr;
                if (distance_index.starts_at_start(net)) {
                    nodl = distance_index.get_bound(net, false, true);
                    nodr = distance_index.get_bound(net, true, false);
                } else {
                    nodl = distance_index.get_bound(net, true, true);
                    nodr = distance_index.get_bound(net, false, false);
                }

                path_trav.add_node_handle(nodl, distance_index);

                // does the chain contain only two nodes
                bool chain_2node = true;
                int child_count = 0;

                distance_index.for_each_child(net, [&](const handlegraph::net_handle_t& child) {
                    ++child_count;
                    if (!distance_index.is_node(child)) {
                        chain_2node = false;
                        return false; // stop early
                    }
                    return true;
                });
                
                if (!(chain_2node && child_count == 2)) {
                    path_trav.add_node(0, false);
                }

                path_trav.add_node_handle(nodr, distance_index);

                // Fail case 
                #ifdef DEBUG_SNARL_DATA_T
                // stoat::LOG_DEBUG();
                assert(distance_index.maximum_length(net) != static_cast<size_t>(INT_MAX) && "Overflow max distance");
                assert(distance_index.minimum_length(net) != static_cast<size_t>(INT_MAX) && "Overflow min distance");
                #endif

                // Add the minimum/maximum lengths of the chain
                path_trav.add_min_allele_len(distance_index.minimum_length(net));
                path_trav.add_max_allele_len(distance_index.maximum_length(net));
            }
        }

        // add to min/max distance the size of the nodes and trivial chains saved above
        // but without including the boundary nodes
        for (size_t i = 1; i < size_node.size()-1; ++i) {
            path_trav.add_min_allele_len(size_node[i]);
            path_trav.add_max_allele_len(size_node[i]);
        }

        // here we flip the path in an attempt to make it consistent with the other path/allele traversal
        // that we will compare later. It's not absolutely necessary but convenient and might speed up the genotyping a bit
        path_trav.check_path_flip();
        
        path_travs.push_back(path_trav);
    }

    return path_travs;
}

sample_hap_t::sample_hap_t(const handlegraph::PathHandleGraph& graph, const handlegraph::path_handle_t& path) {

    std::string path_name = graph.get_path_name(path);

    std::stringstream stream(path_name);
    if (std::getline(stream, sample, '#')) {
        std::getline(stream, haplotype);
    } else {
        haplotype = "";
    }

    #ifdef DEBUG
    std::string test_sample;
    if (graph.get_sense(path) == handlegraph::PathSense::GENERIC) {
        // Generic paths only have a locus, so return whatever that is
        test_sample = graph.get_locus_name(path);
    } else {
        test_sample = graph.get_sample_name(path);
    }
    assert(sample == test_sample);
    #endif
}


} // end stoat namespace

// vg find -x ../snarl_data/fly.gbz -r 5176878:5176884 -c 10 | vg view -dp - | dot -Tsvg -o ../snarl_data/subgraph.svg
