#include "snarl_data_t.hpp"

//#define DEBUG_SNARL_DATA_T

namespace stoat {

// Function to parse the snarl path file
std::unordered_map<std::string, std::vector<Snarl_data_t>> read_snarl_path(const std::string& file_path) {
    
    std::string line, chr, snarl_handle, snarl_id, start_pos_str, end_pos_str, path_list, al_lens, ref, depth;
    std::vector<Snarl_data_t> snarls;
    std::unordered_map<std::string, std::vector<Snarl_data_t>> chr_to_snarls;
    std::ifstream file(file_path);
    std::string save_chr = "";

    std::getline(file, line);

    // --- Parse and validate header ---
    std::istringstream header_stream(line);
    std::vector<std::string> header_fields;
    std::string field;

    while (std::getline(header_stream, field, '\t')) {
        header_fields.push_back(field);
    }

    const std::vector<std::string> expected_header = {
        "CHR", "START_POS", "END_POS", "SNARL_HANDLEGRAPH",
        "SNARL", "PATHS", "PATH_LENGTHS", "REF", "DEPTH"
    };

    if (header_fields != expected_header) {
        std::ostringstream oss;
        oss << "Error: Invalid header format in file: " << file_path << "\n";
        oss << " > Expected: ";
        for (size_t i = 0; i < expected_header.size(); ++i) {
            oss << expected_header[i];
            if (i < expected_header.size() - 1) oss << "\t";
        }
        oss << "\n > Got:      ";
        for (size_t i = 0; i < header_fields.size(); ++i) {
            oss << header_fields[i];
            if (i < header_fields.size() - 1) oss << "\t";
        }
        throw std::runtime_error(oss.str());
    }

    // Process each line
    // JEAN assuming all the lines are sorted by chromosomes
    while (std::getline(file, line)) {
        std::istringstream ss(line);

        std::getline(ss, chr, '\t');   // chr column
        std::getline(ss, start_pos_str, '\t');   // pos column
        std::getline(ss, end_pos_str, '\t');   // pos column
        std::getline(ss, snarl_handle, '\t');   // snarl column
        std::getline(ss, snarl_id, '\t');   // snarl_id column
        std::getline(ss, path_list, '\t'); // paths column
        std::getline(ss, al_lens, '\t');   // al_lens column
        std::getline(ss, ref, '\t');   // ref column
        std::getline(ss, depth, '\t');   // depth column

        // if we're starting a new chromosome, save the current one and reinit
        // JEAN if the lines are not sorted by chromosome, we might be overwriting the previous chunk for that chromosome here
        if (chr != save_chr && !save_chr.empty()) {
            chr_to_snarls[save_chr] = std::move(snarls);
            snarls.clear();
        }
        save_chr = chr;

        // to parse the information from the different fields
        size_t start_pos = std::stoull(start_pos_str);
        size_t end_pos = std::stoull(end_pos_str);

        // create the Snarl object and add it to the list
        Snarl_data_t snarl(handlegraph::as_net_handle(std::stoll(snarl_handle)), snarl_id, path_list, al_lens, start_pos, end_pos, std::stoull(depth));
        snarls.push_back(snarl);
    }

    // last chromosome
    chr_to_snarls[save_chr] = std::move(snarls);
    file.close();

    // --- Print statistics ---
    std::cout << "\nSnarl statistics per chromosome:\n";
    for (const auto& [chr, snarls] : chr_to_snarls) {
        size_t total_paths = 0;
        for (const auto& snarl : snarls) {
            total_paths += snarl.paths.size();
        }
        std::cout << " > " << chr << ": " << snarls.size() << " snarls, " << total_paths << " total paths\n";
    }

    return chr_to_snarls;
}

// Node_traversal_t
Node_traversal_t::Node_traversal_t(const size_t &id, const bool &rev)
        : node_id(id), is_reverse(rev) {}

// Node_traversal_t from a string
Node_traversal_t::Node_traversal_t(const std::string& str) {
    is_reverse = str.at(0) == '<';
    node_id = std::stoi(std::string(str.begin()+1, str.end()));
}

// Convert Node_traversal_t to node + path representation [string]
std::string Node_traversal_t::to_string() const {
    return (is_reverse ? "<" : ">") + std::to_string(node_id);
}

// Getters for Node_traversal_t
size_t Node_traversal_t::get_node_id() const { return node_id; }
bool Node_traversal_t::get_is_reverse() const { return is_reverse; }

bool Node_traversal_t::operator==(const Node_traversal_t& other) const {
    return node_id == other.node_id && is_reverse == other.is_reverse;
}

// Edge_t
Edge_t::Edge_t(const Node_traversal_t &node_traversal_1, 
               const Node_traversal_t &node_traversal_2) :
    edge(std::make_pair(node_traversal_1, node_traversal_2)) {}

// return a flipped version of this edge
Edge_t Edge_t::get_flipped() const {
    Node_traversal_t first_node = Node_traversal_t(edge.first.get_node_id(),
                                                   !edge.first.get_is_reverse());
    Node_traversal_t second_node = Node_traversal_t(edge.second.get_node_id(),
                                                   !edge.second.get_is_reverse());
    return (Edge_t(second_node, first_node));
}
    
// Convert Edge_t to std::pair<size_t, size_t>
std::pair<size_t, size_t> Edge_t::get_node_pair() const {
    return std::make_pair(edge.first.get_node_id(), edge.second.get_node_id());
}

// Convert Edge_t to std::string
std::string Edge_t::to_string() const {
    return edge.first.to_string() + edge.second.to_string();
}

// Accessor to edge, useful for hashing and comparison
const std::pair<Node_traversal_t, Node_traversal_t>& Edge_t::get_edge() const {
    return edge;
}

bool Edge_t::operator==(const Edge_t &other) const {
    return edge == other.edge;
}

// Add a node with known orientation
void PathTraversal::add_node(const size_t& node, bool is_rev) {
    Node_traversal_t node_traversal(node, is_rev); 
    // Add node to path
    this->path.push_back(node_traversal);
}

    
// Add a node handle and extract information using the std::string representation
void PathTraversal::add_node_handle(const handlegraph::net_handle_t& node_h, const bdsg::SnarlDistanceIndex& distance_index) {
    // find the node orientation
    bool is_rev = distance_index.ends_at(node_h) != bdsg::SnarlDistanceIndex::END;
    Node_traversal_t node_traversal(distance_index.node_id(node_h), is_rev); 
    // Add node to path
    this->path.push_back(node_traversal);
}
    
// add a node traversal to the path
void PathTraversal::add_node_traversal_t(const Node_traversal_t &node_trav) {
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

const std::vector<Node_traversal_t>& PathTraversal::get_path() const { 
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
    
std::string vectorPathToString(const std::vector<stoat::PathTraversal>& vec_paths, bool allele_lengths) {
    std::ostringstream oss;
    for (size_t i = 0; i < vec_paths.size(); ++i) {
        if (i > 0) oss << ",";
        if (allele_lengths) {
            oss << vec_paths[i].get_allele_length();
        } else {
            oss << vec_paths[i].to_string();
        }
    }
    return oss.str();
}

// Snarl constructors
Snarl_data_t::Snarl_data_t(bdsg::net_handle_t net_, const bdsg::SnarlDistanceIndex& distance_index) : 
    net(net_), start_position(0), end_position(0), is_on_ref(false), depth(distance_index.get_depth(net_)) {
    ids = std::make_pair(distance_index.node_id(distance_index.get_node_from_sentinel(distance_index.get_bound(net, false, false))),
                               distance_index.node_id(distance_index.get_node_from_sentinel(distance_index.get_bound(net, true, false))));
}

Snarl_data_t::Snarl_data_t(bdsg::net_handle_t net_,
    std::pair<size_t, size_t> ids_,
    std::vector<PathTraversal> paths_,
    const size_t start_position_, const size_t end_position_,
    size_t depth) :
    net(net_),
    ids(ids_),
    paths(std::move(paths_)),
    start_position(start_position_),
    end_position(end_position_),
    depth(depth) {}

Snarl_data_t::Snarl_data_t(bdsg::net_handle_t net_,
    std::string snarl_ids_,
    std::string paths_,
    std::string allele_lengths_,
    const size_t start_position_, const size_t end_position_,
    size_t depth) :
    net(net_),
    start_position(start_position_),
    end_position(end_position_),
    depth(depth) {

    // extract the boundary nodes from the snarl ID
    size_t underscorePos = snarl_ids_.find('_');
    if (underscorePos == std::string::npos) {
        throw std::runtime_error("Input snarl ID " + snarl_ids_ + " does not contain an underscore separator");
    }
    size_t first_id = std::stoull(snarl_ids_.substr(0, underscorePos));
    size_t second_id = std::stoull(snarl_ids_.substr(underscorePos + 1));
    ids = std::make_pair(first_id, second_id);
    
    paths = string_to_path_traversals(paths_, allele_lengths_);
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
                path.add_node_traversal_t(stoat::Node_traversal_t(node_id, is_reverse));
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
    unique_ptr<bdsg::SnarlDistanceIndex>,
    unique_ptr<handlegraph::PathHandleGraph>>
    load_graph_tree(
        const std::string& graph_file, 
        const std::string& dist_file) {

    // Tell the IO library about libvg types.
    if (!stoat::io::register_libvg_io()) {
        throw std::runtime_error("error[stoat vgio]: Could not register libvg types with libvgio");
    }

    // Load the graph and make it a PathPositionHandleGraph
    unique_ptr<handlegraph::PathHandleGraph> graph = vg::io::VPKG::load_one<handlegraph::PathHandleGraph>(graph_file);

    // Load the distance index
    unique_ptr<bdsg::SnarlDistanceIndex> distance_index = std::make_unique<bdsg::SnarlDistanceIndex>();
    distance_index->deserialize(dist_file);

    return std::make_tuple(
        std::move(distance_index),
        std::move(graph)
    );
}

std::unordered_map<std::string, std::vector<Snarl_data_t>> list_all_snarls_with_pos(
        const bdsg::SnarlDistanceIndex& distance_index, 
        handlegraph::PathHandleGraph& graph, 
        std::unordered_set<std::string>& ref_paths) {

    // we'll output a list of snarls for each chromosome
    std::unordered_map<std::string, std::vector<Snarl_data_t>> chr_to_snarls;
    
    // load the position overlay needed to get position on reference paths
    bdsg::PositionOverlay ppo = bdsg::PositionOverlay(&graph);

    // were reference paths/chrs provided
    bool ref_path_provided = !ref_paths.empty();

    // Given a node handle (dist index), return a position (chr, start, end) if on chr reference path
    auto get_node_position = [&](handlegraph::net_handle_t node) -> std::tuple<std::string, size_t, size_t> {
        // path_name, position
        std::tuple<std::string, size_t, size_t> ret_pos;

        auto step_callback = [&](const handlegraph::step_handle_t& step_handle) {
            const auto path_handle = graph.get_path_handle_of_step(step_handle);
            const std::string& chr_path = graph.get_path_name(path_handle);

            // save this position if the path matched the provided chr names,
            // or no chr names provided and it's a reference path
            bool save_position = ref_path_provided
                ? (ref_paths.count(chr_path) > 0)
                : (graph.get_sense(path_handle) == handlegraph::PathSense::REFERENCE);

            if (save_position) {
                const size_t pos = ppo.get_position_of_step(step_handle);

                std::get<0>(ret_pos) = chr_path;
                std::get<1>(ret_pos) = pos;
                std::get<2>(ret_pos) = pos + distance_index.node_length(node);

                return false; // Stop iteration once a position was found
            }

            return true; // Continue iteration
        };

        // look for reference path across "steps" using node handle on the graph
        handlegraph::handle_t node_h = distance_index.get_handle(node, &graph);
        graph.for_each_step_on_handle(node_h, step_callback);
        
        return ret_pos;
    };

    // Given a handle in the dist index, return a position (chr, start, end)
    auto get_net_start_position = [&](handlegraph::net_handle_t net) -> std::tuple<std::string, size_t, size_t> {

        // if it's a node just retun the node position
        if (distance_index.is_node(net)) {
            return get_node_position(net);
        }

        // otherwise find the position of the boundary nodes and return MIN_POS,MAX_POS
        handlegraph::net_handle_t bnode1 = distance_index.get_bound(net, true, false);
        std::tuple<std::string, size_t, size_t> bnode1_p = get_node_position(bnode1);

        handlegraph::net_handle_t bnode2 = distance_index.get_bound(net, false, false); // verify false true ?
        std::tuple<std::string, size_t, size_t> bnode2_p = get_node_position(bnode2);

        // if a boundary doesn't touch a reference path return empty position
        // Check if the std::string part of the pair is empty
        if (std::get<0>(bnode1_p).empty()) return bnode1_p;
        if (std::get<0>(bnode2_p).empty()) return bnode2_p;

#ifdef DEBUG_SNARL_DATA_T
        //LOG_DEBUG(
        // JEAN this should probably be always on, not just for debugging
        assert(std::get<0>(bnode1_p) == std::get<0>(bnode2_p)); // Ensure they are on the same reference path
#endif

        // start pos as the end position of the upstream boundary node
        // end pos as the start position of the downstream boundary node
        if (std::get<1>(bnode1_p) < std::get<1>(bnode2_p)) {
            return make_tuple(std::get<0>(bnode1_p),
                              std::get<2>(bnode1_p),
                              std::get<1>(bnode2_p));
        } else {
            return make_tuple(std::get<0>(bnode1_p),
                              std::get<2>(bnode2_p),
                              std::get<1>(bnode1_p));
        }

    };

    // cache the position of snarls already processed
    // in particular, the parent snarls are processed first, so we can easily get and
    // reuse those position for children that are not on the reference
    // map a snarl handle as string to a position tuple (chr, start, end)
    // JEAN this tuple could be a small struct...
    unordered_map<std::string, std::tuple<std::string, size_t, size_t>> snarl_pos_cache;

    // function to check each element of the tree and fill the vector/map above
    function<void(handlegraph::net_handle_t)> process_tree_element = [&](handlegraph::net_handle_t net) {

        // try to get the position of this element in the snarl tree
        std::tuple<std::string, size_t, size_t> snarl_path_pos = get_net_start_position(net);
        bool is_on_ref = true;
        
        // if we couldn't find a position, use the parent's that we should have
        // found and saved earlier
        if (std::get<0>(snarl_path_pos).empty()) {
            auto par_net = distance_index.get_parent(net);
            snarl_path_pos = snarl_pos_cache[distance_index.net_handle_as_string(par_net)];
            is_on_ref = false;
        }

        // save this position in the cache
        snarl_pos_cache[distance_index.net_handle_as_string(net)] = snarl_path_pos;

        // save snarl if positioned relative to a chromosome of interest
        // if the chr is still empty at this point, it was not on a chr of interest and neither any of its parents
        if (distance_index.is_snarl(net) && !std::get<0>(snarl_path_pos).empty()) {
            Snarl_data_t snarl(net, distance_index);
            snarl.start_position = std::get<1>(snarl_path_pos);
            snarl.end_position = std::get<2>(snarl_path_pos);
            snarl.is_on_ref = is_on_ref;

            // add this snarl to the appropriate chromosome list
            chr_to_snarls[std::get<0>(snarl_path_pos)].emplace_back(std::move(snarl));
        }

        // explore children (if they can have children)
        if (!distance_index.is_node(net) && !distance_index.is_sentinel(net)) {
            distance_index.for_each_child(net, process_tree_element);
        }
    };

    // process the tree starting from the root
    handlegraph::net_handle_t root = distance_index.get_root();
    distance_index.for_each_child(root, process_tree_element);

    // JEAN the path overlay pointer used to be "cleaned up with" .reset() before, not sure if that's still necessary?

    // print the total number of snarls found
    size_t total_snarls = 0;
    for (const auto& [chr, snarls] : chr_to_snarls) {
        total_snarls += snarls.size();
    }
    stoat::LOG_INFO("Total number of snarls : " + std::to_string(total_snarls));
    
    return chr_to_snarls;
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

        // JEAN I think this is just for convenience but should make sure the path/edge orientation is NOT important when parsing the AT
        path_trav.check_path_flip();
        
        path_travs.push_back(path_trav);
    }

    return path_travs;
}

// {chr : matrix(snarl, paths, start_pos, end_pos, type)}
std::unordered_map<std::string, std::vector<Snarl_data_t>> write_snarls_with_paths(
    const bdsg::SnarlDistanceIndex& distance_index,
    std::unordered_map<std::string, std::vector<Snarl_data_t>>& chr_to_snarls,
    handlegraph::PathHandleGraph& graph,
    const std::string& output_dir,
    const size_t& children_threshold,
    const size_t& path_length_threshold,
    const size_t& cycle_threshold) {

    // JEAN would be nice to let the user decide how to name that?
    std::string output_snarl_excluded = output_dir + "/snarl_not_analyse.tsv";
    std::string output_file = output_dir + "/snarl_analyse.tsv";

    // start output files with headers
    std::ofstream out_snarl(output_file);
    std::ofstream out_fail(output_snarl_excluded);
    // write_snarl_data_header(out_snarl);
    out_snarl << "CHR\tSTART_POS\tEND_POS\tSNARL_HANDLEGRAPH\tSNARL\tPATHS\tPATH_LENGTHS\tREF\tDEPTH" << std::endl;
    out_fail << "CHR\tSTART_POS\tEND_POS\tSNARL_HANDLEGRAPH\tSNARL\tFILTER\tCHILDREN" << std::endl;

    // we'll output a list of snarls for each chromosome
    std::unordered_map<std::string, std::vector<Snarl_data_t>> out_chr_to_snarls;
    // JEAN we could maybe just update the input chr_to_snarls but not sure if it works well with omp?

    // metrics for the log
    size_t paths_analyzed = 0;
    size_t snarls_failed = 0;
    size_t paths_failed = 0;

    // decompose the snarls in each chromosome in parallel
    for (auto& [chr, snarls] : chr_to_snarls) {
#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < snarls.size(); ++i) {
            Snarl_data_t snarl = snarls[i];
            std::string snarl_id_str = pairToString(snarl.ids);
            
            // Count children
            size_t children = 0;
            distance_index.for_each_child(snarl.net, [&](const handlegraph::net_handle_t& net) {
                children++;
                return true;
            });

            if (children > children_threshold) {
#pragma omp critical(out_fail)
                out_fail << chr << "\t" 
                          << snarl.start_position << "\t" 
                          << snarl.end_position << "\t" 
                          << handlegraph::as_integer(snarl.net) << "\t" 
                          << snarl_id_str << "\ttoo_many_children\t" << children << "\n";
                snarls_failed++;
                continue;
            }

            // list of paths being explored
            // init with one starting at the first bound
            std::vector<std::vector<handlegraph::net_handle_t>> paths = {
                {distance_index.get_bound(snarl.net, false, true)}
            };

            // once a path reaches the other bound, add to this list of finished paths
            // we use this simpler structure here with vector and handles so avoid doing too much work working preparing a PathTraversal that might be filtered?
            // JEAN maybe simpler to use PathTraversal objects directly here?
            std::vector<std::vector<handlegraph::net_handle_t>> finished_paths;

            // records if we gave up and need to skip this snarl if a path it too long
            // JEAN do we want this? Or should we test the snarl with just the "short" paths? Not sure
            bool skip_snarl = false;

            // extend the paths
            while (!paths.empty()) {
                std::vector<handlegraph::net_handle_t> path = std::move(paths.back());
                paths.pop_back();

                auto add_to_path = [&](const handlegraph::net_handle_t& next_child) {
                    // If this is the bound of the snarl then we're done && next_child is different that the first node
                    if (distance_index.is_sentinel(next_child)) {
                        size_t next_child_node_id = distance_index.node_id(distance_index.get_node_from_sentinel(next_child));
                        size_t first_element_path_node_id = distance_index.node_id(distance_index.get_node_from_sentinel(path[0]));
                        if (next_child_node_id != first_element_path_node_id) {
                            finished_paths.emplace_back(path);
                            finished_paths.back().push_back(next_child);
                        }

                    } else {
                        // stop if cycling too much
                        int next_child_count = std::count(path.begin(), path.end(), next_child);
                        if(next_child_count > cycle_threshold) {
                            return true;
                        }

                        // otherwise extend path
                        paths.emplace_back(path);
                        paths.back().push_back(next_child);
                    }
                    return true;
                };

                // stop if path would be too long
                // JEAN we could do this in the add_to_path function and return false to abort but we would need a newer version of vg and remake all our distance index files, so not doing that yet...
                if (path.size() + 1 > path_length_threshold) {
#pragma omp critical(out_fail)
                    out_fail << chr << "\t" 
                             << snarl.start_position << "\t" 
                             << snarl.end_position << "\t" 
                             << handlegraph::as_integer(snarl.net) << "\t" 
                             << snarl_id_str << "\tpath_too_long\t" << children << "\n";
                    skip_snarl = true;
                    paths_failed++;
                    break;
                }
                
                // Follow edges from the last element in path
                if (!path.empty()) {
                    distance_index.follow_net_edges(path.back(), &graph, false, add_to_path);
                }

            }

            // this snarl should be filtered (some paths too long)
            if (skip_snarl) {
                snarls_failed++;
                continue;
            }

            // stop if only one path
            if (finished_paths.size() < 2) {
                snarls_failed++;
                continue;
            } // avoid special case single path

            // otherwise prepare PathTraversal objects
            auto path_travs = convert_path_traversals(distance_index, graph, finished_paths);
            std::string str_reference = snarl.is_on_ref ? "1" : "0";

            // Output result
#pragma omp critical(out_snarl)
            out_snarl << chr << "\t" 
                      << snarl.start_position << "\t" 
                      << snarl.end_position << "\t" 
                      << handlegraph::as_integer(snarl.net) << "\t" 
                      << snarl_id_str << "\t"
                      << vectorPathToString(path_travs) << "\t"
                      << vectorPathToString(path_travs, true) << "\t"
                      << str_reference << "\t" 
                      << snarl.depth << "\n";

            paths_analyzed += path_travs.size();

            // add the enumerated paths to the Snarl object
            snarl.paths = std::move(path_travs);
#pragma omp critical(out_chr_to_snarls)
            out_chr_to_snarls[chr].emplace_back(std::move(snarl));

        }
    }
    
    stoat::LOG_INFO("Total number of snarl filtered : " + std::to_string(snarls_failed));
    stoat::LOG_INFO("Total number of paths : " + std::to_string(paths_analyzed));
    stoat::LOG_INFO("Total number of paths filtered : " + std::to_string(paths_failed));

    if (paths_analyzed == 0) {
        throw std::runtime_error("Total number of paths = 0. This may indicate that the graph does not contain a flagged reference path. Please use -r/--chr to specify the reference paths.");
    }

    for (const auto& [chr, snarls] : out_chr_to_snarls) {
        stoat::LOG_INFO("chr : " + chr + ", number of snarl : " + std::to_string(snarls.size()));
    }

    return out_chr_to_snarls;
}

} // end stoat namespace

// vg find -x ../snarl_data/fly.gbz -r 5176878:5176884 -c 10 | vg view -dp - | dot -Tsvg -o ../snarl_data/subgraph.svg
