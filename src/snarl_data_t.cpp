#include "snarl_data_t.hpp"

//#define DEBUG_SNARL_DATA_T

namespace stoat {

// Function to parse the snarl path file
std::unordered_map<std::string, std::vector<Snarl_data_t>> read_snarl_path(const std::string& file_path) {
    
    std::string line, chr, snarl, snarl_id, start_pos_str, end_pos_str, path_list, type_var, ref, depth;
    std::vector<Snarl_data_t> snarl_paths;
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
        std::getline(ss, snarl, '\t');   // snarl column
        std::getline(ss, snarl_id, '\t');   // snarl_id column
        std::getline(ss, path_list, '\t'); // paths column
        std::getline(ss, type_var, '\t');   // type_var column
        std::getline(ss, ref, '\t');   // ref column
        std::getline(ss, depth, '\t');   // depth column

        std::istringstream path_stream(path_list);
        std::istringstream type_stream(type_var);
        std::vector<std::string> type;
        size_t start_pos = std::stoull(start_pos_str);
        size_t end_pos = std::stoull(end_pos_str);
        std::string paths_str;
        bool first = true;

        // Reconstruct path string and count paths
        while (std::getline(path_stream, path_list, ',')) {
            if (!first) {
                paths_str += ",";
            }
            paths_str += path_list;
            first = false;
        }

        // create a vector of types
        while (std::getline(type_stream, type_var, ',')) {
            type.push_back(type_var);
        }

        if (chr != save_chr && !save_chr.empty()) {
            chr_to_snarls[save_chr] = std::move(snarl_paths);
            snarl_paths.clear();
        }

        save_chr = chr;

        std::pair<size_t, size_t> snarl_ids = stringToPair(snarl_id);
        std::vector<stoat::PathTraversal> paths = stringToVectorPath(paths_str);
        Snarl_data_t snarl_path(handlegraph::as_net_handle(std::stoll(snarl)), snarl_ids, paths, start_pos, end_pos, type, std::stoull(depth));
        snarl_paths.push_back(snarl_path);
    }

    // last chromosome
    chr_to_snarls[save_chr] = std::move(snarl_paths);
    file.close();

    // --- Print statistics ---
    std::cout << "\nSnarl statistics per chromosome:\n";
    for (const auto& [chr, snarls] : chr_to_snarls) {
        size_t total_paths = 0;
        for (const auto& snarl : snarls) {
            total_paths += snarl.snarl_paths.size();
        }
        std::cout << " > " << chr << ": " << snarls.size() << " snarls, " << total_paths << " total paths\n";
    }

    return chr_to_snarls;
}

void write_snarl_data_header(std::ostream& outstream) {
    outstream << "CHR\tSTART_POS\tEND_POS\tSNARL_HANDLEGRAPH\tSNARL\tPATHS\tPATH_LENGTHS\tREF\tDEPTH" << std::endl;
}

void write_snarl_data_fail_header(std::ostream& outstream) {
    outstream << "SNARL\tREASON" << std::endl;
}

// Node_traversal_t
Node_traversal_t::Node_traversal_t(const size_t &id, const bool &rev)
        : node_id(id), is_reverse(rev) {}

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

// Convert Edge_t to std::pair<size_t, size_t>
std::pair<size_t, size_t> Edge_t::print_pair_edge() const {
    return std::make_pair(edge.first.get_node_id(), edge.second.get_node_id());
}

// Convert Edge_t to std::string
std::string Edge_t::print_string_edge() const {
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
    this->paths.push_back(node_traversal);
}

    
// Add a node handle and extract information using the std::string representation
void PathTraversal::add_node_handle(const handlegraph::net_handle_t& node_h, const bdsg::SnarlDistanceIndex& distance_index) {
    // find the node orientation
    bool is_rev = distance_index.ends_at(node_h) != bdsg::SnarlDistanceIndex::END;
    Node_traversal_t node_traversal(distance_index.node_id(node_h), is_rev); 
    // Add node to path
    this->paths.push_back(node_traversal);
}
    
// add a node traversal to the path
void PathTraversal::add_node_traversal_t(const Node_traversal_t &node_trav) {
    this->paths.push_back(node_trav);
}

void PathTraversal::add_min_allele_len(size_t len){
    min_allele_len += len;
}

void PathTraversal::add_max_allele_len(size_t len){
    max_allele_len += len;
}

    
// TODO : change sum_path to definition using the length of the path including in the boundary nodes
// Matis ans : i don t know how to do it
std::string PathTraversal::get_allele_length() {
    if (paths.size() >= 3) {
        // If there is at least one node other than the boundaries
        if (min_allele_len != max_allele_len) {
            // a "complex" variant with no fixed size because of nested variants
            // return a range of possible lengths
            return std::to_string(min_allele_len) + "/" + std::to_string(max_allele_len);
        } else {
            // one length when simple variant or when nested variants are SNPs or balanced MNPs
            return std::to_string(min_allele_len);
        }        
    } else if (paths.size() == 2) {
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

    if (paths[0].get_node_id() > paths.back().get_node_id()) {
        // flip the path
        path_flip();
    }
}

// Flip the PathTraversal
void PathTraversal::path_flip() {
    std::reverse(paths.begin(), paths.end());

    for (size_t i = 0; i < paths.size(); ++i) {

        paths[i].set_is_reverse(!paths[i].get_is_reverse());    
    }
}

// convert PathTraversal to path representation
std::string PathTraversal::to_string() const {
    std::string result;
    for (const auto& node : paths) {
        result += node.to_string();
    }
    return result;
}

const std::vector<Node_traversal_t>& PathTraversal::get_paths() const { 
    return paths; 
};

// Get the size of the path
size_t PathTraversal::size() const {
    return paths.size();
}
    
std::string pairToString(const std::pair<size_t, size_t>& name) {
    std::ostringstream oss;
    oss << name.first << "_" << name.second;
    return oss.str();
}

std::pair<size_t, size_t> stringToPair(const std::string& str) {
    size_t underscorePos = str.find('_');
    if (underscorePos == std::string::npos) {
        throw std::runtime_error("Input std::string does not contain an underscore separator");
    }

    std::string firstPart = str.substr(0, underscorePos);
    std::string secondPart = str.substr(underscorePos + 1);

    size_t first = std::stoull(firstPart);
    size_t second = std::stoull(secondPart);

    return {first, second};
}

std::string vectorPathToString(const std::vector<stoat::PathTraversal>& vec_paths) {
    std::ostringstream oss;
    for (size_t i = 0; i < vec_paths.size(); ++i) {
        if (i > 0) oss << ",";
        oss << vec_paths[i].to_string();
    }
    return oss.str();
}

std::vector<stoat::PathTraversal> stringToVectorPath(std::string& input) {
    std::vector<stoat::PathTraversal> vec_paths;
    std::istringstream iss(input);
    std::string path_str;

    while (std::getline(iss, path_str, ',')) {
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

        vec_paths.emplace_back(std::move(path));
    }

    return vec_paths;
}

// Add a snarl
Snarl_data_t::Snarl_data_t(bdsg::net_handle_t snarl_, const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index) : 
    snarl(snarl_), start_positions(0), end_positions(0), depth(distance_index.get_depth(snarl_)) {
    snarl_ids = std::make_pair(distance_index.node_id(distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, false))),
                               distance_index.node_id(distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, true, false))));
}

Snarl_data_t::Snarl_data_t(bdsg::net_handle_t snarl_,
    std::pair<size_t, size_t> snarl_ids_,
    std::vector<PathTraversal> snarl_paths_,
    const size_t start_positions_, const size_t end_positions_,
    std::vector<std::string> type_variants_,
    size_t depth) :
    snarl(snarl_),
    snarl_ids(snarl_ids_),
    snarl_paths(std::move(snarl_paths_)),
    start_positions(start_positions_),
    end_positions(end_positions_),
    type_variants(std::move(type_variants_)),
    depth(depth) {}
    
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

    // unique_ptr<bdsg::PositionOverlay> position_overlay = std::make_unique<bdsg::PositionOverlay>(graph.get());

    // // Get root of snarl tree
    // handlegraph::net_handle_t root = distance_index->get_root();

    return std::make_tuple(
        std::move(distance_index),
        std::move(graph)
    );
}

std::vector<std::tuple<handlegraph::net_handle_t, 
    std::string, size_t, size_t, bool>> list_all_snarls_path_pos(
        const bdsg::SnarlDistanceIndex& distance_index, 
        handlegraph::PathHandleGraph& graph, 
        std::unordered_set<std::string>& ref_paths) {

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

    // vector with snarl info: handle, ref path name, start position, end position, on ref?
    std::vector<std::tuple<handlegraph::net_handle_t, std::string, size_t, size_t, bool>> snarls;
    // map: snarl handle (as string) -> position (chr, start, end)
    unordered_map<std::string, std::tuple<std::string, size_t, size_t>> snarls_pos;

    // function to check each element of the tree and fill the vector/map above
    function<void(handlegraph::net_handle_t)> process_tree_element = [&](handlegraph::net_handle_t net) {

        // try to get the position of this element in the snarl tree
        std::tuple<std::string, size_t, size_t> snarl_path_pos = get_net_start_position(net);
        bool is_on_ref = true;

        // if we couldn't find a position, use the parent's that we should have
        // found and saved earlier
        if (std::get<0>(snarl_path_pos).empty()) {
            auto par_net = distance_index.get_parent(net);
            snarl_path_pos = snarls_pos[distance_index.net_handle_as_string(par_net)];
            is_on_ref = false;
        }

        // save this position in the cache
        snarls_pos[distance_index.net_handle_as_string(net)] = snarl_path_pos;

        // save snarl
        if (distance_index.is_snarl(net)) {
            // handlegraph::net_handle_t snarl, chr_ref, pos, is_on_ref_bool
            snarls.push_back(std::make_tuple(net, std::get<0>(snarl_path_pos), std::get<1>(snarl_path_pos), std::get<2>(snarl_path_pos), is_on_ref));
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
    
    stoat::LOG_INFO("Total number of snarls : " + std::to_string(snarls.size()));
    return snarls;
}

std::vector<stoat::PathTraversal> convert_path_traversals(
    const bdsg::SnarlDistanceIndex& distance_index, 
    handlegraph::PathHandleGraph& graph, 
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
std::unordered_map<std::string, std::vector<Snarl_data_t>> loop_over_snarls_write(
    const bdsg::SnarlDistanceIndex& distance_index,
    std::vector<std::tuple<handlegraph::net_handle_t, std::string, size_t, size_t, bool>>& snarls,
    handlegraph::PathHandleGraph& graph,
    const std::string& output_file,
    const std::string& output_snarl_excluded,
    const size_t& children_threshold,
    const size_t& path_length_threshold,
    const size_t& cycle_threshold) {

    // start output files with headers
    std::ofstream out_snarl(output_file);
    std::ofstream out_fail(output_snarl_excluded);
    write_snarl_data_header(out_snarl);
    write_snarl_data_fail_header(out_fail);

    // we'll output a list of snarls for each chromosome
    std::unordered_map<std::string, std::vector<Snarl_data_t>> chr_to_snarls;

    // metrics for the log
    size_t paths_analyzed = 0;
    size_t snarls_failed = 0;
    size_t paths_failed = 0;

    // decompose the snarls in parallel
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < snarls.size(); ++i) {
        const auto& snarl_path_pos = snarls[i];
        handlegraph::net_handle_t snarl = std::get<0>(snarl_path_pos);
        std::pair<size_t, size_t> snarl_id = stoat::find_snarl_id(distance_index, snarl);
        std::string snarl_id_str = pairToString(snarl_id);

        // Count children
        size_t children = 0;
        distance_index.for_each_child(snarl, [&](const handlegraph::net_handle_t& net) {
            children++;
            return true;
        });

        if (children > children_threshold) {
            #pragma omp critical(out_fail)
            out_fail << snarl_id_str << "\ttoo_many_children = " << children << " children\n";
            snarls_failed++;
            continue;
        }

        // list of paths being explored
        // init with one starting at the first bound
        std::vector<std::vector<handlegraph::net_handle_t>> paths = {
            {distance_index.get_bound(snarl, false, true)}
        };

        // once a path reaches the other bound, add to this list of finished paths
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
                    // stop if path would be too long
                    if (path.size() + 1 > path_length_threshold) {
#pragma omp critical(out_fail)
                        out_fail << snarl_id_str << "\titeration_calculation_out = " << children << " children\n";
                        skip_snarl = true;
                        paths_failed++;
                        return true;
                    }

                    // otherwise extend path
                    paths.emplace_back(path);
                    paths.back().push_back(next_child);
                }
                return true;
            };

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

        // otherwise 
        auto path_travs = convert_path_traversals(distance_index, graph, finished_paths);
        // JEAN for now remaking allele_lengths, but I think it should go in PathTraversal
        std::vector<std::string> allele_lengths;
        for(auto path_trav: path_travs){
            allele_lengths.push_back(path_trav.get_allele_length());
        }

        const std::string& chr = std::get<1>(snarl_path_pos);
        if (chr.empty()) {continue;}
        // JEAN that also shouldn't happen and maybe we should raise an error

        size_t start_pos = std::get<2>(snarl_path_pos);
        size_t end_pos = std::get<3>(snarl_path_pos);
        size_t depth = distance_index.get_depth(snarl);
        std::string str_reference = std::get<4>(snarl_path_pos) ? "1" : "0";

        // Output result
        #pragma omp critical(out_snarl)
        out_snarl << chr << "\t" 
                    << start_pos << "\t" 
                    << end_pos << "\t" 
                    << handlegraph::as_integer(snarl) << "\t" 
                    << snarl_id_str << "\t"
                    << vectorPathToString(path_travs) << "\t"
                    << stoat::vectorToString(allele_lengths) << "\t"
                    << str_reference << "\t" 
                    << depth << "\n";

        // JEAN if all this stuff goes in Snarl_data_t anyway, why not put it there from the start and avoid carrying around those tuples?
        Snarl_data_t snarl_path(snarl, snarl_id, path_travs, start_pos, end_pos, allele_lengths, depth);
        
        #pragma omp critical(chr_to_snarls)
        chr_to_snarls[chr].emplace_back(std::move(snarl_path));

        paths_analyzed += path_travs.size();
    }

    stoat::LOG_INFO("Total number of snarl filtered : " + std::to_string(snarls_failed));
    stoat::LOG_INFO("Total number of paths : " + std::to_string(paths_analyzed));
    stoat::LOG_INFO("Total number of paths filtered : " + std::to_string(paths_failed));

    if (paths_analyzed == 0) {
        throw std::runtime_error("Total number of paths = 0. This may indicate that the graph does not contain a flagged reference path. Please use -r/--chr to specify the reference paths.");
    }

    for (const auto& [chr, snarls] : chr_to_snarls) {
        stoat::LOG_INFO("chr : " + chr + ", number of snarl : " + std::to_string(snarls.size()));
    }

    return chr_to_snarls;
}

} // end stoat namespace

// vg find -x ../snarl_data/fly.gbz -r 5176878:5176884 -c 10 | vg view -dp - | dot -Tsvg -o ../snarl_data/subgraph.svg
