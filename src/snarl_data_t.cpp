#include "snarl_data_t.hpp"

//#define DEBUG_SNARL_DATA_T

namespace stoat {

// Function to parse the snarl path file
std::unordered_map<std::string, std::vector<Snarl_data_t>> parse_snarl_path(const std::string& file_path) {
    
    std::string line, chr, snarl, snarl_id, start_pos_str, end_pos_str, path_list, type_var, ref, depth;
    std::vector<Snarl_data_t> snarl_paths;
    std::unordered_map<std::string, std::vector<Snarl_data_t>> chr_snarl_matrix;
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
        "SNARL", "PATHS", "TYPE", "REF", "DEPTH"
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
        size_t start_pos = std::stoi(start_pos_str);
        size_t end_pos = std::stoi(end_pos_str);
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
            chr_snarl_matrix[save_chr] = std::move(snarl_paths);
            snarl_paths.clear();
        }

        save_chr = chr;

        std::pair<size_t, size_t> snarl_ids = stringToPair(snarl_id);
        std::vector<stoat::Path_traversal_t> paths = stringToVectorPath(paths_str);
        Snarl_data_t snarl_path(handlegraph::as_net_handle(std::stoll(snarl)), snarl_ids, paths, start_pos, end_pos, type, std::stoi(depth));
        snarl_paths.push_back(snarl_path);
    }

    // last chromosome
    chr_snarl_matrix[save_chr] = std::move(snarl_paths);
    file.close();

    // --- Print statistics ---
    std::cout << "\nSnarl statistics per chromosome:\n";
    for (const auto& [chr, snarls] : chr_snarl_matrix) {
        size_t total_paths = 0;
        for (const auto& snarl : snarls) {
            total_paths += snarl.snarl_paths.size();
        }
        std::cout << " > " << chr << ": " << snarls.size() << " snarls, " << total_paths << " total paths\n";
    }

    return chr_snarl_matrix;
}

void write_snarl_data_output(std::ostream& outstream) {
    outstream << "CHR\tSTART_POS\tEND_POS\tSNARL_HANDLEGRAPH\tSNARL\tPATHS\tTYPE\tREF\tDEPTH" << std::endl;
}

void write_snarl_data_fail(std::ostream& outstream) {
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

// add a node traversal to the path
void Path_traversal_t::add_node_traversal_t(const Node_traversal_t &node) {
    this->paths.push_back(node);
}

// convert Path_traversal_t to path representation
std::string Path_traversal_t::to_string() const {
    std::string result;
    for (const auto& node : paths) {
        result += node.to_string();
    }
    return result;
}

const std::vector<Node_traversal_t>& Path_traversal_t::get_paths() const { 
    return paths; 
};

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

    size_t first = std::stoul(firstPart);
    size_t second = std::stoul(secondPart);

    return {first, second};
}

std::string vectorPathToString(const std::vector<stoat::Path_traversal_t>& vec_paths) {
    std::ostringstream oss;
    for (size_t i = 0; i < vec_paths.size(); ++i) {
        if (i > 0) oss << ",";
        oss << vec_paths[i].to_string();
    }
    return oss.str();
}

std::vector<stoat::Path_traversal_t> stringToVectorPath(std::string& input) {
    std::vector<stoat::Path_traversal_t> vec_paths;
    std::istringstream iss(input);
    std::string path_str;

    while (std::getline(iss, path_str, ',')) {
        stoat::Path_traversal_t path;
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
    std::vector<Path_traversal_t> snarl_paths_,
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

Path::Path() {}

// Add a node with known orientation
void Path::addNode(const size_t& node, bool orient) {
    nodes.push_back(node);
    orients.push_back(orient);
}

// Add a node handle and extract information using the std::string representation
bool Path::addNodeHandle(const handlegraph::net_handle_t& node_h, const bdsg::SnarlDistanceIndex& distance_index) {

    // Found the orientation
    bool node_o = distance_index.ends_at(node_h) == bdsg::SnarlDistanceIndex::END;

    // Add node to path
    nodes.push_back(distance_index.node_id(node_h));
    orients.push_back(node_o);
    return node_o;
}

// Get the std::string representation of the path
Path_traversal_t Path::print() const {
    Path_traversal_t out_path;
    for (size_t i = 0; i < nodes.size(); ++i) {
        Node_traversal_t node_traversal(nodes[i], !orients[i]); // because is reverse is false for '>' and true for '<'
        out_path.add_node_traversal_t(node_traversal);
    }
    return out_path;
}

// Flip the path orientation
void Path::flip() {
    std::reverse(nodes.begin(), nodes.end());
    std::reverse(orients.begin(), orients.end());
    for (size_t i = 0; i < orients.size(); ++i) {
        if (nodes[i] == 0) {
            continue;
        }
        orients[i] = !orients[i];    
    }
}

// Get the size of the path
size_t Path::size() const {
    return nodes.size();
}

// Count the number of reversed nodes
size_t Path::nreversed() const {
    return std::count(orients.begin(), orients.end(), false);
}

// Function to calculate the type of variant
// tuple<std::string, size_t, size_t, size_t>
//minimum_distance, maximum_distance, the number of nodes in the path (including boundary nodes), sum_path
std::vector<std::string> calcul_pos_type_variant(const std::vector<std::tuple<size_t, size_t, size_t>>& list_length_paths) {
    std::vector<std::string> list_type_variant;

    for (const auto& tuple_info : list_length_paths) {
        size_t min_length = std::get<0>(tuple_info);
        size_t max_length = std::get<1>(tuple_info);
        size_t path_length = std::get<2>(tuple_info); // The number of nodes in the path

        if (path_length >= 3) {
            // If there is at least one node representing this allele
            if (min_length != max_length) { // Case nested
                // If this is a complex variant (includes nested variants), then return a range of possible lengths
                std::string nested = std::to_string(min_length) + "/" + std::to_string(max_length);
                list_type_variant.push_back(nested);
            } else { // Case nodes chain (ex : INS+SNP+...)
                list_type_variant.push_back(std::to_string(min_length));
            }

        } else if (path_length == 2) { // case Deletion
            list_type_variant.push_back("0");
        } else { // Case path_lengths is empty or == 1
            // This should probably never happen right ?
            stoat::LOG_WARN("path_lengths is empty");
        }
    }
    return list_type_variant;
}

std::tuple<
    unique_ptr<bdsg::SnarlDistanceIndex>,
    unique_ptr<handlegraph::PathHandleGraph>,
    handlegraph::net_handle_t,
    unique_ptr<bdsg::PositionOverlay>>
    parse_graph_tree(
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

    unique_ptr<bdsg::PositionOverlay> position_overlay = std::make_unique<bdsg::PositionOverlay>(graph.get());

    // Get root of snarl tree
    handlegraph::net_handle_t root = distance_index->get_root();

    return std::make_tuple(
        std::move(distance_index),
        std::move(graph),
        root,
        std::move(position_overlay)
    );
}

void follow_edges(
                const bdsg::SnarlDistanceIndex& distance_index,
                std::vector<std::vector<handlegraph::net_handle_t>>& finished_paths,
                const std::vector<handlegraph::net_handle_t>& path,
                std::vector<std::vector<handlegraph::net_handle_t>>& paths,
                handlegraph::PathHandleGraph& graph,
                const bool& cycle) {

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

            if (cycle) { // Case where we find a loop
                return false;
            }
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

std::vector<std::tuple<handlegraph::net_handle_t, 
    std::string, size_t, size_t, bool>> save_snarls(
        const bdsg::SnarlDistanceIndex& distance_index, 
        handlegraph::net_handle_t& root,
        handlegraph::PathHandleGraph& graph, 
        std::unordered_set<std::string>& ref_chr,
        const bdsg::PositionOverlay& ppo) {

    std::vector<std::tuple<handlegraph::net_handle_t, std::string, size_t, size_t, bool>> snarls;
    unordered_map<std::string, std::tuple<std::string, size_t, size_t>> snarls_pos;
    bool get_ref = ref_chr.empty() ? true : false;

    // Given a node handle (dist index), return a position if on chr reference path
    auto get_node_position = [&](handlegraph::net_handle_t node) -> std::tuple<std::string, size_t, size_t> { // node : handlegraph::net_handle_t
        handlegraph::handle_t node_h = distance_index.get_handle(node, &graph);

        // path_name, position
        std::tuple<std::string, size_t, size_t> ret_pos;

        auto step_callback = [&](const handlegraph::step_handle_t& step_handle) {
            const auto path_handle = graph.get_path_handle_of_step(step_handle);

            // Determine if this path is a candidate (reference or in ref_chr)
            bool is_candidate = get_ref
                ? (graph.get_sense(path_handle) == handlegraph::PathSense::REFERENCE)
                : (ref_chr.count(graph.get_path_name(path_handle)) > 0);

            if (is_candidate) {
                const std::string& chr_path = graph.get_path_name(path_handle);
                const size_t pos = ppo.get_position_of_step(step_handle);

                std::get<0>(ret_pos) = chr_path;
                std::get<1>(ret_pos) = pos + distance_index.node_length(node);
                std::get<2>(ret_pos) = pos + 1;

                return false; // Stop iteration once a candidate is found
            }

            return true; // Continue iteration
        };

        graph.for_each_step_on_handle(node_h, step_callback);
        return ret_pos;
    };

    auto get_net_start_position = [&](handlegraph::net_handle_t net) -> std::tuple<std::string, size_t, size_t> {

        if (distance_index.is_node(net)) {
            return get_node_position(net);
        }

        handlegraph::net_handle_t bnode1 = distance_index.get_bound(net, true, false);
        std::tuple<std::string, size_t, size_t> bnode1_p = get_node_position(bnode1);

        handlegraph::net_handle_t bnode2 = distance_index.get_bound(net, false, false); // verify false true ?
        std::tuple<std::string, size_t, size_t> bnode2_p = get_node_position(bnode2);

        // Check if the std::string part of the pair is empty
        if (std::get<0>(bnode1_p).empty()) return bnode1_p;
        if (std::get<0>(bnode2_p).empty()) return bnode2_p;

        #ifdef DEBUG_SNARL_DATA_T
            //LOG_DEBUG(
            assert(std::get<0>(bnode1_p) == std::get<0>(bnode2_p)); // Ensure they are on the same reference path
        #endif

        size_t start;
        size_t end;

        // smaller boundary is the start
        // larger boundary is the end
        if (std::get<1>(bnode1_p) < std::get<1>(bnode2_p)) {
            start = std::get<1>(bnode1_p);
            end = std::get<2>(bnode2_p);
        } else {
            start = std::get<1>(bnode2_p);
            end = std::get<2>(bnode1_p);
        }

        // tuple<std::string, size_t, size_t> snarl_start_end;
        return make_tuple(std::get<0>(bnode1_p), start, end);
    };

    function<void(handlegraph::net_handle_t)> save_snarl_tree_node;
    save_snarl_tree_node = [&](handlegraph::net_handle_t net) {

        std::tuple<std::string, size_t, size_t> snarl_path_pos = get_net_start_position(net);
        bool ref = true;

        // if we couldn't find a position, use the parent's that we should have
        // found and saved earlier
        if (std::get<0>(snarl_path_pos).empty()) {
            auto par_net = distance_index.get_parent(net);
            snarl_path_pos = snarls_pos[distance_index.net_handle_as_string(par_net)];
            ref = false;
        }

        // save this position
        snarls_pos[distance_index.net_handle_as_string(net)] = snarl_path_pos;

        // save snarl
        if (distance_index.is_snarl(net)) {
            // handlegraph::net_handle_t snarl, chr_ref, pos, is_on_ref_bool
            snarls.push_back(std::make_tuple(net, std::get<0>(snarl_path_pos), std::get<1>(snarl_path_pos), std::get<2>(snarl_path_pos)-1, ref));
        }

        // explore children
        if (!distance_index.is_node(net) && !distance_index.is_sentinel(net)) {
            distance_index.for_each_child(net, save_snarl_tree_node);
        }
    };

    distance_index.for_each_child(root, save_snarl_tree_node);
    stoat::LOG_INFO("Total number of snarls : " + std::to_string(snarls.size()));
    return snarls;
}

std::tuple<std::vector<stoat::Path_traversal_t>, std::vector<std::string>> fill_pretty_paths(
    const bdsg::SnarlDistanceIndex& distance_index, 
    handlegraph::PathHandleGraph& graph, 
    std::vector<std::vector<handlegraph::net_handle_t>>& finished_paths) {

    // list of paths
    std::vector<stoat::Path_traversal_t> pretty_paths;

    // seq_net, minimum_distance, maximun_distance, size_path, sum_path
    // Used to calculate the type of variant
    std::vector<std::tuple<size_t, size_t, size_t>> seq_net_paths;

    for (const auto& path : finished_paths) {
        Path ppath;
        size_t minimum_distance=0;
        size_t maximun_distance=0;
        std::vector<size_t> size_node;
        size_node.resize(path.size(), 0);

        for (int i=0; i<path.size(); i++) {
            handlegraph::net_handle_t net = path[i];

            if (distance_index.is_sentinel(net)) {
                net = distance_index.get_node_from_sentinel(net);
            }

            // Node case
            if (distance_index.is_node(net)) {
                bool rev = ppath.addNodeHandle(net, distance_index);
                handlegraph::nid_t node_start_id = distance_index.node_id(net);
                handlegraph::handle_t node_handle = graph.get_handle(node_start_id);
                size_node[i] = graph.get_length(node_handle);
            }

            // Trivial chain case
            else if (distance_index.is_trivial_chain(net)) {
                bool rev = ppath.addNodeHandle(net, distance_index);
                auto stn_start = distance_index.starts_at_start(net) ? distance_index.get_bound(net, false, true) : distance_index.get_bound(net, true, true);
                handlegraph::nid_t node_start_id = distance_index.node_id(stn_start);
                handlegraph::handle_t net_trivial_chain = graph.get_handle(node_start_id);
                size_node[i] = graph.get_length(net_trivial_chain);
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

                ppath.addNodeHandle(nodl, distance_index);

                bool chain_2node = true;
                int child_count = 0;
                size_t sum_node = 0;

                distance_index.for_each_child(net, [&](const handlegraph::net_handle_t& child) {
                    ++child_count;
                    if (!distance_index.is_node(child)) {
                        chain_2node = false;
                        return false; // stop early
                    } else {
                        sum_node += graph.get_length(graph.get_handle(distance_index.node_id(child)));
                    }
                    return true;
                });
                
                if (!(chain_2node && child_count == 2)) {
                    ppath.addNode(0, true);
                } else {
                    size_node[i] = sum_node;
                }
                ppath.addNodeHandle(nodr, distance_index);

                // Fail case 
                #ifdef DEBUG_SNARL_DATA_T
                // stoat::LOG_DEBUG();
                assert(distance_index.maximum_length(net) != static_cast<size_t>(INT_MAX) && "Overflow max distance");
                assert(distance_index.minimum_length(net) != static_cast<size_t>(INT_MAX) && "Overflow min distance");
                #endif

                // Add the minimum/maximum lengths of the chain
                minimum_distance += distance_index.minimum_length(net);
                maximun_distance += distance_index.maximum_length(net);
            }
        }

        if (ppath.nreversed() > ppath.size() / 2) {
            ppath.flip();
        }

        for (size_t i = 1; i < size_node.size()-1; ++i) {
            maximun_distance += size_node[i];
            minimum_distance += size_node[i];
        }

        pretty_paths.push_back(ppath.print());
        // The number of nodes (may be chains) in the path, including boundary nodes
        size_t size_path = ppath.size();
        seq_net_paths.push_back(std::make_tuple(minimum_distance, maximun_distance, size_path));
    }

    std::vector<std::string> type_variants = calcul_pos_type_variant(seq_net_paths);
    return std::make_tuple(pretty_paths, type_variants);
}

// {chr : matrix(snarl, paths, start_pos, end_pos, type)}
std::unordered_map<std::string, std::vector<Snarl_data_t>> loop_over_snarls_write(
    const bdsg::SnarlDistanceIndex& distance_index,
    std::vector<std::tuple<handlegraph::net_handle_t, std::string, size_t, size_t, bool>>& snarls,
    handlegraph::PathHandleGraph& graph,
    const std::string& output_file,
    const std::string& output_snarl_not_analyse,
    const size_t& children_threshold,
    const size_t& path_length_threshold,
    const size_t& cycle_threshold) {

    std::ofstream out_snarl(output_file);
    std::ofstream out_fail(output_snarl_not_analyse);

    write_snarl_data_output(out_snarl);
    write_snarl_data_fail(out_fail);

    std::unordered_map<std::string, std::vector<Snarl_data_t>> chr_snarl_matrix;
    size_t paths_number_analysis = 0;
    size_t snarl_fail = 0;
    size_t paths_fail = 0;

    // Parallel loop
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
            paths_fail += children;
            snarl_fail++;
            continue;
        }

        // Path exploration
        std::vector<std::vector<handlegraph::net_handle_t>> paths = {
            {distance_index.get_bound(snarl, false, true)}
        };

        std::vector<std::vector<handlegraph::net_handle_t>> finished_paths;

        size_t itr = 0;
        bool break_snarl = false;

        while (!paths.empty()) {
            std::vector<handlegraph::net_handle_t> path = std::move(paths.back());
            paths.pop_back();

            std::unordered_map<handlegraph::net_handle_t, size_t> dict_path_occ;
            bool cycle = false;

            for (const auto& net : path) {
                if (++dict_path_occ[net] > cycle_threshold + 1) {
                    cycle = true;
                    break;
                }
            }

            if (itr++ > path_length_threshold) {
                #pragma omp critical(out_fail)
                out_fail << snarl_id_str << "\titeration_calculation_out = " << children << " children\n";
                break_snarl = true;
                paths_fail++;
                break;
            }

            follow_edges(distance_index, finished_paths, path, paths, graph, cycle);
        }

        if (break_snarl) {continue;}

        auto [pretty_paths, type_variants] = fill_pretty_paths(distance_index, graph, finished_paths);
        
        if (pretty_paths.size() < 2) {
            snarl_fail++;
            continue;
        } // avoid special case single path

        const std::string& chr = std::get<1>(snarl_path_pos);
        if (chr.empty()) {continue;}

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
                    << vectorPathToString(pretty_paths) << "\t"
                    << stoat::vectorToString(type_variants) << "\t"
                    << str_reference << "\t" 
                    << depth << "\n";
        Snarl_data_t snarl_path(snarl, snarl_id, pretty_paths, start_pos, end_pos, type_variants, depth);
        
        #pragma omp critical(chr_snarl_matrix)
        chr_snarl_matrix[chr].emplace_back(std::move(snarl_path));

        paths_number_analysis += pretty_paths.size();
    }

    stoat::LOG_INFO("Total number of snarl filtered : " + std::to_string(snarl_fail));
    stoat::LOG_INFO("Total number of paths : " + std::to_string(paths_number_analysis));
    stoat::LOG_INFO("Total number of paths filtered : " + std::to_string(paths_fail));

    if (paths_number_analysis == 0) {
        throw std::runtime_error("Total number of paths = 0. This may indicate that the graph does not contain a flagged reference path. Please use -r/--chr to specify the reference paths.");
    }

    for (const auto& [chr, snarls] : chr_snarl_matrix) {
        stoat::LOG_INFO("chr : " + chr + ", number of snarl : " + std::to_string(snarls.size()));
    }

    return chr_snarl_matrix;
}

} // end stoat namespace

// vg find -x ../snarl_data/fly.gbz -r 5176878:5176884 -c 10 | vg view -dp - | dot -Tsvg -o ../snarl_data/subgraph.svg
