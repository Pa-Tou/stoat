#include "path_partitioner.hpp"
#include "log.hpp"
#include <fstream>

#define DEBUG_PATH_PARTITIONER

using namespace stoat;
namespace stoat_graph {


// This is supposed to partition the paths in the snarl by the walks they take through the netgraph.
// Instead of explicitly enumerating the paths, it actually finds the sets of edges that each path takes.
// But since paths may loop, it actually finds the order and number of outgoing edges from each node.
// I think this is equivalent to partitioning by the actual sets of unique walks.
std::vector<size_t> partition_embedded_paths_in_snarl(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                          const net_handle_t& snarl,
                          const std::vector<stoat::sample_hap_t>& all_sample_haplotypes) {

    #ifdef DEBUG_PATH_PARTITIONER
    std::cerr <<  "Get walk sets of " << distance_index.net_handle_as_string(snarl) << std::endl;;
    #endif

    // Map each sample(plus haplotype) to its index in all_sample_haplotypes
    std::map<stoat::sample_hap_t, std::size_t> sample_to_index;
    for (size_t i = 0 ; i < all_sample_haplotypes.size() ; i++) {
        sample_to_index[all_sample_haplotypes[i]] = i;
    }

    // For each path (along all_sample_haplotypes), the index of the set it is currently in
    // Everything starts out in the same set, 0, which represents not being in the snarl
    std::vector<size_t> old_sets (all_sample_haplotypes.size(), 0);
    size_t old_set_count = 1;

    // TODO: I think this is unnecessary
    //// First, check for paths that don't go through this snarl
    //std::vector<handlegraph::PathSense> senses = {handlegraph::PathSense::GENERIC,
    //                                              handlegraph::PathSense::REFERENCE,
    //                                              handlegraph::PathSense::HAPLOTYPE};
    //
    ////TODO: This could also use steps_of_handle() but it doesn't seem to work 
    //for (const auto& sense : senses) {
    //    graph.for_each_step_of_sense(distance_index.get_handle(distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, true)), &graph),
    //        sense, [&](const handlegraph::step_handle_t& step) {
    //        old_sets[sample_to_index[stoat::get_sample_and_haplotype(graph, graph.get_path_handle_of_step(step))]] = 1;
    //        old_set_count = 2;
    //    });
    //}


    // A struct representing a path's edge going out from child
    // Represents the next node and orientation and the offset along the path of the edge
    // Use (0,0,0,false) as a default value to indicate that the path didn't traverse this child
    // This becomes a linked list when a path traverses the child multiple times (this doesn't happen often) 
    struct path_edge_t {
        size_t offset;           // The offset along the path
        size_t additional_edge;  // If this path traverses the child multiple times, index into additional_steps for others. max() for end of linked list 
        handlegraph::nid_t id;   // The id of the next node
        bool rev;                // Is the edge going to the right side of the next node?

        path_edge_t() : offset(0), additional_edge(std::numeric_limits<size_t>::max()), id(0), rev(false){};
        path_edge_t(size_t offset, size_t edge, handlegraph::nid_t id, bool rev) : 
            offset(offset), additional_edge(edge), id(id), rev(rev){};
    };


    // Helper function to check the edges leaving child
    // This gets called in for_each_child, and also from the start node going in
    // Goes through the edges leaving a child and checks with paths take which edges.
    // old_sets and old_set_count are updated to split up sets that take different edges at this node 
    // If require_node is true (used for the bounds going in), then any path that doesn't go through this node
    // is considered to not be in the snarl and its set number is reset to 0.
    // Second required node is true if this is the boundary node checked, in which case the path should be skipped
    // if it didn't pass through the first boundary node 
    // For this to work, check_outgoing_edges must be called on the two bounds last
    auto check_outgoing_edges = [&] (const handlegraph::net_handle_t& child, bool go_left, bool require_node, bool second_required_node) {

        #ifdef DEBUG_PATH_PARTITIONER
        std::cerr << "At snarl child " << distance_index.net_handle_as_string(child) << " going " << (go_left ? "left" : "right") << std::endl;
        #endif
        
        std::vector<path_edge_t> next_steps (all_sample_haplotypes.size());
        // For paths that go through multiple steps 
        std::vector<path_edge_t> additional_steps;
        
        // Get a handle to the end of the child
        handlegraph::handle_t handle;
        if ( distance_index.is_trivial_chain(child)) {
            handle =  (go_left ? graph.flip(distance_index.get_handle(child, &graph)) 
                               : distance_index.get_handle(child, &graph));
        } else if (distance_index.is_sentinel(child)) {
            handle = distance_index.get_handle(child, &graph);
        } else {
            handle = distance_index.get_handle(distance_index.get_bound(child, go_left, false), &graph);
        }
        
        std::vector<handlegraph::PathSense> senses = {handlegraph::PathSense::GENERIC,
                                                      handlegraph::PathSense::REFERENCE,
                                                      handlegraph::PathSense::HAPLOTYPE};

        #ifdef DEBUG_PATH_PARTITIONER
        std::cerr << "\tgraph handle " << graph.get_id(handle) << " going " << (graph.get_is_reverse(handle) ? "left" : "right") << std::endl;
        #endif

        for (const auto& sense : senses) {
            graph.for_each_step_of_sense(handle, sense, [&](const handlegraph::step_handle_t& step) {
                // For each step on the node handle, keep track of which paths take different steps

                #ifdef DEBUG_PATH_PARTITIONER
                std::cerr << "\ton path " << graph.get_path_name(graph.get_path_handle_of_step(step)) << std::endl;
                #endif

                sample_hap_t step_sample_haplotype (graph, graph.get_path_handle_of_step(step));

                // If this is not a sample of interest, skip it
                if (!sample_to_index.count(step_sample_haplotype)) {
                    return true;
                }
        
                //Do we go forwards in the path? We need to check the direction of the handle in the path
                bool go_forwards = graph.get_is_reverse(handle) == graph.get_is_reverse(graph.get_handle_of_step(step));
        
                //In the case where a path doesn't go all the way through the snarl, stop when the path stops
                if ((go_forwards && !graph.has_next_step(step)) || (!go_forwards && !graph.has_previous_step(step))){
                    return true;
                }
        
                //Get the next step and make an edge
                handlegraph::step_handle_t next_step = go_forwards ? graph.get_next_step(step) : graph.get_previous_step(step);
                handlegraph::handle_t next_handle = graph.get_handle_of_step(next_step);
        
                #ifdef DEBUG_PATH_PARTITIONER
                std::cerr << "" << "\t\tgoing to " << graph.get_id(next_handle) << std::endl;
                #endif
        
                path_edge_t edge (graph.get_position_of_step(step), 
                                  std::numeric_limits<size_t>::max(),
                                  graph.get_id(next_handle), 
                                  graph.get_is_reverse(next_handle));
        
                size_t sample_num = sample_to_index[step_sample_haplotype];
        
                if (next_steps[sample_num].id == 0) {
                    // If this path hasn't been seen before
                    next_steps[sample_num] = std::move(edge);
                } else if (next_steps[sample_num].offset > edge.offset) {
                    // If the new edge comes before the edge stored in the vector, replace it
                    edge.additional_edge = additional_steps.size();
                    additional_steps.emplace_back(std::move(next_steps[sample_num]));
                    next_steps[sample_num] = std::move(edge);
                } else {
                    // If the new edge comes after something in additional_steps, walk through the linked list to find its place

                    // The index into the linked list. Need a bool to know if it's an index into next_steps or additional_steps
                    bool old_edge_first = true;
                    size_t old_edge_index = sample_num;
                    size_t old_next_edge = next_steps[old_edge_index].additional_edge;
                    while (old_next_edge != std::numeric_limits<size_t>::max() &&
                           additional_steps[old_next_edge].offset < edge.offset) {
                        // Step through the linked list
                        old_edge_first = false;
                        old_edge_index = old_next_edge;
                        old_next_edge = additional_steps[old_edge_index].additional_edge;
                    }
                    //old_edge_index now points to the item just smaller than edge
                    if (old_edge_first) {
                        next_steps[old_edge_index].additional_edge = additional_steps.size();
                    } else {
                        additional_steps[old_edge_index].additional_edge = additional_steps.size();
                    }
                    edge.additional_edge = old_next_edge;
                    additional_steps.emplace_back(std::move(edge));
                }
        
                return true;
            });
        }
        
        // We now have the edges for all paths going out of the node in one direction
        // Now we want get an "intermediate set" for each path based on the edge(s) it took from this node/direction
        // Equality in this case is that the edges (as node id and orientation) are in the same order along the path
        // This is later used to split the paths into sets for each pair of old_set and intermediate_set

        #ifdef DEBUG_PATH_PARTITIONER
        // To check that we properly used all additional edges, count how many we follow
        size_t additional_edge_count = 0;
        #endif
        
        //This maps each edge list for this node/direction to the intermediate set index
        std::map<std::vector<std::pair<handlegraph::nid_t, bool>>, size_t> edge_to_intermediate_set;

        //Everything starts in the same set, representing not going through this node
        std::vector<size_t> intermediate_sets (old_sets.size(), 0);
        size_t intermediate_set_count = 1;
        for (size_t path_i = 0 ; path_i < next_steps.size() ; path_i++) {
            const path_edge_t& edge = next_steps[path_i];
            if (edge.id != 0) {
                //If this is a real edge
                std::vector<std::pair<handlegraph::nid_t, bool>> edge_list;
                edge_list.emplace_back(edge.id, edge.rev);
                size_t next_i = edge.additional_edge;
                while (next_i != std::numeric_limits<size_t>::max()) {
                    edge_list.emplace_back(additional_steps[next_i].id, additional_steps[next_i].rev);
                    auto& new_edge = additional_steps[next_i]; 
                    next_i = new_edge.additional_edge;
                    #ifdef DEBUG_PATH_PARTITIONER
                    additional_edge_count++;
                    #endif
                }
                if (edge_to_intermediate_set.count(edge_list) == 0) {
                    edge_to_intermediate_set[edge_list] = intermediate_set_count;
                    intermediate_set_count++;
                }
                intermediate_sets[path_i] = edge_to_intermediate_set[edge_list];
            }
        }
        #ifdef DEBUG_PATH_PARTITIONER
        std::cerr << "Intermediate sets: " << std::endl;
        for (size_t i = 0 ; i < all_sample_haplotypes.size() ; i++) {
            std::cerr << "" << "\t" << all_sample_haplotypes[i] << ": " << intermediate_sets[i] << std::endl;
        } 
        assert(additional_edge_count == additional_steps.size());
        #endif
        
        // We now have an old set and an intermediate set for each path
        // Assign the path to a new set. Everything gets a new set
        std::vector<size_t> new_sets (intermediate_sets.size(), std::numeric_limits<size_t>::max());

        // Map pairs of <old_set, intermediate_set> to new set number
        std::map<std::pair<size_t, size_t>, size_t>  old_to_new_set;

        // Reserve 0 for paths that didn't go through this snarl
        old_to_new_set[std::make_pair(0,0)] = 0;
        size_t new_set_count = 1;

        for (size_t path_i = 0 ; path_i < new_sets.size() ; path_i++) {
            if (require_node && (intermediate_sets[path_i] == 0 || (second_required_node && old_sets[path_i] == 0 ))) {
                // If this is a boundary node and the path never traverses it, then add it to the set
                // that doesn't go through the snarl
                new_sets[path_i] = 0;
            } else {
                std::pair<size_t, size_t> old_set (old_sets[path_i], intermediate_sets[path_i]);
                if (old_to_new_set.count(old_set) == 0) {
                    old_to_new_set[old_set] = new_set_count;
                    new_set_count++;
                }
                new_sets[path_i] = old_to_new_set[old_set];
            }
        }
        
        old_sets = std::move(new_sets);
        old_set_count = new_set_count;
        
        #ifdef DEBUG_PATH_PARTITIONER
        std::cerr << "New sets: " << std::endl;
        for (size_t i = 0 ; i < all_sample_haplotypes.size() ; i++) {
            std::cerr << "" << "\t" << all_sample_haplotypes[i] << ": " << old_sets[i] << std::endl;
        } 
        #endif
    };

    // Now do the work of going through the edges for the start bound and each child in both directions


    if (!distance_index.is_regular_snarl(snarl, true, &graph)) {
        // Go through each child of the snarl and check the paths on outgoing edges.
        // Split up sets if the paths have different edges leaving this child
        // TODO: This is doubling the work because each edges is looked at twice
        distance_index.for_each_child(snarl, [&] (const handlegraph::net_handle_t& child) {
            for (bool go_left : {true, false}) {
                check_outgoing_edges(child, go_left, false, false);
            }
            return true;
        });// end for_each_child of the snarl
    }

    check_outgoing_edges(distance_index.get_bound(snarl, false, true), false, true, false);
    check_outgoing_edges(distance_index.get_bound(snarl, true, true), false, true, true);

    // Now go through the sets and change 0 to inf, and decrement all others
    for (size_t i = 0 ; i < old_sets.size() ; i++) {
        if (old_sets[i] == 0) {
            old_sets[i] = std::numeric_limits<size_t>::max();
        } else {
            --old_sets[i];
        }
    }



    #ifdef DEBUG_PATH_PARTITIONER
    std::vector<std::set<sample_hap_t>> sample_sets_by_allele(old_set_count-1);
    std::cerr << "Check for " << all_sample_haplotypes.size() << " sample" << std::endl;
    for (size_t i = 0 ; i < all_sample_haplotypes.size() ; i++) {
        if (old_sets[i] != std::numeric_limits<size_t>::max()) {
            sample_sets_by_allele[old_sets[i]].emplace(all_sample_haplotypes[i]);
        }
    }
    std::cerr << "Found walk sets " << std::endl;
    for (const auto& s : sample_sets_by_allele) {
        std::cerr <<  "Set" << std::endl;
        for (const auto& x : s) {
            std::cerr << "" << "\t" << x << std::endl;
        }
    }
    #endif

    // old_sets is now a vector that assigns each sample/haplotype in all_sample_haplotypes to a allele number
    // Return this
    return old_sets;
}



// Get the traversals through the snarl from the gbwt
// This is heavily based on vg/haplotype_extracter.cpp
size_t get_gbwt_traversals(const handlegraph::PathPositionHandleGraph& graph, const gbwt::GBWT& gbwt, const bdsg::SnarlDistanceIndex& distance_index,     
                           const handlegraph::net_handle_t& snarl,
                           std::vector<std::tuple<std::vector<handlegraph::net_handle_t>, gbwt::SearchState, size_t>>& finished_paths) {
    #ifdef DEBUG_PATH_PARTITIONER
    std::cerr << "Get threads through snarl " << distance_index.net_handle_as_string(snarl) << std::endl;
    #endif

    size_t path_count = 0;


    // Get the bounds of the snarl: start facing in and end facing out. It doesn't matter which one is which since we get all start-end traversals
    handlegraph::net_handle_t start_net = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, true));
    handlegraph::net_handle_t end_net = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, true, false));

    handlegraph::handle_t start_in = distance_index.get_handle(start_net, &graph);


    // The bounds leaving the snarl
    handlegraph::net_handle_t snarl_start = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, false));
    handlegraph::net_handle_t snarl_end = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, true, false));

    // Look up the start node in GBWT and start a path
    gbwt::node_type start_node = gbwt::Node::encode(graph.get_id(start_in), graph.get_is_reverse(start_in));
    std::vector<handlegraph::net_handle_t> first_path = {start_net};
    gbwt::SearchState first_state = gbwt.find(start_node);

    // The list of intermediate paths and the gbwt::SearchState they end on
    // The search state encompasses the haplotypes that followed the path
    // The size_t is an identifier for the path, since some paths may be identical
    std::vector<std::tuple<std::vector<handlegraph::net_handle_t>, gbwt::SearchState, size_t>> intermediate_paths;

    // The GBWT traversals may split up in nested chains, meaning that we'd get a separate SearchState for the same
    // path through the netgraph. Since these separate threads might then take different paths, we need to make sure that
    // the distinct walks get the same id.
    // Keep track of the path id and the walk that it takes

    // For a path id and its length, has a next step been found?
    // This is used to determine if any new next step should be given a new path id
    // One path can keep the same id as the current path but each  other branching path needs a new one. 
    // The length is the number of net_handle_t's in the path before taking the next step
    std::unordered_set<std::pair<size_t, size_t>> path_step_was_continued; 

    // For a path as its id, length not including the next step, and the next step, what is the path id of 
    // the path plus next step?
    std::unordered_map<std::tuple<size_t, size_t, handlegraph::net_handle_t>, size_t> path_next_step_to_id; 

    //TODO: These two could be combined into an unordered_map<std::pair<size_t, size_t>, unordered_map<net_handle_t, size_t>> but I think that's a bad idea

    // Helper function to get the path id of the next step out from an intermediate path
    auto get_next_path_id = [&] (size_t current_path_id, size_t current_path_length, const handlegraph::net_handle_t& next_net)  {

        // Check if we've already branched from this point. If yes, then get the correct path id. Otherwise, record the new branch
        size_t path_id = current_path_id;
        bool already_branched = path_step_was_continued.count(std::make_pair(path_id, current_path_length));
        if (already_branched && path_next_step_to_id.count(std::make_tuple(path_id, current_path_length, next_net))) {
            path_id = path_next_step_to_id.at(std::make_tuple(path_id, current_path_length, next_net));
        } else {
            size_t new_path_id = already_branched ? path_count++ : path_id;
            path_next_step_to_id[std::make_tuple(path_id, current_path_length, next_net)] = new_path_id;
            path_id = new_path_id;
            path_step_was_continued.emplace(path_id, current_path_length);
        }
        return path_id;
    };



#ifdef DEBUG_PATH_PARTITIONER
    std::cerr << "Start with state " << first_state << " for node " << gbwt::Node::id(start_node)  << ":"
         << gbwt::Node::is_reverse(start_node) << std::endl;
#endif

    if (!first_path.empty()) {
        intermediate_paths.emplace_back(first_path, first_state, path_count++);
    }

    while (!intermediate_paths.empty()) {
        // For one intermediate path, add the next step


        std::tuple<std::vector<handlegraph::net_handle_t>, gbwt::SearchState, size_t> current_path = std::move(intermediate_paths.back()); 
        intermediate_paths.pop_back();
        #ifdef DEBUG_PATH_PARTITIONER
        std::cerr << "Continue path " << std::get<2>(current_path) << ":\t";
        for (const auto& net : std::get<0>(current_path)) {
            std::cerr << distance_index.net_handle_as_string(net) << ",";
        }
        std::cerr << std::endl;
        #endif


        // The next steps out from the current path, as a handle, a node in the gbwt, and a search state
        std::vector<std::tuple<handle_t, gbwt::node_type, gbwt::SearchState>> next_steps;

        // Get the net handle to the last node in the path. If it was a chain, make sure to get the boundary node leaving the chain in the right direction
        handlegraph::net_handle_t last_net = distance_index.is_node(std::get<0>(current_path).back()) 
                                          ? std::get<0>(current_path).back()
                                          : (distance_index.ends_at(std::get<0>(current_path).back()) == handlegraph::SnarlDecomposition::END 
                                               ? distance_index.get_bound(std::get<0>(current_path).back(), true, false)
                                               : distance_index.get_bound(std::get<0>(current_path).back(), false, false));

        #ifdef DEBUG_PATH_PARTITIONER
            std::cerr << "\tFrom last net " << distance_index.net_handle_as_string(last_net) << std::endl;
        #endif

        graph.follow_edges(distance_index.get_handle(last_net, &graph), false, [&](const handle_t& next) {
            // extend the last node of the thread using gbwt
            auto gbwt_next = gbwt::Node::encode(graph.get_id(next), graph.get_is_reverse(next));
            auto new_state = gbwt.extend(std::get<1>(current_path), gbwt_next);
            if (!new_state.empty()) {
                next_steps.push_back(std::make_tuple(next, gbwt_next, new_state));
            }
        });

    
        for (auto& next_step : next_steps) {
            // For each of the next steps, extend the search state and path and add it to the intermediate paths

            // Does this next step count as a branch? 
            // When the path continues in a nested chain, it doesn't count as a branch, only when there is a new node or 
            // child chain and this isn't the first new node or child chain.
            bool branch = false;

            const handlegraph::handle_t& next = std::get<0>(next_step);
            gbwt::node_type& gbwt_next = std::get<1>(next_step);
            gbwt::SearchState& new_state = std::get<2>(next_step);

            size_t current_path_length = std::get<0>(current_path).size();

            std::vector<handlegraph::net_handle_t> updated_path;
            if (&next_step == &next_steps.back()) {
                // avoid a copy by re-using the vector for the last thread. this way simple cases
                // like scanning along one path don't blow up to n^2
                updated_path = std::move(std::get<0>(current_path));
            } else {
                updated_path = std::get<0>(current_path);
            }

            // Add this node to the path, keeping track of if it is a node or a nested chain
            // For nested chains, we will add the chain itself to the path, then the node in the traversal.
            // As the path is traversed, the last thing in the path is popped and replaced with the current node
            // until the end of the chain is reached, in which case the path will end with the chain

            handlegraph::net_handle_t next_net = distance_index.get_net(next, &graph);
            // The parent is always going to be a chain, trivial or not
            handlegraph::net_handle_t next_net_parent = distance_index.get_parent(next_net);
            // if this is a node we're interested in, then the grandparent is the current snarl
            handlegraph::net_handle_t next_net_grandparent = distance_index.get_parent(next_net_parent);
            #ifdef DEBUG_PATH_PARTITIONER
                std::cerr << "\tReached next net " << distance_index.net_handle_as_string(next_net) << std::endl;
            #endif


            if (next_net == snarl_start || next_net == snarl_end) {
                // If this handle is leaving the snarl, then add the completed path to the list of completed paths
                if (next_net == end_net) {

                    updated_path.push_back(next_net); 
                    finished_paths.emplace_back(std::move(updated_path), new_state, get_next_path_id(std::get<2>(current_path), current_path_length, next_net));
                    #ifdef DEBUG_PATH_PARTITIONER
                        std::cerr << "\tFinished_path num " << std::get<2>(finished_paths.back()) << ":\t";
                        for (const auto& net : std::get<0>(finished_paths.back())) {
                            std::cerr << distance_index.net_handle_as_string(net) << ",";
                        }
                        std::cerr << std::endl;
                    #endif
                }
                branch = true;

            } else {
                // If this handle isn't leaving the snarl, then add the path back to the list of intermediate paths to continue it
                if (distance_index.start_end_traversal_of(next_net_grandparent) == distance_index.start_end_traversal_of(snarl)) {
                    // If the grandparent is the snarl we're traversing
                    if (distance_index.is_trivial_chain(next_net_parent)) {
                        // If this is a trivial chain whose parent is the snarl, add it as the node
                        branch = true;
                        #ifdef DEBUG_PATH_PARTITIONER
                            std::cerr << "\t\tnew path with node child " << distance_index.net_handle_as_string(next_net) << std::endl;
                        #endif
                        updated_path.push_back(next_net);
                        // TODO: Could also get the sequence here
                    } else {
                        // Otherwise, this is a node in a child of the chain
                        // If this is the start, then add the chain then the node
                        // If this is the end, then take out the last thing in the path so that the last thing becomes the chain
                        // Otherwise, do as for other nested nodes and replace the last thing with the node
                        handlegraph::net_handle_t chain_start = distance_index.get_bound(next_net_parent, false, true);
                        handlegraph::net_handle_t chain_end = distance_index.get_bound(next_net_parent, true, true);
                        if (next_net == chain_start || next_net == chain_end) {
                            branch = true;
                            #ifdef DEBUG_PATH_PARTITIONER
                                std::cerr << "\t\tnew path with chain child" << std::endl;
                            #endif
                            updated_path.push_back(next_net_parent);
                            updated_path.push_back(next_net);
                        } else if (next_net == distance_index.flip(chain_start) || next_net == distance_index.flip(chain_end)) {
                            #ifdef DEBUG_PATH_PARTITIONER
                                std::cerr << "\t\t finish chain child" << std::endl;
                            #endif
                            updated_path.pop_back();
                        } else {
                            #ifdef DEBUG_PATH_PARTITIONER
                                std::cerr << "\t\t continue chain child" << std::endl;
                            #endif
                            updated_path.pop_back();
                            updated_path.push_back(next_net);
                        }

                    }
                } else {
                    // Otherwise, this is nested and we need to replace the last thing in the path with this node
                    #ifdef DEBUG_PATH_PARTITIONER
                        std::cerr << "\t\tnested node, replace last node in path" << std::endl;
                    #endif
                    updated_path.pop_back();
                    updated_path.push_back(next_net);
                }
                intermediate_paths.emplace_back(std::move(updated_path), new_state, 
                                                branch ? get_next_path_id(std::get<2>(current_path), current_path_length, next_net) : std::get<2>(current_path));
            }
        }
    } // End while loop going through intermediate paths
    #ifdef DEBUG_PATH_PARTITIONER
    std::cerr << "Found " << finished_paths.size() << " threads through " << distance_index.net_handle_as_string(snarl) << std::endl;
    std::cerr << "\tthere are " << path_count << " distinct walks" << std::endl;
    for (const auto& current_path : finished_paths) {
        std::cerr << "\t" << std::get<2>(current_path) << ":";
        for (const auto& net : std::get<0>(current_path)) {
            std::cerr << distance_index.net_handle_as_string(net) << ",";
        }
        std::cerr << std::endl;
    }
    #endif

    return path_count;
}

std::vector<size_t> partition_embedded_paths_in_snarl_with_gbwt(const handlegraph::PathPositionHandleGraph& graph, 
                          const gbwt::GBWT& gbwt, const bdsg::SnarlDistanceIndex& distance_index,
                          const net_handle_t& snarl,
                          const std::vector<stoat::sample_hap_t>& all_sample_haplotypes,
                          std::vector<PathTraversal>& paths_per_allele) {



    // The final boundary-to-boundary paths
    // The size_t is an identifier for the path, since some paths may be identical
    std::vector<std::tuple<std::vector<handlegraph::net_handle_t>, gbwt::SearchState, size_t>> finished_paths;

    // Get all start-end or start-start traversals and put them in  finished_paths 
    // Path count will be an upper bound of the number of distinct paths
    size_t path_count = get_gbwt_traversals(graph, gbwt, distance_index, snarl, finished_paths);

    // At this point, finished paths holds a path for each distinct, haplotype-supported walk through th esnarl netgraph
    // Since a sample may go through the snarl multiple times, we now want to partition the samples based on how many
    // of each walk they take.

    // Map each sample(plus haplotype) to its index in all_sample_haplotypes
    std::map<stoat::sample_hap_t, std::size_t> sample_to_index;
    for (size_t i = 0 ; i < all_sample_haplotypes.size() ; i++) {
        sample_to_index[all_sample_haplotypes[i]] = i;
    }

    // For each sample, remember which paths through the snarl it takes. Since it may take multiple, save this as a linked list.
    // Since it doesn't happen very often, there will be one vector of all the heads of the linked list, and a separate vector of
    // the subsequent nodes/repeats.

    // A struct representing the allele/path taken by a sample.
    // This becomes a linked list when a sample traverses the snarl multiple times (this doesn't happen often) 
    // A value of {max(), max()} indicates that no path was taken
    struct path_link_t {
        size_t path_id;           // The path
        size_t additional_edge;  // Index into a list of additional edges

        path_link_t() : path_id(std::numeric_limits<size_t>::max()), additional_edge(std::numeric_limits<size_t>::max()){};
        path_link_t(size_t path_id) : 
            path_id(path_id), additional_edge(std::numeric_limits<size_t>::max()){};
    };
    std::vector<path_link_t> path_by_sample (all_sample_haplotypes.size());
    std::vector<path_link_t> additional_edges;

    auto add_path_to_sample = [&](size_t sample_id, size_t path_id) {
        if (path_by_sample.at(sample_id).path_id == std::numeric_limits<size_t>::max()) {
            path_by_sample.at(sample_id).path_id = path_id;
        } else {
            // Follow the linked list to the end and then add this
            size_t next_index = path_by_sample.at(sample_id).additional_edge;
            if (next_index == std::numeric_limits<size_t>::max()) {
                path_by_sample.at(sample_id).additional_edge = additional_edges.size();
                additional_edges.emplace_back(path_id);
            } else {
                while (additional_edges.at(next_index).additional_edge != std::numeric_limits<size_t>::max()) {
                    next_index = additional_edges.at(next_index).additional_edge;
                }
                additional_edges.at(next_index).additional_edge = additional_edges.size();
                additional_edges.emplace_back(path_id);
            }
        }
    };

    // Remember the paths themselves, once per path_id
    std::vector<std::vector<handlegraph::net_handle_t>> paths_per_path_id (path_count);

    // For each path (along all_sample_haplotypes), the index of the set it is currently in
    // Everything starts out in the same set, 0, which represents not being in the snarl
    std::vector<size_t> old_sets (all_sample_haplotypes.size(), 0);
    size_t old_set_count = 1;

    // This will be the count of the current allele for each haplotype
    std::vector<size_t> intermediate_sets (all_sample_haplotypes.size(), 0);

    std::vector<size_t> new_sets (all_sample_haplotypes.size(), 0);
    size_t new_set_count = 1;

    // Sort the finished paths so by path id so that identical paths are consecutive in the vector.
    // Since actually sorting the vector would be slow, make a vector of indices and sort that
    std::vector<size_t> sort_order(finished_paths.size(), 0);
    for (size_t i = 0 ; i < sort_order.size() ; i++) {
        sort_order[i] = i;
    }
    std::sort(sort_order.begin(), sort_order.end(), [&](size_t a, size_t b) { return std::get<2>(finished_paths.at(a)) < std::get<2>(finished_paths.at(b)); });

    auto gbwt_reference_samples = gbwtgraph::parse_reference_samples_tag(gbwt);


    // Go through all finished paths by path id
    for (size_t sorted_i = 0 ; sorted_i < sort_order.size() ; sorted_i++) {
        const std::tuple<std::vector<handlegraph::net_handle_t>, gbwt::SearchState, size_t>& current_state = finished_paths.at(sort_order.at(sorted_i));
        #ifdef DEBUG_PATH_PARTITIONER
            std::cerr << "At path id " << std::get<2>(current_state) << ":\t";
            for (const auto& net : std::get<0>(current_state)) {
                std::cerr << distance_index.net_handle_as_string(net) << ",";
            }
            std::cerr << std::endl;
        #endif


        //locate() finds the path identifiers for the search state
        std::vector<gbwt::size_type> path_ids = gbwt.locate(std::get<1>(current_state));
        std::cerr << "Found " << path_ids.size() << " path ids for this path " << std::endl;
        for (const gbwt::size_type id : path_ids) {
            gbwt::size_type gbwt_path_id = gbwt::Path::id(id);

            handlegraph::PathSense sense = gbwtgraph::get_path_sense(gbwt, gbwt_path_id, gbwt_reference_samples);

            std::string path_name = handlegraph::PathMetadata::create_path_name(sense,
                                        gbwtgraph::get_path_sample_name(gbwt, gbwt_path_id, sense),
                                        gbwtgraph::get_path_locus_name(gbwt, gbwt_path_id, sense),
                                        gbwtgraph::get_path_haplotype(gbwt, gbwt_path_id, sense),
                                        gbwtgraph::get_path_phase_block(gbwt, gbwt_path_id, sense),
                                        gbwtgraph::get_path_subrange(gbwt, gbwt_path_id, sense)); 
            #ifdef DEBUG_PATH_PARTITIONER
                std::cerr << "\tpath " << path_name << " takes this path" << std::endl;
            #endif

            intermediate_sets.at(sample_to_index.at(sample_hap_t(path_name)))++;

            add_path_to_sample(sample_to_index.at(sample_hap_t(path_name)), std::get<2>(current_state));
        }

        if (sorted_i == sort_order.size()-1 || std::get<2>(current_state) != std::get<2>(finished_paths.at(sort_order.at(sorted_i+1)))) {
            // If we finished going through the last set of paths for this path id
            std::map<std::pair<size_t, size_t>, size_t> old_to_new_set;
            old_to_new_set[std::make_pair(0,0)] = 0;

            for (size_t sample_i = 0 ; sample_i < old_sets.size() ; sample_i++) {
                if (old_to_new_set.count(std::make_pair(old_sets.at(sample_i), intermediate_sets.at(sample_i))) == 0) {
                    new_sets[sample_i] = new_set_count;
                    old_to_new_set[std::make_pair(old_sets.at(sample_i), intermediate_sets.at(sample_i))] = new_set_count++; 
                } else {
                    new_sets.at(sample_i) = old_to_new_set.at(std::make_pair(old_sets.at(sample_i), intermediate_sets.at(sample_i))); 
                }
            }
            old_sets = std::move(new_sets);
            old_set_count = new_set_count;
            new_sets.assign(all_sample_haplotypes.size(), 0);
            intermediate_sets.assign(all_sample_haplotypes.size(), 0);
            new_set_count = 1;

            // Now deal with the walks
            paths_per_path_id.at(std::get<2>(current_state)) = std::move(std::get<0>(current_state));
        }
    }

    // The old_set_count originally counted the set with nothing, so decrement it
    old_set_count--;

    // For each allele, get one sample with this allele. Used to get the paths
    std::vector<size_t> sample_per_allele (old_set_count, std::numeric_limits<size_t>::max());
    paths_per_allele.clear();

    // Now go through the sets and change 0 to inf, and decrement all others
    for (size_t i = 0 ; i < old_sets.size() ; i++) {
        if (old_sets[i] == 0) {
            old_sets[i] = std::numeric_limits<size_t>::max();
        } else {
            --old_sets[i];
            if (sample_per_allele.at(old_sets[i]) == std::numeric_limits<size_t>::max()) {
                sample_per_allele.at(old_sets[i]) = i;
            }
        }
    }

    for (size_t allele_num = 0 ; allele_num < old_set_count ; allele_num++) {

        paths_per_allele.emplace_back();

        size_t sample_num = sample_per_allele.at(allele_num);
        size_t path_id = path_by_sample.at(sample_num).path_id;
        size_t next_edge = path_by_sample.at(sample_num).additional_edge;
        assert (path_id != std::numeric_limits<size_t>::max());

        const std::vector<handlegraph::net_handle_t>& path_as_net_handles = paths_per_path_id.at(path_id);
        
        // Add the first path through the snarl
        for (size_t i = 0 ; i < path_as_net_handles.size() ; i++) {
            const handlegraph::net_handle_t& net = path_as_net_handles.at(i);
            paths_per_allele.back().add_net_handle(net, distance_index, i == 0 || i == path_as_net_handles.size());
        }

        // If there are additional paths through the snarl, add them too
        while (next_edge != std::numeric_limits<size_t>::max()) {
            paths_per_allele.back().add_out_of_snarl_walk();

            path_id = additional_edges.at(next_edge).path_id;
            next_edge = additional_edges.at(next_edge).additional_edge;

            const std::vector<handlegraph::net_handle_t>& next_path_as_net_handles = paths_per_path_id.at(path_id);
            for (size_t i = 0 ; i < next_path_as_net_handles.size() ; i++) {
                const handlegraph::net_handle_t& net = next_path_as_net_handles.at(i);
                paths_per_allele.back().add_net_handle(net, distance_index, i == 0 || i == next_path_as_net_handles.size());
            }
        }
    }

    return old_sets;
}

}
