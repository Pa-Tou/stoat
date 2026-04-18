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
    for (size_t i = 0 ; i < all_sample_haplotypes.size() ; i++) {
        if (old_sets[i] != 0) {
            sample_sets_by_allele[old_sets[i]-1].emplace(all_sample_haplotypes[i]);
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


void get_gbwt_traversals(const handlegraph::PathPositionHandleGraph& graph, const gbwt::GBWT& gbwt, const bdsg::SnarlDistanceIndex& distance_index,     
                         const handlegraph::net_handle_t& snarl,
                         std::vector<std::tuple<std::vector<handlegraph::net_handle_t>, gbwt::SearchState, size_t>>& finished_paths, handlegraph::net_handle_t start_net, 
                         size_t& path_count, bool only_loops) {
    #ifdef DEBUG_PATH_PARTITIONER
    std::cerr << "Get threads through snarl " << distance_index.net_handle_as_string(snarl) << std::endl;
    #endif


    // Get the traversals through the snarl from the gbwt
    // This is heavily based on vg/haplotype_extracter.cpp

    handlegraph::net_handle_t end_net = distance_index.flip(start_net);
    handlegraph::handle_t start_in = distance_index.get_handle(start_net, &graph);

    // Look up the start node in GBWT and start a path
    gbwt::node_type start_node = gbwt::Node::encode(graph.get_id(start_in), graph.get_is_reverse(start_in));
    std::vector<handlegraph::net_handle_t> first_path = {distance_index.get_node_from_sentinel(start_net)};
    gbwt::SearchState first_state = gbwt.find(start_node);

    // The list of intermediate paths and the gbwt::SearchState they end on
    // The search state encompasses the haplotypes that followed the path
    // The size_t is an identifier for the path, since some paths may be identical
    std::vector<std::tuple<std::vector<handlegraph::net_handle_t>, gbwt::SearchState, size_t>> intermediate_paths;

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
        std::cerr << "\tContinue with net handle " << distance_index.net_handle_as_string(std::get<0>(current_path).back()) << std::endl;
        #endif


        // The next steps out from the current path, as a handle, a node in the gbwt, and a search state
        std::vector<std::tuple<handle_t, gbwt::node_type, gbwt::SearchState>> next_steps;
        // Get the net handle to the last node in the path. If it was a chain, make sure to get the boundary node leaving the chain in the right direction
        // TODO: I'm actually not sure if we can skip through the gbwt like this
        handlegraph::net_handle_t last_net = distance_index.is_node(std::get<0>(current_path).back()) 
                                          ? std::get<0>(current_path).back()
                                          : (distance_index.ends_at(std::get<0>(current_path).back()) == handlegraph::SnarlDecomposition::END 
                                               ? distance_index.get_bound(std::get<0>(current_path).back(), true, false)
                                               : distance_index.get_bound(std::get<0>(current_path).back(), false, false));

        std::cerr << "From last net " << distance_index.net_handle_as_string(last_net) << std::endl;

        graph.follow_edges(distance_index.get_handle(last_net, &graph), false, [&](const handle_t& next) {
            // extend the last node of the thread using gbwt
            auto gbwt_next = gbwt::Node::encode(graph.get_id(next), graph.get_is_reverse(next));
            auto new_state = gbwt.extend(std::get<1>(current_path), gbwt_next);
#ifdef DEBUG_PATH_PARTITIONER
              std::cerr << "Extend state " << std::get<1>(current_path) << " to " << new_state << " with " << gbwt::Node::id(gbwt_next) << std::endl;
#endif
            if (!new_state.empty()) {
                next_steps.push_back(std::make_tuple(next, gbwt_next, new_state));
            }
        });

    
        // We need to assign an id to each new path. One path can keep the same id as the current path but each
        // other branching path needs a new one. When the path continues in a nested chain, it doesn't count as 
        // a branch so it can keep the same id. So only use a new id when there is a new node or child chain
        // and this isn't the first new node or child chain.
        bool first_branch = true;
        for (auto& next_step : next_steps) {
            // For each of the next steps, extend the search state and path and add it to the intermediate paths
            bool branched = false;

            const handlegraph::handle_t& next = std::get<0>(next_step);
            gbwt::node_type& gbwt_next = std::get<1>(next_step);
            gbwt::SearchState& new_state = std::get<2>(next_step);

            std::vector<handlegraph::net_handle_t> updated_path;
            if (&next_step == &next_steps.back()) {
                // avoid a copy by re-using the vector for the last thread. this way simple cases
                // like scanning along one path don't blow up to n^2
                updated_path = std::move(std::get<0>(current_path));
            } else {
                // TODO: idk why we don't move in both cases 
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


            // The bounds leaving the snarl
            handlegraph::net_handle_t snarl_start = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, false));
            handlegraph::net_handle_t snarl_end = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, true, false));
             if (next_net == snarl_start || next_net == snarl_end) {
                // If this handle is leaving the snarl, then add the completed path to the list of completed paths
                // If this handle isn't leaving the snarl, then add the path back to the list of intermediate paths to continue it
                // If this is leaving the snarl
                if (!only_loops || next_net == end_net) {
                    updated_path.push_back(next_net);
                    finished_paths.emplace_back(std::move(updated_path), new_state, std::get<2>(current_path));
                    #ifdef DEBUG_PATH_PARTITIONER
                    std::cerr << "Finished_path num " << std::get<2>(current_path) << ":\t";
                    for (const auto& net : std::get<0>(finished_paths.back())) {
                        std::cerr << distance_index.net_handle_as_string(net) << ",";
                    }
                    std::cerr << std::endl;
                    #endif
                }

            } else {
                if (distance_index.start_end_traversal_of(next_net_grandparent) == distance_index.start_end_traversal_of(snarl)) {
                    // If the grandparent is the snarl we're traversing
                    if (distance_index.is_trivial_chain(next_net_parent)) {
                        // If this is a trivial chain whose parent is the snarl, add it as the node
                        std::cerr << "\t\tnew path with node child " << distance_index.net_handle_as_string(next_net) << std::endl;
                        updated_path.push_back(next_net);
                        branched = true;
                        // TODO: Could also get the sequence here
                    } else {
                        // Otherwise, this is a node in a child of the chain
                        // If this is the start, then add the chain then the node
                        // If this is the end, then take out the last thing in the path so that the last thing becomes the chain
                        // Otherwise, do as for other nested nodes and replace the last thing with the node
                        handlegraph::net_handle_t chain_start = distance_index.get_bound(next_net_parent, false, true);
                        handlegraph::net_handle_t chain_end = distance_index.get_bound(next_net_parent, true, true);
                        if (next_net == chain_start || next_net == chain_end) {
                        std::cerr << "\t\tnew path with chain child" << std::endl;
                            updated_path.push_back(next_net_parent);
                            updated_path.push_back(next_net);
                            branched = true;
                        } else if (next_net == distance_index.flip(chain_start) || next_net == distance_index.flip(chain_end)) {
                        std::cerr << "\t\t finish chain child" << std::endl;
                            updated_path.pop_back();
                        } else {
                        std::cerr << "\t\t continue chain child" << std::endl;
                            updated_path.pop_back();
                            updated_path.push_back(next_net);
                        }

                    }
                } else {
                    // Otherwise, this is nested and we need to replace the last thing in the path with this node
                    std::cerr << "\t\tnested node, replace last node in path" << std::endl;
                    updated_path.pop_back();
                    updated_path.push_back(next_net);
                }
                if (!first_branch && branched) {
                    // If this is a new branch and not the first branch, then we use a new id for the path
                    intermediate_paths.emplace_back(std::move(updated_path), new_state, path_count++);
                    #ifdef DEBUG_PATH_PARTITIONER
                        std::cerr << "\t\tNew path num" << (path_count-1) << std::endl;
                    #endif
                } else {
                    // Otherwise, use the same id
                    intermediate_paths.emplace_back(std::move(updated_path), new_state, std::get<2>(current_path));
                }
                if (branched) {
                    first_branch = false;
                }
            }
        }
    } // End while loop going through intermediate paths
    #ifdef DEBUG_PATH_PARTITIONER
    std::cerr << "Found " << finished_paths.size() << " threads through " << distance_index.net_handle_as_string(snarl) << std::endl;
    #endif

}

std::vector<size_t> partition_embedded_paths_in_snarl_with_gbwt(const handlegraph::PathPositionHandleGraph& graph, 
                          const gbwt::GBWT& gbwt, const bdsg::SnarlDistanceIndex& distance_index,
                          const net_handle_t& snarl,
                          const std::vector<stoat::sample_hap_t>& all_sample_haplotypes,
                          std::vector<PathTraversal>& paths_per_allele) {

    // Get the bounds of the snarl
    handlegraph::net_handle_t start_in = distance_index.get_bound(snarl, false, true);
    handlegraph::net_handle_t end_in = distance_index.get_bound(snarl, true, true);


    // The final boundary-to-boundary paths
    // The size_t is an identifier for the path, since some paths may be identical
    std::vector<std::tuple<std::vector<handlegraph::net_handle_t>, gbwt::SearchState, size_t>> finished_paths;


    // Get all traversals and put them in  finished_paths 
    // Path count will be the number of distinct paths
    size_t path_count = 0;
    get_gbwt_traversals(graph, gbwt, distance_index, snarl, finished_paths, start_in, path_count, false);

    get_gbwt_traversals(graph, gbwt, distance_index, snarl, finished_paths, end_in, path_count, true);



    // At this point, finished_paths holds a path for each distinct, haplotype-supported path going through the snarl
    // Now, we need to assign haplotypes to each of these paths
    std::unordered_map<stoat::sample_hap_t, size_t> sample_to_allele;

    auto gbwt_reference_samples = gbwtgraph::parse_reference_samples_tag(gbwt);

    // For each distinct path, have we found it yet?
    // When we find it, move into unique_paths
    std::vector<bool> found_path (path_count, false);
    std::cerr << "Found " << path_count << " unique paths with " << finished_paths.size() << " total paths" << std::endl;
    std::vector<std::vector<handlegraph::net_handle_t>> unique_paths(path_count);

    for (const std::tuple<std::vector<handlegraph::net_handle_t>, gbwt::SearchState, size_t>& state : finished_paths) {

        size_t path_id = std::get<2>(state);
        std::cerr << "At path id " << path_id << ":\t";
        for (const auto& net : std::get<0>(state)) {
            std::cerr << distance_index.net_handle_as_string(net) << ",";
        }
        std::cerr << std::endl;

        //locate() finds the path identifiers for the search state
        std::vector<gbwt::size_type> path_ids = gbwt.locate(std::get<1>(state));
        for (const gbwt::size_type id : path_ids) {
            gbwt::size_type path_id = gbwt::Path::id(id);

            handlegraph::PathSense sense = gbwtgraph::get_path_sense(gbwt, path_id, gbwt_reference_samples);

            std::string path_name = handlegraph::PathMetadata::create_path_name(sense,
                                        gbwtgraph::get_path_sample_name(gbwt, path_id, sense),
                                        gbwtgraph::get_path_locus_name(gbwt, path_id, sense),
                                        gbwtgraph::get_path_haplotype(gbwt, path_id, sense),
                                        gbwtgraph::get_path_phase_block(gbwt, path_id, sense),
                                        gbwtgraph::get_path_subrange(gbwt, path_id, sense)); 

            sample_to_allele.emplace(sample_hap_t(path_name), 
                                     path_id);
        }

        //TODO: MOving these paths around isn't very efficient. Would probably be better to make found_path be a vector of indices into finished_paths,
        // then build the vector PathTraversals one at a time in order
        if (!found_path.at(path_id)) {
            found_path.at(path_id) = true;
            unique_paths.at(path_id) = std::get<0>(state);
        }
    }
    

    // Now get the allele assignments to return
    std::vector<size_t> allele_assignments;
    allele_assignments.reserve(all_sample_haplotypes.size());
    for (const stoat::sample_hap_t& samp : all_sample_haplotypes) {
        if (sample_to_allele.count(samp) == 0) {
            allele_assignments.emplace_back(std::numeric_limits<size_t>::max());

        } else {
            allele_assignments.emplace_back(sample_to_allele.at(samp));
        }
    }

    // Get the paths in the right format
    paths_per_allele = stoat::convert_path_traversals(distance_index, graph, unique_paths);
    return allele_assignments;
}

}
