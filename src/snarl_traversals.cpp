#include "snarl_traversals.hpp"
#include <fstream>
#include <filesystem>
#include "utils.hpp"

//#define DEBUG_SNARL_TRAVERSALS

namespace stoat {

void get_all_walks_through_snarl(
        const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
        const net_handle_t& snarl, std::vector<stoat::PathTraversal>& walks, 
        size_t walk_cycle_limit, size_t walk_steps_limit) {

#ifdef DEBUG_SNARL_TRAVERSALS
#pragma omp critical(cerr)
{
    std::cerr << "Get all possible walks through snarl " << distance_index.net_handle_as_string(snarl) << std::endl;
}
#endif

    // Path exploration
    std::vector<std::vector<handlegraph::net_handle_t>> paths = {
        {distance_index.get_bound(snarl, false, true)}
    };
    
    std::vector<std::vector<handlegraph::net_handle_t>> walks_as_net_handles;
    
    bool break_snarl = false;
    
    // For each incomplete path in paths, walk out from the end and add a copy of the path plus each next step to paths
    // Do this until the path reaches the end
    while (!paths.empty()) {
        std::vector<handlegraph::net_handle_t> path = std::move(paths.back());
        paths.pop_back();
    
        std::unordered_map<handlegraph::net_handle_t, size_t> dict_path_occ;
        bool cycle = false;
    
        // TODO: Put this back
        for (const auto& net : path) {
            if (++dict_path_occ[net] > walk_cycle_limit + 1) {
                cycle = true;
                break;
            }
        }
    
        // TODO: Add out_fail
        // TODO: Get the child count properly
        if (path.size() > walk_steps_limit) {
           // #pragma omp critical(out_fail)
           // out_fail << distance_index.node_id(distance_index.get_bound(snarl, false, true)) << "_" << distance_index.node_id(distance_index.get_bound(snarl, true, true)) 
           //          << "\titeration_calculation_out = " << std::endl;// << children << " children\n";
           break_snarl = true;
           break;
        }
	
        // Follow edges from the last element in path
        if (!path.empty()) {
            distance_index.follow_net_edges(path.back(), &graph, false, [&](const handlegraph::net_handle_t& next_child) {
                // If this is the bound of the snarl then we're done
                if (distance_index.is_sentinel(next_child)) {

                    size_t next_child_node_id = distance_index.node_id(distance_index.get_node_from_sentinel(next_child));
                    size_t first_element_path_node_id = distance_index.node_id(distance_index.get_node_from_sentinel(path[0]));

                    // Only keep the walk if it entered and exited the snarl at opposite sides
                    if (next_child_node_id != first_element_path_node_id) {
                        walks_as_net_handles.emplace_back(path);
                        walks_as_net_handles.back().push_back(next_child);
                    }
            
                } else {
                    //TODO: Look for the cycle sooner
            
                    if (cycle) { // Case where we find a loop
                        return false;
                    }
                    paths.emplace_back(path);
                    paths.back().push_back(next_child);
                }
                return true;
            });
        }
    }
    
    if (break_snarl) {
        walks_as_net_handles.clear();
    }

#ifdef DEBUG_SNARL_TRAVERSALS
    // Validate paths
    std::set<std::vector<handlegraph::net_handle_t>> found_walks;
    for (const auto& walk : walks_as_net_handles) {
        for (size_t i = 0 ; i < walk.size() - 1 ; i++) {
            assert(distance_index.distance_in_parent(snarl, walk[i], distance_index.flip(walk[i+1])) == 0);
        }
        assert(found_walks.count(walk) == 0);
        assert(walk.front() == distance_index.get_bound(snarl, false, true));
        assert(walk.back() == distance_index.get_bound(snarl, true, false));
        found_walks.insert(walk);
    }
#endif

    walks = stoat::net_handles_to_path_traversals(distance_index, graph, walks_as_net_handles);  
 
#ifdef DEBUG_SNARL_TRAVERSALS
#pragma omp critical(cerr)
{
    std::cerr << "Found " << walks.size() << " paths through the snarl" << std::endl;
    for (const auto& walk : walks) {
        
        std::cerr << "\t" << walk.to_string() << std::endl;
    }
}
#endif

    return;
}

void get_haplotype_walks_through_snarl(
        const handlegraph::PathPositionHandleGraph& graph, const gbwt::GBWT& gbwt,
        const gbwt::FastLocate& r_index,const bdsg::SnarlDistanceIndex& distance_index,
        const net_handle_t& snarl, std::vector<stoat::PathTraversal>& walks) {

    // Get all distinct walks with some possible duplicates 
    std::vector<gbwt_path_t> gbwt_walks;
    size_t max_path_count = get_gbwt_traversals(graph, gbwt, r_index, distance_index, snarl, gbwt_walks);

    // Copy the traversals into walks, making sure not to duplicate 
    std::vector<bool> got_path (max_path_count, false);
    walks.clear();
    for (gbwt_path_t& gbwt_walk : gbwt_walks) {
        if (!got_path[gbwt_walk.identifier]) {
            got_path[gbwt_walk.identifier] = true;
            walks.emplace_back(net_handles_to_path_traversal(distance_index, graph, gbwt_walk.path));
        }
    }
}

// Get the traversals through the snarl from the gbwt, filling in finished_paths
// This is heavily based on vg/haplotype_extracter.cpp
size_t get_gbwt_traversals(const handlegraph::PathPositionHandleGraph& graph, const gbwt::GBWT& gbwt,
                           const gbwt::FastLocate& r_index,
                           const bdsg::SnarlDistanceIndex& distance_index,     
                           const handlegraph::net_handle_t& snarl,
                           std::vector<gbwt_path_t>& finished_paths) {
#ifdef DEBUG_SNARL_TRAVERSALS
#pragma omp critical(cerr)
{
    std::cerr << "Get threads through snarl " << distance_index.net_handle_as_string(snarl) << std::endl;
}
#endif

    size_t path_count = 0;


    // Get the bounds of the snarl: start facing in and end facing out. It doesn't matter which one is which since we get all start-end traversals
    handlegraph::net_handle_t start_net = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, true));
    handlegraph::net_handle_t end_net = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, true, false));

    handlegraph::handle_t start_in = distance_index.get_handle(start_net, &graph);


    // The bounds leaving the snarl
    handlegraph::net_handle_t snarl_start = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, false));
    handlegraph::net_handle_t snarl_end = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, true, false));

    // The list of intermediate paths and the gbwt::SearchState they end on
    // The search state encompasses the haplotypes that followed the path
    std::vector<gbwt_path_t> intermediate_paths;

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


    // Define a struct representing one step in a walk through the gbwt
    // The same as gbwt_path_t but with one step instead of a path
    struct gbwt_step_t {
        handle_t node_handle;
        gbwt::SearchState search_state;
        gbwt::size_type sa_index;

        // destructively construct the struct by move()ing the search state
        gbwt_step_t(handle_t node_handle, gbwt::SearchState search_state, gbwt::size_type sa_index) :
             node_handle(node_handle), search_state(std::move(search_state)), sa_index(sa_index) {}
    };





    // Look up the start node in GBWT and start a path
    gbwt::node_type start_node = gbwt::Node::encode(graph.get_id(start_in), graph.get_is_reverse(start_in));
    std::vector<handlegraph::net_handle_t> first_path = {start_net};

    // When using the r-index, we keep track of the first occurrence in the suffix array of our range to get back to the location in the original text
    gbwt::size_type first_sa_index;
    gbwt::SearchState first_state = r_index.find(start_node, first_sa_index);

    intermediate_paths.emplace_back(first_path, first_state, path_count++, first_sa_index);

#ifdef DEBUG_SNARL_TRAVERSALS
#pragma omp critical(cerr)
{
    std::cerr << "Start with state " << first_state << " for node " << gbwt::Node::id(start_node)  << ":"
         << gbwt::Node::is_reverse(start_node) << " with sa index " << first_sa_index << std::endl;
}
#endif

    while (!intermediate_paths.empty()) {
        // For one intermediate path, add the next step


        gbwt_path_t current_path = std::move(intermediate_paths.back()); 
        intermediate_paths.pop_back();
#ifdef DEBUG_SNARL_TRAVERSALS
#pragma omp critical(cerr)
{
        std::cerr << "Continue path " << current_path.identifier << ":\t";
        for (const auto& net : current_path.path) {
            std::cerr << distance_index.net_handle_as_string(net) << ",";
        }
        std::cerr << std::endl;
}
#endif


        // The next steps out from the current path, as a handle, a node in the gbwt, a search state, and the suffix array offset for the range in the search state
        std::vector<gbwt_step_t> next_steps;

        // Get the net handle to the last node in the path. If it was a chain, make sure to get the boundary node leaving the chain in the right direction
        handlegraph::net_handle_t last_net = distance_index.is_node(current_path.path.back()) 
                                          ? current_path.path.back()
                                          : (distance_index.ends_at(current_path.path.back()) == handlegraph::SnarlDecomposition::END 
                                               ? distance_index.get_bound(current_path.path.back(), true, false)
                                               : distance_index.get_bound(current_path.path.back(), false, false));

#ifdef DEBUG_SNARL_TRAVERSALS
#pragma omp critical(cerr)
{
            std::cerr << "\tFrom last net " << distance_index.net_handle_as_string(last_net) << std::endl;
}
#endif

        graph.follow_edges(distance_index.get_handle(last_net, &graph), false, [&](const handle_t& next) {
            // extend the last node of the thread using gbwt
            auto gbwt_next = gbwt::Node::encode(graph.get_id(next), graph.get_is_reverse(next));
            gbwt::size_type next_sa_index = current_path.sa_index;
            auto new_state = r_index.extend(current_path.search_state, gbwt_next, next_sa_index);
            if (!new_state.empty()) {
                next_steps.emplace_back(next, new_state, next_sa_index);
            }
        });

    
        for (auto& next_step : next_steps) {
            // For each of the next steps, extend the search state and path and add it to the intermediate paths

            // Does this next step count as a branch? 
            // When the path continues in a nested chain, it doesn't count as a branch, only when there is a new node or 
            // child chain and this isn't the first new node or child chain.
            bool branch = false;

            size_t current_path_length = current_path.path.size();

            std::vector<handlegraph::net_handle_t> updated_path;
            if (&next_step == &next_steps.back()) {
                // avoid a copy by re-using the vector for the last thread. this way simple cases
                // like scanning along one path don't blow up to n^2
                updated_path = std::move(current_path.path);
            } else {
                updated_path = current_path.path;
            }

            // Add this node to the path, keeping track of if it is a node or a nested chain
            // For nested chains, we will add the chain itself to the path, then the node in the traversal.
            // As the path is traversed, the last thing in the path is popped and replaced with the current node
            // until the end of the chain is reached, in which case the path will end with the chain

            handlegraph::net_handle_t next_net = distance_index.get_net(next_step.node_handle, &graph);
            // The parent is always going to be a chain, trivial or not
            handlegraph::net_handle_t next_net_parent = distance_index.get_parent(next_net);
            // if this is a node we're interested in, then the grandparent is the current snarl
            handlegraph::net_handle_t next_net_grandparent = distance_index.get_parent(next_net_parent);
#ifdef DEBUG_SNARL_TRAVERSALS
#pragma omp critical(cerr)
{
                std::cerr << "\tReached next net " << distance_index.net_handle_as_string(next_net) << std::endl;
}
#endif


            if (next_net == snarl_start || next_net == snarl_end) {
                // If this handle is leaving the snarl, then add the completed path to the list of completed paths
                if (next_net == end_net) {

                    updated_path.push_back(next_net); 
                    finished_paths.emplace_back(updated_path, next_step.search_state, 
                                                get_next_path_id(current_path.identifier, current_path_length, next_net),
                                                next_step.sa_index);
#ifdef DEBUG_SNARL_TRAVERSALS
#pragma omp critical(cerr)
{
                        std::cerr << "\tFinished_path num " << finished_paths.back().identifier << ":\t";
                        for (const auto& net : finished_paths.back().path) {
                            std::cerr << distance_index.net_handle_as_string(net) << ",";
                        }
                        std::cerr << std::endl;
}
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
#ifdef DEBUG_SNARL_TRAVERSALS
#pragma omp critical(cerr)
{
                            std::cerr << "\t\tnew path with node child " << distance_index.net_handle_as_string(next_net) << std::endl;
}
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
                            // If this is going into the child chain, then add the chain and then this node
                            branch = true;
                            updated_path.push_back(next_net_parent);
                            updated_path.push_back(next_net);
                        } else if (next_net == distance_index.flip(chain_start) || next_net == distance_index.flip(chain_end)) {
                            // If this is leaving the child chain, then just pop the extra node so that the path finishes on the chain
                            updated_path.pop_back();
                        } else {
                            // Otherwise, this is continuing the walk in the chain so just replace the last traversal in the path
                            updated_path.pop_back();
                            updated_path.push_back(next_net);
                        }

                    }
                } else {
                    // Otherwise, this is nested and we need to replace the last thing in the path with this node
                    updated_path.pop_back();
                    updated_path.push_back(next_net);
                }
                intermediate_paths.emplace_back(std::move(updated_path), next_step.search_state, 
                                                branch ? get_next_path_id(current_path.identifier, current_path_length, next_net) : current_path.identifier,
                                                next_step.sa_index);
            }
        }
    } // End while loop going through intermediate paths
#ifdef DEBUG_SNARL_TRAVERSALS
#pragma omp critical(cerr)
{
    std::cerr << "Found " << finished_paths.size() << " threads through " << distance_index.net_handle_as_string(snarl) << std::endl;
    std::cerr << "\tthere are " << path_count << " distinct walks" << std::endl;
    for (const auto& current_path : finished_paths) {
        std::cerr << "\t" << current_path.identifier << ":";
        for (const auto& net : current_path.path) {
            std::cerr << distance_index.net_handle_as_string(net) << ",";
        }
        std::cerr << std::endl;
    }
}
#endif

    return path_count;
}


}// end namespace stoat


