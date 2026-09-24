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
    std::cerr << "Get all possible walks through snarl " << distance_index.net_handle_as_string(snarl) << std::endl;
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
    // Validate paths
    std::cerr << "Found " << walks.size() << " paths through the snarl" << std::endl;
    for (const auto& walk : walks) {
        
        std::cerr << "\t" << walk.to_string() << std::endl;
    }
#endif

    return;
}
}// end namespace stoat
