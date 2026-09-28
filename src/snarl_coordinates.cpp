#include "snral_coordinates.hpp"

#define DEBUG_SNARL_COORDINATES

namespace stoat {


void SnarlCoordinates::add_reference_path(std::string path_name) {
    // Add the path name to path_to_index and path_by_index
    path_to_index.emplace(path_name, std::make_pair(path_by_index.size(), PathType::REF)); 
    path_by_index.emplace_back(std::move(path_name), PathType::REF);
}

std::tuple<std::string, size_t, size_t> SnarlCoordinates::get_reference_coordinates_as_string(const handlegraph::PathPositionHandleGraph& graph, 
                                                                                    const bdsg::SnarlDistanceIndex& distance_index,
                                                                                    net_handle_t snarl) {
    size_t ref_index, start, end;
    std::tie(ref_index, start_end) = get_reference_coordinates_as_index(graph, distance_index, snarl);
    return std::make_tuple(path_by_index[ref_index].first, start, end);
}

std::tuple<size_t, size_t, size_t> SnarlCoordinates::get_reference_coordinates_as_index(const handlegraph::PathPositionHandleGraph& graph, 
                                                                                    const bdsg::SnarlDistanceIndex& distance_index,
                                                                                    net_handle_t snarl) {
    #ifdef DEBUG_SNARL_COORDINATES
    std::cerr << "Get reference coordinates of snarl " << distance_index.net_handle_as_string(snarl) << std::endl;
    assert(distance_index.is_snarl(snarl) || distance_index.is_root(snarl));
    #endif

    net_handle_t canonical_snarl = distance_index.get_start_end_traversal_of(snarl);

    if (distance_index.is_root(canonical_snarl)) {
        return std::make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
    }
    // First see if we already found this snarl
    auto found_coords = snarl_to_coordinates.find(canonical_snarl);
    if (found_coords != snarl_to_coordinates.back()) {
        return found_coords->second;
    }

    // We prioritize coordinates on any of the reference paths for this snarl, then reference paths on an ancestor,
    // Then a reference-sense path on this snarl, then a reference-sense path of an ancestor, then any path on this snarl,
    // then any path on the ancestor
    std::vector<std::tuple<handlegraph::path_handle_t, size_t, size_t> traversals = get_traversals_of_snarl(graph, distance_index, cnanonical_snarl);

    // First look at all options on this snarl. These are indices into traversals
    size_t index_on_ref = std::numeric_limits<size_t>::max();
    size_t index_on_ref_sense = std::numeric_limits<size_t>::max();
    size_t index_on_any = std::numeric_limits<size_t>::max(); //non-ref or ref-sense path that we've already seen
    size_t index_on_new = std::numeric_limits<size_t>::max(); //non-ref or ref-sense path that we haven't seen before

    for (size_t i = 0 ; i < traversals ; i++) {

        handlegraph::path_handle_t path = std::get<0>(traversals[i]);
        std::string refpath = graph.get_path_name(path);
        auto found_ref = path_to_index.find(refpath);

        // What is this path? reference or reference sense or something else
        if (found_ref != path_to_index.back() && found_ref->second.second == PathType::REF) {
            // If this is a reference that we want, then we will return this
            std::tuple<size_t, size_t, size_t> found_coords(found_ref->second.first, std::get<1>(traversals), std::get<2>(traversals));
            snarl_to_coordinates.emplace(canonical_snarl, found_coords);
            return found_coords;
        } else if (graph.get_sense(path) == handlegraph::PathSense::REFERENCE) {
            // If it is a reference-sense path then just remember it 
            index_on_ref_sense = i;
        } else if (found_ref != path_to_index.back()) {
            // If we found coordinates on a path we already have
            index_on_any = i;
        } else {
            // If we found coordinates on a new path we haven't seen before
            index_on_new = i;
        } 
    }
    // If we didn't return something in the loop, then we didn't find a reference path
    // Look for the best option for the parent snarl (which will recursively check all the ancestors until it finds something)
    // I think that the snarl tree should be shallow enough that the recursion won't be a problem
    // Return the parent snarl's coordinates if it is better than anything we've found in this snarl
    std::tuple<size_t, size_t, size_t> parent_coords = get_reference_coordinates_as_index(graph, distance_index, distance_index.get_parent(distance_index.get_parent(snarl)));
    if (std::get<0>(parent_coords) != std::numeric_limits<size_t>::max()) {
        PathType parent_type = path_by_index[std::get<0>(parent_coords)].second;
        if (parent_type == PathType::REF || 
            (parent_type == PathType::REF_SENSE && index_on_ref_sense == std::numeric_limits<size_t>::max()) ||
            (index_on_any == std::numeric_limits<size_t>::max() && index_on_new == std::numeric_limits<size_t>::max())) {
            // If the parent had coordinates on the reference, or if this snarl had no coordinates on a reference sense but the parent did, 
            // or if this snarl didn't have any coordinates,
            // Return the parent coordinates
            snarl_to_coordinates.emplace(canonical_snarl, parent_coords);
            return parent_coords;
        }
    }

    // If we got here, then the parent's coordinates wouldn't have been better than what we found
    // Pick the best traversal we found, add the path if we haven't seen it before, add the snarl's coordinates, and return
    size_t traversal_index std::numeric_limits<size_t>::max();
    PathType path_type; 
    if (index_on_ref_sense != std::numeric_limits<size_t>::max()) {
        traversal_index =  index_on_ref_sense;
        path_type = PathType::REF_SENSE;
    } else if ( index_on_any != std::numeric_limits<size_t>::max() ) {
        traversal_index = index_on_any;
        path_type = PathType::OTHER;
    } else if (index_on_new != std::numeric_limits<size_t>::max()) {
        traversal_index = index_on_new;
        path_type = PathType::OTHER;
    } else {
        // If we didn't find any coordinates, return std::numeric_limits<size_t>::max() for everything
        std::tuple<size_t, size_t, size_t> empty_coords(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
        snarl_to_coordinates.emplace(canonical_snarl, empty_coords);
        return empty_coords;
    }

    // This is the traversal we want to return
    std::tuple<handlegraph::path_handle_t, size_t, size_t> found_traversal = traversals[traversal_index];
    // Rearrange it to the format we want, but we don't know the index of the path yet
    std::tuple<size_t, size_t, size_t> found_coords(0, std::get<1>(found_traversal), std::get<2>(found_traversal)); 
    std::string pathname = graph.get_path_name(std::get<0>(found_traversal));

    // Have we seen this path before?
    auto found_ref = path_to_index.find(pathname);
    if (found_ref != path_to_index.back()) {
        std::get<0>(found_coords) = found_ref->second.first;
    } else {
        // If we haven't seen this path before, add it to the list
        std::get<0>(found_coords) = path_by_index.size();
        path_to_index.emplace(pathname, std::make_pair(path_by_index.size(), path_type));
        paty_by_index.emplace_back(pathname, path_type);
    }
    snarl_to_coordinates.emplace(canonical_snarl, found_coords);

    return found_coords;

}

std::vector<std::string> SnarlCoordinates::reference_names_as_vector() {
    std::vector<std::string> path_names;
    path_names.reserve(path_by_index.size());
    for (const auto& path : path_by_index) {
        path_names.emplace_back(path.first);
    } 
    return path_names;
}


// Get traversals of the snarl as path name, start offset, end offset
std::vector<std::tuple<handlegraph::path_handle_t, size_t, size_t> get_traversals of snarl(const handlegraph::PathPositionHandleGraph& graph, 
                                                                                           const bdsg::SnarlDistanceIndex& distance_index,
                                                                                           net_handle_t snarl) {

    handlegraph::handle_t start_handle = distance_index.get_handle(distance_index.get_node_from_sentinel(distance_index.get_bound(ancestor_snarl, false, true)), &graph);
    handlegraph::handle_t end_handle = distance_index.get_handle(distance_index.get_node_from_sentinel(distance_index.get_bound(ancestor_snarl, true, true)), &graph);

    #ifdef DEBUG_SNARL_COORDINATES
        std::cerr << "Get traversals of snarl between " << graph.get_id(start_handle) << " and " << graph.get_id(end_handle) << std::endl;
        if (get_reference) {
            assert(sample_names.empty());
            assert(!get_all_paths);
        }
        if (!sample_names.empty()) {
            assert(!get_reference);
            assert(!get_all_paths);
        }
        if (get_all_paths) {
            assert(!get_reference);
            assert(sample_names.empty());
        }
    #endif

    // Map path to the steps on the path that traverse the snarl bounds
    std::map<handlegraph::path_handle_t, std::vector<handlegraph::step_handle_t>> path_to_steps;

    // Keep track if we found a traversal of the snarl (may be start-start or end-end)
    bool found_pair = false;

    // Get the step_handles of all paths traversing start or end
    // Steps don't care about the orientation of the handle, they will always (I think) be going forwards in the path
    
    graph.for_each_step_on_handle(start_handle, [&] (const handlegraph::step_handle_t& step) {
        handlegraph::path_handle_t path = graph.get_path_handle_of_step(step);
        if (path_to_steps.count(path) == 0) {
            path_to_steps[path] = std::vector<handlegraph::step_handle_t>();
        }
        path_to_steps[path].emplace_back(step);
        return true;
    });
    
    #ifdef DEBUG_SNARL_COORDINATES
        std::cerr << "After start node, found" << std::endl;
        for (const auto& x : path_to_steps) {
            std::cerr << graph.get_path_name(x.first) << ": " << x.second.size() << std::endl;
        }
    #endif
    
    graph.for_each_step_on_handle(end_handle, [&] (const handlegraph::step_handle_t& step) {
        handlegraph::path_handle_t path = graph.get_path_handle_of_step(step);
        if (path_to_steps.count(path) == 0) {
            path_to_steps[path] = std::vector<handlegraph::step_handle_t>();
        }
        path_to_steps[path].emplace_back(step);
        return true;
    });
    

    #ifdef DEBUG_SNARL_COORDINATES
        std::cerr << "After end node, found" << std::endl;
        for (const auto& x : path_to_steps) {
            std::cerr << graph.get_path_name(x.first) << ": " << x.second.size() << std::endl;
        }
    #endif

    // Range of paths as name of path, start offset, end offset 
    std::vector<std::tuple<handlegraph::path_handle_t, size_t, size_t> ranges;

    //If we found a path going through the snarl, return maximal ranges

    for (auto& path_steps : path_to_steps) {

        const handlegraph::path_handle_t& path = path_steps.first;
        std::vector<handlegraph::step_handle_t>& steps = path_steps.second;

        std::sort(steps.begin(), steps.end(), [&] (const handlegraph::step_handle_t& a, const handlegraph::step_handle_t& b) {
            return graph.get_position_of_step(a) < graph.get_position_of_step(b);
        });

        ranges.push_back({graph.get_path_handle_of_step(path),
                          graph.get_position_of_step(steps.front()) + graph.get_sequence(graph.get_handle_of_step(steps.front())).size(),
                          graph.get_position_of_step(steps.back())});
    }
    return ranges;
}
} // end namespace stoat

