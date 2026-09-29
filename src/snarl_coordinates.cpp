#include "snarl_coordinates.hpp"

//#define DEBUG_SNARL_COORDINATES

namespace stoat {


void SnarlCoordinates::add_reference_path(std::string path_name) {
    // Add the path name to path_to_index and path_by_index
    #pragma omp critical (SC_references) 
    {
    path_to_index.emplace(path_name, std::make_pair(path_by_index.size(), PathType::REF)); 
    path_by_index.emplace_back(std::move(path_name), PathType::REF);
    }
}

std::string SnarlCoordinates::get_path_name_from_index(size_t i) const {
    if (i == std::numeric_limits<size_t>::max()) {
        return "NA";
    } else {
        std::string pathname; 
        #pragma omp critical (SC_references) 
        {
        pathname = path_by_index[i].first;
        }
        return pathname;
    }
}

void SnarlCoordinates::clear() {
    #pragma omp critical (SC_references) 
    {
    path_by_index.clear();
    //TODO: I have no idea why this could segfault
    path_to_index.clear();
    }
    #pragma omp critical(SC_snarls) 
    {
    snarl_to_coordinates.clear();
    }
}

std::tuple<std::string, size_t, size_t> SnarlCoordinates::get_reference_coordinates_as_string(const handlegraph::PathPositionHandleGraph& graph, 
                                                                                    const bdsg::SnarlDistanceIndex& distance_index,
                                                                                    net_handle_t snarl) {
    size_t ref_index, start, end;
    std::tie(ref_index, start,end) = get_reference_coordinates_as_index(graph, distance_index, snarl);
    if (ref_index == std::numeric_limits<size_t>::max()) {
        return std::make_tuple("NA", start, end);
    }
    return std::make_tuple(get_path_name_from_index(ref_index), start, end);
}

std::tuple<size_t, size_t, size_t> SnarlCoordinates::get_reference_coordinates_as_index(const handlegraph::PathPositionHandleGraph& graph, 
                                                                                    const bdsg::SnarlDistanceIndex& distance_index,
                                                                                    net_handle_t snarl) {
    #ifdef DEBUG_SNARL_COORDINATES
    std::cerr << "Get reference coordinates of snarl " << distance_index.net_handle_as_string(snarl) << std::endl;
    assert(distance_index.is_snarl(snarl) || distance_index.is_root(snarl));
    #endif

    net_handle_t canonical_snarl = distance_index.start_end_traversal_of(snarl);

    if (distance_index.is_root(canonical_snarl)) {
        return std::make_tuple(std::numeric_limits<size_t>::max(), 0, 0);
    }
    // First see if we already found this snarl
    bool already_found = false;
    std::tuple<size_t, size_t, size_t> already_found_coords;
    #pragma omp critical(SC_snarls) 
    {
    auto found_coords = snarl_to_coordinates.find(canonical_snarl);
    if (found_coords != snarl_to_coordinates.end()) {
        already_found = true;
        already_found_coords = found_coords->second;
    }
    }
    if (already_found) {
        return already_found_coords;
    }

    // We prioritize coordinates on any of the reference paths for this snarl, then reference paths on an ancestor,
    // Then a reference-sense path on this snarl, then a reference-sense path of an ancestor, then any path on this snarl,
    // then any path on the ancestor
    std::tuple<handlegraph::path_handle_t, size_t, size_t> traversal = get_traversal_of_snarl(graph, distance_index, canonical_snarl);
    bool found_traversal = std::get<1>(traversal) != std::numeric_limits<size_t>::max();

    // Decide what type of path we just found
    bool traversal_is_ref = false;
    bool traversal_is_ref_sense = false;
    bool traversal_is_any = false; //non-ref or ref-sense path that we've already seen
    bool traversal_is_new = false; //non-ref or ref-sense path that we haven't seen before


    handlegraph::path_handle_t path = std::get<0>(traversal);
    std::string refpath = found_traversal ? graph.get_path_name(path) : "NA";
    std::tuple<size_t, size_t, size_t> traversal_coords(std::numeric_limits<size_t>::max(), std::get<1>(traversal), std::get<2>(traversal)); 
    #pragma omp critical (SC_references) 
    {
    if (found_traversal) {
        auto found_ref = path_to_index.find(refpath);
        // If we did find something for this snarl

        // What is this path? reference or reference sense or something else
        if (found_ref != path_to_index.end() && found_ref->second.second == PathType::REF) {
            traversal_coords = std::make_tuple(found_ref->second.first, std::get<1>(traversal), std::get<2>(traversal));
            // If this is a reference that we want, then we will return this
            #pragma omp critical(SC_snarls) 
            {
            snarl_to_coordinates.emplace(canonical_snarl, traversal_coords);
            }
            traversal_is_ref = true;
        } else if (graph.get_sense(path) == handlegraph::PathSense::REFERENCE) {
            // If it is a reference-sense path then just remember it 
            traversal_is_ref_sense = true;
        } else if (found_ref != path_to_index.end()) {
            // If we found coordinates on a path we already have
            traversal_is_any = true;
        } else {
            // If we found coordinates on a new path we haven't seen before
            traversal_is_new = true;
        } 

    }
    } //end omp critical

    // If we already found a reference path through the snarl, return it
    if (traversal_is_ref) {
        return traversal_coords;
    }
    
    // Look for the best option for the parent snarl (which will recursively check all the ancestors until it finds something)
    // I think that the snarl tree should be shallow enough that the recursion won't be a problem
    // Return the parent snarl's coordinates if it is better than anything we've found in this snarl
    std::tuple<size_t, size_t, size_t> parent_coords = get_reference_coordinates_as_index(graph, distance_index, distance_index.get_parent(distance_index.get_parent(snarl)));
    if (std::get<0>(parent_coords) != std::numeric_limits<size_t>::max()) {
        PathType parent_type;
        #pragma omp critical (SC_references) 
        {
        parent_type = path_by_index[std::get<0>(parent_coords)].second;
        }
        if (!found_traversal || parent_type == PathType::REF || 
            (parent_type == PathType::REF_SENSE && !traversal_is_ref_sense) ||
            (!traversal_is_ref_sense && !traversal_is_any && !traversal_is_new)) {
            // If the parent had coordinates on the reference, or if this snarl had no coordinates on a reference sense but the parent did, 
            // or if this snarl didn't have any coordinates,
            // Return the parent coordinates
            #pragma omp critical(SC_snarls) 
            {
            snarl_to_coordinates.emplace(canonical_snarl, parent_coords);
            }
            return parent_coords;
        }
    }
    if (!found_traversal) {
        return std::make_tuple(std::numeric_limits<size_t>::max(), 0, 0);
    }

    // If we got here, then the parent's coordinates wouldn't have been better than what we found
    // Add the path if we haven't seen it before, add the snarl's coordinates, and return
    PathType path_type; 
    if (traversal_is_ref_sense) {
        path_type = PathType::REF_SENSE;
    } else if ( traversal_is_any ) {
        path_type = PathType::OTHER;
    } else if (traversal_is_new) {
        path_type = PathType::OTHER;
    } else {
        // If we didn't find any coordinates, return std::numeric_limits<size_t>::max() for everything
        std::tuple<size_t, size_t, size_t> empty_coords(std::numeric_limits<size_t>::max(), 0, 0);
        #pragma omp critical(SC_snarls) 
        {
        snarl_to_coordinates.emplace(canonical_snarl, empty_coords);
        }
        return empty_coords;
    }

    // This is the traversal we want to return
    // Rearrange it to the format we want, but we don't know the index of the path yet
    std::string pathname = graph.get_path_name(std::get<0>(traversal));

    #pragma omp critical (SC_references) 
    {
    // Have we seen this path before?
    auto found_ref = path_to_index.find(pathname);
    if (found_ref != path_to_index.end()) {
        std::get<0>(traversal_coords) = found_ref->second.first;
    } else {
        // If we haven't seen this path before, add it to the list
        std::get<0>(traversal_coords) = path_by_index.size();
        path_to_index.emplace(pathname, std::make_pair(path_by_index.size(), path_type));
        path_by_index.emplace_back(pathname, path_type);
    }
    }// end omp critical
    #pragma omp critical(SC_snarls) 
    {
    snarl_to_coordinates.emplace(canonical_snarl, traversal_coords);
    }


    return traversal_coords;

}

std::vector<std::string> SnarlCoordinates::reference_names_as_vector() const {
    std::vector<std::string> path_names;
    #pragma omp critical (SC_references) 
    {
    path_names.reserve(path_by_index.size());
    for (const auto& path : path_by_index) {
        path_names.emplace_back(path.first);
    } 
    }
    return path_names;
}


// Get a traversal of the snarl as path name, start offset, end offset
// Should be the best we can find, prioritizing reference, then reference-sense, then any path that we've already seen
std::tuple<handlegraph::path_handle_t, size_t, size_t> SnarlCoordinates::get_traversal_of_snarl(const handlegraph::PathPositionHandleGraph& graph, 
                                                                                           const bdsg::SnarlDistanceIndex& distance_index,
                                                                                           net_handle_t snarl) {


    handlegraph::handle_t start_handle = distance_index.get_handle(distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, true)), &graph);
    handlegraph::handle_t end_handle = distance_index.get_handle(distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, true, true)), &graph);

    #ifdef DEBUG_SNARL_COORDINATES
        std::cerr << "Get traversals of snarl between " << graph.get_id(start_handle) << " and " << graph.get_id(end_handle) << std::endl;
    #endif

    // Map path to the steps on the path that traverse the snarl bounds
    std::map<handlegraph::path_handle_t, std::vector<handlegraph::step_handle_t>> path_to_steps;

    // Keep track if we found a traversal of the snarl (may be start-start or end-end)
    bool found_pair = false;

    // Get the step_handles of all paths traversing start or end
    // Steps don't care about the orientation of the handle, they will always (I think) be going forwards in the path

    
    // Going through all the path and finding their offsets is super slow so try to limit the number of paths we keep.
    // If we find a reference path, don't look for anything else. If we find a reference-sense path, don't look for non-ref paths
    bool only_ref = false;
    bool only_ref_sense = false;
    bool only_found = false;
    auto want_ref_path = [&](const handlegraph::path_handle_t& path ) {
        bool is_ref = false;
        bool is_ref_sense = false;
        bool is_found = false;
        #pragma omp critical (SC_references) 
        {
            auto found_ref = path_to_index.find(graph.get_path_name(path));

            // What is this path? reference or reference sense or something else
            if (found_ref != path_to_index.end() && found_ref->second.second == PathType::REF) {
                is_ref = true;
                only_ref = true;
            } else if (graph.get_sense(path) == handlegraph::PathSense::REFERENCE) {
                is_ref_sense = true;
                only_ref_sense = true;
            } else if (found_ref != path_to_index.end()) {
                is_found = true;
                only_found = true;
            } 
        } //end omp critical
        bool keep_this_path = false;

        if ((only_ref && is_ref) || 
            (!only_ref && only_ref_sense && is_ref_sense) ||
            (!only_ref && !only_ref_sense && only_found && is_found) ||
            (!only_ref && !only_ref_sense && !only_found)) {
            return true;
        } else {
            return false;
        }
    };

    graph.for_each_step_on_handle(start_handle, [&] (const handlegraph::step_handle_t& step) {
        handlegraph::path_handle_t path = graph.get_path_handle_of_step(step);

        if (want_ref_path(path)) {

            if (path_to_steps.count(path) == 0) {
                path_to_steps[path] = std::vector<handlegraph::step_handle_t>();
            }
            path_to_steps[path].emplace_back(step);
        }
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
        if (want_ref_path(path)) {
            if (path_to_steps.count(path) == 0) {
                path_to_steps[path] = std::vector<handlegraph::step_handle_t>();
            }
            path_to_steps[path].emplace_back(step);
        }
        return true;
    });
    

    #ifdef DEBUG_SNARL_COORDINATES
        std::cerr << "After end node, found" << std::endl;
        for (const auto& x : path_to_steps) {
            std::cerr << graph.get_path_name(x.first) << ": " << x.second.size() << std::endl;
        }
    #endif


    //If we found a path going through the snarl, return maximal ranges

    for (auto& path_steps : path_to_steps) {

        const handlegraph::path_handle_t& path = path_steps.first;
        if (want_ref_path(path)) {
            std::vector<handlegraph::step_handle_t>& steps = path_steps.second;

            std::sort(steps.begin(), steps.end(), [&] (const handlegraph::step_handle_t& a, const handlegraph::step_handle_t& b) {
                return graph.get_position_of_step(a) < graph.get_position_of_step(b);
            });
            //TODO: If there is just one traversal, need to decide if we need to add the node offset or not (depending on if it is going into or out of the snarl
            size_t start_offset = graph.get_position_of_step(steps.front()) + graph.get_sequence(graph.get_handle_of_step(steps.front())).size();
            size_t end_offset = steps.size() == 1 ? start_offset : graph.get_position_of_step(steps.back());

            // want_ref_path() will tell us if this was the best path, so return it immediately
            #ifdef DEBUG_SNARL_COORDINATES
            std::cerr << "Found best range for snarl " << distance_index.net_handle_as_string(snarl) << ": " << graph.get_path_name(path) << ":" << start_offset << "-" << end_offset << std::endl;
            #endif
            return std::make_tuple(path, start_offset, end_offset);
        }
    }
    #ifdef DEBUG_SNARL_COORDINATES
    std::cerr << "Couldn't find a traversal of this snarl" << std::endl;
    #endif
    return std::make_tuple(handlegraph::path_handle_t(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
}
} // end namespace stoat

