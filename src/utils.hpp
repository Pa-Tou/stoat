#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <string>
#include <tuple>
#include <unordered_set>

#include <bdsg/snarl_distance_index.hpp>
#include <bdsg/overlays/packed_path_position_overlay.hpp>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>

#include <boost/multiprecision/cpp_dec_float.hpp>

#include "types_and_structs.hpp"

#include "log.hpp"


namespace stoat {

    // std::string format_group_paths(const std::vector<size_t>& g0, const std::vector<size_t>& g1);
std::string set_precision(const double& value);

bool is_na(const std::string& s);
double string_to_pvalue(const std::string& p1);

std::pair<double, size_t> adjusted_hochberg(const std::vector<double>& p_values);

template <typename T>
std::vector<T> stringToVector(const std::string& str);

// Given a path, return its sample name
std::string get_sample_name_from_path(const handlegraph::PathHandleGraph& graph, const handlegraph::path_handle_t& path);

/// Given a snarl, return a vector of path_ranges of that snarl (the boundary nodes).
/// Since a path can traverse a snarl multiple times, this returns each start-to-end (or end-to-start) range
/// of step_handle_t's, ordered according to the order of the path.
/// If get_reference is true, return a reference path and its coordinates.
/// This will first try to find a path with the sample name, if not empty, then a reference-sense path, then with any path traversing the snarl.
/// If get_reference is false, try to find coordinates on a path containing the given sample name, or if it fails, with any path.
/// If get_reference is false and sample_name is empty and get_all_paths is true, return all coordinates for all paths
/// For finding a specific path or reference, if the snarl is not on the desired path, then walk up the snarl tree to find an ancestor snarl on the path
/// An ancestor also takes priority over a different path, so if there is a reference-sense path of the snarl, but sample_name on the ancestor, return sample_name on the ancestor
std::vector<stoat::path_range_t> get_coordinates_of_snarl(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                                                   const handlegraph::net_handle_t& snarl, bool get_reference, const std::unordered_set<std::string>& sample_names, bool get_all_paths);

/// The function that gets called by get_coordinates_of_snarl
/// This either looks for a particular sample, or a reference-sense path, or all paths
/// It finds paths between two nodes, which should be the start and end bounds of a snarl, pointing into each other
std::vector<stoat::path_range_t> get_coordinates_between_nodes(const handlegraph::PathPositionHandleGraph& graph, const handlegraph::handle_t& start_handle,
                                                          const handlegraph::handle_t& end_handle, bool get_reference, const std::unordered_set<std::string>& sample_names, bool get_all_paths);

/// Given a path_range_t representing a path going through a snarl (with the start and end step_handle_t's representing the boundary nodes)
/// Return the path name and range in the path of the snarl, not including the boundary nodes
std::tuple<std::string, size_t, size_t> get_name_and_offsets_of_snarl_path_range(const handlegraph::PathPositionHandleGraph& graph, 
                                                                                 const stoat::path_range_t& range);

/// Print ids of all nodes present in a snarl to stderr, one per line
/// Useful for debugging with `vg find -N`
void print_nodes_in_snarl(const bdsg::SnarlDistanceIndex& distance_index, const handlegraph::net_handle_t& snarl);


// equality within a given epsilon
template<typename T>
bool is_equal(T a, T b, T e = std::numeric_limits<T>::epsilon()) {
    return std::fabs(a-b) <= e;
};

} // namespace stoat


#endif
