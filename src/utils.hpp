#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <string>
#include <tuple>
#include <unordered_set>
#include <string>
#include <charconv>
#include <fstream>
#include <stdexcept>

#include <bdsg/snarl_distance_index.hpp>
#include <bdsg/overlays/packed_path_position_overlay.hpp>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>

#include "types_and_structs.hpp"

#include "log.hpp"


namespace stoat {

// Helper function to trim whitespace from a string 
std::string trim(const std::string& value);

// Parse a haplotype count from a string. Throws an exception if the string is not a valid number.
size_t parse_count(const std::string& count_text, const std::string& context);

std::string set_precision(const double& value);

bool is_na(const std::string& s);
double string_to_pvalue(const std::string& p1);

std::pair<double, size_t> adjusted_hochberg(const std::vector<double>& p_values);

template <typename T>
std::vector<T> stringToVector(const std::string& str);

// Given a path, return its sample name
std::string get_sample_name_from_path(const handlegraph::PathHandleGraph& graph, const handlegraph::path_handle_t& path);

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
