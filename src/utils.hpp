#ifndef UTILS_HPP
#define UTILS_HPP

#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <tuple>
#include <unordered_set>
#include <iomanip>
#include <Eigen/Dense>
#include <fstream>

#include <bdsg/hash_graph.hpp>
#include <bdsg/packed_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include <bdsg/overlays/packed_path_position_overlay.hpp>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include "log.hpp"

using namespace std;

namespace stoat {

std::string format_group_paths(const std::vector<size_t>& g0, const std::vector<size_t>& g1);
std::string set_precision(const double& value);
std::string set_precision_float_50(const boost::multiprecision::cpp_dec_float_50& value);

bool is_na(const std::string& s);
double string_to_pvalue(const std::string& p1);

bool isPValueSignificant(const double& pvalue_threshold, const std::string& pvalue_str);
void retain_indices(std::vector<double>& vec, const std::unordered_set<size_t>& indices_to_keep);
std::vector<double> adjusted_holm(const std::vector<double>& p_values);

template <typename T>
std::string vectorToString(const std::vector<T>& vec);

template <typename T>
std::vector<T> stringToVector(const std::string& str);

// Given a path, return its sample name
std::string get_sample_name_from_path(const handlegraph::PathPositionHandleGraph& graph, const handlegraph::path_handle_t& path);

// This stores a sample name and haplotype number
struct sample_hap_t {
    std::string sample;
    std::size_t haplotype;

    const inline bool operator==(const sample_hap_t& other) const {
        return (sample==other.sample && haplotype==other.haplotype);
    }
    const inline bool operator<(const sample_hap_t& other) const {
        if (sample == other.sample) {
            return haplotype < other.haplotype;
        } else {
            return sample < other.sample;
        }
    }
};

inline std::ostream& operator<<(std::ostream& out, const sample_hap_t& sample) {
    return out << sample.sample << "#" << sample.haplotype;
}

// Given a path, return its sample name and haplotype 
sample_hap_t get_sample_and_haplotype(const handlegraph::PathPositionHandleGraph& graph, const handlegraph::path_handle_t& path);

// A struct for holding a range along the path
struct path_range_t {
    handlegraph::step_handle_t start;
    handlegraph::step_handle_t end;
};

/// Given a snarl, return a vector of path_ranges of that snarl (the boundary nodes).
/// Since a path can traverse a snarl multiple times, this returns each start-to-end (or end-to-start) range
/// of step_handle_t's, ordered according to the order of the path.
/// If the path leaves by the same bound (for example start-> start<- start-> end->), then the range will include
/// the outermost start->end range.
/// If get_reference is true, return a reference path and its coordinates.
/// This will first try to find a path with the sample name, if not empty, then a reference-sense path, then with any path traversing the snarl.
/// If get_reference is false, try to find coordinates on a path containing the given sample name, or if it fails, with any path.
/// If get_reference is false and sample_name is empty and get_all_paths is true, return all coordinates for all paths
/// For finding a specific path or reference, if the snarl is not on the desired path, then walk up the snarl tree to find an ancestor snarl on the path
std::vector<path_range_t> get_coordinates_of_snarl(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                                                   const handlegraph::net_handle_t& snarl, bool get_reference, std::string sample_name, bool get_all_paths);

/// The function that gets called by get_coordinates_of_snarl
/// This either looks for a particular sample, or a reference-sense path, or all paths
std::vector<path_range_t> get_coordinates_of_snarl_helper(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                                                          const handlegraph::net_handle_t& snarl, bool get_reference, std::string sample_name, bool get_all_paths);

/// Given a path_range_t representing a path going through a snarl (with the start and end step_handle_t's representing the boundary nodes)
/// Return the path name and range in the path of the snarl, not including the boundary nodes
std::tuple<std::string, size_t, size_t> get_name_and_offsets_of_snarl_path_range(const handlegraph::PathPositionHandleGraph& graph, 
                                                                                 const bdsg::SnarlDistanceIndex& distance_index, const path_range_t& range);

/// Function to find snarl ID- the start and end ids as a pair of size_t's
std::pair<size_t, size_t> find_snarl_id(const bdsg::SnarlDistanceIndex& stree, const handlegraph::net_handle_t& snarl);

/// Print ids of all nodes present in a snarl to stderr, one per line
/// Useful for debugging with `vg find -N`
void print_nodes_in_snarl(const bdsg::SnarlDistanceIndex& distance_index, const handlegraph::net_handle_t& snarl);


// equality within a given epsilon
template<typename T>
bool is_equal(T a, T b, T e = std::numeric_limits<T>::epsilon()) {
    return std::fabs(a-b) <= e;
};

enum phenotype_type_t { BINARY = 1, QUANTITATIVE, EQTL };

} // namespace stoat

#endif
