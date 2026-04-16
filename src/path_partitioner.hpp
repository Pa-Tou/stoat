#ifndef STOAT_PATH_PARTITIONER_HPP_INCLUDED
#define STOAT_PATH_PARTITIONER_HPP_INCLUDED

#include <iostream>
#include <handlegraph/path_position_handle_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include "utils.hpp"
#include "snarl_data_collection.hpp"

using namespace stoat;

namespace stoat_graph {


/// Given a snarl, assign each path in the snarl to an allele. Each allele corresponds to a distinct walk through the snarl taken by a path.
/// This walk may leave the snarl and come back
/// If a walk doesn't traverse both bounds of the snarl then it is excluded
/// The returned vector contains a allele assignment for all sample/haplotypes in all_sample_haplotypes. 
/// std::numeric_limits<size_t>::max() for a path that didn't traverse this snarl
/// All paths in the snarl must correspond to exactly one sample_hap_t in all_sample_haplotypes.
/// This can be used by a SnarlDataCollection with a little adjustment
std::vector<size_t> partition_embedded_paths_in_snarl(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                          const net_handle_t& snarl,
                          const std::vector<stoat::sample_hap_t>& all_sample_haplotypes);


/// The same as partition_embedded_paths_in_snarl, except using a GBWT. 
/// Fills in allele assignments, paths_per_allele, and, optionally, sequences_per_allele
/// TODO: This finds all steps along the path including those going through nested snarls. Could do all nested snarls at the same time
void partition_embedded_paths_in_snarl_with_gbwt(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                          const net_handle_t& snarl,
                          const std::vector<stoat::sample_hap_t>& all_sample_haplotypes,
                          std::vector<size_t>& allele_assignments,
                          std::vector<PathTraversal>& paths_per_allele);

/// Helper function for partition_embedded_paths_in_snarl_with_gbwt
/// Starting from first_path and first_state, find all paths traversing the snarl and add them to finished_paths and finished_search_states
/// If only_loops is true, then only look for paths that loop back to the start node (first_net, which points into the snarl)
void get_gbwt_traversals(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                         const net_handle_t& snarl, std::vector<std::vector<handlegraph::net_handle_t>>& finished_paths,
                         std::vector<gbwt::SearchState>& finished_search_states, vector<gbwt::node_type> first_path,
                         gbwt::SearchState first_state, const handlegraph::net_handle_t& first_net, bool only_loops);

#endif
