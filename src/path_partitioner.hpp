#ifndef STOAT_PATH_PARTITIONER_HPP_INCLUDED
#define STOAT_PATH_PARTITIONER_HPP_INCLUDED

#include <iostream>
#include <handlegraph/path_position_handle_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include "utils.hpp"
#include "snarl_data_collection.hpp"
#include <gbwt/gbwt.h>
#include <gbwtgraph/gbwtgraph.h>

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
/// This only finds start-end/end-start paths
/// TODO: This finds all steps along the path including those going through nested snarls. Could do all nested snarls at the same time
std::vector<size_t> partition_embedded_paths_in_snarl_with_gbwt(const handlegraph::PathPositionHandleGraph& graph, const gbwt::GBWT& gbwt, 
                                                                const bdsg::SnarlDistanceIndex& distance_index,
                                                                const net_handle_t& snarl,
                                                                const std::vector<stoat::sample_hap_t>& all_sample_haplotypes,
                                                                std::vector<PathTraversal>& paths_per_allele);

/// A struct for holding the growing path through the gbwt
struct gbwt_path_t {

    std::vector<handlegraph::net_handle_t> path;

    // a search state points to one node and a set of haplotypes passing through the node represented 
    // as a range in the bwt of the node. The search state will maintain a range corresponding to the
    // set of haplotypes that take the same path
    gbwt::SearchState search_state;

    // a unique identifier for the path. 
    size_t identifier;

    //the suffix array index for the first occurrence in the range (see https://github.com/jltsiren/gbwt/wiki/Fast-Locate)
    gbwt::size_type sa_index;

    // Destructively construct the struct by move()ing components
    gbwt_path_t (std::vector<handlegraph::net_handle_t> path, gbwt::SearchState search_state, size_t identifier, gbwt::size_type sa_index) :
                path(std::move(path)), search_state(std::move(search_state)), identifier(identifier), sa_index(sa_index) {}
};


/// Helper function for partition_embedded_paths_in_snarl_with_gbwt
/// Starting from first_path and first_state, find all paths traversing the snarl and add them to finished_paths and finished_search_states
/// The paths returned will contain chains instead of the full paths through them, so there may be identical chains in finished_paths
/// Return the number of distinct paths through the netgraph that were found.
/// Note that not all paths found may be returned (if there are fragmented paths that don't traverse both bounds), then so the return 
/// value is only useful as the maximum path identifier plus one.
size_t get_gbwt_traversals(const handlegraph::PathPositionHandleGraph& graph, const gbwt::GBWT& gbwt, 
                           const bdsg::SnarlDistanceIndex& distance_index,
                           const net_handle_t& snarl,
                           std::vector<gbwt_path_t>& finished_paths);
}

#endif
