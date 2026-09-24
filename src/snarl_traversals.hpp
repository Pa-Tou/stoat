#ifndef STOAT_SNARL_TRAVERSALS_HPP_INCLUDED
#define STOAT_SNARL_TRAVERSALS_HPP_INCLUDED

#include <handlegraph/path_position_handle_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include "types_and_structs.hpp"
#include <gbwt/gbwt.h>
#include <gbwt/fast_locate.h>
#include <gbwtgraph/gbwtgraph.h>


namespace stoat {

/// Find all possible walks through the snarl's netgraph. Fills in walks
/// Walks include the start and end bounds of the snarl
/// If a path cycles more than walk_cycle_limit times, stop looking for more cycles
/// If a path takes more than walk_steps_limit steps, stop extending this path 
void get_all_walks_through_snarl(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index, 
                                 const net_handle_t& snarl, std::vector<stoat::PathTraversal>& walks,
				                 size_t walk_cycle_limit = 1, size_t walk_steps_limit = 50);

/// Find all walks through the snarl netgraph that are taken by haplotypes in the gbwt and fill in walks
void get_haplotype_walks_through_snarl(const handlegraph::PathPositionHandleGraph& graph, const gbwt::GBWT& gbwt, 
                           const gbwt::FastLocate& r_index,const bdsg::SnarlDistanceIndex& distance_index, 
                                 const net_handle_t& snarl, std::vector<stoat::PathTraversal>& walks);

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
                           const gbwt::FastLocate& r_index,
                           const bdsg::SnarlDistanceIndex& distance_index,
                           const net_handle_t& snarl,
                           std::vector<gbwt_path_t>& finished_paths);
} //end stoat namespace
#endif
