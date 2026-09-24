#ifndef STOAT_SNARL_TRAVERSALS_HPP_INCLUDED
#define STOAT_SNARL_TRAVERSALS_HPP_INCLUDED

#include <handlegraph/path_position_handle_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include "types_and_structs.hpp"

namespace stoat {

/// Find all possible walks through the snarl's netgraph. Fills in walks
/// Walks include the start and end bounds of the snarl
/// If a path cycles more than walk_cycle_limit times, stop looking for more cycles
/// If a path takes more than walk_steps_limit steps, stop extending this path 
void get_all_walks_through_snarl(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index, 
                                 const net_handle_t& snarl, std::vector<stoat::PathTraversal>& walks,
				                 size_t walk_cycle_limit = 1, size_t walk_steps_limit = 50);
}

#endif
