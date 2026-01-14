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


}

#endif
