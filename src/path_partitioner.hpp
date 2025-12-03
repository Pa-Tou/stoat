#ifndef STOAT_PATH_PARTITIONER_HPP_INCLUDED
#define STOAT_PATH_PARTITIONER_HPP_INCLUDED

#include <iostream>
#include <handlegraph/path_position_handle_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include "utils.hpp"
#include "snarl_data_collection.hpp"

using namespace std;
using namespace stoat;

namespace stoat_graph {


/// Given a snarl, fill in sample_sets_by_allele. Each allele corresponds to a distinct walk through the snarl taken by a path.
/// This walk may leave the snarl and come back
/// If a walk doesn't traverse both bounds of the snarl then it is excluded
/// Fill in sample_sets_by_allele with the sample/haplotype taking each allele
/// This can be used by a SnarlDataCollection with a little adjustment
void partition_embedded_paths_in_snarl(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                          const net_handle_t& snarl,
                          const std::set<stoat::sample_hap_t>& all_sample_haplotypes,
                          std::vector<std::set<sample_hap_t>>& sample_sets_by_allele);


}

#endif
