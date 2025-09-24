#ifndef STOAT_PARTITIONER_HPP_INCLUDED
#define STOAT_PARTITIONER_HPP_INCLUDED

#include <iostream>
#include <handlegraph/path_position_handle_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include "utils.hpp"

using namespace std;
using namespace stoat;

namespace stoat_graph {

/***
    General template class for finding partitions of samples in a snarl and serializing the sample partitions.
***/
class Partitioner {

    public:
        Partitioner(std::set<stoat::sample_hap_t> all_sample_haplotypes) :
            all_sample_haplotypes(std::move(all_sample_haplotypes)) {}

        /// The main function of this class
        /// Given a snarl, return a partition of samples that will be used to determine association with samples of interest.
        /// This is a template function that must be implemented by inherited classes
        /// If save_partitions (member variable of class) is true, then this should add the result to snarl_to_partitions 
        virtual std::vector<std::set<std::string>> partition_samples_in_snarl(const handlegraph::PathPositionHandleGraph& graph, 
                                                                              const bdsg::SnarlDistanceIndex& distance_index, 
                                                                              const handlegraph::net_handle_t& snarl) const = 0;
        
        void serialize(const std::string& filename);
        void deserialize(const std::string& filename);

    protected:

        /// A set of all samples+haplotypes in the graph
        std::set<stoat::sample_hap_t> all_sample_haplotypes;

        /// If this is true, then partition_samples_in_snarl should also add its results 
        bool save_partitions;
        // If we want to serialize the snarls we found

        //TODO: Use indices instead of actual sample names 
        /// A data structure to store all information about a snarl
        /// - size (as the "maximum" length of the snarl)
        /// - partitions of samples
        struct snarl_sample_partition_t {
            size_t max_length;
            std::vector<std::set<sample_hap_t>> partitions;
        };
        
        std::unordered_map<std::pair<size_t, size_t>, snarl_sample_partition_t> snarl_to_partitions;


};

class PathPartitioner : public Partitioner {
    public:
        
        PathPartitioner(std::set<stoat::sample_hap_t> all_sample_haplotypes) : Partitioner(std::move(all_sample_haplotypes)) {}

        std::vector<std::set<std::string>> partition_samples_in_snarl(const handlegraph::PathPositionHandleGraph& graph, 
                                                                      const bdsg::SnarlDistanceIndex& distance_index, 
                                                                      const handlegraph::net_handle_t& snarl) const;


        /// Given a snarl, partition the paths going through the snarl based on the walks they take in the netgraph.
        /// Unlike get_start_edge_sets, any path that doesn't traverse the snarl will not be included in any set
        /// Returns sets of samples + haplotypes
        std::vector<std::set<stoat::sample_hap_t>> get_walk_sets(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index, const bdsg::net_handle_t& snarl, bool only_bound) const;


};

}

#endif
