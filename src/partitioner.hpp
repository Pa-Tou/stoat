#ifndef STOAT_PARTITIONER_HPP_INCLUDED
#define STOAT_PARTITIONER_HPP_INCLUDED

#include <iostream>
#include <handlegraph/path_position_handle_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include "utils.hpp"
#include "snarl_data_t.hpp"

using namespace std;
using namespace stoat;

namespace stoat_graph {

/***
    General template class for finding partitions of samples in a snarl and serializing the sample partitions.
    This class can find partitions of snarls from the distance index and graph, and optionally serialize those partitions,
    or it can read serialized partitions without the distance index.

***/
class SnarlTraverserAndPartitioner {

    public:
        SnarlTraverserAndPartitioner(const std::set<stoat::sample_hap_t>& all_sample_haplotypes, const bdsg::SnarlDistanceIndex* distance_index, const std::string& reference_sample, bool save_partitions) :
            all_sample_haplotypes(all_sample_haplotypes),
            distance_index(distance_index),
            reference_sample(reference_sample),
            save_partitions(save_partitions),
            check_distances(distance_index == nullptr ? false : distance_index->has_distances()) {}


        /// This runs the function iteratee for every snarl, either using the distance index or, if given, the serialized snarls
        /// If using the distance index, check if the snarl is eligible before computing snarl_info.
        /// If using the distance index and save_partitions is true, then also serialize the snarl
        void for_each_snarl_partition(const handlegraph::PathPositionHandleGraph& graph, 
                            const std::function<bool(const handlegraph::net_handle_t& net)>& snarl_is_eligible, 
                            const std::function<void(const snarl_partition_t& snarl_info)>& iteratee);

        
        /// Serialize all the snarl_partition_t's
        void serialize(const std::string& filename);

        /// Desrialize all the snarl_partition_t's
        void deserialize(const std::string& filename);


    protected:

        ///////////////////////////////////////////////////////// Member functions

        /// If this is nullptr, then we are only loading the snarls
        const bdsg::SnarlDistanceIndex* distance_index;

        // The reference sample that we try to get coordinates against
        const std::string& reference_sample;

        /// A set of all samples+haplotypes in the graph
        const std::set<stoat::sample_hap_t>& all_sample_haplotypes;

        // Does the distance index contain distances?
        bool check_distances;

        /// If this is true, then we want to serialize the snarls we found
        bool save_partitions;

        

        /// This holds all snarl info.
        /// If we are not loading or saving the snarls, then this doesn't get filled in
        std::vector<snarl_partition_t> snarl_partitions;

        // The names of the references we use to output snarl coordinates
        // This is only filled in if we're going to save the snarls
        std::set<std::string> all_references;

        /// This goes at the beginning of a serialized file to ensure that it is the right file type and version
        inline const static std::string file_header = "#SNARL_PARTITIONS_v1.0";


        //////////////////////////////////////////////////////////// protected functions
        

        /// Given a snarl, return a partition of samples that will be used to determine association with samples of interest.
        /// This function uses the distance index to fill in snarl_info
        /// This is a template function that must be implemented by inherited classes
        /// If save_partitions (member variable of class) is true, then this should add the result to snarl_to_partitions 
        virtual std::vector<std::set<sample_hap_t>> partition_samples_in_snarl(const handlegraph::PathPositionHandleGraph& graph, 
                                                                              const handlegraph::net_handle_t& snarl) const = 0;


};

class SnarlTraverserAndPathPartitioner : public SnarlTraverserAndPartitioner {
    public:
        
        SnarlTraverserAndPathPartitioner(const std::set<stoat::sample_hap_t>& all_sample_haplotypes, const bdsg::SnarlDistanceIndex* distance_index, const std::string& reference_sample, bool save_partitions)
            : SnarlTraverserAndPartitioner(all_sample_haplotypes, distance_index, reference_sample, save_partitions) {}

    protected:

        std::vector<std::set<sample_hap_t>> partition_samples_in_snarl(const handlegraph::PathPositionHandleGraph& graph, 
                                                                      const handlegraph::net_handle_t& snarl) const;


        /// Given a snarl, partition the paths going through the snarl based on the walks they take in the netgraph.
        /// Unlike get_start_edge_sets, any path that doesn't traverse the snarl will not be included in any set
        /// Returns sets of samples + haplotypes
        std::vector<std::set<stoat::sample_hap_t>> get_walk_sets(const handlegraph::PathPositionHandleGraph& graph, const bdsg::net_handle_t& snarl, bool only_bound) const;


};

}

#endif
