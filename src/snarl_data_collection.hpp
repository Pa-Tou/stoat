#ifndef STOAT_SNARL_DATA_COLLECTION_HPP_INCLUDED
#define STOAT_SNARL_DATA_COLLECTION_HPP_INCLUDED

#include <iostream>
#include <handlegraph/path_position_handle_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include "utils.hpp"
#include "snarl_data_t.hpp"

using namespace std;
using namespace stoat;

namespace stoat {

// This holds snarl information from the distance index and graph: id of the snarl, reference coordinates, walks, 
// the sequences of the walks and sets samples/haplotypes taking each walk
// The latter two may be empty
struct snarl_info_t {
    public:

        // Constructor from elements
        snarl_info_t(stoat::Node_traversal_t start, stoat::Node_traversal_t end_node, std::string ref_path, 
                     size_t start_position, size_t end_position, size_t depth, const std::string& variant_type, 
                     const std::vector<Path_traversal_t>& snarl_walks, const std::vector<std::set<sample_hap_t>>& partitions, 
                     const std::vector<std::string>& sequences) :

                     start_node(start_node), end_node(end_node), 
                     ref_path(ref_path), start_position(start_position), end_position(end_position), depth(depth),
                     variant_type(variant_type), snarl_walks(snarl_walks),  partitions(partitions), sequences(sequences)
                     {};

        // Start and end nodes, both pointing into the snarl
        stoat::Node_traversal_t start_node;
        stoat::Node_traversal_t end_node;

        std::string ref_path; //TODO: I think this could get pretty big, might want to save it as an index into a list of reference paths
        // Start and end offset along the reference path
        size_t start_position;
        size_t end_position;
        size_t depth;

        // The "variant type" of the snarl, which represents the min/max length (or 0 for a deletion) of each walk in snarl_walks
        const std::string& variant_type;

        // All possible walks through the snarl
        const std::vector<Path_traversal_t>& snarl_walks;

        // For each walk in snarl_walks, which sample/haplotypes take that walk
        const std::vector<std::set<sample_hap_t>>& partitions;

        // For each walk in snarl_walks, what is its sequence?
        const std::vector<std::string>& sequences;

};

/// A class for holding per-snarl data for a collection of snarls
/// Publicly, this allow access to snarl_info_t's for holding basic information about the snarl, the walks through the snarl 
///    the groups of haplotypes following each walk, and the sequence of each walk.
/// The latter two are both optional and may be empty if the SnarlDataCollection was built without them.
/// Internally, the basic snarl information, the walks, haplotype partitions, and sequences are each stored separately. 
/// Build the collection either from the distance index or by loading from a file
/// Access snarls by chromosome name or all at once

class SnarlDataCollection {



    public:
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////// Public functions for building and interrogating snarl info

        SnarlDataCollection(size_t allele_size_limit, size_t snarl_child_limit, size_t walk_cycle_limit, size_t walk_steps_limit);

        /// Fill in the SnarlDataCollection for all snarls in the distance index
        /// If partition_requested is true, then call find_sample_partitions and save the output partitions.
        /// If sequence_requested is true, then find the sequence of each walk 
        /// If reference_sample is not empty, get coordinates on this reference path
        /// Write failed snarls to out_fail
        void fill_in_snarl_info(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                            bool partition_requested,
                            const std::function<std::vector<std::set<sample_hap_t>>(const net_handle_t& snarl)>& find_sample_partitions,
                            bool sequence_requested, const string& reference_sample, bool check_distances,
                            std::ostream& out_fail);

        
        /// If the snarl partitions (which assigns sample/haplotypes to each snarl_walk) were not found during construction, go through
        /// all snarls and call find_sample_partitions to add partitions
        /// find_sample_partitions takes a snarl_info_t (which is not guaranteed to contain partitions or sequences) and returns the partitions
        /// corresponding to the snarl_walks in snarl_data. This should return the partitions field of the snarl_info_t, but since the snarl_info_t
        /// doesn't exist in the SnarlDataCollection it is immutable and the partitions must be returned and saved separately 
        /// Note that this will overwrite any existing partitions. 
        /// If chr is specified, run this only on snarls on the given chromosome (as reference path name)
        //TODO:  I think this should be fine to run multithreaded as long as the list of snarl data doesn't move around, and I think the object
        //        stores a reference to the vector somewhere else in memory
        void add_snarl_partitions(const std::function<std::vector<std::set<sample_hap_t>>(const snarl_info_t& snarl_data,
                                                                                          const std::vector<Path_traversal_t>& walks)>& find_sample_partitions,
                                  std::string chr="");

        /// Run iteratee for all snarls
        void for_each_snarl(const std::function<void(const snarl_info_t& snarl_info)>& iteratee) const;

        /// Write the collection of snarls to the given file
        void write_snarl_data_collection(std::ostream& outstream) const;
        
        /// Load the collection of snarls from the given file
        /// Warn if the allele_size_limit or snarl_child_limit of the file are less permissive than this SnarlDataCollection
        void load_snarl_data_collection(std::istream& instream); 

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////// Private data members

    private:

        /// This stores the basic information from the snarl_info_t
        struct snarl_info_internal_t {
    
            // Start and end nodes, both pointing into the snarl
            stoat::Node_traversal_t start_node = stoat::Node_traversal_t(0, false);
            stoat::Node_traversal_t end_node= stoat::Node_traversal_t(0, false);

            // Index into reference_names to get the string representation of the reference path
            size_t reference_index;
    
            // Start and end offset along the reference path
            size_t start_position;
            size_t end_position;

            // Depth of the snarl in the snarl decomposition
            size_t depth;

            // The "variant type" of the snarl, which represents the min/max length (or 0 for a deletion) of each walk in snarl_walks
            std::string variant_type;

        };

        //////////////////////////// The stuff holding the data 

        /// This holds the snarl data as a map from the chromosome name to the data
        //TODO: Make sure that this gets the chr name the way Matis did it
        //TODO: idk if I want it to be a map from chr to vector of snarl data or just a vector and then check the chromosome for the per-chromosome calls 
        //std::unordered_map<std::string, std::vector<snarl_info_internal_t>> chr_to_snarl_data;
        std::vector<snarl_info_internal_t> all_snarl_data;


        /// Map snarl (as the start node, which uniquely identifies the snarl) to the walks through the snarl.
        std::unordered_map<stoat::Node_traversal_t, std::vector<Path_traversal_t>> snarl_to_walks;

        /// Map snarl (as the start node, which uniquely identifies the snarl) to the partitions.
        /// The vector follows the vector of walks in snarl_to_walks, meaning that each set in snarl_to_partitions
        ///  represents the sample/haplotype (as an index into sample_haplotypes) that take the walk at the same
        ///  index in the snarl's entry in snarl_to_walks. 
        std::unordered_map<stoat::Node_traversal_t, std::vector<std::set<size_t>>> snarl_to_partitions;

        /// Map snarl (as the start node, which uniquely identifies the snarl) to the sequences of the walks.
        /// Each string in the vector is the sequence for the walk in the corresponding vector in snarl_to_walks. 
        std::unordered_map<stoat::Node_traversal_t, std::vector<std::string>> snarl_to_sequences;


        ///////////////// Lists of strings and stuff that are stored as indexes in the real data structures instead of duplicating them a bunch

        /// The string representations of reference paths
        /// Used by the snarl_info_internal_t's reference_index
        std::vector<std::string> reference_names;

        /// This stores all sample_hap_t's that are stored as indexes by snarl_to_partitions
        vector<sample_hap_t> sample_haplotypes;


        //////////////////////////// Extra housekeeping stuff

        /// This goes at the beginning of the file to ensure that it is the right file type and version
        inline const static std::string file_header = "#SNARL_DATA_v1.0";

        /// Skip snarls if their maximum length is smaller than this
        size_t allele_size_limit;

        /// Skip snarls if they have more children than this
        size_t snarl_child_limit;

        /// Don't include walks that have more than this many cycles
        size_t walk_cycle_limit;

        /// Don't include snarls if enumerating all its walks takes more than this many steps
        size_t walk_steps_limit;


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////// Private functions
    private:

        // Given a snarl, enumerate the walks going through the snarl and get the variant_type string 
        // If a snarl is skipped because the walks can't be found, remove it and write it to out_fail
        std::tuple<std::vector<stoat::Path_traversal_t>, std::vector<std::string>> get_walks_through_snarl(const handlegraph::PathPositionHandleGraph& graph, 
                const bdsg::SnarlDistanceIndex& distance_index, const net_handle_t& snarl, bool check_distances, std::ostream& out_fail) const;

        // Given the partitions of haplotype paths in a snarl, find the walks as Path_traversal_t's, the variant_type (walk lengths as length or min/max per walk), 
        //and the sequences corresponding to each partition.
        std::tuple<std::vector<stoat::Path_traversal_t>, std::vector<std::string>> get_walks_from_partitions(
                const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index, const net_handle_t& snarl,
                const std::vector<std::set<sample_hap_t>>& sample_partitions, bool check_distances) const; 
    
        // Given the walks through the snarl, find the sequence. The sequence will just be a concatination of sequences of nodes, ignoring anything else
        std::vector<std::string> get_sequences_from_walks(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                const std::vector<stoat::Path_traversal_t>& paths) const; 
    

        // Do we want to analyze this snarl, based on the various limits we were given?
        bool snarl_is_eligible( const bdsg::SnarlDistanceIndex& distance_index, const handlegraph::net_handle_t& snarl, bool check_distances) const; 
};
}

#endif
