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
// The latter three may be empty if they were not filled in by the SnarlDataCollection
struct snarl_info_t {
    public:

        // Constructor from elements
        snarl_info_t(stoat::Node_traversal_t start_node, stoat::Node_traversal_t end_node, std::string ref_path, 
                     size_t start_position, size_t end_position, size_t depth, const std::string& variant_type, 
                     const std::vector<Path_traversal_t>& walks_by_allele, const std::vector<std::set<sample_hap_t>>& sample_sets_by_allele, 
                     const std::vector<std::string>& sequences_by_allele) :

                     start_node(start_node), end_node(end_node), 
                     ref_path(ref_path), start_position(start_position), end_position(end_position), depth(depth),
                     variant_type(variant_type), walks_by_allele(walks_by_allele), sample_sets_by_allele(sample_sets_by_allele), sequences_by_allele(sequences_by_allele)
                     {};

        // Start and end nodes, both pointing into the snarl
        stoat::Node_traversal_t start_node;
        stoat::Node_traversal_t end_node;

        // The reference chromosome/path
        std::string ref_path; 

        // Start and end offset along the reference path
        size_t start_position;
        size_t end_position;
        size_t depth;

        // The "variant type" of the snarl, which represents the min/max length (or 0 for a deletion) of each walk in walks_by_allele
        const std::string& variant_type;

        // For each allele, the walk through the snarl
        // The walk includes the boundary nodes of the snarl
        // Nested chains are included in the walk as the boundary node of the chain going into the chain, an empty node 0 going forward, 
        //    and the other bound of the chain going out.
        // When a walk leaves the snarl and comes back, it will include the boundary node of the snarl leaving it, an empty node 0 going backward, and the 
        //    boundary of the snarl going back in
        const std::vector<Path_traversal_t>& walks_by_allele;

        // For each allele, the set of sample/haplotypes that contain that allele
        const std::vector<std::set<sample_hap_t>>& sample_sets_by_allele;

        // For each allele, what is its sequence?
        // The sequences don't include sequences of the boundary nodes
        const std::vector<std::string>& sequences_by_allele;

};

/// A class for holding per-snarl data for a collection of snarls
/// Publicly, this allow access to snarl_info_t's for holding basic information about the snarl, the walks through the snarl 
///    the groups of haplotypes following each walk, and the sequence of each walk.
/// The latter three optional and may be empty if the SnarlDataCollection was built without them.
/// Internally, the basic snarl information, the walks, haplotype samples sets, and sequences are each stored separately. 
/// Build the collection either from the distance index or by loading from a file
/// Access snarls by chromosome name or all at once

class SnarlDataCollection {



    public:
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////// Public functions for building and interrogating snarl info

        /// Make a SnarlDataCollection with limits on which snarls to include
        /// Ignore snarls whose maximum length is less than allele_size_limit
        /// Ignore snarls with more children than snarl_child_limit
        /// Ignore snarls if traversing the paths takes more than walk_steps_limit steps
        SnarlDataCollection(size_t allele_size_limit, size_t snarl_child_limit, size_t walk_steps_limit);

        /// Fill in the SnarlDataCollection for all snarls in the distance index
        /// If sample_set_requested is true, then call find_sample_sample_sets and save the output sample sets.
        /// If walks_requested is true, then find the walks using find_walks.
        /// If sequence_requested is true, then find the sequence of each walk.
        /// The walks, sample sets, and sequences must all match each other, and the walks may be dependent on the sample sets or vice versa 
        /// find_sample_sets_first is true if we want to find the sample sets first and then walks based on the sample sets, and false to do the opposite.
        ///    If only one or the other of sample sets and walks is requested then find_sample_sets_first doesn't matter.
        /// find_sample_sample_sets and find_walks must be check the walks or sample sets accordingly to make sure that they match.
        /// Sequences are always found from the walks and cannot be found if walks_requested is false.
        /// The SnarlDataCollection provides default implementations get_all_walks_through_snarl and get_walks_from_sample_sets that may be used for find_walks
        /// If reference_sample is not empty, get coordinates on this reference path
        /// Since the distance index may not contain distances, use check_distances=false to skip distance checking
        /// Write failed snarls to out_fail
        void fill_in_snarl_info(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                                bool find_sample_sets_first,
                                bool walks_requested,
                                const std::function<void(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index, 
                                                         const net_handle_t& snarl, const snarl_info_t& snarl_data, 
                                                         std::vector<Path_traversal_t>& walks, 
                                                         std::vector<std::string>& walk_lengths)>& find_walks,
                                bool sample_set_requested,
                                const std::function<void(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                                                          const net_handle_t& snarl, const snarl_info_t& snarl_data, 
                                                          std::vector<std::set<sample_hap_t>>& sample_sets_by_allele)>& find_sample_sets,
                                bool sequence_requested, 
                                const string& reference_sample, bool check_distances);

        
        /// If the snarl sample_sets_by_allele (which assigns sample/haplotypes to each snarl_walk) were not found during construction, go through
        /// all snarls and call find_sample_sample_sets to add sample_sets_by_allele
        /// find_sample_sample_sets takes a snarl_info_t (which is not guaranteed to contain sample_sets_by_allele or sequences) and returns the sample_sets_by_allele
        /// corresponding to the walks_by_allele in snarl_data. This should return the sample_sets_by_allele field of the snarl_info_t, but since the snarl_info_t
        /// doesn't exist in the SnarlDataCollection it is immutable and the sample_sets_by_allele must be returned and saved separately 
        /// Note that this will overwrite any existing sample_sets_by_allele. 
        /// If chr is not empty, run this only on snarls on the given chromosome (as reference path name)
        //TODO:  I think this should be fine to run multithreaded as long as the list of snarl data doesn't move around, and I think the object
        //        stores a reference to the vector somewhere else in memory
        void add_snarl_sample_sets(std::unordered_map<stoat::sample_hap_t, size_t>& sample_haplotype_to_index, 
                                  const std::function<std::vector<std::set<sample_hap_t>>(const snarl_info_t& snarl_data)>& find_sample_sets,
                                  std::string chr="");

        /// Run iteratee for all snarls
        void for_each_snarl(const std::function<void(const snarl_info_t& snarl_info)>& iteratee) const;

        /// Write the collection of snarls to the given file
        void write_snarl_data_collection(std::ostream& outstream) const;
        
        /// Load the collection of snarls from the given file
        /// Warn if the allele_size_limit or snarl_child_limit of the file are less permissive than this SnarlDataCollection
        void load_snarl_data_collection(std::istream& instream); 


        ///////////////////////////// Helper functions for use in add_snarl_sample_sets

        /// TODO: The walks must include the start and end bounds of the snarl. Could be taken out 
        /// TODO: This should probably not be a member function but I'm leaving it here for now instead of changing Matis's code to use my data types
        /// Helper function for finding all possible walks through the snarl. Fills in walks
        /// snarl_data will not have the sample_sets_by_allele filled in
        /// If a path cycles more than walk_cycle_limit times, stop looking for more cycles
        static void get_all_walks_through_snarl(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index, 
                           const net_handle_t& snarl, const snarl_info_t& snarl_data, std::vector<stoat::Path_traversal_t>& walks,
                           std::vector<std::string>& walk_lengths, size_t walk_cycle_limit = 1);

        /// Helper function for finding all possible walks through the snarl. Fills in walks
        /// snarl_data is assumed to have the sample_sets_by_allele filled in and walks must be filled in to match the sample_sets_by_allele
        static void get_walks_from_sample_sets(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index, 
                           const net_handle_t& snarl, const snarl_info_t& snarl_data, std::vector<stoat::Path_traversal_t>& walks,
                           std::vector<std::string>& walk_lengths);


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

            // The "variant type" of the snarl, which represents the min/max length (or 0 for a deletion) of each walk in walks_by_allele
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

        /// Map snarl (as the start node, which uniquely identifies the snarl) to the sample_sets_by_allele.
        /// The vector follows the vector of walks in snarl_to_walks, meaning that each set in snarl_to_sample_sets
        ///  represents the sample/haplotype (as an index into sample_haplotypes) that take the walk at the same
        ///  index in the snarl's entry in snarl_to_walks. 
        std::unordered_map<stoat::Node_traversal_t, std::vector<std::set<size_t>>> snarl_to_sample_sets;

        /// Map snarl (as the start node, which uniquely identifies the snarl) to the sequences of the walks.
        /// Each string in the vector is the sequence for the walk in the corresponding vector in snarl_to_walks. 
        std::unordered_map<stoat::Node_traversal_t, std::vector<std::string>> snarl_to_sequences;


        ///////////////// Lists of strings and stuff that are stored as indexes in the real data structures instead of duplicating them a bunch

        /// The string representations of reference paths
        /// Used by the snarl_info_internal_t's reference_index
        std::vector<std::string> reference_names;

        /// This stores all sample_hap_t's that are stored as indexes by snarl_to_sample_sets
        vector<sample_hap_t> sample_haplotypes;


        //////////////////////////// Extra housekeeping stuff

        /// This goes at the beginning of the file to ensure that it is the right file type and version
        inline const static std::string file_header = "#SNARL_DATA_v1.0";

        /// Skip snarls if their maximum length is smaller than this
        size_t allele_size_limit;

        /// Skip snarls if they have more children than this
        size_t snarl_child_limit;

        /// Don't include snarls if enumerating all its walks takes more than this many steps
        size_t walk_steps_limit;


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////// Private functions
    private:
    
        // Given the walks through the snarl, find the sequence. The sequence will just be a concatination of sequences of nodes, ignoring anything else
        std::vector<std::string> get_sequences_from_walks(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                const std::vector<stoat::Path_traversal_t>& paths) const; 
    

        // Do we want to analyze this snarl, based on the various limits we were given?
        bool snarl_is_eligible( const bdsg::SnarlDistanceIndex& distance_index, const handlegraph::net_handle_t& snarl, bool check_distances) const; 

        // Helper function to add the given sample sample_sets_by_allele to the collection as integers.
        // Thread safe with other functions modifying the collection
        void add_sample_sets_to_collection(std::unordered_map<stoat::sample_hap_t, size_t>& sample_haplotype_to_index, const std::vector<std::set<sample_hap_t>>& sample_sets, const Node_traversal_t& snarl_start);
};
}

#endif
