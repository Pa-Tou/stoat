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

/// A class for holding per-snarl data for a collection of snarls
/// Publicly, this allow access to snarl_data_t's for holding basic information about the snarl, the paths through the snarl 
///    the groups of haplotypes following each path, and the sequence of each path.
/// The latter two are both optional and may be empty if the SnarlDataCollection was built without them.
/// Internally, the basic snarl information, the paths, haplotype partitions, and sequences are each stored separately. 
/// Build the collection either from the distance index or by loading from a file
/// Access snarls by chromosome name or all at once

class SnarlDataCollection {



    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////// Public struct used for storing snarl information
    public:

    // This holds snarl information from the distance index and graph: id of the snarl, reference coordinates, walks, 
    // the sequences of the walks and sets samples/haplotypes taking each walk
    // The latter two may be empty
    struct snarl_data_t {
        public:
    
            // Constructor from elements
            snarl_data_t(stoat::Node_traversal_t start_node, stoat::Node_traversal_t end_node, std::string ref_path, 
                         size_t start_position, size_t end_position, size_t depth, const std::string& variant_type, 
                         const std::vector<Path_traversal_t>& snarl_paths, const std::vector<std::set<sample_hap_t>>& partitions, 
                         const std::vector<std::string>& sequences) :

                         start_node(start_node), end_node(end_node), 
                         ref_path(ref_path), start_position(start_position), end_position(end_position), depth(depth),
                         variant_type(variant_type), snarl_paths(snarl_paths),  partitions(partitions)
                         {};
    
            // Start and end nodes, both pointing into the snarl
            stoat::Node_traversal_t start_node;
            stoat::Node_traversal_t end_node;
    
            std::string ref_path; //TODO: I think this could get pretty big, might want to save it as an index into a list of reference paths
            // Start and end offset along the reference path
            size_t start_position;
            size_t end_position;
            size_t depth;

            // The "variant type" of the snarl, which represents the min/max length (or 0 for a deletion) of each path in snarl_paths
            const std::string& variant_type;
 
            // All possible walks through the snarl
            const std::vector<Path_traversal_t>& snarl_paths;

            // For each walk in snarl_paths, which sample/haplotypes take that path
            const std::vector<std::set<sample_hap_t>>& partitions;

            // For each walk in snarl_paths, what is its sequence?
            const std::vector<std::string>& sequences;

    };

    public:
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////// Public functions for building and interrogating snarl info

        /// Fill in the SnarlDataCollection for all snarls in the distance index
        /// If partition_requested is true, then call find_sample_partitions and save the output partitions.
        /// If sequence_requested is true, then find the sequence of each walk 
        SnarlDataCollection(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                            size_t allele_size_limit, size_t snarl_child_limit, 
                            bool partition_requested,
                            const std::function<std::vector<std::set<sample_hap_t>>(const net_handle_t& snarl, const std::vector<Path_traversal_t>& paths)>& find_sample_partitions,
                            bool sequence_requested);

        
        /// If the snarl partitions (which assigns sample/haplotypes to each snarl_path) were not found during construction, go through
        /// all snarls and call find_sample_partitions to add partitions
        /// find_sample_partitions takes a snarl_data_t (which is not guaranteed to contain partitions or sequences) and returns the partitions
        /// corresponding to the snarl_paths in snarl_data. This should return the partitions field of the snarl_data_t, but since the snarl_data_t
        /// doesn't exist in the SnarlDataCollection it is immutable and the partitions must be returned and saved separately 
        /// Note that this will overwrite any existing partitions. 
        /// If chr is specified, run this only on snarls on the given chromosome (as reference path name)
        //TODO:  I think this should be fine to run multithreaded as long as the list of snarl data doesn't move around, and I think the object
        //        stores a reference to the vector somewhere else in memory
        void add_snarl_partitions(const std::function<std::vector<std::set<sample_hap_t>>(const snarl_data_t& snarl_data,
                                                                                          const std::vector<Path_traversal_t>& paths)>& find_sample_partitions,
                                  std::string chr="");

        /// Run iteratee for all snarls
        void for_each_snarl(const std::function<void(const snarl_data_t& snarl_info)>& iteratee) const;

        /// Write the collection of snarls to the given file
        void write_snarl_data_collection(const string& filename) const;
        
        /// Load the collection of snarls from the given file
        /// Warn if the allele_size_limit or snarl_child_limit of the file are less permissive than this SnarlDataCollection
        void load_snarl_data_collection(const string& filename); 

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////// Private data members

    private:

        /// This stores the basic information from the snarl_data_t
        struct snarl_data_internal_t {
    
            // Start and end nodes, both pointing into the snarl
            stoat::Node_traversal_t start_node;
            stoat::Node_traversal_t end_node;

            // Index into reference_names to get the string representation of the reference path
            size_t reference_index;
    
            // Start and end offset along the reference path
            size_t start_position;
            size_t end_position;

            // Depth of the snarl in the snarl decomposition
            size_t depth;

            // The "variant type" of the snarl, which represents the min/max length (or 0 for a deletion) of each path in snarl_paths
            std::string variant_type;

        };

        //////////////////////////// The stuff holding the data 

        /// This holds the snarl data as a map from the chromosome name to the data
        //TODO: Make sure that this gets the chr name the way Matis did it
        //TODO: idk if I want it to be a map from chr to vector of snarl data or just a vector and then check the chromosome for the per-chromosome calls 
        //std::unordered_map<std::string, std::vector<snarl_data_internal_t>> chr_to_snarl_data;
        std::vector<snarl_data_internal_t> all_snarl_data;


        /// Map snarl (as the start node, which uniquely identifies the snarl) to the paths through the snarl.
        std::unordered_map<stoat::Node_traversal_t, std::vector<Path_traversal_t>> snarl_to_paths;

        /// Map snarl (as the start node, which uniquely identifies the snarl) to the partitions.
        /// The vector follows the vector of paths in snarl_to_paths, meaning that each set in snarl_to_partitions
        ///  represents the sample/haplotype (as an index into sample_haplotypes) that take the path at the same
        ///  index in the snarl's entry in snarl_to_paths. 
        std::unordered_map<stoat::Node_traversal_t, std::vector<std::set<size_t>>> snarl_to_partitions;

        /// Map snarl (as the start node, which uniquely identifies the snarl) to the sequences of the paths.
        /// Each string in the vector is the sequence for the path in the corresponding vector in snarl_to_paths. 
        std::unordered_map<stoat::Node_traversal_t, std::vector<std::string>> snarl_to_sequences;


        ///////////////// Lists of strings and stuff that are stored as indexes in the real data structures instead of duplicating them a bunch

        /// The string representations of reference paths
        /// Used by the snarl_data_internal_t's reference_index
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


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////// Private functions
    private:
        std::vector<stoat::Path_traversal_t> get_snarl_paths(graph, distance_index, snarl) const;
    
};
}

#endif
