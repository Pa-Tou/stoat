#ifndef STOAT_SNARL_DATA_COLLECTION_HPP_INCLUDED
#define STOAT_SNARL_DATA_COLLECTION_HPP_INCLUDED

#include <iostream>
#include <handlegraph/path_position_handle_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include "utils.hpp"
#include "snarl_data_t.hpp"
#include "feature_tables.hpp"

using namespace std;
using namespace stoat;

namespace stoat {

// Used to store the allele assignment for a vector of samples. 
// This struct only makes sense when it is paired with a vector of samples or sample_hap_t's (as a member of a snarl_info_t). 
// Each entry in alleles is the identifier of an allele
// The allele identifiers are also used as indices into other per-allele vectors, so they should start from 0 and go up to the number of alleles-1.
// allele_count is the number of non-inf alleles, or 1 + max value in alleles_by_sample. This is not enforced though
struct allele_by_sample_t{
    size_t allele_count;
    std::vector<size_t> alleles;

    allele_by_sample_t(size_t count, std::vector<size_t> alleles) :
        allele_count(count),
        alleles(std::move(alleles)) {}
    allele_by_sample_t() :
        allele_count(0), alleles(0) {}
};

// This holds snarl information from the distance index and graph: id and reference coordinates of the snarl, genotypes, reference coordinates, walks, 
// the sequences of the walks and sets samples/haplotypes taking each walk
// The latter four may be empty if they were not filled in by the SnarlDataCollection
struct snarl_info_t {
    public:

        // Constructor from elements
        snarl_info_t(stoat::Node_traversal_t start_node, stoat::Node_traversal_t end_node, std::string ref_path, 
                     size_t start_position, size_t end_position, size_t depth,
                     const GenotypeTable& genotypes, 
                     const std::vector<stoat::sample_hap_t>& all_sample_haplotypes,
                     const allele_by_sample_t& alleles_by_sample,
                     const std::vector<PathTraversal>& walks_by_allele, 
                     const std::vector<std::string>& sequences_by_allele) :

                     start_node(start_node), end_node(end_node), 
                     ref_path(ref_path), start_position(start_position), end_position(end_position), depth(depth),
                     genotypes(genotypes), all_sample_haplotypes(all_sample_haplotypes), alleles_by_sample(alleles_by_sample),
                     walks_by_allele(walks_by_allele), sequences_by_allele(sequences_by_allele)
                     {};

        // Start and end nodes, both pointing into the snarl
        stoat::Node_traversal_t start_node;
        stoat::Node_traversal_t end_node;

        // The reference chromosome/path
        std::string ref_path; 

        // Start and end offset along the reference path
        size_t start_position;
        size_t end_position;

        // The depth of the snarl in the snarl tree
        size_t depth;

        // A genotype table of the counts of each allele for each sample (see feature_table.hpp).
        // Alleles have the same numbering as walks_by_allele and sequences_by_allele
        const GenotypeTable& genotypes;

        // This stores all the sample/haplotypes 
        const std::vector<stoat::sample_hap_t>& all_sample_haplotypes;

        // For each haplotype in all_sample_haplotypes, an assignment to an allele number
        const allele_by_sample_t& alleles_by_sample;

        // For each allele, the walk through the snarl
        // The walk includes the boundary nodes of the snarl
        // Nested chains are included in the walk as the boundary node of the chain going into the chain, an empty node 0 going forward, 
        //    and the other bound of the chain going out.
        // When a walk leaves the snarl and comes back, it will include the boundary node of the snarl leaving it, an empty node 0 going backward, and the 
        //    boundary of the snarl going back in
        const std::vector<PathTraversal>& walks_by_allele;

        // For each allele, what is its sequence?
        // The sequences don't include sequences of the boundary nodes
        const std::vector<std::string>& sequences_by_allele;

};

/// A class for holding per-snarl data for a collection of snarls
/// Publicly, this allow access to snarl_info_t's for holding basic information about the snarl, the genotypes, the walks through the snarl 
///    the groups of haplotypes following each walk, and the sequence of each walk.
/// The latter four are optional and may be empty if the SnarlDataCollection was built without them.
/// Internally, the basic snarl information, the walks, haplotype samples sets, and sequences are each stored separately. 
/// Build the collection either from the distance index or by loading from a file
/// Access snarls by chromosome name or all at once
///
/// The constructor only fills in constraints about which snarls are included in the collection.
/// The collection can be filled in from a file using load_snarl_data_collection()
/// Otherwise, it can be filled in by calling fill_in_snarl_info(). This function has different options for which information is filled in:
/// walks_requested, alleles_requested, sequence_requested 
/// stoat graph fills in everything at the same time, finding walks and sequences based on the alleles
/// stoat vcf fills in the walks first, saves the intermediate file, then fills in the alleles. Sequences are generally not included.
/// This can be done by calling fill_in_snarl_info() with alleles_requested and sequence_requested set to false.
/// Alleles can be filled in later by calling add_alleles_by_sample()
///
/// Once the collection is filled in with snarls, the data can be accessed by calling for_each_snarl()

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
        /// sample_haplotypes gets copied and kept around as all_sample_haplotypes. Fill in sample_to_index based on sample_haplotypes 
        /// If alleles_requested is true, then call find_alleles_by_sample to assign each sample/haplotype in all_sample_haplotypes to an allele
        /// find_alleles_by_sample must return a vector of length all_sample_haplotypes with the allele of each sample_hap_t (std::numeric_limits<size_t>::max() if not present).
        /// If walks_requested is true, then find the walks using find_walks.
        /// If sequence_requested is true, then find the sequence of each walk.
        /// The walks, sample sets, and sequences must all match each other, and the walks may be dependent on the sample sets or vice versa 
        /// find_alleles_first is true if we want to find the sample sets first and then walks based on the sample sets, and false to do the opposite.
        ///    If only one or the other of sample sets and walks is requested then find_alleles_first doesn't matter.
        /// find_alleles_by_sample and find_walks must check the walks or sample sets accordingly to make sure that they match.
        /// Sequences are always found from the walks and cannot be found if walks_requested is false.
        /// The SnarlDataCollection provides default implementations get_all_walks_through_snarl and get_walks_from_alleles that may be used for find_walks
        /// If reference_samples is not empty, get coordinates on one of these reference path. If it is empty then the coordinates will be on any path
        /// Since the distance index may not contain distances, use check_distances=false to skip distance checking
        /// Write failed snarls to out_fail
        void fill_in_snarl_info(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                                const std::vector<stoat::sample_hap_t>& sample_haplotypes,
                                bool find_alleles_first,
                                bool walks_requested,
                                const std::function<void(const net_handle_t& snarl, const snarl_info_t& snarl_data, 
                                                         std::vector<PathTraversal>& walks)>& find_walks,
                                bool alleles_requested,
                                const std::function<std::vector<size_t>(const net_handle_t& snarl, const snarl_info_t& snarl_data, 
                                                                        const std::vector<stoat::sample_hap_t>& all_sample_haplotypes)>& find_alleles_by_sample,
                                bool sequence_requested, 
                                const std::unordered_set<string>& reference_samples, bool check_distances);

        
        /// Use if the snarl allele_by_sample (which assigns sample/haplotypes to each snarl_walk) were not found during construction. Go through
        /// all snarls and call find_alleles_by_sample, which returns a vector of allele assignments for each sample_hap_t in all_sample_haplotypes.
        ///  find_alleles_by_sample() should return the allele_by_sample field of the snarl_info_t, but since the snarl_info_t
        /// doesn't exist in the SnarlDataCollection it is immutable and the allele_by_sample must be returned and saved separately 
        /// Note that this will overwrite any existing allele_by_sample. 
        /// If chr is not empty, run this only on snarls on the given chromosome (as reference path name)
        //TODO:  I think this should be fine to run multithreaded as long as the list of snarl data doesn't move around, and I think the object
        //        stores a reference to the vector somewhere else in memory
        void add_alleles_by_sample(const std::function<std::vector<size_t>(const snarl_info_t& snarl_data, 
                                                                           const std::vector<stoat::sample_hap_t>& all_sample_haplotypes)>& find_alleles_by_sample,
                                   std::string chr);

        /// Run iteratee for all snarls
        void for_each_snarl(const std::function<void(const snarl_info_t& snarl_info)>& iteratee) const;

        /// Write the collection of snarls to the given file
        void write_snarl_data_collection(std::ostream& outstream, const bool output_samples = true) const;
        
        /// Load the collection of snarls from the given file
        /// Warn if the allele_size_limit or snarl_child_limit of the file are less permissive than this SnarlDataCollection
        void load_snarl_data_collection(std::string& filename); 

        size_t size() const {return all_snarl_data.size();}

        // Get a reference to the reference path names that are stored in the collection
        const std::vector<std::string>& get_reference_names() const {
            return reference_names;
        };

        ///////////////////////////// Helper functions for use in add_alleles_by_sample

        /// TODO: The walks must include the start and end bounds of the snarl. Could be taken out 
        /// TODO: This should probably not be a member function but I'm leaving it here for now instead of changing Matis's code to use my data types (from write_snarls_with_paths() in snarl_data_t.hpp)
        /// Helper function for finding all possible walks through the snarl. Fills in walks
        /// snarl_data will not have the allele_by_sample filled in
        /// If a path cycles more than walk_cycle_limit times, stop looking for more cycles
        static void get_all_walks_through_snarl(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index, 
                           const net_handle_t& snarl, const snarl_info_t& snarl_data, std::vector<stoat::PathTraversal>& walks,
                           size_t walk_cycle_limit = 1);

        /// Helper function for finding walks through the snarl. Fills in walks
        /// the collection is assumed to have the allele_by_sample filled in and walks must be filled in to match the allele_by_sample
        /// This requires a PathPositionHandleGraph to make sure that multiple traversals of the snarl are properly ordered
        static void get_walks_from_alleles(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index, 
                           const net_handle_t& snarl,  const snarl_info_t& snarl_data, std::vector<stoat::PathTraversal>& walks);


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

        };

        //////////////////////////// The stuff holding the data 


        /// This holds the snarl data as a map from the chromosome name to the data
        //TODO: Make sure that this gets the chr name the way Matis did it
        //TODO: idk if I want it to be a map from chr to vector of snarl data or just a vector and then check the chromosome for the per-chromosome calls 
        //std::unordered_map<std::string, std::vector<snarl_info_internal_t>> chr_to_snarl_data;
        std::vector<snarl_info_internal_t> all_snarl_data;


        /// Map snarl (as the start node, which uniquely identifies the snarl) to the walks through the snarl.
        /// Each entry in the vector is considered to be one allele.
        std::unordered_map<stoat::Node_traversal_t, std::vector<PathTraversal>> snarl_to_walks;


        /// Map snarl (as the start node, which uniquely identifies the snarl) to the alleles_by_sample vector.
        /// So each entry in the alleles_by_sample vector is the allele number of one sample_hap_t. The allele number
        /// corresponds to the index of an allele in snarl_to_walks vectors, and snarl_to_sequences vectors
        std::unordered_map<stoat::Node_traversal_t, allele_by_sample_t> snarl_to_alleles_by_sample;

        /// Map snarl (as the start node, which uniquely identifies the snarl) to the sequences of the walks.
        /// Each string in the vector is the sequence for the walk in the corresponding vector in snarl_to_walks. 
        std::unordered_map<stoat::Node_traversal_t, std::vector<std::string>> snarl_to_sequences;


        ///////////////// Lists of strings and stuff that are stored as indexes in the real data structures instead of duplicating them a bunch

        /// The string representations of reference paths
        /// Used by the snarl_info_internal_t's reference_index
        std::vector<std::string> reference_names;

        /// This stores all sample_hap_t's that are stored as indexes by snarl_to_alleles_by_sample
        /// It is given by fill_in_snarl_info or loaded from the file
        vector<stoat::sample_hap_t> all_sample_haplotypes;

        // This maps all sample names from all_sample_haplotypes to a unique identifier, (which is used as an index into a vector 
        // so it must start from 0 and go up to the number of samples-1).
        // Note that since there may be more sample_hap_t's than samples, the indices don't necessarily match
        // This gets made based on all_sample_halpotypes
        std::unordered_map<std::string, size_t> sample_to_index;


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
                const std::vector<stoat::PathTraversal>& paths) const; 
    

        // Do we want to analyze this snarl, based on the various limits we were given?
        bool snarl_is_eligible( const bdsg::SnarlDistanceIndex& distance_index, const handlegraph::net_handle_t& snarl, bool check_distances) const; 

    public:
        // Return true if the two SnarlDataCollections contain the same information (but possibly in a different order)
        // Otherwise, throw an error
        // This is very inefficient and should only be used for testing small examples
        static bool is_equivalent (const SnarlDataCollection& collection1, const SnarlDataCollection& collection2); 

};
}

#endif
