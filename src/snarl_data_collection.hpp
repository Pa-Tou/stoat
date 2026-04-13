#ifndef STOAT_SNARL_DATA_COLLECTION_HPP_INCLUDED
#define STOAT_SNARL_DATA_COLLECTION_HPP_INCLUDED

#include <handlegraph/path_position_handle_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include "types_and_structs.hpp"
#include "writer.hpp"
#include "vcf_parser.hpp"

namespace stoat {

/// A class for holding per-snarl data for a collection of snarls
/// Publicly, this allow access to snarl_info_t's for holding basic information about the snarl, the genotypes, the walks through the snarl 
/// the groups of haplotypes following each walk, and the sequence of each walk.
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

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////// Public functions for building and interrogating snarl info
    public:

        /// Make a SnarlDataCollection with limits on which snarls to include
        /// Ignore snarls whose maximum length is less than allele_size_limit
        /// Ignore snarls with more children than snarl_child_limit
        /// Ignore snarls if traversing the paths takes more than walk_steps_limit steps
        SnarlDataCollection(size_t allele_size_limit, size_t snarl_child_limit, size_t walk_steps_limit);

        /// Construct using already computed snarl_data_collection (used for removed sample case / phenotype not present)
        SnarlDataCollection(size_t allele_size_limit, size_t snarl_child_limit, size_t walk_steps_limit, std::unordered_map<std::string, size_t> sample_to_index);

        /// Fill in the SnarlDataCollection for all snarls in the distance index
        /// sample_haplotypes gets copied and kept around as all_sample_haplotypes. Fills in sample_to_index based on sample_haplotypes 
        /// If alleles_requested is true, then call find_alleles_by_sample to assign each sample/haplotype in all_sample_haplotypes to an allele
        /// find_alleles_by_sample must return a vector of length all_sample_haplotypes with the allele of each sample_hap_t (std::numeric_limits<size_t>::max() if not present).
        /// If walks_requested is true, then find the walks using find_walks.
        /// If sequence_requested is true, then find the sequence of each walk.
        /// The walks, sample sets, and sequences must all match each other, and the walks may be dependent on the sample sets or vice versa 
        /// find_alleles_first is true if we want to find the sample sets first and then walks based on the sample sets, and false to do the opposite.
        /// If only one or the other of sample sets and walks is requested then find_alleles_first doesn't matter.
        /// find_alleles_by_sample and find_walks must check the walks or sample sets accordingly to make sure that they match.
        /// Sequences are always found from the walks and cannot be found if walks_requested is false.
        /// The SnarlDataCollection provides default implementations get_all_walks_through_snarl and get_walks_from_alleles that may be used for find_walks
        /// If reference_samples is not empty, get coordinates on one of these reference path. If it is empty then the coordinates will be on any path
        /// Since the distance index may not contain distances, use check_distances=false to skip distance checking
        /// If out_filename is not empty, then write the snarl collection to a file as well. If keep_snarls is False, then delete the snarls as they are found
        /// instead of keeping them in the collection. This saves memory if writing the snarls. If out_filename is empty, then keep_snarls should be True

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
                                const std::unordered_set<std::string>& reference_samples, bool check_distances,
                                Writer& out_writer, bool keep_snarls);

        
        /// Use if the snarl allele_by_sample (which assigns sample/haplotypes to each snarl_walk) were not found during construction. Go through
        /// all snarls and call find_alleles_by_sample, which returns a vector of allele assignments for each sample_hap_t in all_sample_haplotypes.
        ///  find_alleles_by_sample() should return the allele_by_sample field of the snarl_info_t, but since the snarl_info_t
        /// doesn't exist in the SnarlDataCollection it is immutable and the allele_by_sample must be returned and saved separately 
        /// Note that this will overwrite any existing allele_by_sample. 
        /// If chr is not empty, run this only on snarls on the given chromosome (as reference path name)
        // TODO:  I think this should be fine to run multithreaded as long as the list of snarl data doesn't move around, and I think the object
        //       stores a reference to the vector somewhere else in memory
        void add_alleles_by_sample(const std::function<std::vector<size_t>(const snarl_info_t& snarl_data, 
                                                                           const std::vector<stoat::sample_hap_t>& all_sample_haplotypes)>& find_alleles_by_sample,
                                   std::string chr);

        /// Run iteratee for all snarls
        void for_each_snarl(const std::function<void(snarl_info_t& snarl_info)>& iteratee) const;

        /// Starting from an empty SnarlDataCollection, run iteratee for all snarls, but one line at a time from a file instead of loading the whole thing into memory
        /// This will load the header but not keep any of the snarls in the SnarlDataCollection
        void for_each_snarl_in_file(stoat::Reader& in_reader, 
            const std::function<void(snarl_info_t& snarl_info)>& iteratee);

        /// Write the collection of snarls to the given file
        void write_snarl_data_collection(Writer& out_writer) const;
        
        /// Load the collection of snarls from the given file
        /// Warn if the allele_size_limit or snarl_child_limit of the file are less permissive than this SnarlDataCollection
        /// also a mode to load just the header used to reuse the same sample_to_index for other objects and then run snarl file line by line (although it means loading the sample_to_index map twice technically)
        void load_snarl_data_collection(stoat::Reader& in_reader, const bool header_only = false); 

        std::unordered_map<std::string, size_t> get_sample_to_index_copy() const;
        std::unordered_map<std::string, size_t>& get_sample_to_index_reference();

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
    
    // fill in the genotypes of the snarls based on the edge matrix built on a VCF stream
    // the VCF is read and parsed by chromosome
    // The vcf parser is assumed to have loaded the header and be pointing to the start of the actual records
    void genotype_snarls_by_chr_from_vcf(std::vector<std::string>& sample_names, stoat_vcf::VCFParser& vcf_parser);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////// Private data members
    private:

        /// This stores the basic information from the snarl_info_t
        struct snarl_info_internal_t {
    
            // Start and end nodes, both pointing into the snarl
            stoat::node_traversal_t start_node = stoat::node_traversal_t(0, false);
            stoat::node_traversal_t end_node = stoat::node_traversal_t(0, false);

            // Index into reference_names to get the string representation of the reference path
            size_t reference_index;
    
            // Start and end offset along the reference path
            size_t start_position;
            size_t end_position;

            // Depth of the snarl in the snarl decomposition
            size_t depth;

        };

        //////////////////////////// The stuff holding the data for each snarl, indexed by snarl start node (which uniquely identifies the snarl)

        /// This holds the snarl data as a map from the chromosome name to the data
        // TODO: Make sure that this gets the chr name the way Matis did it
        // TODO: idk if I want it to be a map from chr to vector of snarl data or just a vector and then check the chromosome for the per-chromosome calls 
        //std::unordered_map<std::string, std::vector<snarl_info_internal_t>> chr_to_snarl_data;
        std::vector<snarl_info_internal_t> all_snarl_data;

        /// Map snarl (as the start node, which uniquely identifies the snarl) to the walks through the snarl.
        /// Each entry in the vector is considered to be one allele.
        std::unordered_map<stoat::node_traversal_t, std::vector<PathTraversal>> snarl_to_walks;

        /// Map snarl (as the start node, which uniquely identifies the snarl) to the alleles_by_sample vector.
        /// So each entry in the alleles_by_sample vector is the allele number of one sample_hap_t. The allele number
        /// corresponds to the index of an allele in snarl_to_walks vectors, and snarl_to_sequences vectors
        std::unordered_map<stoat::node_traversal_t, allele_by_sample_t> snarl_to_alleles_by_sample;

        /// Map snarl (as the start node, which uniquely identifies the snarl) to the sequences of the walks.
        /// Each string in the vector is the sequence for the walk in the corresponding vector in snarl_to_walks. 
        std::unordered_map<stoat::node_traversal_t, std::vector<std::string>> snarl_to_sequences;

        ///////////////// Lists of strings and stuff that are stored as indexes in the real data structures instead of duplicating them a bunch

        /// The string representations of reference paths
        /// Used by the snarl_info_internal_t's reference_index
        std::vector<std::string> reference_names;

        /// This stores all sample_hap_t's that are stored as indexes by snarl_to_alleles_by_sample
        /// It is given by fill_in_snarl_info or loaded from the file
        std::vector<stoat::sample_hap_t> all_sample_haplotypes;

        // This maps all sample names from all_sample_haplotypes to a unique identifier, (which is used as an index into a vector 
        // so it must start from 0 and go up to the number of samples-1).
        // Note that since there may be more sample_hap_t's than samples, the indices don't necessarily match
        // This gets made based on all_sample_halpotypes
        std::unordered_map<std::string, size_t> sample_to_index;

        std::set<size_t> removed_sample_idx;

        //////////////////////////// Extra housekeeping stuff

        /// This goes at the beginning of the file to ensure that it is the right file type and version
        inline const static std::string file_header = "#SNARL_DATA_v1.0";

        /// Skip snarls if their maximum length is smaller than this
        size_t allele_size_limit;

        /// Skip snarls if they have more children than this
        size_t snarl_child_limit;

        /// Don't include snarls if enumerating all its walks takes more than this many steps
        size_t walk_steps_limit;

        /// These are just for logging purposes to keep track of snarls info
        size_t number_snarl_limit_distance;
        size_t number_snarl_limit_children;
        size_t number_snarl_analyzed;
        size_t number_paths_analyzed;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////// Private functions
    private:
    
        // Given the walks through the snarl, find the sequence. The sequence will just be a concatination of sequences of nodes, ignoring anything else
        std::vector<std::string> get_sequences_from_walks(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                const std::vector<stoat::PathTraversal>& paths) const; 
    
        // Do we want to analyze this snarl, based on the various limits we were given?
        bool snarl_is_eligible(const bdsg::SnarlDistanceIndex& distance_index, const handlegraph::net_handle_t& snarl, bool check_distances); 

        //////////////////// Helper functions for writing and loading stuff from files

        /// Write just the header of the collection of snarls to the given file
        void write_snarl_data_collection_header(Writer& out_writer) const;

        /// Write just one snarl from the collection
        void write_snarl_data_line(Writer& out_writer, const snarl_info_internal_t& snarl_info) const;
        void write_snarl_data_line(Writer& out_writer, const snarl_info_internal_t& snarl_info, const std::vector<stoat::PathTraversal>* walks_by_allele, 
                                    const std::vector<std::string>* snarl_sequences, const allele_by_sample_t* alleles_by_sample) const;

        /// Given a stream to the start of the file, load just the header. The stream will be advanced to point to the beginning of the snarl records
        void load_snarl_data_collection_header(stoat::Reader& in_reader);

        /// Given a string representing a line in the file, load one snarl_info_internal_t
        /// This assumes that load_snarl_collection_header() has already been called
        snarl_info_internal_t load_snarl_data_line(std::string& line);

        // Helper function for for_each_snarl() to go from internal snarl info to running iteratee on the snarl_info_t
        // Note that this can't be a function to return the snarl_info_t because it has references to tables and stuff that would go out of scope
        void run_iteratee_on_one_snarl(const snarl_info_internal_t& internal_snarl_info, const std::function<void(snarl_info_t& snarl_info)>& iteratee) const;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////// Public functions
    public:
        // Return true if the two SnarlDataCollections contain the same information (but possibly in a different order)
        // Otherwise, throw an error
        // This is very inefficient and should only be used for testing small examples
        static bool is_equivalent(const SnarlDataCollection& collection1, const SnarlDataCollection& collection2); 

};
}

#endif
