#ifndef NODE_DATA_COLLECTION_HPP
#define NODE_DATA_COLLECTION_HPP

#include "snarl_data_collection.hpp"
#include "types_and_structs.hpp"
#include "writer.hpp"
#include "vcf_parser.hpp"

namespace stoat {

class NodeDataCollection {

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////// Public functions for building and interrogating node info
    public:

        /// Make a NodeDataCollection with limits on which nodes to include
        /// Ignore nodes whose maximum length is less than allele_size_limit
        /// Ignore nodes with more children than node_child_limit
        /// Ignore nodes if traversing the paths takes more than walk_steps_limit steps
        NodeDataCollection();
        NodeDataCollection(std::unordered_map<std::string, size_t> sample_to_index);

        void fill_in_node_info(const handlegraph::PathPositionHandleGraph& graph, 
                        const bdsg::SnarlDistanceIndex& distance_index,
                        const std::vector<stoat::sample_hap_t>& sample_haplotypes,
                        bool find_alleles_first,
                        bool walks_requested,
                        const std::function<void(const net_handle_t& node, const node_info_t& node_data, 
                                                    std::vector<PathTraversal>& walks)>& find_walks,
                        bool alleles_requested,
                        const std::function<std::vector<size_t>(const net_handle_t& node, const node_info_t& node_data, 
                                                                const std::vector<stoat::sample_hap_t>& all_sample_haplotypes)>& find_alleles_by_sample,
                        bool sequence_requested, 
                        const std::unordered_set<std::string>& reference_samples, bool check_distances,
                        Writer& out_writer, bool keep_nodes);

        /// Use if the node allele_by_sample (which assigns sample/haplotypes to each node_walk) were not found during construction. Go through
        /// all nodes and call find_alleles_by_sample, which returns a vector of allele assignments for each sample_hap_t in all_sample_haplotypes.
        ///  find_alleles_by_sample() should return the allele_by_sample field of the node_info_t, but since the node_info_t
        /// doesn't exist in the NodeDataCollection it is immutable and the allele_by_sample must be returned and saved separately 
        /// Note that this will overwrite any existing allele_by_sample. 
        /// If chr is not empty, run this only on nodes on the given chromosome (as reference path name)
        // TODO:  I think this should be fine to run multithreaded as long as the list of node data doesn't move around, and I think the object
        //       stores a reference to the vector somewhere else in memory
        void add_alleles_by_sample(const std::function<std::vector<size_t>(const node_info_t& node_data, 
                                                                           const std::vector<stoat::sample_hap_t>& all_sample_haplotypes)>& find_alleles_by_sample,
                                   std::string chr);

        /// Run iteratee for all nodes
        void for_each_node(const std::function<void(node_info_t& node_info)>& iteratee) const;

        /// Starting from an empty NodeDataCollection, run iteratee for all nodes, but one line at a time from a file instead of loading the whole thing into memory
        /// This will load the header but not keep any of the nodes in the NodeDataCollection
        void for_each_node_in_file(stoat::Reader& in_reader, 
            const std::function<void(node_info_t& node_info)>& iteratee);

        /// Load the collection of nodes from the given file
        /// Warn if the allele_size_limit or node_child_limit of the file are less permissive than this NodeDataCollection
        /// also a mode to load just the header used to reuse the same sample_to_index for other objects and then run node file line by line (although it means loading the sample_to_index map twice technically)
        void load_node_data_collection(stoat::Reader& in_reader, const bool header_only = false); 

        std::unordered_map<std::string, size_t> get_sample_to_index_copy() const;
        std::unordered_map<std::string, size_t>& get_sample_to_index_reference();

        size_t size() const {return all_node_data.size();}

        // Get a reference to the reference path names that are stored in the collection
        const std::vector<std::string>& get_reference_names() const {
            return reference_names;
        };

    ///////////////////////////// Helper functions for use in add_alleles_by_sample

    static void get_all_walks_through_node(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index, 
                        const net_handle_t& node, const node_info_t& node_data, std::vector<stoat::PathTraversal>& walks,
                        size_t walk_cycle_limit = 1);

    /// Helper function for finding walks through the node. Fills in walks
    /// the collection is assumed to have the allele_by_sample filled in and walks must be filled in to match the allele_by_sample
    /// This requires a PathPositionHandleGraph to make sure that multiple traversals of the node are properly ordered
    // static void get_walks_from_alleles(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index, 
    //                     const net_handle_t& node,  const node_info_t& node_data, std::vector<stoat::PathTraversal>& walks);

    // fill in the genotypes of the nodes based on the edge matrix built on a VCF stream
    // the VCF is read and parsed by chromosome
    // The vcf parser is assumed to have loaded the header and be pointing to the start of the actual records
    void genotype_nodes_by_chr_from_vcf(std::vector<std::string>& sample_names, stoat_vcf::VCFParser& vcf_parser);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////// Private data members
    private:

        /// This stores the basic information from the node_info_t
        struct node_info_internal_t {

            // node info
            stoat::node_traversal_t node = stoat::node_traversal_t(0, false);

            // Index into reference_names to get the string representation of the reference path
            size_t reference_index;
    
            // Start and end offset along the reference path
            size_t position = 0;

            // Depth of the node in the node decomposition
            size_t depth = 0;

        };

        //////////////////////////// The stuff holding the data for each node, indexed by node start node (which uniquely identifies the node)

        std::vector<node_info_internal_t> all_node_data;

        /// Map node (as the start node, which uniquely identifies the node) to the alleles_by_sample vector.
        /// So each entry in the alleles_by_sample vector is the allele number of one sample_hap_t. The allele number
        /// corresponds to the index of an allele in node_to_walks vectors, and node_to_sequences vectors
        std::unordered_map<stoat::node_traversal_t, allele_by_sample_t> node_to_alleles_by_sample;

        std::unordered_map<stoat::node_traversal_t, std::string> node_to_sequences;

        ///////////////// Lists of strings and stuff that are stored as indexes in the real data structures instead of duplicating them a bunch

        /// The string representations of reference paths
        /// Used by the node_info_internal_t's reference_index
        std::vector<std::string> reference_names;

        /// This stores all sample_hap_t's that are stored as indexes by node_to_alleles_by_sample
        /// It is given by fill_in_node_info or loaded from the file
        std::vector<stoat::sample_hap_t> all_sample_haplotypes;

        // This maps all sample names from all_sample_haplotypes to a unique identifier, (which is used as an index into a vector 
        // so it must start from 0 and go up to the number of samples-1).
        // Note that since there may be more sample_hap_t's than samples, the indices don't necessarily match
        // This gets made based on all_sample_halpotypes
        std::unordered_map<std::string, size_t> sample_to_index;

        //////////////////////////// Extra housekeeping stuff

        /// This goes at the beginning of the file to ensure that it is the right file type and version
        inline const static std::string file_header = "#NODE_DATA_v1.0";
        size_t number_node_analyzed = 0;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////// Private functions
    private:
    
        //////////////////// Helper functions for writing and loading stuff from files

        /// Write just the header of the collection of nodes to the given file
        void write_node_data_collection_header(Writer& out_writer) const;

        /// Write one line for one node to the given file. This assumes that write_node_data_collection_header() has already been called
        void write_node_data_collection_line(Writer& out_writer, const node_info_internal_t& node_info) const;

        /// Given a stream to the start of the file, load just the header. The stream will be advanced to point to the beginning of the node records
        void load_node_data_collection_header(stoat::Reader& in_reader);

        /// Given a string representing a line in the file, load one node_info_internal_t
        /// This assumes that load_node_collection_header() has already been called
        node_info_internal_t load_node_data_line(std::string& line);

        // Helper function for for_each_node() to go from internal node info to running iteratee on the node_info_t
        // Note that this can't be a function to return the node_info_t because it has references to tables and stuff that would go out of scope
        void run_iteratee_on_one_node(const node_info_internal_t& internal_node_info, const std::function<void(node_info_t& node_info)>& iteratee) const;

};
}

#endif
