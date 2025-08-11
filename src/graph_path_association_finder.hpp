#ifndef STOAT_GRAPH_PATH_ASSOCIATION_FINDER_HPP_INCLUDED
#define STOAT_GRAPH_PATH_ASSOCIATION_FINDER_HPP_INCLUDED

#include <iostream>
#include <handlegraph/path_position_handle_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include "partitioner.hpp"

using namespace std;
using namespace stoat;

namespace stoat_graph{

/***
    General template class for finding associations in a graph.
    This will implement the helper functions needed to traverse the graph, filter snarls that 
    are too small, write the output, etc.
    Inherited classes must implement is_snarl_associated() to test each snarl.
***/
class AssociationFinder {

    protected:
        // At a minimum, an AssociationFinder must have a graph (with path information for printing reference coordinates),
        // a distance index, the names of the samples we are interested in, and, optionally, the name of the reference
        const handlegraph::PathPositionHandleGraph& graph;
        const bdsg::SnarlDistanceIndex& distance_index; 
        const std::set<std::string>& samples_of_interest;
        const std::string& reference_sample;
        const std::string& test_method;
        const std::string& output_format;
        size_t total_sample_count;
        size_t allele_size_limit;
        std::ostream& out_associated = std::cout;
        std::ostream& out_unassociated = std::cout;
        bool check_distances;


        // object for finding partitions of samples in a snarl
        std::shared_ptr<Partitioner> partitioner;


    public:

        /// Create an association finder with the graph and distance index, a set of samples of interest for which we want associated
        /// variants, a string of the reference sample name (may be empty), the output format (tsv or fasta), 
        /// filenames for writing associated alleles and unassociated alleles, and a size limit for the minimum length of snarl reported,
        /// measured as the "maximum" length of a snarl
        AssociationFinder(const handlegraph::PathPositionHandleGraph& graph, 
                          const bdsg::SnarlDistanceIndex& distance_index, 
                          std::shared_ptr<Partitioner> partitioner,
                          const std::set<std::string>& samples_of_interest, 
                          const std::string& reference_sample,
                          const std::string& test_method,
                          const std::string& output_format,
                          size_t allele_size_limit,
                          std::ostream& out_associated,
                          std::ostream& out_unassociated);

        
        /// Main function that gets called to go through all snarls in the graph, check if they are eligible with snarl_is_eligible(),
        /// use the Tester::test_snarl on each eligible snarl, then use the Writer to write snarls that pass the Tester/MultipleTester
        void test_snarls() const;

    protected:

       //////////////////////////////// Other helper functions

        // Do we care about this snarl? Based on allele_size_limit
        bool snarl_is_eligible(const handlegraph::net_handle_t& snarl) const;

};


}

#endif
