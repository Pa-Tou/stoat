#ifndef STOAT_VCF_UNTANGLER_HPP_INCLUDED
#define STOAT_VCF_UNTANGLER_HPP_INCLUDED

#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <vector>
#include <unordered_map>
#include <string>
#include "types_and_structs.hpp"



using namespace stoat;

namespace stoat_vcf {

/// This class is used to store information needed to resolve nested snarls in a VCF, since the current output of vg call can produce calls that are inconsistent. 
/// For example, if there are multiple SNPs nested inside an insertion, the caller may pick the best version of the insertion, which is wrong for one SNP. 
/// The SNP would be called correctly in the nested snarl, but it would be inconsistent with the outer snarl. It is also possible for the outer snarl to be called
/// as a deletion, but still have a call for the SNP.
///
/// This class stores the information needed to resolve these contradictions. For the first issue, it stores the bounds of nested snarls so that they can be 
/// skipped when going through the paths through parent snarls. For the nested snarls, it stores the genotypes from the parent snarls to determine how many copies
/// it should actually have.
/// TODO: There could still be an issue if the outer snarl says it is het for the SNP, but the SNP also says it it het
/// TODO: This assumes that all samples are always in the same order, which I guess they are
class VCFUntangler {



    public:

    // Initialize the untangler by opening the vcf and making it point to the first records
    VCFUntangler(const std::string& vcf_filename);

    /// Given the bound of a snarl and a sample (as the index of the sample and haplotype), return true if that sample goes through the snarl, according to the parent snarl.
    /// TODO: This assumes that the samples are always in a consistent order, which is the order in the VCF. This is done to match the make_edge_matrix() in snarl_analyzer.cpp
    bool does_sample_have_snarl(size_t sample_hap_index, stoat::node_traversal_t snarl_bound);

    /// Given a string representing a traversal of a node in a path through a snarl, check if the node is
    /// entering another snarl. If yes, return the string traversal of the end node of the snarl, to skip it in the parent. 
    /// Otherwise return the input snarl_bound
    stoat::node_traversal_t get_opposite_snarl_bound(stoat::node_traversal_t snarl_bound); 

    /// Clear the previous contents and get all the snarls for the given chromosome. Advance the file pointers
    void fill_in_snarls_for_chromosome(const std::string& chr);

    protected:

    /////////////////////////////////// Private data members

    // How many haplotypes total? eg, 2* sample count for diploid
    size_t sample_hap_count;
    // How many snarls do we have
    size_t snarl_count;
    
    // Stuff for keeping track of the vcf file, from htslib
    // Keep two sets, since we need two passes through the vcf for each chromosome for fill_in_nested_snarl_bounds and fill_in_nested_genotypes

    // I think this is a pointer into the vcf
    htsFile* ptr_vcf_bounds; 
    // idk what this is but the htslib stuff needs it
    bcf_hdr_t* hdr_bounds; 
    // Struct representing one line of the vcf
    bcf1_t* rec_bounds;

    htsFile* ptr_vcf_genotypes; 
    bcf_hdr_t* hdr_genotypes; 
    bcf1_t* rec_genotypes;


    /// For each nested snarl (everything except top-level snarls), map the start bound going in to the end bound going out, and to an index for genotypes.
    /// also stores the reverse to find the snarl from the end bound. 
    std::unordered_map<stoat::node_traversal_t, std::pair<stoat::node_traversal_t, size_t>> snarl_in_to_out;

    /// Genotype for each nested snarl. Indexed by snarl and sample/hap
    /// The index for the snarl is in snarl_bounds, which must be filled in before trying to use this
    /// There are sample_hap_count bools per snarl, and an unknown number of snarls
    std::vector<bool> genotypes;


    /////////////////// Helper functions

    // Get the index into the genotypes vector
    size_t genotype_index(size_t sample_hap_index, const stoat::node_traversal_t& snarl_bound);


    ///////////////////////////// Functions to fill in the snarls

    /// Do this before fill_in_nested_genotypes. This fills in snarl_bounds for all LV>=1
    void fill_in_nested_snarl_bounds(const std::string& chr);

    /// Go through the vcf a second time and for each snarl, check if it has nested snarls and fill in genotypes
    void fill_in_nested_genotypes(const std::string& chr);


};
}

#endif
