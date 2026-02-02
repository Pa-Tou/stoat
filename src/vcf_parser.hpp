#ifndef STOAT_VCF_PARSER_HPP_INCLUDED
#define STOAT_VCF_PARSER_HPP_INCLUDED

#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <vector>
#include <unordered_map>
#include <string>
#include "types_and_structs.hpp"



using namespace stoat;

namespace stoat_vcf {


/// A struct for holding one line of the VCF
/// Everything is a reference so it can only be used in the context of for_each_record_on_chromosome unless its contents are copied 
/// Note: this currently only holds stuff that we need. To add more, it must be retrieved in for_each_record_on_chromosome()
struct vcf_info_t {
    const size_t lv; //LV field of the VCF
    //std::string snarl; // The string representation of the snarl, from the ID field
    const std::vector<int>& genotype; // The GT field TODO: If the VCF parser contains a VCFUntangler, this removes a hom call if the parent was het for this snarl
    const std::vector<std::vector<node_traversal_t>>& paths; // The paths through the snarl, from the AT field. TODO: If the VCF parser contains a VCFUntangler, this skips nested snarls 
};

/// This is used to walk through a VCF file with a VCFUntangler and keep everything synchronized
/// A workflow would be 
/// 1. Make a VCFParser and call initialize_parser() to parse the header
/// 2. Call get_next_chromosome to see what chromosome this is currently pointing at
/// 3. If it is a chromosome we want, call for_each_record_on_chromosome() to do whatever on everything on this chromosome
///    Otherwise, call skip_to_next_chromosome() to continue reading until the next chromosome
/// 4. Repeat 2 and 3 until get_next_chromosome returns an empty string 
/// 5. Call close_vcf()
class VCFParser {

    public:

    /// This does nothing. initialize_parser() must be called to actually fill stuff in from a file.
    VCFParser() {};

    /// Parse the header and return the sample names
    std::vector<std::string> initialize_parser(const std::string& vcf_path); 

    /// From wherever in the VCF file we currently are, what is the chromosome?
    /// This doesn't advance the VCF file so it can be called multiple times pointing to the same thing
    /// Return an empty string if we've finished the VCF
    std::string get_next_chromosome_name();

    /// Given a chromosome name, run iteratee on every consecutive VCF record on this chromosome, starting with the record currently being pointed at.
    /// This advances through the VCF until it is pointing at the first thing not on this chromosome (or the end of the file)
    /// If the VCFParser is not pointing to a record on this chromosome, do nothing.
    void for_each_record_on_chromosome(const std::string& chr, const std::function<void(const vcf_info_t& vcf_info)>& iteratee);

    /// Read through the VCF file until we find a new chromosome that is not chr.
    /// This is equivalent to running for_each_record_on_chromosome and doing nothing in iteratee
    void skip_to_next_chromosome(const std::string& chr);

    /// This should be called after running through the vcf to close the file
    void close_vcf(){
        bcf_destroy(rec);
        bcf_hdr_destroy(hdr);
        bcf_close(ptr_vcf);

        bcf_destroy(rec_bounds);
        bcf_hdr_destroy(hdr_bounds);
        bcf_close(ptr_vcf_bounds);

        bcf_destroy(rec_genotypes);
        bcf_hdr_destroy(hdr_genotypes);
        bcf_close(ptr_vcf_genotypes);
    }

    // How many haplotypes ? Number of samples * ploidy
    size_t hap_count;

    protected:

    // I think this is a file handle for the vcf
    htsFile* ptr_vcf;
    // idk what this is but the htslib stuff needs it
    bcf_hdr_t* hdr;
    // Struct representing one line of the vcf
    bcf1_t* rec;

    // What is the return value of the last call to bcf_read?
    // 0 if we've finished reading the file
    int read_status;



    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////// Functions and data for resolving inconsistent nested calls

    protected:

    /// The parser also stores information needed to resolve nested snarls in a VCF, since the current output of vg call can produce calls that are inconsistent. 
    /// For example, if there are multiple SNPs nested inside an insertion, the caller may pick the best version of the insertion, which is wrong for one SNP. 
    /// The SNP would be called correctly in the nested snarl, but it would be inconsistent with the outer snarl. It is also possible for the outer snarl to be called
    /// as a deletion, but still have a call for the SNP.
    ///
    /// For the first issue, store the bounds of nested snarls so that they can be skipped when going through the paths through parent snarls. 
    /// For the nested snarls, store the genotypes from the parent snarls to determine how many copies it should actually have.
    /// TODO: There could still be an issue if the outer snarl says it is het for the SNP, but the SNP also says it it het
    /// TODO: This assumes that all samples are always in the same order, which I guess they are



    /// Given the bound of a snarl and a sample (as the index of the sample and haplotype), return true if that sample goes through the snarl, according to the parent snarl.
    /// TODO: This assumes that the samples are always in a consistent order, which is the order in the VCF. This is done to match the make_edge_matrix() in snarl_analyzer.cpp
    bool does_sample_have_snarl(size_t sample_hap_index, stoat::node_traversal_t snarl_bound);

    /// Given a string representing a traversal of a node in a path through a snarl, check if the node is
    /// entering another snarl. If yes, return the string traversal of the end node of the snarl, to skip it in the parent. 
    /// Otherwise return the input snarl_bound
    stoat::node_traversal_t get_opposite_snarl_bound(stoat::node_traversal_t snarl_bound); 

    /// Clear the previous contents and get all the snarls for the given chromosome. Advance the file pointers
    void fill_in_snarls_for_chromosome(const std::string& chr);

    /////////////////////////////////// Private data members

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
    /// There are hap_count bools per snarl, and an unknown number of snarls
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
