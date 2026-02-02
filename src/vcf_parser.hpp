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
/// 1. Make a VCFParser and call initialize_parser
/// 2. Call get_next_chromosome to see what chromosome this is currently pointing at
/// 3. If it is a chromosome we want, call for_each_record_on_chromosome to do whatever on everything on this chromosome
/// 4. Repeat 2 and 3 until get_next_chromosome returns an empty string 
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

    /// This should be called after running through the vcf to close the file
    void close_vcf(){
        bcf_destroy(rec);
        bcf_hdr_destroy(hdr);
        bcf_close(ptr_vcf);
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



};
}

#endif
