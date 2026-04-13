#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>

#include "types_and_structs.hpp"
#include "vcf_parser.hpp"

using namespace std;

namespace stoat_vcf {

// A class to store a 2d bit-matrix
// Rows represent keys and the index of each key can be found from the row_header
// Columns represent samples/haplotypes
template<class ValueType> 
class AlleleBySampleMatrix {
public:
    AlleleBySampleMatrix(const std::vector<std::string>& sample_names, size_t n_keys);
    ~AlleleBySampleMatrix()=default;

    // Operator to get the value
    bool operator()(size_t row, size_t col) const;
    bool get_key(size_t row, size_t col) const;
    
    // Add this key to the matrix
    void add_sample_key(const ValueType& key, size_t col_index);

    // Set this value to true
    void set_key(size_t row, size_t col);

    // Double the size of the matrix
    void expand_matrix();

    // Shrink to use the minimum amount of memory possible allowing the current number of rows
    void shrink();

    // Clear the memory and re-initialize
    void clear_keys(size_t new_n_keys);
    
    std::string get_sample_name(size_t sample_idx) const;

    // query the matrix: find sample-haplotypes that have all keys along the queried path
    std::vector<size_t> get_samples_on_path(const stoat::PathTraversal &path_trav) const;

    std::vector<size_t> get_samples_on_node(const stoat::node_traversal_t &node_trav) const;

    // clear the matrix and load keys in a VCF for a chunk (corresponding to a chromosome)
    // the returned pointers are the VCF file stream after reading that chunk
    void load_vcf_chunk(stoat_vcf::VCFParser& vcf_parser, std::string &chr);

protected:
    size_t n_samp_haps;
    size_t max_keys;
    std::vector<bool> matrix_1D;
    std::unordered_map<ValueType, size_t> row_header;
    std::vector<std::string> sample_names;
    size_t number_node_not_found = 0;

};

} // end namespace stoat

#endif
