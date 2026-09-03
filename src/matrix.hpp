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

struct pending_edge_t {
    stoat::edge_t edge;
    size_t hap_num;
};

// A class to store a 2d bit-matrix
// Rows represent edges and the index of each edge can be found from the row_header
// Columns represent samples/haplotypes
class EdgeBySampleMatrix {
public:
    EdgeBySampleMatrix(const std::vector<std::string>& sample_names, size_t n_edges);
    ~EdgeBySampleMatrix()=default;

    // Operator to get the value
    bool operator()(size_t row, size_t col) const;
    bool get_edge(size_t row, size_t col) const;
    
    // Add this edge to the matrix
    void add_sample_edge(const stoat::edge_t& edge, size_t col_index);

    // Set this value to true
    void set_edge(size_t row, size_t col);

    // Double the size of the matrix
    void expand_matrix();

    // Shrink to use the minimum amount of memory possible allowing the current number of rows
    void shrink();

    // Clear the memory and re-initialize
    void clear_edges(size_t new_n_edges);
    
    std::string get_sample_name(size_t sample_idx) const;

    // query the matrix: find sample-haplotypes that have all edges along the queried path
    std::vector<size_t> get_samples_on_path(const stoat::PathTraversal &path_trav) const;

    // clear the matrix and load edges in a VCF for a chunk (corresponding to a chromosome)
    // the returned pointers are the VCF file stream after reading that chunk
    void load_vcf_chunk(stoat_vcf::VCFParser& vcf_parser, std::string &chr);

protected:
    size_t n_samp_haps;
    size_t max_edges;
    std::vector<bool> matrix_1D;
    std::unordered_map<stoat::edge_t, size_t> row_header;
    std::vector<std::string> sample_names;

};



} // end namespace stoat

#endif
