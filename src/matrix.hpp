#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <cstdint>

#include "snarl_data_t.hpp"

using namespace std;

namespace stoat_vcf {

// A class to store a 2d bit-matrix
// Rows represent edges and the index of each edge can be found from the row_header
// Columns represent samples/haplotypes
class EdgeBySampleMatrix {
public:
    EdgeBySampleMatrix(const std::vector<std::string>& sample_names, size_t n_edges, size_t n_samp_haps);
    ~EdgeBySampleMatrix()=default;

    // Operator to get the value
    bool operator()(size_t row, size_t col) const;
    bool get_edge(size_t row, size_t col) const;
    
    // Add this edge to the matrix
    void add_sample_edge(const stoat::Edge_t& edge, size_t col_index);

    // Set this value to true
    void set_edge(size_t row, size_t col);

    // Double the size of the matrix
    void expand_matrix();

    // Shrink to use the minimum amount of memory possible allowing the current number of rows
    void shrink();

    // Clear the memory and re-initialize
    void reset(size_t new_n_edges, size_t new_n_samp_haps);
    
    std::string get_sample_name(size_t sample_idx) const;

    std::vector<size_t> get_samples_on_path(const std::vector<stoat::Edge_t> &path) const;

protected:
    size_t n_samp_haps;
    size_t max_edges;
    std::vector<bool> matrix_1D;
    // size_t cols_;
    // size_t MaxElement;
    // std::vector<uint8_t> matrix_1D;
    std::unordered_map<stoat::Edge_t, size_t> row_header;
    std::vector<std::string> sample_names;

};

} // end namespace stoat

#endif
