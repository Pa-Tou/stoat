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
    EdgeBySampleMatrix(const std::vector<std::string>& sampleNames, size_t rows, size_t cols);
    ~EdgeBySampleMatrix()=default;

    // Operator to get the value
    bool operator()(size_t row, size_t col) const;

    // Add this edge to the matrix
    void push_matrix(const stoat::Edge_t& EdgePath, size_t indexColumn);

    // Set this value to true
    void set(size_t row, size_t col);

    // Get the maximum index into the vector representing the matrix
    size_t getMaxElement() const;
    
    // Double the size of the matrix
    void expandMatrix();

    // Shrink to use the minimum amount of memory possible allowing the current number of rows
    void shrink();

    // Clear the memory and re-initialize
    void reset(const std::vector<std::string>& newSampleNames, size_t rows, size_t cols);
    
    // Return the index of the edge in row_header, std::numeric_limits<size_t>::max() if the edge does not exist
    size_t find_edge(const stoat::Edge_t& edge) const;

    // Retrieve the index of `edge` if it exists. Otherwise, add it and return the new index.
    size_t getOrAddIndex(const stoat::Edge_t& key, const size_t& size_edge_index_dict);

protected:
    size_t cols_;
    size_t MaxElement;
    std::vector<uint8_t> matrix_1D;
    std::unordered_map<stoat::Edge_t, size_t> row_header;

// TODO: This shouldn't be public
public:
    std::vector<std::string> sampleNames;
};

} // end namespace stoat

#endif
