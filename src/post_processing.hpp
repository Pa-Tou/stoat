#ifndef POST_PROCESSING_HPP
#define POST_PROCESSING_HPP

#include <vector>
#include <tuple>
#include <string>
#include "types_and_structs.hpp"

#include <handlegraph/path_position_handle_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>


namespace stoat{

// Given a vector of <p-value, 1.0, line index from the input file>, fill in the vector with the adjusted p-value
// and sort the vector by adjusted p-value
void adjust_pvalues_with_BH(std::vector<std::tuple<double, double, size_t>>& data);

// Read a tsv from input_file, collect the p-values from the correct column (depending on phenotype_type), 
// and write the same file.
void add_BH_adjusted_column(
    const std::string& input_file,
    const stoat::phenotype_type_t& phenotype_type);

// The same, but specify the column number (0-indexed) of the p-value
// If p_col is std::numeric_limits<size_t>::max(), then check the header for the column number.
// If p_col is given, then still check the header to make sure that the column label is some sort of P-value
// Add the adjusted p-value column after the p-value column
void add_BH_adjusted_column(
    const std::string& input_file,
    size_t p_col);


/// Copy the input file (which is the output of stoat) and re-write it with reference coordinates on any path starting with reference_prefix
/// write result to stdout
/// This deals with its own file handling
void change_reference(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index, 
                      const std::string& input_file, const std::unordered_set<std::string>& reference_names);

} // namespace stoat_vcf

#endif // ADJUSTED_PVALUE_HPP
