#ifndef WRITER_INCLUDED
#define WRITER_INCLUDED

#include <iostream>
#include <handlegraph/path_position_handle_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include "utils.hpp"
#include "snarl_data_t.hpp"

using namespace std;
namespace stoat {

// Write headers
void write_binary_header(std::ostream& outstream);
void write_binary_covar_header(std::ostream& outstream);
void write_quantitative_header(std::ostream& outstream);
void write_eqtl_header(std::ostream& outstream);

// Write lines
void write_binary(std::ostream& outstream, const std::string& chr, const Snarl_data_t& snarl_data_s, const std::string& type_var_str,
                    const std::string& fastfisher_p_value, const std::string& chi2_p_value, const std::string& p_value_adjusted, const std::string& group_paths);

void write_binary_covar(std::ostream& outstream, const std::string& chr, const Snarl_data_t& snarl_data_s, const std::string& type_var_str,
                        const std::string& p_value, const std::string& p_value_adjusted, const std::string& r2,
                        const std::string& beta, const std::string& se, const std::vector<size_t>& allele_paths);

void write_quantitative(std::ostream& outstream, const std::string& chr, const Snarl_data_t& snarl_data_s, const std::string& type_var_str,
                        const std::string& p_value, const std::string& p_value_adjusted, const std::string& r2,
                        const std::string& beta, const std::string& se, const std::vector<size_t>& allele_paths);

void write_eqtl(std::ostream& outstream, const std::string& chr, const Snarl_data_t& snarl_data_s, const std::string& type_var_str,
                    const std::string& gene_name, const std::string& p_value, const std::string& p_value_adjusted, const std::string& r2,
                    const std::string& beta, const std::string& se, const std::vector<size_t>& allele_paths);

// Write the fasta for paths in a snarl. If samples is given, only write the fast for samples present to outstream_associated if the sample maps to true, and
// to outstream_unassociated if it maps to false. If samples is not given, write all samples.
void write_fasta(std::ostream& outstream_associated, std::ostream& outstream_unassociated, const handlegraph::PathPositionHandleGraph& graph, 
                 const bdsg::SnarlDistanceIndex& distance_index, const handlegraph::net_handle_t& snarl, 
                 const std::unordered_map<std::string, bool>& samples, const string& reference_name);

void writeSignificantTableToTSV(
    const std::vector<std::vector<double>>& table,
    const std::vector<std::string>& list_snarl,
    const std::vector<std::string>& list_samples,
    const std::string& filename);

} //end namespace

#endif
