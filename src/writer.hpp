#ifndef WRITER_INCLUDED
#define WRITER_INCLUDED

#include <iostream>
#include <handlegraph/path_position_handle_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include "utils.hpp"
#include "snarl_data_t.hpp"
#include "snarl_data_collection.hpp"

using namespace std;
namespace stoat {

// Write headers
void write_stoat_output_header(std::ostream& outstream, stoat::phenotype_type_t phenotype_type);
void write_binary_header(std::ostream& outstream);

// Write lines
void write_binary(std::ostream& outstream, const std::string& chr, const Snarl_data_t& snarl_data_s, const std::string& type_var_str,
                    const std::string& fastfisher_p_value, const std::string& chi2_p_value, const std::string& group_paths);

void write_binary(std::ostream& outstream, const snarl_info_t& snarl_data,
                    const std::string& fastfisher_p_value, const std::string& chi2_p_value, const std::string& group_paths);

void write_binary_covar(std::ostream& outstream, const std::string& chr, const Snarl_data_t& snarl_data_s, const std::string& type_var_str,
                        const std::string& p_value,
                        const std::string& beta, const std::string& se, const std::vector<size_t>& allele_paths);

void write_quantitative(std::ostream& outstream, const std::string& chr, const Snarl_data_t& snarl_data_s, const std::string& type_var_str,
                        const std::string& p_value, const std::string& r2, const std::vector<size_t>& allele_paths);

void write_eqtl(std::ostream& outstream, const std::string& chr, const Snarl_data_t& snarl_data_s, const std::string& type_var_str,
                    const std::string& gene_name, const std::string& p_value, const std::string& r2, const std::vector<size_t>& allele_paths);

void write_vcf(std::ostream& outstream, const std::string& chr, const size_t& pos, const std::string& id,
               const std::string& ref, const std::string& alt, const std::string& paths, 
               const std::vector<std::vector<char>>& genotype);

void write_fasta(std::ostream& outstream, const handlegraph::PathPositionHandleGraph& graph, 
                 const bdsg::SnarlDistanceIndex& distance_index, const snarl_info_t& snarl_info);

void writeSignificantTableToTSV(
    const std::vector<std::vector<double>>& table,
    const std::vector<std::string>& list_snarl,
    const std::vector<std::string>& list_samples,
    const std::string& filename);

} //end namespace

#endif
