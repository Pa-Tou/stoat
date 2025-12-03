#include "writer.hpp"

//#define DEBUG_WRITER

namespace stoat {

void write_stoat_output_header(std::ostream& outstream, stoat::phenotype_type_t phenotype_type) {
    if (phenotype_type == stoat::BINARY) {
        outstream << "#CHR\tSTART_POS\tEND_POS\tSNARL\tPATH_LENGTHS\tP_FISHER\tP_CHI2\tGROUP_PATHS\tDEPTH" << std::endl;
    } else if (phenotype_type == stoat::BINARY_COVAR) {
        outstream << "#CHR\tSTART_POS\tEND_POS\tSNARL\tPATH_LENGTHS\tP\tBETA\tSE\tALLELE_PATHS\tDEPTH" << std::endl;
    } else if (phenotype_type == stoat::QUANTITATIVE) {
        outstream << "#CHR\tSTART_POS\tEND_POS\tSNARL\tPATH_LENGTHS\tP\tRSQUARE\tALLELE_PATHS\tDEPTH" << std::endl;
    } else if (phenotype_type == stoat::EQTL) {
        outstream <<  "#CHR\tSTART_POS\tEND_POS\tSNARL\tPATH_LENGTHS\tGENE\tP\tRSQUARE\tALLELE_PATHS\tDEPTH" << std::endl;
    } 
}

void write_binary_header(std::ostream& outstream) {
    outstream << "#CHR\tSTART_POS\tEND_POS\tSNARL\tPATH_LENGTHS\tP_FISHER\tP_CHI2\tGROUP_PATHS\tDEPTH" << std::endl;
}

void write_binary(std::ostream& outstream, const std::string& chr, const Snarl_data_t& snarl_data_s, const std::string& type_var_str,
                        const std::string& fastfisher_p_value, const std::string& chi2_p_value, const std::string& group_paths) {

    outstream << chr << "\t"
              << snarl_data_s.start_position << "\t"
              << snarl_data_s.end_position << "\t"
              << stoat::pairToString(snarl_data_s.ids) << "\t"
              << type_var_str << "\t"
              << fastfisher_p_value << "\t"
              << chi2_p_value << "\t"
              << group_paths << "\t"
              << snarl_data_s.depth << endl;
}

void write_binary(std::ostream& outstream, const snarl_info_t& snarl_data,
                  const std::string& fastfisher_p_value, const std::string& chi2_p_value, const std::string& group_paths) {

    outstream << snarl_data.ref_path << "\t"
              << snarl_data.start_position << "\t"
              << snarl_data.end_position << "\t"
              << stoat::pairToString(std::make_pair(snarl_data.start_node.get_node_id(), snarl_data.end_node.get_node_id())) << "\t"
              << ((snarl_data.walks_by_allele.empty()) 
                    ? "." 
                    : stoat::vectorPathToString(snarl_data.walks_by_allele, true)) << "\t"
              << fastfisher_p_value << "\t"
              << chi2_p_value << "\t"
              << group_paths << "\t"
              << snarl_data.depth << endl;
}

void write_binary_covar(std::ostream& outstream, const std::string& chr, const Snarl_data_t& snarl_data_s, const std::string& type_var_str,
                        const std::string& p_value, const std::string& beta, const std::string& se, const std::vector<size_t>& allele_paths) {

    outstream << chr << "\t"
              << snarl_data_s.start_position << "\t"
              << snarl_data_s.end_position << "\t"
              << stoat::pairToString(snarl_data_s.ids) << "\t"
              << type_var_str << "\t"
              << p_value << "\t"
              << beta << "\t"
              << se << "\t"
              << stoat::vectorToString(allele_paths) << "\t"
              << snarl_data_s.depth << endl;
}

void write_quantitative(std::ostream& outstream, const std::string& chr, const Snarl_data_t& snarl_data_s, const std::string& type_var_str,
                        const std::string& p_value,  const std::string& r2, const std::vector<size_t>& allele_paths) {

    outstream << chr << "\t"
              << snarl_data_s.start_position << "\t"
              << snarl_data_s.end_position << "\t"
              << stoat::pairToString(snarl_data_s.ids) << "\t"
              << type_var_str << "\t"
              << p_value  << "\t"
              << r2 << "\t"
              << stoat::vectorToString(allele_paths) << "\t"
              << snarl_data_s.depth << "\n";

}

void write_eqtl(std::ostream& outstream, const std::string& chr, const Snarl_data_t& snarl_data_s, const std::string& type_var_str,
                   const std::string& gene_name, const std::string& p_value,  const std::string& r2, const std::vector<size_t>& allele_paths) {

    outstream << chr << "\t"
              << snarl_data_s.start_position << "\t"
              << snarl_data_s.end_position << "\t"
              << stoat::pairToString(snarl_data_s.ids) << "\t"
              << type_var_str << "\t"
              << gene_name << "\t"
              << p_value  << "\t"
              << r2 << "\t"
              << stoat::vectorToString(allele_paths) << "\t"
              << snarl_data_s.depth << endl;

}

void write_vcf(std::ostream& outstream,
               const std::string& chr,
               const size_t& pos,
               const std::string& id,
               const std::string& ref,
               const std::string& alt,
               const std::string& paths,
               const std::vector<std::vector<char>>& genotype) {
    
    // Build genotype string for all samples (e.g., "0/1", "1/1", "0/0")
    std::ostringstream gt_stream;
    for (size_t i = 0; i < genotype.size(); ++i) {
        const auto& sample_gt = genotype[i];
        // Join alleles with '/'
        for (size_t j = 0; j < sample_gt.size(); ++j) {
            gt_stream << sample_gt[j];
            if (j + 1 < sample_gt.size()) {
                gt_stream << '/';
            }
        }
        if (i + 1 < genotype.size()) {
            gt_stream << '\t'; // separate multiple samples by tab
        }
    }

    outstream << chr << '\t'
              << pos << '\t'
              << id << '\t'
              << ref << '\t'
              << alt << '\t'
              << '.' << '\t'             // QUAL placeholder
              << '.' << '\t'             // FILTER placeholder
              << "AT=" << paths << '\t'  // INFO field
              << "GT" << '\t'            // FORMAT field
              << gt_stream.str()          // Sample genotypes
              << '\n';
}

void write_fasta(std::ostream& outstream, const handlegraph::PathPositionHandleGraph& graph,
                 const snarl_partition_t& snarl_info, const std::unordered_map<std::string, bool>& samples, const string& reference_name) {
    
    string ref_coordinates = snarl_info.ref_path + ":" + std::to_string(snarl_info.start_position) + "-" + std::to_string(snarl_info.end_position);

    
    // Now go through each path that goes through the snarl and print the sequence
    std::vector<stoat::path_range_t> path_ranges = get_coordinates_between_nodes(graph, snarl_info.start_handle, snarl_info.end_handle, false, "", true);
    for (const stoat::path_range_t& path_range : path_ranges) {
        handlegraph::path_handle_t path = graph.get_path_handle_of_step(path_range.start);
        string sample_name = stoat::get_sample_name_from_path(graph, path);
        if (samples.empty() || samples.count(sample_name) != 0) {
            //If we aren't checking samples, or if this is a sample we want
    
            std::tuple<std::string, size_t, size_t> range_coordinates = get_name_and_offsets_of_snarl_path_range(graph, path_range);
            // Print the header
            outstream << ">snarl:" << graph.get_id(snarl_info.start_handle) << "-" << graph.get_id(snarl_info.end_handle) << "|"
                << ref_coordinates << "|"
                << std::get<0>(range_coordinates) << ":"
                << std::get<1>(range_coordinates) << "-"    
                << std::get<2>(range_coordinates) << endl;
    
            // Now print the sequence in 80bp chunks.
            // Keep a buffer to print 80 bp at a time
            std::string sequence_buffer = "";
            handlegraph::step_handle_t next_step = graph.get_next_step(path_range.start);
            while (next_step != path_range.end) {
                std::string node_seq = graph.get_sequence(graph.get_handle_of_step(next_step));
                while (node_seq.size() != 0) {
    
                    // Fill in sequence_buffer up to 80 characters
                    size_t to_add = 80 - sequence_buffer.size();
                    sequence_buffer += node_seq.substr(0, to_add);
                    node_seq.erase(0, to_add);
    
                    // If the buffer is full, write it and clear it
                    if (sequence_buffer.size() == 80) {
                        outstream << sequence_buffer << endl;
                        sequence_buffer.clear();
                    }
                }
                handlegraph::step_handle_t step = next_step;
                if (!graph.has_next_step(step)) {
                    break;
                }
                next_step = graph.get_next_step(step);
            }
            if (!sequence_buffer.empty()) {
                outstream << sequence_buffer << endl;
            }
        }
    }
}

void write_fasta(std::ostream& outstream, const handlegraph::PathPositionHandleGraph& graph,
                 const bdsg::SnarlDistanceIndex& distance_index, const snarl_info_t& snarl_info) {
    
    // The reference coordinates of this snarl as a string
    string ref_coordinates = snarl_info.ref_path + ":" + std::to_string(snarl_info.start_position) + "-" + std::to_string(snarl_info.end_position);


    // Get the coordinate range of every sample going through this snarl
    // We will match each sample to one of these ranges later
    handlegraph::handle_t start_node = graph.get_handle(snarl_info.start_node.get_node_id(), snarl_info.start_node.get_is_reverse());
    handlegraph::handle_t end_node = graph.get_handle(snarl_info.end_node.get_node_id(), snarl_info.end_node.get_is_reverse());
    std::vector<stoat::path_range_t> path_ranges = get_coordinates_between_nodes(graph, start_node, end_node, false, "", true);

    for (size_t allele_i = 0 ; allele_i < snarl_info.walks_by_allele.size() ; allele_i++) {
        // For each allele

        // These sequence of the allele
        std::string sequence = snarl_info.sequences_by_allele[allele_i];

        std::tuple<std::string, size_t, size_t> range_coordinates ("", std::numeric_limits<size_t>::max(), 0);

        // For any one haplotype with this allele, find the coordinates from path_ranges
        for (const stoat::path_range_t& path_range : path_ranges) {

            handlegraph::path_handle_t path = graph.get_path_handle_of_step(path_range.start);

            if (snarl_info.sample_sets_by_allele[allele_i].count(stoat::sample_hap_t(graph, path))) {
                // If the sample/haplotype of this path is in the set for this allele

                auto this_range_coordinates = get_name_and_offsets_of_snarl_path_range(graph, path_range);
                std::get<0>(range_coordinates) = std::get<0>(this_range_coordinates);
                std::get<1>(range_coordinates) = std::min(std::get<1>(this_range_coordinates), std::get<1>(range_coordinates));
                std::get<2>(range_coordinates) = std::max(std::get<2>(this_range_coordinates), std::get<2>(range_coordinates));
            }
        }

        //////// Now print the fasta 

        // Print the header
        outstream << ">snarl:" << stoat::pairToString(std::make_pair(snarl_info.start_node.get_node_id(), snarl_info.end_node.get_node_id())) << "|"
            << ref_coordinates << "|"
            << std::get<0>(range_coordinates) << ":"
            << std::get<1>(range_coordinates) << "-"    
            << std::get<2>(range_coordinates) << endl;

        // Now print the sequence in 80bp chunks.
        // Keep a buffer to print 80 bp at a time
        std::string sequence_buffer = "";
        while (sequence.size() != 0) { 
            // Fill in sequence_buffer up to 80 characters
            size_t to_add = 80 - sequence_buffer.size();
            sequence_buffer += sequence.substr(0, to_add);
            sequence.erase(0, to_add);
            
            // If the buffer is full, write it and clear it
            if (sequence_buffer.size() == 80) {
                outstream << sequence_buffer << endl;
                sequence_buffer.clear();
            }

        }
        if (!sequence_buffer.empty()) {
            outstream << sequence_buffer << endl;
        }
    }
}

// Write the table to a TSV file
void writeSignificantTableToTSV(
    const std::vector<std::vector<double>>& table,
    const std::vector<std::string>& list_snarl,
    const std::vector<std::string>& list_samples,
    const std::string& filename) {

    std::ofstream outFile(filename);

    // Write header
    outFile << "sample_name";
    for (const auto& snarl_name : list_snarl) {
        outFile << "\t" << snarl_name;
    }
    outFile << "\n";

    // Write each sample's data
    size_t itr = 0;
    for (const auto& allele_vector : table) {
        outFile << list_samples[itr];

        for (size_t i=0; i < allele_vector.size(); ++i) {
            outFile << "\t" << allele_vector[i];
        }
        outFile << "\n";
        ++itr;
    }
    outFile.close();
}

}//end namespace

