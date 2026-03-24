#include "post_processing.hpp"
#include "utils.hpp"
#include "writer.hpp"

namespace stoat {

// Adjust p-values using the Benjamini-Hochberg procedure
void adjust_pvalues_with_BH(std::vector<std::tuple<double, double, size_t>>& data) {
    size_t n = data.size();
    if (n == 0) return;

    // Sort by raw p-value
    std::sort(data.begin(), data.end(), [](const auto& a, const auto& b) {
        return std::get<0>(a) < std::get<0>(b);
    });

    std::vector<double> adjusted(n, 0.0);

    // Compute adjusted p-values
    for (size_t i = 0; i < n; ++i) {
        double p = std::get<0>(data[i]);
        adjusted[i] = p * n / (i + 1);
    }

    // Ensure monotonicity (adjusted[i - 1] <= adjusted[i])
    for (size_t i = n - 1; i > 0; --i) {
        adjusted[i - 1] = std::min(adjusted[i - 1], adjusted[i]);
    }

    // Clamp to [0, 1] and assign adjusted p-values back
    for (size_t i = 0; i < n; ++i) {
        std::get<1>(data[i]) = std::min(1.0, adjusted[i]);
    }

    // Restore original order by index
    std::sort(data.begin(), data.end(), [](const auto& a, const auto& b) {
        return std::get<2>(a) < std::get<2>(b);
    });
}

// Main Function
void add_BH_adjusted_column(
    const std::string& input_file,
    const std::string& output_dir,
    const std::string& output_file_significant,
    const stoat::phenotype_type_t& phenotype_type) {

    size_t adjusted_col_index;

    if (phenotype_type == stoat::BINARY || phenotype_type == stoat::EQTL) {
        adjusted_col_index = 7;
    } else if (phenotype_type == stoat::QUANTITATIVE || phenotype_type == stoat::BINARY_COVAR) {
        adjusted_col_index = 6;
    }

    add_BH_adjusted_column(input_file, output_dir, output_file_significant, adjusted_col_index-1);
}

// Main Function
void add_BH_adjusted_column(
    const std::string& input_file,
    const std::string& output_dir,
    const std::string& output_file_significant,
    size_t p_col_index) {

    // Get the reader
    std::shared_ptr<stoat::Reader> reader;
    if ((input_file.compare(input_file.length()-3, 3, ".gz") == 0) ||
        (input_file.compare(input_file.length()-4, 4, ".bgz") == 0)) {
        reader.reset(new BgzReader(input_file));
    } else {
        reader.reset(new StdReader(input_file));
    }


    std::string col;

    // First pass: Collect p-values
    std::vector<std::tuple<double, double, size_t>> pvalues;
    std::string line;
    size_t line_index = 0;

    // Read the header line
    std::string header_line;
    reader->getline(header_line);
    std::stringstream header_ss(header_line);
    std::vector<std::string> headers;
    while (std::getline(header_ss, col, '\t')) {
        headers.push_back(col);
    }

    // Get the column number of the p-value or check that the given column is a p-value
    bool find_column = true;
    if (p_col_index != std::numeric_limits<size_t>::max() && 
        (headers[p_col_index] == "P" || headers[p_col_index] == "P_FISHER" || headers[p_col_index] == "P_CHI2")) {
        find_column = false;
    } else if (p_col_index != std::numeric_limits<size_t>::max()) {
        std::cerr << "warning [stoat BHcorrect]: given column with header " << headers[p_col_index] << " is not a valid p-value. Checking header for a valid column" << std::endl;
        p_col_index = std::numeric_limits<size_t>::max();
    }
    // Look for a column with a P-value. Always choose the "P" column. Prioritize "P-CHI2" over "P_FISHER"
    if (find_column) {
        for (size_t i = 0 ; i < headers.size() ; i++) {
            if ( headers[i] == "P") {
                // Always pick column "P"
                p_col_index = i;
                break;
            } else if (headers[i] == "P_CHI2") {
                // Always pick column "P_CHI2"
                p_col_index = i;
                break;
            } else if (headers[i] == "P_FISHER" && p_col_index == std::numeric_limits<size_t>::max()) {
                // Pick column "P_FISHER" if we don't have anything else but let it be overridden if there's something better
                p_col_index = i;
            }

        }
    }

    while (reader->getline(line)) {
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> columns;

        while (std::getline(ss, token, '\t')) {
            columns.push_back(token);
        }


        double pval =stoat::string_to_pvalue(columns[p_col_index]);
        //if (phenotype_type ==stoat::BINARY) {
        //    // combine both p-value
        //    //pval = stoat::set_precision_float_50(columns[4], columns[5]);
        //}
        pvalues.emplace_back(pval, 1.0, line_index++);
    }
    reader->close();

    // Apply BH correction
    adjust_pvalues_with_BH(pvalues);

    // Second pass: rewrite with BH-adjusted values

    // Open the reader again
    std::shared_ptr<stoat::Reader> new_reader;
    std::shared_ptr<stoat::Writer> writer;
    std::shared_ptr<stoat::Writer> writer_significant;
    const std::string output_temp_file = output_dir + "/temp_output.tsv";
    if ((input_file.compare(input_file.length()-3, 3, ".gz") == 0) ||
        (input_file.compare(input_file.length()-4, 4, ".bgz") == 0)) {
        reader.reset(new BgzReader(input_file));
        writer.reset(new BgzWriter(output_temp_file + ".gz"));
    } else {
        reader.reset(new StdReader(input_file));
        writer.reset(new StdWriter(output_temp_file));
    }
    if ((output_file_significant.compare(output_file_significant.length()-3, 3, ".gz") == 0) ||
        (output_file_significant.compare(output_file_significant.length()-4, 4, ".bgz") == 0)) {
        writer_significant.reset(new BgzWriter(output_file_significant));
    } else {
        writer_significant.reset(new StdWriter(output_file_significant));
    }
    std::ofstream outfile_significant(output_file_significant);

    // Write headers
    for (size_t i = 0 ; i < headers.size() ; i++ ) {
        writer->write(headers[i]);
        if (i == p_col_index) {
            writer->write("\tP_ADJUSTED");
            writer_significant->write("\tP_ADJUSTED");
        } 

        if (i == headers.size()-1) {
            writer->write("\n");
            writer_significant->write("\n");
        } else {
            writer->write("\t");
            writer_significant->write("\t");
        }
    }

    new_reader->getline(line); // Skip header again
    line_index = 0;

    while (new_reader->getline(line)) {
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> columns;

        while (std::getline(ss, token, '\t')) {
            columns.push_back(token);
        }

        double adjusted_p = std::get<1>(pvalues[line_index]);
        std::string adj_str = stoat::set_precision(adjusted_p);

        // Write updated line
        for (size_t i = 0; i < columns.size(); ++i) {
            writer->write(columns[i]);
            if (i == p_col_index) {
                writer->write("\t");
                writer->write(adj_str);
            }
            if (i != columns.size() - 1) writer->write("\t");
        }

        writer->write("\n");

        if (adjusted_p < 1e-5) {
            for (size_t i = 0; i < columns.size(); ++i) {
                writer_significant->write(columns[i]);
                if (i != columns.size() - 1) writer_significant->write("\t");
            }
            writer_significant->write("\n");
        }
        ++line_index;
    }

    new_reader->close();
    writer->close();
    writer_significant->close();

    // Replace original file
    std::remove(input_file.c_str());
    std::rename(output_temp_file.c_str(), input_file.c_str());
}

void change_reference(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index, 
                      const std::string& input_file, const std::unordered_set<std::string>& reference_names) {
    std::ifstream instream;
    instream.open(input_file);

    // Read the header line
    std::string header_line;
    std::getline(instream, header_line);
    std::stringstream header_ss(header_line);
    std::vector<std::string> headers;
    std::string col;
    while (std::getline(header_ss, col, '\t')) {
        headers.push_back(col);
    }

    if (headers[0] != "#CHR" || headers[1] != "START_OFFSET" || headers[2] != "END_OFFSET" || headers[3] != "START_NODE" || headers[4] != "END_NODE") {
        throw std::runtime_error("[stoat::change_reference]: error: input tsv has the wrong format");
    }

    // Write the header line
    std::cout << header_line << std::endl;

    std::string line;
    while (std::getline(instream, line)) {

        // Parse the line into a vector of strings
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> columns;

        while (std::getline(ss, token, '\t')) {
            columns.push_back(token);
        }

        // Get the snarl bounds
        stoat::node_traversal_t start_node(columns[3]);
        stoat::node_traversal_t end_node(columns[4]);

        // The snarl is actually just id1_id2, so we don't know the direction needed to get it
        // This could be done just with the graph but I prefer to use the function we already have with the snarls
        // to be consistent 
        handlegraph::net_handle_t node_net = distance_index.get_node_net_handle(start_node.get_node_id(), start_node.get_is_reverse());
        handlegraph::net_handle_t snarl_net = distance_index.get_root();

        // Get the snarl, which should be next thing in the chain
        distance_index.follow_net_edges(node_net, &graph, false, [&](const handlegraph::net_handle_t& next) {
            snarl_net = next;
            return true;
        });


        // Get the offsets of the start and end nodes along the reference
        std::string new_reference_name;
        size_t new_reference_start;
        size_t new_reference_end;
        std::vector<stoat::path_range_t> ranges = stoat::get_coordinates_of_snarl(graph, distance_index, snarl_net, true, reference_names, false);
        if (ranges.size() != 0) {
            // Check if we have already seen the reference path and if not add it
        
            auto reference_range = get_name_and_offsets_of_snarl_path_range(graph, ranges.front());
            new_reference_name = std::get<0>(reference_range);
            new_reference_start = std::get<1>(reference_range);
            new_reference_end = std::get<2>(reference_range);
        
        } else {
            new_reference_start = 0;
            new_reference_end = 0;
            new_reference_name = "NA";
        }

        // Now write the new reference coordinates and the rest of the line
        // If we didn't find reference coordinates, write the old coordinates
        if (new_reference_name == "NA" && new_reference_start == 0 && new_reference_end == 0) { 
            std::cout << columns[0] << "\t" << columns[1] << "\t" << columns[2]; 
        } else {
            std::cout << new_reference_name << "\t" << new_reference_start << "\t" << new_reference_end; 
        }

        for (size_t i = 3 ; i < columns.size() ; i++ ) {
            std::cout << "\t" << columns[i];
        }
        std::cout << std::endl;
    }

    instream.close();
}

} // namespace stoat_vcf
