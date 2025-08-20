#include "post_processing.hpp"
#include "utils.hpp"

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

    add_BH_adjusted_column(input_file, output_dir, output_file_significant, adjusted_col_index-1, adjusted_col_index);
}

// Main Function
void add_BH_adjusted_column(
    const std::string& input_file,
    const std::string& output_dir,
    const std::string& output_file_significant,
    size_t p_col_index, size_t adjusted_col_index) {

    std::ifstream infile(input_file);
    std::string col;

    // First pass: Collect p-values
    std::vector<std::tuple<double, double, size_t>> pvalues;
    std::string line;
    size_t line_index = 0;

    // Read the header line
    std::string header_line;
    std::getline(infile, header_line);
    std::stringstream header_ss(header_line);
    std::vector<std::string> headers;
    while (std::getline(header_ss, col, '\t')) {
        headers.push_back(col);
    }

    while (std::getline(infile, line)) {
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
    infile.close();

    // Apply BH correction
    adjust_pvalues_with_BH(pvalues);

    // Second pass: rewrite with BH-adjusted values
    infile.open(input_file);
    const std::string output_temp_file = output_dir + "/temp_output.tsv";
    std::ofstream outfile(output_temp_file);
    std::ofstream outfile_significant(output_file_significant);

    // Write headers
    outfile << header_line << '\n';
    outfile_significant << header_line << '\n';

    std::getline(infile, line); // Skip header again
    line_index = 0;

    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> columns;

        while (std::getline(ss, token, '\t')) {
            columns.push_back(token);
        }

        double adjusted_p = std::get<1>(pvalues[line_index]);
        std::string adj_str = stoat::set_precision(adjusted_p);
        columns[adjusted_col_index] = adj_str;

        // Write updated line
        for (size_t i = 0; i < columns.size(); ++i) {
            outfile << columns[i];
            if (i != columns.size() - 1) outfile << '\t';
        }

        outfile << '\n';

        if (adjusted_p < 1e-5) {
            for (size_t i = 0; i < columns.size(); ++i) {
                outfile_significant << columns[i];
                if (i != columns.size() - 1) outfile_significant << '\t';
            }
            outfile_significant << '\n';
        }
        ++line_index;
    }

    infile.close();
    outfile.close();
    outfile_significant.close();

    // Replace original file
    std::remove(input_file.c_str());
    std::rename(output_temp_file.c_str(), input_file.c_str());
}

} // namespace stoat_vcf
