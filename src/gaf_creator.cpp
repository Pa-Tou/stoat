#include "gaf_creator.hpp"

namespace stoat_vcf {

std::pair<double, double> calcul_proportion_signi(size_t number_ind_group0, size_t number_ind_group1, double p_value) {

    // Step 1: Calculate initial proportions based on a total of 60
    size_t total_ind = number_ind_group0 + number_ind_group1;
    if (total_ind == 0) {
        return {0, 0};
    }

    double proportion_group0 = (static_cast<double>(number_ind_group1) / total_ind) * 60.0;
    double proportion_group1 = 60.0 - proportion_group0;

    // Step 2: Calculate the adjustment factor based on a logarithmic scale of p_value
    constexpr double epsilon = 1e-20;  // Small value to avoid log(0)
    double adjustment_factor = -std::log(std::max(p_value, epsilon));

    // Step 3: Apply the adjustment to the group with the higher initial proportion
    double adjusted_group0, adjusted_group1;
    if (proportion_group0 > proportion_group1) {
        adjusted_group0 = proportion_group0 + adjustment_factor;
        adjusted_group1 = proportion_group1 - adjustment_factor;
    } else {
        adjusted_group0 = proportion_group0 - adjustment_factor;
        adjusted_group1 = proportion_group1 + adjustment_factor;
    }

    // Step 4: Ensure values remain within bounds [0, 60] and maintain a sum of 60
    adjusted_group0 = std::clamp(adjusted_group0, 0.0, 60.0);
    adjusted_group1 = std::clamp(adjusted_group1, 0.0, 60.0);

    // Rescale if needed to ensure the total sum is exactly 60
    double total = adjusted_group0 + adjusted_group1;
    if (total != 60.0) {
        double scale_factor = 60.0 / total;
        adjusted_group0 *= scale_factor;
        adjusted_group1 *= scale_factor;
    }

    return {adjusted_group0, adjusted_group1};
}

std::string addSuffixToFilename(const std::string& filename, const std::string& suffix) {
    size_t lastDotPos = filename.find_last_of(".");
    if (lastDotPos == std::string::npos) {
        // No extension found, return the filename with suffix appended
        return filename + suffix;
    }
    
    std::string base = filename.substr(0, lastDotPos);
    std::string ext = filename.substr(lastDotPos);
    
    return base + suffix + ext;
}

void writeGafLines(const std::string& sequenceName, const std::string& path, 
                    int length, int proportion, std::ofstream& outFile) {
    if (!outFile.is_open()) {
        throw std::runtime_error("Error: Output file stream is not open.");
    }
    
    outFile << sequenceName << "\t" << length << "\t0\t" << length << "\t+\t"
            << path << "\t" << length << "\t0\t" << length << "\t"
            << length << "\t" << length << "\t" << proportion 
            << "\tcs:Z::" << length << "\n";
}

// Adds a suffix to a filename before the file extension
string add_suffix_to_filename(const std::string& filename, const std::string& suffix) {
    size_t dotPos = filename.find_last_of(".");
    if (dotPos == std::string::npos) {
        return filename + suffix;
    }
    return filename.substr(0, dotPos) + suffix + filename.substr(dotPos);
}

// TODO : change this to use already implemented function in utils.hpp
std::vector<size_t> decompose_snarl(const std::string& snarl) {
    std::vector<size_t> snarl_node;
    regex re("\\d+");
    sregex_iterator begin(snarl.begin(), snarl.end(), re), end;
    
    for (auto it = begin; it != end; ++it) {
        snarl_node.push_back(stoi(it->str()));
    }
    return snarl_node;
}

int calcul_path_length(handlegraph::PathHandleGraph& graph, const std::string& snarl) {
    std::vector<size_t> snarl_nodes = decompose_snarl(snarl);
    int length_node = 0;

    for (size_t node : snarl_nodes) {
        handlegraph::handle_t handle = graph.get_handle(node);
        length_node += graph.get_length(handle);
    }
    return length_node;
}

// Writes a formatted line to the output file
void write_gaf_lines(const std::string& sequence_name, const std::string& path, int length, double prop, ofstream& outfile) {
    outfile << sequence_name << "\t" << path << "\t" << length << "\t" << prop << "\n";
}

// Parses the input file and processes data into two output files
void gaf_creation(const std::string& input_file, 
    std::unordered_map<std::string, std::vector<stoat::Snarl_data_t>>& snarl_chr,
    handlegraph::PathHandleGraph& graph, 
    const std::string& output_file) {

    std::string output_file_1 = add_suffix_to_filename(output_file, "_0");
    std::string output_file_2 = add_suffix_to_filename(output_file, "_1");

    ifstream infile(input_file);
    ofstream outfile1(output_file_1), outfile2(output_file_2);
    if (!infile || !outfile1 || !outfile2) {
        throw std::runtime_error("Error opening files");
    }

    std::string line;
    size_t count_line = 0;
    getline(infile, line); // Skip header

    while (getline(infile, line)) {
        stringstream ss(line);
        std::vector<std::string> columns;
        std::string token;
        while (ss >> token) {
            columns.push_back(token);
        }

        if (columns.size() < 14) continue;

        const std::string& chr = columns[0];
        std::string snarl_list = columns[2];
        double pfisher = stoat::string_to_pvalue(columns[6]);
        double pchi = stoat::string_to_pvalue(columns[7]);
        std::string group_paths = columns[11];
        auto it = snarl_chr.find(chr);
        auto& data = it->second;  
        const std::vector<std::string>& list_path = stoat::stringToVector<std::string>(stoat::vectorPathToString(data[count_line].snarl_paths));

        // Split group paths by comma
        std::vector<std::string> decomposed_group_paths;
        stringstream gp_ss(group_paths);
        while (getline(gp_ss, token, ',')) {
            decomposed_group_paths.push_back(token);
        }

        std::vector<std::string> group_0, group_1, sequence_name_g0, sequence_name_g1;
        for (const auto& path_group : decomposed_group_paths) {
            size_t pos = path_group.find(':');
            if (pos == std::string::npos) continue;
            group_0.push_back(path_group.substr(0, pos));
            group_1.push_back(path_group.substr(pos + 1));
            sequence_name_g0.push_back(snarl_list + "_G0_" + group_0.back() + "_F" + to_string(pfisher) + "_C" + to_string(pchi));
            sequence_name_g1.push_back(snarl_list + "_G1_" + group_1.back() + "_F" + to_string(pfisher) + "_C" + to_string(pchi));
        }

        for (size_t idx = 0; idx < list_path.size(); ++idx) {
            const std::string& path = list_path[idx];

            // Case where "*" is in path
            // TODO : update to the snarl_data_t.hpp implementation
            if (path.find("0") != std::string::npos) {
                // If the path contains '*', split it into two sub-paths
                size_t star_pos = path.find("0");
                std::string star_path_1 = path.substr(0, star_pos - 1);
                std::string star_path_2 = path.substr(star_pos + 1);

                int len1 = calcul_path_length(graph, star_path_1);
                int len2 = calcul_path_length(graph, star_path_2);
                auto [prop_g0, prop_g1] = calcul_proportion_signi(stoi(group_0[idx]), stoi(group_1[idx]), pfisher);

                // Write lines for start_path_1
                write_gaf_lines(sequence_name_g0[idx], star_path_1, len1, prop_g0, outfile1);
                write_gaf_lines(sequence_name_g1[idx], star_path_1, len1, prop_g1, outfile2);
                
                // Write lines for start_path_2
                write_gaf_lines(sequence_name_g0[idx], star_path_2, len2, prop_g0, outfile1);
                write_gaf_lines(sequence_name_g1[idx], star_path_2, len2, prop_g1, outfile2);
            } else {
                // Case where "*" is NOT in path
                int len = calcul_path_length(graph, path);
                auto [prop_g0, prop_g1] = calcul_proportion_signi(stoi(group_0[idx]), stoi(group_1[idx]), pfisher);
                write_gaf_lines(sequence_name_g0[idx], path, len, prop_g0, outfile1);
                write_gaf_lines(sequence_name_g1[idx], path, len, prop_g1, outfile2);
            }
        }
        count_line += 1;
    }
}

} //end stoat namespace
