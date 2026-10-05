#include "utils.hpp"

//#include DEBUG

namespace stoat {

std::string trim(const std::string& value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

size_t parse_count(const std::string& count_text, const std::string& context) {
    size_t count = 0;
    const char* first = count_text.data();
    const char* last = first + count_text.size();
    const auto [ptr, ec] = std::from_chars(first, last, count);
    if (count_text.empty() || ec != std::errc{} || ptr != last || count == 0) {
        throw std::invalid_argument("Haplotype count must be a positive integer" + context);
    }
    return count;
}

std::string set_precision(const double& value) {

    if (std::isnan(value)) {
        return "NA";
    }

    std::ostringstream oss;
    if (std::abs(value) < 1e-1 && value != 0.0) {
        oss << std::scientific << std::setprecision(4);
    } else {
        oss << std::defaultfloat << std::setprecision(4);
    }

    oss << value;
    return oss.str();
}

bool is_na(const std::string& s) {
    return s.empty() || s == "NA";
}

double string_to_pvalue(const std::string& p1) {
    bool na1 = is_na(p1);

    if (!na1) {
        return std::stod(p1);
    } else {
        return 1.0;
    }
}

// Adjust p-values using Hochberg correction
std::pair<double, size_t> adjusted_hochberg(const std::vector<double>& p_values) {
    size_t m = p_values.size();

    // Pair each p-value with its original index
    std::vector<std::pair<double, size_t>> indexed;
    indexed.reserve(m);
    for (size_t i = 0; i < m; ++i) {
        indexed.emplace_back(p_values[i], i);
    }

    // Sort by ASCENDING p-value (Hochberg requires ascending sort)
    std::sort(indexed.begin(), indexed.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });

    // Apply Hochberg step-down correction
    std::vector<double> adjusted(m);
    for (int i = m - 1; i >= 0; --i) {
        size_t rank = m - i;
        adjusted[i] = indexed[i].first * rank;
        if (i < static_cast<int>(m - 1)) {
            adjusted[i] = std::min(adjusted[i], adjusted[i + 1]);
        }
        adjusted[i] = std::min(adjusted[i], 1.0);
    }

    // Reorder adjusted p-values back to original order
    std::vector<double> reordered(m);
    for (size_t i = 0; i < m; ++i) {
        reordered[indexed[i].second] = adjusted[i];
    }

    // Return minimum adjusted p-value and its index in original order
    auto min_iter = std::min_element(reordered.begin(), reordered.end());
    size_t min_index = std::distance(reordered.begin(), min_iter);

    return {*min_iter, min_index};
}

template std::vector<std::string> stringToVector(const std::string& vec);
template std::vector<size_t> stringToVector(const std::string& vec);

template <typename T>
std::vector<T> stringToVector(const std::string& str) {
    std::vector<T> result;
    std::istringstream iss(str);
    std::string token;

    while (std::getline(iss, token, ',')) {
        std::istringstream tokenStream(token);
        T value;
        tokenStream >> value;
        if (tokenStream.fail()) {
            throw std::runtime_error("Failed to parse token: " + token);
        }
        result.push_back(value);
    }

    return result;
}

std::string get_sample_name_from_path(const handlegraph::PathHandleGraph& graph, const handlegraph::path_handle_t& path) {

    if (graph.get_sense(path) == handlegraph::PathSense::GENERIC) {
        // Generic paths only have a locus, so return whatever that is
        return graph.get_locus_name(path);
    } else {
        return graph.get_sample_name(path);
    }

}


void print_nodes_in_snarl(const bdsg::SnarlDistanceIndex& distance_index, const handlegraph::net_handle_t& snarl) {
    std::vector<handlegraph::net_handle_t> to_print;
    to_print.emplace_back(snarl);
    while (!to_print.empty()) {
        handlegraph::net_handle_t net = std::move(to_print.back());
        to_print.pop_back();

        if (distance_index.is_node(net)) {
            stoat::LOG_INFO(std::to_string(distance_index.node_id(net)));
        } else {
            distance_index.for_each_child(net, [&](const handlegraph::net_handle_t& child) {
                to_print.emplace_back(child);
                return true;
            });
        }
    }
}


} // end namespace stoat

