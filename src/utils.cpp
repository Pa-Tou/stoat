#include "utils.hpp"

namespace stoat {

std::string set_precision(const double& value) {
    std::ostringstream oss;
    if (std::abs(value) < 1e-1 && value != 0.0) {
        oss << std::scientific << std::setprecision(4);
    } else {
        oss << std::defaultfloat << std::setprecision(4);
    }

    oss << value;
    return oss.str();
}

std::string set_precision_float_50(const boost::multiprecision::cpp_dec_float_50& value) {
    std::ostringstream oss;

    // Use scientific notation only for very small values
    if (value != 0 && abs(value) < boost::multiprecision::cpp_dec_float_50("1e-1")) {
        oss << std::scientific << std::setprecision(4);
    } else {
        oss << std::defaultfloat << std::setprecision(4);  // Avoid trailing zeros
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

// Function to check significance from a std::string
bool isPValueSignificant(const double& pvalue_threshold, const std::string& pvalue_str) {
    double pvalue;
    try {
        if (pvalue_str == "NA") {
            return false; // Treat "NA" as not significant
        } else {
            pvalue = std::stod(pvalue_str);
        }
    } catch (const std::exception& e) {
        throw std::runtime_error("Error parsing pvalue std::string : " + pvalue_str + " " + e.what());
    }
    return pvalue < pvalue_threshold;
}

// Adjust p-values using Hochberg correction
std::pair<double, size_t> adjusted_hochberg(const std::vector<double>& p_values) {
    int m = p_values.size();

    // Pair each p-value with its original index
    std::vector<std::pair<double, size_t>> indexed;
    indexed.reserve(m);
    for (size_t i = 0; i < m; ++i) {
        indexed.emplace_back(p_values[i], i);
    }

    // Sort by descending p-value
    std::sort(indexed.begin(), indexed.end(), [](const auto& a, const auto& b) {
        return a.first > b.first;
    });

    // Hochberg step-down adjustment
    std::vector<double> adjusted(m);
    for (size_t i = 0; i < m; ++i) {
        adjusted[i] = std::min(1.0, indexed[i].first * (m - i));
    }

    // Enforce monotonicity (non-increasing) from top to bottom
    for (int i = m - 2; i >= 0; --i) {
        adjusted[i] = std::min(adjusted[i], adjusted[i + 1]);
    }

    // Reorder to original order
    std::vector<double> reordered(m);
    for (int i = 0; i < m; ++i) {
        reordered[indexed[i].second] = adjusted[i];
    }

    // Find minimum adjusted p-value and its index in original order
    auto min_iter = std::min_element(reordered.begin(), reordered.end());
    size_t min_index = std::distance(reordered.begin(), min_iter);
    
    return {reordered[min_index], min_index};
}

void retain_indices(std::vector<double>& vec, const std::set<size_t>& indices_to_keep) {
    size_t write_idx = 0;
    for (size_t read_idx = 0; read_idx < vec.size(); ++read_idx) {
        if (indices_to_keep.count(read_idx)) {
            vec[write_idx++] = vec[read_idx];
        }
    }
    vec.resize(write_idx);
}

template std::string vectorToString(const std::vector<std::string>& vec);
template std::string vectorToString(const std::vector<size_t>& vec);

template<typename T>
std::string vectorToString(const std::vector<T>& vec) {
    std::ostringstream oss;
    for (size_t i = 0; i < vec.size(); ++i) {
        if (i > 0) oss << ",";
        oss << vec[i];
    }
    return oss.str();
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

std::string get_sample_name_from_path(const handlegraph::PathPositionHandleGraph& graph, const handlegraph::path_handle_t& path) {

    if (graph.get_sense(path) == handlegraph::PathSense::GENERIC) {
        // Generic paths only have a locus, so return whatever that is
        return graph.get_locus_name(path);
    } else {
        return graph.get_sample_name(path);
    }

}

sample_hap_t get_sample_and_haplotype(const handlegraph::PathPositionHandleGraph& graph, const handlegraph::path_handle_t& path) {
    sample_hap_t result;

    if (graph.get_sense(path) == handlegraph::PathSense::GENERIC) {
        // Generic paths only have a locus, so return whatever that is
        result.sample = graph.get_locus_name(path);
    } else {
        result.sample = graph.get_sample_name(path);
    }
    result.haplotype = graph.get_haplotype(path);

    return result;
}

std::vector<path_range_t> get_coordinates_of_snarl(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                                                   const handlegraph::net_handle_t& snarl, bool get_reference, std::string sample_name, bool get_all_paths) {
    std::vector<path_range_t> ranges;
    // If a sample name is given, then always look for that first
    if (!sample_name.empty() || get_reference) {
        handlegraph::net_handle_t ancestor_snarl = snarl;
        while (!distance_index.is_root(ancestor_snarl)) {
            if (!sample_name.empty()) {
                ranges = get_coordinates_of_snarl_helper(graph, distance_index, ancestor_snarl, false, sample_name, false);
            }
            if (ranges.empty() && get_reference) {
                ranges = get_coordinates_of_snarl_helper(graph, distance_index, ancestor_snarl, true, "", false);
            }
            if (!ranges.empty()) {
                return ranges;
            }
            // If this snarl isn't on the path we want, go up the snarl tree until we find something
            ancestor_snarl = distance_index.get_parent(distance_index.get_parent(ancestor_snarl));
        }
    }
    if (get_all_paths) {
        //If we just want all paths, return that
        ranges = get_coordinates_of_snarl_helper(graph, distance_index, snarl, false, "", true);
        return ranges;

    } else {
        //Try with any path
        handlegraph::net_handle_t start_net = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, true));
        handlegraph::net_handle_t end_net = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, true, true));

        //Get all paths on the start node
        std::set<handlegraph::path_handle_t> start_paths;
        graph.for_each_step_on_handle(distance_index.get_handle(start_net, &graph), [&](const handlegraph::step_handle_t& step) {
            start_paths.insert(graph.get_path_handle_of_step(step));
            return true;
        });

        std::string new_sample_name = "";
        //Now go through the end node and find a path that also goes through the start node
        graph.for_each_step_on_handle(distance_index.get_handle(end_net, &graph), [&](const handlegraph::step_handle_t& step) {
            if (start_paths.count(graph.get_path_handle_of_step(step)) != 0) {
                new_sample_name = graph.get_path_name(graph.get_path_handle_of_step(step));
                return false;
            }
            return true;
        });

        if (new_sample_name.empty()) {
            return ranges;
        } else {
            ranges = get_coordinates_of_snarl_helper(graph, distance_index, snarl, false, new_sample_name, false);
            return ranges;
        }

    }
}

std::vector<path_range_t> get_coordinates_of_snarl_helper(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                                                          const handlegraph::net_handle_t& snarl, bool get_reference, std::string sample_name, bool get_all_paths) {
    #ifdef DEBUG
        cerr << "Get coordinates of " << distance_index.net_handle_as_string(snarl) << endl;
        if (get_reference) {
            assert(sample_name.empty());
            assert(!get_all_paths);
        }
        if (!sample_name.empty()) {
            assert(!get_reference);
            assert(!get_all_paths);
        }
        if (get_all_paths) {
            assert(!get_reference);
            assert(sample_name.empty());
        }
    #endif
    // Bound nodes going into of the snarl
    handlegraph::net_handle_t start_net = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, true));
    handlegraph::net_handle_t end_net = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, true, true));

    // Map path to the steps on the path that traverse the snarl bounds
    std::map<handlegraph::path_handle_t, std::vector<handlegraph::step_handle_t>> path_to_steps;

    // Keep track if we found a traversal of the snarl (may be start-start or end-end)
    bool found_pair = false;

    // Get the step_handles, filtering for the paths we are interested in
    // Steps don't care about the orientation of the handle, they will always (I think) be going forwards in the path
    graph.for_each_step_on_handle(distance_index.get_handle(start_net, &graph), [&] (const handlegraph::step_handle_t& step) {
        handlegraph::path_handle_t path = graph.get_path_handle_of_step(step);
        if ((get_reference && (graph.get_sense(path) == handlegraph::PathSense::REFERENCE)) ||
            (!sample_name.empty() && graph.get_path_name(path).find(sample_name) != std::string::npos) ||
            (get_all_paths)) {
            // If we are looking for a reference path and this is a reference path
            // or if this is the path we want or if we want all paths
            if (path_to_steps.count(path) == 0) {
                path_to_steps[path] = std::vector<handlegraph::step_handle_t>();
            } else {
                found_pair = true;
            }
            path_to_steps[path].emplace_back(step);
        }
        return true;
    });
    #ifdef DEBUG
        cerr << "After start node, found" << endl;
        for (const auto& x : path_to_steps) {
            cerr << graph.get_path_name(x.first) << ": " << x.second.size() << endl;
        }
    #endif
    graph.for_each_step_on_handle(distance_index.get_handle(end_net, &graph), [&] (const handlegraph::step_handle_t& step) {
        handlegraph::path_handle_t path = graph.get_path_handle_of_step(step);
        if ((get_reference && (graph.get_sense(path) == handlegraph::PathSense::REFERENCE)) ||
            (!sample_name.empty() && graph.get_path_name(path).find(sample_name) != std::string::npos) ||
            (get_all_paths)) {
            // If we are looking for a reference path and this is a reference path
            // or if this is the path we are looking for, or we want all paths 
            if (path_to_steps.count(path) == 0) {
                path_to_steps[path] = std::vector<handlegraph::step_handle_t>();
            } else {
                found_pair = true;
            }
            path_to_steps[path].emplace_back(step);
        }
        return true;
    });

    #ifdef DEBUG
        cerr << "After end node, found" << endl;
        for (const auto& x : path_to_steps) {
            cerr << graph.get_path_name(x.first) << ": " << x.second.size() << endl;
        }
    #endif

    vector<path_range_t> ranges;

    if (found_pair) {
        //If we found a path going through the snarl, return the pairs

        for (auto& path_steps : path_to_steps) {
            const handlegraph::path_handle_t& path = path_steps.first;
            std::vector<handlegraph::step_handle_t>& steps = path_steps.second;

            if (steps.size() < 2) {
                continue;
            }
            std::sort(steps.begin(), steps.end(), [&] (const handlegraph::step_handle_t& a, const handlegraph::step_handle_t& b) {
                return graph.get_position_of_step(a) < graph.get_position_of_step(b);
            });

            #ifdef DEBUG
                for (size_t step_i = 0 ; step_i < steps.size() ; step_i++) {
                    if (step_i % 2 == 0) {
                        // If this is an even number, then the path should go into the snarl
                        assert(graph.get_handle_of_step(steps[step_i]) == distance_index.get_handle(start_net, &graph) ||
                            graph.get_handle_of_step(steps[step_i]) == distance_index.get_handle(end_net, &graph));
                    } else {
                        //If this is an odd number, it should go out of the snarl
                        assert(graph.get_handle_of_step(steps[step_i]) == graph.flip(distance_index.get_handle(start_net, &graph)) ||
                            graph.get_handle_of_step(steps[step_i]) == graph.flip(distance_index.get_handle(end_net, &graph)));
                    }
                }
            #endif
            for (size_t i = 0 ; i < steps.size() ; i += 2) {
                ranges.push_back({steps[i], steps[i+1]});
            }
        }
        return ranges;

    } else {
        // If we didn't find anything
        return ranges;
    }

}
std::tuple<std::string, size_t, size_t> get_name_and_offsets_of_snarl_path_range(const handlegraph::PathPositionHandleGraph& graph, 
                                                                                 const bdsg::SnarlDistanceIndex& distance_index, 
                                                                                 const path_range_t& range) {
    return {graph.get_path_name(graph.get_path_handle_of_step(range.start)),
            graph.get_position_of_step(range.start) + graph.get_sequence(graph.get_handle_of_step(range.start)).size(),
            graph.get_position_of_step(range.end)};
}


std::pair<size_t, size_t> find_snarl_id(const bdsg::SnarlDistanceIndex& stree, const handlegraph::net_handle_t& snarl) {

    // Get start and end boundary nodes for the snarl
    auto sstart = stree.get_bound(snarl, false, true);  // False for the left boundary
    auto send = stree.get_bound(snarl, true, true);     // True for the right boundary

    // Convert the sentinels into nodes
    auto start_node = stree.get_node_from_sentinel(sstart);
    auto end_node = stree.get_node_from_sentinel(send);

    // Get the node IDs from bdsg::SnarlDistanceIndex
    // handlegraph::nid_t
    auto start_node_id = stree.node_id(start_node);
    auto end_node_id = stree.node_id(end_node);

    // Convert to size_t
    size_t start_node_id_size_t = static_cast<size_t>(start_node_id);
    size_t end_node_id_size_t = static_cast<size_t>(end_node_id);

    // Construct the snarl ID
    std::pair<size_t, size_t> snarl_id(start_node_id_size_t, end_node_id_size_t);

    return snarl_id;  // Return the generated snarl ID as a std::string
}

void print_nodes_in_snarl(const bdsg::SnarlDistanceIndex& distance_index, const handlegraph::net_handle_t& snarl) {
    std::vector<handlegraph::net_handle_t> to_print;
    to_print.emplace_back(snarl);
    while (!to_print.empty()) {
        handlegraph::net_handle_t net = std::move(to_print.back());
        to_print.pop_back();

        if (distance_index.is_node(net)) {
            cerr << distance_index.node_id(net) << endl;
        } else {
            distance_index.for_each_child(net, [&](const handlegraph::net_handle_t& child) {
                to_print.emplace_back(child);
                return true;
            });
        }
    }
}

} // end namespace stoat
