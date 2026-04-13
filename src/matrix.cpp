#include "matrix.hpp"

namespace stoat_vcf {

// Constructor implementation
template<class ValueType>
AlleleBySampleMatrix<ValueType>::AlleleBySampleMatrix(const std::vector<std::string>& sample_names, 
    size_t n_keys) : 
    sample_names(sample_names), 
    max_keys(n_keys) {

    // assuming (at max) two copies per individuals
    n_samp_haps = sample_names.size() * 2;

    // nothing to do if no keys or no samples
    if (n_keys == 0 || n_samp_haps == 0) {
        return;
    }

    // Initialize with "zeros"
    matrix_1D.resize(n_keys * n_samp_haps, false); 
    row_header.rehash(n_keys); // JEAN not sure if that's useful?
}

template<class ValueType>
std::string AlleleBySampleMatrix<ValueType>::get_sample_name(size_t sample_idx) const {
    return (sample_names[sample_idx]);
}

    
// Add True to the matrix if key is found
template<class ValueType>
void AlleleBySampleMatrix<ValueType>::add_sample_key(const ValueType& key, size_t col_index) {
    // look for that key in the matrix
    auto it = row_header.find(key);
    if (it != row_header.end()) {
        // it's present, set that key for the specified sample-hap
        set_key(it->second, col_index);
    } else {
        // add a new key to the matrix
        size_t new_idx = row_header.size();
        row_header[key] = new_idx;
        // expand the matrix if we used up all the allocated rows
        if (new_idx >= max_keys) {
            expand_matrix();
        }
        // set that key for the specified sample-hap
        set_key(new_idx, col_index);
    }
}

// Double the number of elements in the matrix
template<class ValueType>
void AlleleBySampleMatrix<ValueType>::expand_matrix() {
    max_keys *= 2;  
    matrix_1D.resize(matrix_1D.size() * 2, false); // Initialize new memory with zeros
}

// Overloaded operator() to access elements as matrix(row, col)
template<class ValueType>
bool AlleleBySampleMatrix<ValueType>::operator()(size_t row, size_t col) const {
    return (matrix_1D.at(row * n_samp_haps + col));
}

template<class ValueType>
bool AlleleBySampleMatrix<ValueType>::get_key(size_t row, size_t col) const {
    return (matrix_1D.at(row * n_samp_haps + col));
}

    
// Function to set a specific element (row, col) to true
template<class ValueType>
void AlleleBySampleMatrix<ValueType>::set_key(size_t row, size_t col) {
    matrix_1D.at(row * n_samp_haps + col) = true;
}

template<class ValueType>
void AlleleBySampleMatrix<ValueType>::shrink() {
    matrix_1D.resize(row_header.size() * n_samp_haps); // Resize
    matrix_1D.shrink_to_fit(); // Free unused capacity
}

template<class ValueType>
void AlleleBySampleMatrix<ValueType>::clear_keys(size_t new_n_keys) { 
    matrix_1D.clear();
    row_header.clear();

    if (new_n_keys == 0 || n_samp_haps == 0) {
        return;
    }

    max_keys = new_n_keys;
    // Initialize with "zeros"
    matrix_1D.resize(new_n_keys * n_samp_haps, false); 
    row_header.rehash(new_n_keys); // JEAN not sure if that's useful?
}

// Function to identify the path in the key matrix
template<>
std::vector<size_t> AlleleBySampleMatrix<stoat::edge_t>::get_samples_on_path(const stoat::PathTraversal &path_trav) const {

    const std::vector<stoat::node_traversal_t> &path = path_trav.get_path();
    
    // get the subset of rows to check for that path and its flipped version
    // we'll check N keys, N being the number of nodes in the path - 1
    size_t path_len = path.size();
    std::vector<size_t> rows_to_check;
    rows_to_check.reserve(path_len - 1);
    std::vector<size_t> rows_to_check_flipped;
    rows_to_check_flipped.reserve(path_len - 1);

    // look for each key in the matrix
    bool skip_path = false;
    bool skip_flipped_path = false;
    for (size_t i = 0; i < path_len - 1; ++i) {
        // Skip if snarl contains '*' (here * == 0) aka nested element
        if (path[i].get_node_id() == 0 || path[i+1].get_node_id() == 0) {
            continue;
        }

        // make an key
        ValueType key(path[i], path[i + 1]);

        // if we can find that key, save its index in the key matrix
        auto itr = row_header.find(key);
        if (itr == row_header.end()) {
            // if at least one key not found, abort early and skip this path below
            skip_path = true; 
        } else {
            rows_to_check.push_back(itr->second);
        }
        // check the flipped key
        ValueType key_flipped = key.get_flipped();
        auto itr_flipped = row_header.find(key_flipped);
        if (itr_flipped == row_header.end()) {
            // if at least one key not found, abort early and skip this path below
            skip_flipped_path = true; 
        } else {
            rows_to_check_flipped.push_back(itr_flipped->second);
        }
    }

    // now let's check the presence of the keys in each sample-hap
    std::vector<size_t> idx_samp_hap;
    idx_samp_hap.reserve(n_samp_haps);

    // loop by columns first (better cache locality in the matrix)
    bool all_ones;
    for (size_t col = 0; col < n_samp_haps; ++col) {
        // check path
        if (!skip_path) {
            all_ones = true;
            for (size_t row : rows_to_check) {
                if (!get_key(row, col)) {
                    all_ones = false;
                    break;
                }
            }
            if (all_ones) {
                // JEAN why this static cast?
                idx_samp_hap.push_back(static_cast<int>(col));
                continue;
            }
        }
        // if we haven't found it, we should check the flipped path too
        if (!skip_flipped_path) {
            all_ones = true;
            for (size_t row : rows_to_check_flipped) {
                if (!get_key(row, col)) {
                    all_ones = false;
                    break;
                }
            }
            if (all_ones) {
                // JEAN why this static cast?
                idx_samp_hap.push_back(static_cast<int>(col));
            }
        }
    }
    return idx_samp_hap;
}

// Function to identify the path in the key matrix
template<>
std::vector<size_t> AlleleBySampleMatrix<stoat::node_traversal_t>::get_samples_on_node(const stoat::node_traversal_t &node_trav) const {

    size_t node_row = 0;
    auto itr = row_header.find(node_trav);
    if (itr == row_header.end()) {
        stoat::LOG_WARN("node traversal : " + node_trav.to_string() + "not found", number_node_not_found++);
    };

    // now let's check the presence of the keys in each sample-hap
    std::vector<size_t> idx_samp_hap;
    idx_samp_hap.reserve(n_samp_haps);

    // loop by columns first (better cache locality in the matrix)
    bool all_ones;
    for (size_t col = 0; col < n_samp_haps; ++col) {
        all_ones = true;

        if (!get_key(node_row, col)) {
            all_ones = false;
            break;
        }

        if (all_ones) {
            idx_samp_hap.push_back(static_cast<int>(col));
            continue;
        }
    }

    return idx_samp_hap;
}

//     // const std::vector<stoat::node_traversal_t> &path = path_trav.get_path();
    
//     // get the subset of rows to check for that path and its flipped version
//     // we'll check N keys, N being the number of nodes in the path - 1
//     size_t path_len = path.size();
//     std::vector<size_t> rows_to_check;
//     rows_to_check.reserve(path_len - 1);
//     std::vector<size_t> rows_to_check_flipped;
//     rows_to_check_flipped.reserve(path_len - 1);

//     // look for each key in the matrix
//     bool skip_path = false;
//     bool skip_flipped_path = false;
//     for (size_t i = 0; i < path_len - 1; ++i) {
//         // Skip if snarl contains '*' (here * == 0) aka nested element
//         if (path[i].get_node_id() == 0 || path[i+1].get_node_id() == 0) {
//             continue;
//         }

//         // make an key
//         ValueType key(path[i], path[i + 1]);

//         // if we can find that key, save its index in the key matrix
//         auto itr = row_header.find(key);
//         if (itr == row_header.end()) {
//             // if at least one key not found, abort early and skip this path below
//             skip_path = true; 
//         } else {
//             rows_to_check.push_back(itr->second);
//         }
//         // check the flipped key
//         ValueType key_flipped = key.get_flipped();
//         auto itr_flipped = row_header.find(key_flipped);
//         if (itr_flipped == row_header.end()) {
//             // if at least one key not found, abort early and skip this path below
//             skip_flipped_path = true; 
//         } else {
//             rows_to_check_flipped.push_back(itr_flipped->second);
//         }
//     }

//     // now let's check the presence of the keys in each sample-hap
//     std::vector<size_t> idx_samp_hap;
//     idx_samp_hap.reserve(n_samp_haps);

//     // loop by columns first (better cache locality in the matrix)
//     bool all_ones;
//     for (size_t col = 0; col < n_samp_haps; ++col) {
//         // check path
//         if (!skip_path) {
//             all_ones = true;
//             for (size_t row : rows_to_check) {
//                 if (!get_key(row, col)) {
//                     all_ones = false;
//                     break;
//                 }
//             }
//             if (all_ones) {
//                 // JEAN why this static cast?
//                 idx_samp_hap.push_back(static_cast<int>(col));
//                 continue;
//             }
//         }
//         // if we haven't found it, we should check the flipped path too
//         if (!skip_flipped_path) {
//             all_ones = true;
//             for (size_t row : rows_to_check_flipped) {
//                 if (!get_key(row, col)) {
//                     all_ones = false;
//                     break;
//                 }
//             }
//             if (all_ones) {
//                 // JEAN why this static cast?
//                 idx_samp_hap.push_back(static_cast<int>(col));
//             }
//         }
//     }
//     return idx_samp_hap;
// }

template<class ValueType>
void AlleleBySampleMatrix<ValueType>::load_vcf_chunk(stoat_vcf::VCFParser& vcf_parser, std::string &chr) {
    // init the key matrix, allocating about 10 000 keys?
    clear_keys(10000);

    vcf_parser.for_each_record_on_chromosome(chr, [&](const stoat_vcf::vcf_info_t& vcf_info) {
        if (vcf_info.lv == 0 || vcf_parser.resolve_nested_calls ) {
            // If we resolve snarls, then keep keys in all snarls. If we aren't resolving snarls, then only keep lv=0 snarls
            for (size_t hap_num = 0 ; hap_num < vcf_parser.hap_count ; hap_num++) {
                int allele_num = vcf_info.genotype[hap_num];
                if (allele_num != -1) { // if the genotype wasn't .
                    const std::vector<stoat::node_traversal_t>& path = vcf_info.paths[allele_num];
                    for (size_t node_i = 0 ; node_i < path.size()-1 ; node_i++) {
                        // Go through each key (as pair of nodes) and add it to the key matrix
                        // ignoring any 0 nodes which indicate the inside of a snarl
                        if (path[node_i].get_node_id() != 0 && path[node_i+1].get_node_id() != 0) {
                            add_sample_key(ValueType(path[node_i], path[node_i+1]), hap_num);
                        }
                    }
                }
            }
        }
    });

    // we're done, we can shrink the matrix to rows set
    shrink();

}
    
} // end namespace stoat

// Apparently these definitions are supposed to be done here, after all the members are defined
template class stoat_vcf::AlleleBySampleMatrix<stoat::edge_t>;
template class stoat_vcf::AlleleBySampleMatrix<stoat::node_traversal_t>;
