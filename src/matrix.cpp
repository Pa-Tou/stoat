#include "matrix.hpp"

namespace stoat_vcf {

// Constructor implementation
EdgeBySampleMatrix::EdgeBySampleMatrix(const std::vector<std::string>& sample_names, size_t n_edges) : sample_names(sample_names), max_edges(n_edges) {

    // assuming (at max) two copies per individuals
    n_samp_haps = sample_names.size() * 2;

    // nothing to do if no edges or no samples
    if (n_edges == 0 || n_samp_haps == 0) {
        return;
    }

    // Initialize with "zeros"
    matrix_1D.resize(n_edges * n_samp_haps, false); 
    row_header.rehash(n_edges); // JEAN not sure if that's useful?
}

std::string EdgeBySampleMatrix::get_sample_name(size_t sample_idx) const {
    return (sample_names[sample_idx]);
}

    
// Add True to the matrix if edge is found
void EdgeBySampleMatrix::add_sample_edge(const stoat::edge_t& edge, size_t col_index) {
    // look for that edge in the matrix
    auto it = row_header.find(edge);
    if (it != row_header.end()) {
        // it's present, set that edge for the specified sample-hap
        set_edge(it->second, col_index);
    } else {
        // add a new edge to the matrix
        size_t new_idx = row_header.size();
        row_header[edge] = new_idx;
        // expand the matrix if we used up all the allocated rows
        if (new_idx >= max_edges) {
            expand_matrix();
        }
        // set that edge for the specified sample-hap
        set_edge(new_idx, col_index);
    }
}

// Double the number of elements in the matrix
void EdgeBySampleMatrix::expand_matrix() {
    max_edges *= 2;  
    matrix_1D.resize(matrix_1D.size() * 2, false); // Initialize new memory with zeros
}

// Overloaded operator() to access elements as matrix(row, col)
bool EdgeBySampleMatrix::operator()(size_t row, size_t col) const {
    return (matrix_1D.at(row * n_samp_haps + col));
}

bool EdgeBySampleMatrix::get_edge(size_t row, size_t col) const {
    return (matrix_1D.at(row * n_samp_haps + col));
}

    
// Function to set a specific element (row, col) to true
void EdgeBySampleMatrix::set_edge(size_t row, size_t col) {
    matrix_1D.at(row * n_samp_haps + col) = true;
}

void EdgeBySampleMatrix::shrink() {
    matrix_1D.resize(row_header.size() * n_samp_haps); // Resize
    matrix_1D.shrink_to_fit(); // Free unused capacity
}

void EdgeBySampleMatrix::clear_edges(size_t new_n_edges) { 
    matrix_1D.clear();
    row_header.clear();

    if (new_n_edges == 0 || n_samp_haps == 0) {
        return;
    }

    max_edges = new_n_edges;
    // Initialize with "zeros"
    matrix_1D.resize(new_n_edges * n_samp_haps, false); 
    row_header.rehash(new_n_edges); // JEAN not sure if that's useful?
}

// Function to identify the path in the edge matrix
std::vector<size_t> EdgeBySampleMatrix::get_samples_on_path(const stoat::PathTraversal &path_trav) const {

    const std::vector<stoat::node_traversal_t> &path = path_trav.get_path();
    
    // get the subset of rows to check for that path and its flipped version
    // we'll check N edges, N being the number of nodes in the path - 1
    size_t path_len = path.size();
    std::vector<size_t> rows_to_check;
    rows_to_check.reserve(path_len - 1);
    std::vector<size_t> rows_to_check_flipped;
    rows_to_check_flipped.reserve(path_len - 1);

    // look for each edge in the matrix
    bool skip_path = false;
    bool skip_flipped_path = false;
    for (size_t i = 0; i < path_len - 1; ++i) {
        // Skip if snarl contains '*' (here * == 0) aka nested element
        if (path[i].get_node_id() == 0 || path[i+1].get_node_id() == 0) {
            continue;
        }

        // make an edge
        stoat::edge_t edge(path[i], path[i + 1]);

        // if we can find that edge, save its index in the edge matrix
        auto itr = row_header.find(edge);
        if (itr == row_header.end()) {
            // if at least one edge not found, abort early and skip this path below
            skip_path = true; 
        } else {
            rows_to_check.push_back(itr->second);
        }
        // check the flipped edge
        stoat::edge_t edge_flipped = edge.get_flipped();
        auto itr_flipped = row_header.find(edge_flipped);
        if (itr_flipped == row_header.end()) {
            // if at least one edge not found, abort early and skip this path below
            skip_flipped_path = true; 
        } else {
            rows_to_check_flipped.push_back(itr_flipped->second);
        }
    }

    // now let's check the presence of the edges in each sample-hap
    std::vector<size_t> idx_samp_hap;
    idx_samp_hap.reserve(n_samp_haps);

    // loop by columns first (better cache locality in the matrix)
    bool all_ones;
    for (size_t col = 0; col < n_samp_haps; ++col) {
        // check path
        if (!skip_path) {
            all_ones = true;
            for (size_t row : rows_to_check) {
                if (!get_edge(row, col)) {
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
                if (!get_edge(row, col)) {
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

void EdgeBySampleMatrix::load_vcf_chunk(stoat_vcf::VCFParser& vcf_parser, std::string &chr) {
    // init the edge matrix, allocating about 10 000 edges?
    clear_edges(10000);

    vcf_parser.for_each_record_on_chromosome(chr, [&](const stoat_vcf::vcf_info_t& vcf_info) {
        if (vcf_info.lv == 0 || vcf_parser.resolve_nested_calls ) {
            // If we resolve snarls, then keep edges in all snarls. If we aren't resolving snarls, then only keep lv=0 snarls
            for (size_t hap_num = 0 ; hap_num < vcf_parser.hap_count ; hap_num++) {
                int allele_num = vcf_info.genotype[hap_num];
                if (allele_num != -1) { // if the genotype wasn't .
                    const std::vector<stoat::node_traversal_t>& path = vcf_info.paths[allele_num];
                    for (size_t node_i = 0 ; node_i < path.size()-1 ; node_i++) {
                        // Go through each edge (as pair of nodes) and add it to the edge matrix
                        // ignoring any 0 nodes which indicate the inside of a snarl
                        if (path[node_i].get_node_id() != 0 && path[node_i+1].get_node_id() != 0) {
                            add_sample_edge(stoat::edge_t(path[node_i], path[node_i+1]), hap_num);
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
