#include "load_tables.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>

bool binary_table_values_t::operator==(const binary_table_values_t& other) const {
    if (chr != other.chr ||
        start_pos != other.start_pos ||
        end_pos != other.end_pos ||
        depth != other.depth) {
        return false;
    }

    // For the snarl, the bounds can be in either order
    // TODO: This should maybe take into account the orientation
    if (!((snarl_ids.first == other.snarl_ids.first && snarl_ids.second == other.snarl_ids.second) ||
        (snarl_ids.first == other.snarl_ids.second && snarl_ids.second == other.snarl_ids.first))) {
        return false;
    }

    // Check the p-values. TODO: Do this better
    if (p_fishers != other.p_fishers ||
        p_chi2 != other.p_chi2) {
        return false;
    }

    // The path_lengths and group_paths 
    if (path_lengths.size() != other.path_lengths.size() ||
        group_paths.size() != other.group_paths.size()) {
        return false;
    }

    // Put the corresponding path length and group_paths together in a set
    std::unordered_set<std::string> path_lengths_and_groups;
    for (size_t i = 0 ; i < path_lengths.size() ; i++) {
        path_lengths_and_groups.emplace(path_lengths[i] + " " + group_paths[i]);
    }
    std::unordered_set<std::string> other_path_lengths_and_groups;
    for (size_t i = 0 ; i < other.path_lengths.size() ; i++) {
        other_path_lengths_and_groups.emplace(other.path_lengths[i] + " " + other.group_paths[i]);
    }
    // the sets must be the same
    if (path_lengths_and_groups != other_path_lengths_and_groups) {
        return false;
    }

    return true;
}


//CHR\tSTART_POS\tEND_POS\tSNARL\tPATH_LENGTHS\tP_FISHER\tP_CHI2\tGROUP_PATHS\tDEPTH
binary_table_values_t load_binary_snarl_line(std::string line) {


    binary_table_values_t vals;
    // Components:
    //std::string chr;
    //size_t start_pos;
    //size_t end_pos;
    //std::pair<size_t, size_t> snarl_ids;
    //std::vector<std::pair<size_t, size_t>> path_lengths;
    //std::string p_fishers;
    //std::string p_chi2;
    //std::vector<std::pair<size_t, size_t>> group_paths;
    //size_t depth;

    std::stringstream linestream(line);
    std::string part;
    // chromosome/sample name
    std::getline(linestream, part, '\t');
    vals.chr = part;

    // start pos
    std::getline(linestream, part, '\t');
    vals.start_pos = std::stoull(part);

    // end pos
    std::getline(linestream, part, '\t');
    vals.end_pos = std::stoull(part);

    // start id
    std::getline(linestream, part, '_');
    vals.snarl_ids.first = std::stoull(part);

    // end id
    std::getline(linestream, part, '\t');
    vals.snarl_ids.second = std::stoull(part);

    // path lengths, comma separated by allele, /-separated for min/max in one allele
    std::getline(linestream, part, '\t');

    std::stringstream pathlength_stream(part);
    std::string lengths;
    while (std::getline(pathlength_stream, lengths, ',')){
        vals.path_lengths.emplace_back(lengths);
    }

    // Fishers p-value
    std::getline(linestream, part, '\t');
    vals.p_fishers = part;

    // chi2 p-value
    std::getline(linestream, part, '\t');
    vals.p_chi2 = part;

    // group paths, comma separated by allele, :-separated for the counts per phenotype group
    std::getline(linestream, part, '\t');

    std::stringstream grouppath_stream(part);
    std::string counts;
    while (std::getline(grouppath_stream, counts, ',')){
        vals.group_paths.emplace_back(counts);
    }

    // depth
    std::getline(linestream, part, '\t');
    vals.depth = std::stoull(part);

    return vals;
}
