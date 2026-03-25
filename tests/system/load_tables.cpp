#include "load_tables.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>

bool snarl_genotype_values_t::operator==(const snarl_genotype_values_t& other) const {
    if (!((start_node == other.start_node && end_node == other.end_node) || (end_node == other.start_node && start_node == other.end_node)) ||
        reference != other.reference ||
        start_offset != other.start_offset ||
        end_offset != other.end_offset ||
        depth != other.depth ||
        allele_lengths.size() != other.allele_lengths.size() ||
        walks.size() != other.walks.size() ||
        sequences.size() != other.sequences.size() ||
        sample_to_allele.size() != other.sample_to_allele.size()) {
        return false;
    }
    // Check that the alleles are all the same.
    // Alleles may be in a different order, so try to map from this object's alleles to the other's alleles
    std::unordered_map<size_t, size_t> allele_conversion;


    // Go through all i walks and all j other walks and see what matches, and make sure that everything has a match
    if (walks.size() != 0) {
        std::vector<bool> found_j(other.walks.size() , false);
        for (size_t i = 0 ; i < walks.size() ; i++) {
            bool found_i = false;
            for (size_t j = 0 ; j < other.walks.size() ; j++) {
                if (walks[i] == other.walks[j]) {
                    found_i = true;
                    if (found_j[j]) {
                        std::cerr << "Duplicate walk? " << other.walks[j] << std::endl;
                    }
                    found_j[j] = true;
                    allele_conversion[i] = j;
                }
            }
            if (!found_i) {
                std::cerr << "Walk not found in other snarl genotype: " << walks[i] << std::endl;
                return false;
            }
        }
        for (size_t j= 0 ; j < other.walks.size() ; j++) {
            if (!found_j[j]) {
                std::cerr << "Walk  not found in this snarl genotype: " << other.walks[j] << std::endl;
                return false;
            }
        }
    }

    // Check the sequences
    if (allele_conversion.size() == 0) {
        // If we didn't already get the conversion, try again with the sequences
        std::vector<bool> found_j(other.sequences.size() , false);
        for (size_t i = 0 ; i < sequences.size() ; i++) {
            bool found_i = false;
            for (size_t j = 0 ; j < other.sequences.size() ; j++) {
                if (sequences[i] == other.sequences[j]) {
                    found_i = true;
                    if (found_j[j]) {
                        std::cerr << "Duplicate sequences? " << other.sequences[j] << std::endl;
                    }
                    found_j[j] = true;
                    allele_conversion[i] = j;
                }
            }
            if (!found_i) {
                std::cerr << "Sequence not found in other snarl genotype: " << sequences[i] << std::endl;
                return false;
            }
        }
        for (size_t j= 0 ; j < other.sequences.size() ; j++) {
            if (!found_j[j]) {
                std::cerr << "Sequence  not found in this snarl genotype: " << other.sequences[j] << std::endl;
                return false;
            }
        }
    } else {
        for (size_t i = 0 ; i < sequences.size() ; i++) {
            if (sequences[i] != other.sequences[allele_conversion[i]]) {
                return false;
            }
        }
    }

    // Check the lengths
    if (allele_conversion.size() == 0) {

        // If we didn't get the conversion already then sort the lengths and check if they are the same
        std::vector<std::string> sorted_lengths = allele_lengths;
        std::sort(sorted_lengths.begin(), sorted_lengths.end());
        std::vector<std::string> other_sorted_lengths = other.allele_lengths;
        std::sort(other_sorted_lengths.begin(), other_sorted_lengths.end());

        for (size_t i = 0 ; i < sorted_lengths.size() ; i++ ){
            if (sorted_lengths[i] != other_sorted_lengths[i]) {
                return false;
            }
        }

    } else {

        for (size_t i = 0 ; i < allele_lengths.size() ; i++) {
            if (allele_lengths[i] != other.allele_lengths[allele_conversion[i]]) {
                return false;
            }
        }
    }
    for (const auto& x : sample_to_allele) {
        if (other.sample_to_allele.count(x.first) == 0) {
            std::cerr << "Other missing sample " << x.first << std::endl;
        }
        if (x.second == std::numeric_limits<size_t>::max()) {
            if (other.sample_to_allele.at(x.first) != std::numeric_limits<size_t>::max()) {
                return false;
            }
        } else {
            if (other.sample_to_allele.at(x.first) == std::numeric_limits<size_t>::max()) {
                return false;
            } else {
                if (allele_conversion.count(other.sample_to_allele.at(x.first)) == 0) {
                    allele_conversion[other.sample_to_allele.at(x.first)] = x.second;
                } else if (x.second != allele_conversion[other.sample_to_allele.at(x.first)]) {
                    return false;
                }
            }
        }
    } 

    
    return true;
}



std::vector<snarl_genotype_values_t> load_genotype_file(const std::string& infile) {
    std::ifstream instream;
    instream.open(infile);

    std::string line;

    // First header line should be the version
    std::getline(instream, line);
    if (line != "#SNARL_DATA_v1.0" ) {
        std::cerr << "Warning: genotype file is version " << line << std::endl;
    }
    // allele size limit
    std::getline(instream, line);
    // snarl child limit
    std::getline(instream, line);
    // walk steps limit 
    std::getline(instream, line);
    // Start of refs
    std::getline(instream, line);

    std::vector<std::string> references;
    while (line != "#SNARLS") {
        std::getline(instream, line);
        std::string ref (line.begin()+1, line.end());
        references.emplace_back(ref); 
    }


    //Next line is the header with the samples
    std::getline(instream, line);
    std::vector<std::string> sample_names;

    std::stringstream linestream(line);
    std::string sample_name;
    // First go through the first 9 things and ignore them
    for (size_t i = 0 ; i < 9 ; i++) {
        std::getline(linestream, sample_name, '\t');
    }
    while (std::getline(linestream, sample_name, '\t')) {
        sample_names.emplace_back(sample_name);
    }

    std::vector<snarl_genotype_values_t> snarls;

    // Now go through the actual snarls
    while (std::getline(instream, line)) {
        snarl_genotype_values_t snarl;

        std::stringstream snarlstream(line);

        std::string part;
        //Snarl start node traversal
        std::getline(snarlstream, part, '\t');
        snarl.start_node = part;
        
        //Snarl start node traversal
        std::getline(snarlstream, part, '\t');
        snarl.end_node = part;
        
        // Index of reference
        std::getline(snarlstream, part, '\t');
        snarl.reference = part == "18446744073709551615" ? "." : references.at(std::stoull(part));
        
        // Start offset along reference
        std::getline(snarlstream, part, '\t');
        snarl.start_offset = std::stoull(part);
        
        // End offset along reference
        std::getline(snarlstream, part, '\t');
        snarl.end_offset = std::stoull(part);
        
        // Snarl depth
        std::getline(snarlstream, part, '\t');
        snarl.depth = std::stoull(part);

        // Lengths of paths
        std::getline(snarlstream, part, '\t');
        std::vector<std::string> lengths;
        if (part != ".") {
            std::string length_str;
            std::stringstream lengthstream(part);
            while (std::getline(lengthstream, length_str, ',')) {
                lengths.emplace_back(length_str);
            }
        }
        snarl.allele_lengths = std::move(lengths);

        // Walks
        std::vector<std::string> walks;
        std::getline(snarlstream, part, '\t');
        if (part != ".") {
            std::string path_str;
            std::stringstream pathstream(part);
            while (std::getline(pathstream, path_str, ',')) {
                walks.emplace_back(std::move(path_str));
            }
        }
        snarl.walks = std::move(walks);


        // Sequences
        std::getline(snarlstream, part, '\t');
        if (part != ".") {
            std::vector<std::string> sequences;
            std::stringstream seqstream(part);
            std::string seq;
            bool last_empty_sequence = part.size() > 0 && part.at(part.size()-1) == ',';
            while (std::getline(seqstream, seq, ',')) {
                sequences.emplace_back(seq);
            }
        
            if (last_empty_sequence) {
                sequences.emplace_back("");
            }
        
            snarl.sequences = std::move(sequences);
        }

        std::unordered_map<std::string, size_t> sample_to_allele;
        size_t sample_num = 0;
        while (std::getline(snarlstream, part, '\t')) {
            sample_to_allele[sample_names.at(sample_num)] = (part == "." ? std::numeric_limits<size_t>::max() : std::stoull(part));
            sample_num ++;
        };

        snarl.sample_to_allele = std::move(sample_to_allele);
    

        snarls.emplace_back(std::move(snarl));
    }

    instream.close();
    return snarls;
}
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

    // The group paths may be na
    bool no_group = group_paths[0] == "NA";
    if (no_group && other.group_paths[0] != "NA") {
            return false;
    }

    // Put the corresponding path length and group_paths together in a set
    std::unordered_set<std::string> path_lengths_and_groups;
    for (size_t i = 0 ; i < path_lengths.size() ; i++) {
        path_lengths_and_groups.emplace(path_lengths[i] + " " + (no_group ? "" : group_paths[i]));
    }
    std::unordered_set<std::string> other_path_lengths_and_groups;
    for (size_t i = 0 ; i < other.path_lengths.size() ; i++) {
        other_path_lengths_and_groups.emplace(other.path_lengths[i] + " " + (no_group ? "" : other.group_paths[i]));
    }
    // the sets must be the same
    if (path_lengths_and_groups != other_path_lengths_and_groups) {
        return false;
    }

    return true;
}


//CHR\tSTART_POS\tEND_POS\tSTART_NODE\tEND_NODE\tALLELE_LENGTHS\tP_FISHER\tP_CHI2\tALLELE_COUNT_PER_PHENO\tDEPTH
binary_table_values_t load_binary_snarl_line(std::string line) {


    binary_table_values_t vals;
    // Components:
    //std::string chr;
    //size_t start_pos;
    //size_t end_pos;
    //std::pair<size_t, size_t> snarl_ids;
    //std::vector<std::string> path_lengths;
    //std::string p_fishers;
    //std::string p_chi2;
    //std::vector<std::string> group_paths;
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
    std::getline(linestream, part, '\t');
    vals.snarl_ids.first = part;

    // end id
    std::getline(linestream, part, '\t');
    vals.snarl_ids.second = part;

    // path lengths, comma separated by allele, /-separated for min/max in one allele
    std::getline(linestream, part, '\t');

    std::stringstream pathlength_stream(part);
    std::string lengths;
    while (std::getline(pathlength_stream, lengths, ',')){
        if (lengths != "NA") {
            vals.path_lengths.emplace_back(lengths);
        }
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
