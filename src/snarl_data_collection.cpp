#include "snarl_data_collection.hpp"
#include <fstream>

//#define DEBUG_SNARL_DATA_COLLECTION

using namespace std;
namespace stoat {

/////////////////////////////////////////// Writing/reading the snarl data
/*
This needs to hold snarl_data_internal_t's which contain:

- Node_traversal_t start_node;
- Node_traversal_t end_node;
- size_t reference_index (which points to a string representing the reference name);
- size_t start_position;
- size_t end_position;
- size_t depth;
- string variant_type;

Each snarl will also contain a vector of Path_traversal_t's, (optionally) a vector of haplotypes that take each Path_traversal_t, 
and (optionally) the sequences of each Path_traversal_t


Start with the header (SnarlDataCollection::file_header)

Next the limits for testing snarls, in case we skipped small snarls for example:
allele_size_limit:N
snarl_child_limit:N

Next a list of reference path names. These get stored in reference_names and the order will be the index for reference_index
This section starts with #REFS

Next all sample/haplotypes. These get stored in sample_haplotypes and the order will be the index for sample/haplotypes used when storing partitions
This section starts with #SAMPLES


Each snarl_partition_t is then stored, one per line, as a tab-separated vector of integers/strings.
the first 7 items are the contents of the snarl_data_internal_t.
The 8th item is all of the paths, comma separated.
The haplotype partitions are stored next (if present). Each partition starts with the number of sample/haplotypes for that partition, followed with
that number of entries for the sample_hap_t. There should be N partitions.
The next N items (if present) are the sequences of the Path_traversal_t's. These can be distinguished from the haplotype partitions because they can only contain A/C/G/T/N, or may be empty
This section starts with #SNARLS

*/
void SnarlDataCollection::save(const std::string& filename) {
    ofstream outstream;
    outstream.open(filename);
    if (!outstream.good()) {
        throw std::runtime_error("stoat: could not open file " + filename + " for writing");
    }
    // Write the header
    outstream << file_header << endl;

    outstream << "allele_size_limit:" << allele_size_limit << endl;
    outstream << "snarl_child_limit:" << snarl_child_limit << endl;

    // Next will be a list of reference path names.
    outstream << "#REFS" << endl;
    for (const auto& ref : reference_names) {
        outstream << ref << endl;
    }
    
    // Next will be a list of sample/haplotypes
    outstream << "#SAMPLES" << endl;
    for (const auto& samp : sample_haplotypes) {
        outstream << sample.sample << "\t" << sample.haplotype << endl;
    }

    //Finally the snarls
    outstream << "#SNARLS" << endl;
    for (const snarl_data_internal_t& snarl_data : all_snarl_data) {

        // Start with just the contents of the snarl_data_internal_t
        outstream << snarl_data.start_node.to_string() << "\t"
                  << snarl_data.end_node.to_string() << "\t"
                  << snarl_data.reference_index << "\t"
                  << snarl_data.start_position << "\t"
                  << snarl_data.end_position << "\t"
                  << snarl_data.depth << "\t"
                  << snarl_data.variant_type << "\t";

        // Next, always include the paths as a single comma-separated string
        stoat::vectorPathToString(snarl_to_paths[snarl_data.start_node]);

        // Next add the partitions, if there are any
        // Each partition is saved as the number of items in it then the items
        if (!snarl_to_partitions.empty()) {
            #ifdef DEBUG_SNARL_DATA_COLLECTION
            assert(snarl_to_paths[snarl_data.start_node].size() == snarl_to_partitions[snarl_data.start_node].size());
            #endif
            for (const std::set<size_t>& partition : snarl_to_partitions[snarl_data.start_node]){
                cerr << partition.size() << "\t";
                for (const size_t& x : partition) {
                    cerr << x << "\t";
                }
            }
        }

        // Add the sequences, if any, as comma separated strings
        if (!snarl_to_sequences.empty()) {
            #ifdef DEBUG_SNARL_DATA_COLLECTION
            assert(snarl_to_paths[snarl_data.start_node].size() == snarl_to_sequences[snarl_data.start_node].size());
            #endif
            for (const std::string& seq : snarl_to_sequences[snarl_data.start_node]){
                cerr << seq << ",";
            }
        }
        
        outstream << endl;
    }


    outstream.close();
    return;
}

void SnarlDataCollection::load(const std::string& filename, const handlegraph::PathPositionHandleGraph& graph) {
    ifstream instream;
    instream.open(filename);

    string line;

    // Read the first line, which must match the header
    std::getline(instream, line);
    if (line != file_header) {
        throw std::runtime_error("stoat: Snarl data file " +filename+ " contains the wrong header: " + line);
    }

    // Next should be the allele size limit
    std::getline(instream, line);
    std::stringstream linestream(line);
    string limit_str;
    std::getline(linestream, limit_str, ':');
    #ifdef DEBUG_SNARL_DATA_COLLECTION
    assert(limit_str == "allele_size_limit");
    #endif
    std::getline(linestream, limit_str, ':');
    if (allele_size_limit < std::stoull(limit_str)) {
        cerr << "warning [stoat]: The allele_size_limit of the saved snarls file is larger than the given allele_size_limit. Some snarls may be missed" << endl;;
    }

    // And the snarl child limit
    std::getline(instream, line);
    std::stringstream newlinestream(line);
    std::getline(newlinestream, limit_str, ':');
    #ifdef DEBUG_SNARL_DATA_COLLECTION
    assert(limit_str == "snarl_child_limit");
    #endif
    std::getline(newlinestream, limit_str, ':');
    if (snarl_child_limit < std::stoull(limit_str)) {
        cerr << "warning [stoat]: The snarl_child_limit of the saved snarls file is larger than the given snarl_child_limit. Some snarls may be missed" << endl;;
    }

    // Next are the references
    std::getline(instream, line);
    if (line != "#REFS") {
        throw std::runtime_error("stoat: Snarl partitions file " +filename+ " is not formatted correctly");
    }
    // Get the reference path names from the file
    std::getline(instream, line);
    while (line != "#SAMPLES") {
        reference_names.emplace_back(std::move(line));
        std::getline(instream, line);
    }

    // Get the sample/haplotype indexes from the file
    std::getline(instream, line);
    while (line != "#SNARLS") {
        std::stringstream linestream(line);
        string sample_name;
        std::getline(linestream, sample_name, '\t');
        string hap_str;
        std::getline(linestream, hap_str, '\t');
        sample_haplotypes.emplace_back(sample_name, std::stoull(hap_str));
        std::getline(instream, line);
    }


    // Get the snarls
    while (std::getline(instream,line)) {

        all_snarl_data.emplace_back();

        std::stringstream linestream(line);
        string part;
        //Snarl start node traversal
        std::getline(linestream, part, '\t');
        all_snarl_data.back().start_node = stoat::Node_traversal_t(part);

        //Snarl start node traversal
        std::getline(linestream, part, '\t');
        all_snarl_data.back().end_node = stoat::Node_traversal_t(part);

        // Index of reference
        std::getline(linestream, part, '\t');
        all_snarl_data.back().reference_index = std::stoll(part);

        // Start offset along reference
        std::getline(linestream, part, '\t');
        all_snarl_data.back().start_position = std::stoll(part);

        // End offset along reference
        std::getline(linestream, part, '\t');
        all_snarl_data.back().end_position = std::stoll(part);

        // Snarl depth
        std::getline(linestream, part, '\t');
        all_snarl_data.back().depth = std::stoll(part);

        // variant type
        std::getline(linestream, part, '\t');
        all_snarl_data.back().variant_type = part;

        // Next are the paths, all as one comma-separated string
        std::getline(linestream, part, '\t');
        snarl_to_paths[all_snarl_data.back().start_node] = stoat::stringToVectorPath(part);
        size_t path_count = snarl_to_paths[all_snarl_data.back().start_node].size();

        // This may be the end, or there may be partitions and/or sequences
        bool keep_going = std::getline(linestream, part, '\t');
        if (!keep_going) {
            continue;
        }

        if (!(part == "," || part[0] == "A" || part[0] == "C" || part[0] == "G" || part[0] == "T")) {
            // If the next thing isn't a sequence then we need to get the partitions next
            std::vector<std::set<size_t>> partitions;

            for (size_t i = 0 ; i < path_count ; i++) {

                partitions.emplace_back();
                
                // Since we had to check for the end of the line, part is currently the count of samples
                // and we update it at the end of the loop 
                size_t sample_count = std::stoull(part);
                for (size_t i = 0 ; i < sample_count ; i++) {
                    std::getline(linestream, part, '\t');
                    partitions.back().emplace(samples[std::stoull(part)]);
                }

                // Check if this is the end of the line
                keep_going = std::getline(linestream, part, '\t');
            }

            #ifdef DEBUG_SNARL_DATA_COLLECTION
            assert(partitions.size() == path_count);
            #endif

            snarl_to_partitions[all_snarl_data.back().start_node] = std::move(partitions);
        }

        // The last thing we got is now either the end of the line or the first sequence
        if (keep_going) {
            std::vector<std::string> sequences;
            std::stringstream seqstream(part);
            std::string seq;
            while (std::getline(seqstream, seq, ',')) {
                sequences.emplace_back(seq);
            }
            #ifdef DEBUG_SNARL_DATA_COLLECTION
            assert(sequences.size() == path_count);
            #endif

            snarl_to_sequences[all_snarl_data.back().start_node] = std::move(sequences);
        }
    }


    instream.close();
    return;
}
}
