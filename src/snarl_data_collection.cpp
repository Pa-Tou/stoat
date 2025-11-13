#include "snarl_data_collection.hpp"
#include <fstream>

//#define DEBUG_SNARL_DATA_COLLECTION

using namespace std;
namespace stoat {


// Constructor
// This goes through all the snarls and fills in the data
SnarlDataCollection::SnarlDataCollection(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                    size_t allele_size_limit, size_t snarl_child_limit, size_t walk_cycle_limit, size_t walk_steps_limit,
                    bool partition_requested,
                    const std::function<std::vector<std::set<sample_hap_t>>(const handlegraph::PathPositionHandleGraph& graph, 
                                                                            const bdsg::SnarlDistanceIndex& distance_index, 
                                                                            const net_handle_t& snarl)>& find_sample_partitions,
                    bool sequence_requested) :
                    allele_size_limit(allele_size_limit),
                    snarl_child_limit(snarl_child_limit),
                    walk_cycle_limit(walk_cycle_limit),
                    walk_steps_limit(walk_steps_limit) {
    
    check_distances = distance_index == nullptr ? false : distance_index->has_distances()

    // Get a list of all chains in root
    std::vector<handlegraph::net_handle_t> chains;
    chains.reserve(50);
    handlegraph::net_handle_t root = distance_index->get_root();
    distance_index->for_each_child(root, [&] (handlegraph::net_handle_t chain) {
        chains.emplace_back(chain);
        return true;
    });

    // Count the number of chains that we added and the number of chains that we actually process,
    // as a way of debugging the parallelization. Not in ifdef because it needs to go in the omp parallel shared
    size_t chains_added = chains.size();
    size_t chains_processed = 0;

    bool keep_going = !chains.empty();

    // Keep track of which references we've seen and their index in reference_names
    std::unordered_map<std::string, size_t> reference_name_to_index;
    
    // Go through the contents of chains in parallel
    // Everything touching chains needs to be in an omp critical block so they don't collide. 
    #pragma omp parallel shared(chains, keep_going, chains_added, chains_processed, all_snarl_data, snarl_to_walks, snarl_to_partitions, snarl_to_sequences, reference_names, sample_haplotypes)
    {
        // The actual while loop is run on a single thread
        #pragma omp single
        {
            while (keep_going) {
                handlegraph::net_handle_t chain;
                #pragma omp critical(chain_loop)
                {
                chain = chains.back();
                chains.pop_back();
                #ifdef DEBUG_SNARL_DATA_COLLECTION
                chains_processed++;
                #endif
                }
    
                // Everything in here is parallelized
                #pragma omp task
                {
                    
                    distance_index->for_each_child(chain, [&] (handlegraph::net_handle_t snarl) {
                    
                        //TODO: Actually use is_eligible
                        //TODO: For now it's fine to check is_eligible here because it's only checking size and we don't want to look at small chains anyway
                        if (distance_index->is_snarl(snarl) && snarl_is_eligible(snarl)) {
    
                            #ifdef DEBUG_SNARL_DATA_COLLECTION
                            cerr << "Test snarl " << distance_index->net_handle_as_string(snarl) << endl;
                            #endif

                            // Make the snarl_data_internal_t to fill in. Since it's multithreaded its better to move() it instead of adding it here
                            snarl_data_internal_t snarl_data;

                            // Get the start and end nodes
                            // Do it through the graph because it's a pain to get the orientation from the distance index
                            handlegraph::handle_t start_in = distance_index.get_handle(distance_index.get_bound(snarl, false, true), &graph);
                            handlegraph::handle_t end_in = distance_index.get_handle(distance_index.get_bound(snarl, true, true), &graph);

                            snarl_data.start_node = stoat::Node_traversal_t(graph.get_id(start_in),
                                                                            graph.get_is_reverse(start_in));
                            snarl_data.end_node = stoat::Node_traversal_t(graph.get_id(end_in),
                                                                            graph.get_is_reverse(end_in));

                            // Add the depth of the snarl
                            snarl_data.depth = distance_index.get_depth(snarl);
    
    
                            // Get the offsets of the start and end nodes along the reference
                            std::vector<stoat::path_range_t> ranges = stoat::get_coordinates_of_snarl(graph, *distance_index, snarl, true, reference_sample, false);
                            if (ranges.size() != 0) {
                                // Check if we have already seen the reference path and if not add it
                                size_t ref_index;

                                auto reference_range = get_name_and_offsets_of_snarl_path_range(graph, ranges.front());
                                snarl_data.start_position = std::get<1>(reference_range);
                                snarl_data.end_position = std::get<2>(reference_range);
                                #pragma omp critical(chain_loop)
                                {
                                    if (reference_name_to_index.count(std::get<0>(reference_range)) == 0) {
                                        ref_index = reference_names.size();
                                        reference_name_to_index[std::get<0>(reference_range)] = ref_index;
                                        reference_names.emplace_back(std::move(std::get<0>(reference_range)));
                                    } else {
                                        ref_index = reference_name_to_index[std::get<0>(reference_range)];
                                    }
                                }

                            } else {
                                snarl_data.start_position = 0;
                                snarl_data.end_position = 0;
                                snarl_data.ref_index = std::numeric_limits<size_t>::max();
                            }

                            // There will always be the set of walks through the snarl
                            // The partitions and sequences may optionally be filled in.
                            // All three of these vectors must have the same number of entries because they correspond to the same walks. 
                            // If partitions are requested, then build the other two around the partitions that are found.
                            // Otherwise, find the walks and get the sequences from those
                            std::vector<stoat::Path_traversal_t> snarl_walks; 
                            std::vector<std::string> walk_lengths;
                            std::vector<std::set<sample_hap_t>> sample_partitions;
                            std::vector<std::string> snarl_sequences;

                            if (partition_requested) {
                                // If we want partitions, then find them and then only enumerate walks and sequences in the partition

                                sample_partitions = find_sample_partitions(graph, distance_index, snarl);

                                std::tie(snarl_walks, walk_lengths, snarl_sequences) = get_walks_and_sequences_from_partitions(graph, distance_index, snarl, 
                                                                                                                               sample_partitions, sequence_requested);

                            } else {
                                // If we don't want partitions, then enumerate all walks then find the sequences
                                std::tie(snarl_walks, walk_lengths) = get_walks_through_snarl(graph, distance_index, snarl);
                            }

                            // Get the sequences for each walk
    
                            #pragma omp critical(chain_loop)
                            {
   
                                // Add the child chains to the stack
                                distance_index->for_each_child(snarl, [&] (handlegraph::net_handle_t child) {
    
                                    chains.emplace_back(child);
                                    #ifdef DEBUG_SNARL_DATA_COLLECTION
                                    chains_added++;
                                    #endif
                                    return true;
                                });
                            }
    
                        }
                        return true;
                    }); //end for each child
                }// end omp task
                #pragma omp critical(chain_loop)
                {
                    keep_going = !chains.empty();
                }
    
                if (!keep_going) {
                    // Wait for tasks to complete
                    #pragma omp taskwait
    
                    #pragma omp critical(chain_loop)
                    {
                        // Check again if we're done or not
                        keep_going = !chains.empty();
                    }
                }
            }// end while loop
        }// End omp single
    }//end omp shared
    
    #ifdef DEBUG_SNARL_DATA_COLLECTION
    cerr << "Added " << chains_added << " chains and processed " << chains_processed << endl;
    assert(chains_added == chains_processed);
    #endif
}

std::tuple<std::vector<stoat::Path_traversal_t>, std::vector<std::string>> SnarlDataCollection::get_walks_through_snarl(
        const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index, const net_handle_t& snarl) const {


    // Path exploration
    std::vector<std::vector<handlegraph::net_handle_t>> paths = {
        {distance_index.get_bound(snarl, false, true)}
    };
    
    std::vector<std::vector<handlegraph::net_handle_t>> finished_paths;
    
    // How many steps have we taken trying to enumerate paths? Includes all all paths
    size_t steps_taken = 0;

    bool break_snarl = false;
    
    // For each incomplete path in paths, walk out from the end and add a copy of the path plus each next step to paths
    // Do this until the path reaches the end
    while (!paths.empty()) {
        std::vector<handlegraph::net_handle_t> path = std::move(paths.back());
        paths.pop_back();
    
        std::unordered_map<handlegraph::net_handle_t, size_t> dict_path_occ;
        bool cycle = false;
    
        for (const auto& net : path) {
            if (++dict_path_occ[net] > walk_cycle_limit + 1) {
                cycle = true;
                break;
            }
        }
    
        // TODO: Add out_fail
        if (steps_taken++ > walk_steps_limit) {
            #pragma omp critical(out_fail)
            out_fail << snarl_id_str << "\titeration_calculation_out = " << children << " children\n";
            break_snarl = true;
            break;
        }

        // Follow edges from the last element in path
        if (!path.empty()) {
            distance_index.follow_net_edges(path.back(), &graph, false, [&](const handlegraph::net_handle_t& next_child) {
                // If this is the bound of the snarl then we're done
                if (distance_index.is_sentinel(next_child)) {

                    size_t next_child_node_id = distance_index.node_id(distance_index.get_node_from_sentinel(next_child));
                    size_t first_element_path_node_id = distance_index.node_id(distance_index.get_node_from_sentinel(path[0]));

                    // Only keep the walk if it entered and exited the snarl at opposite sides
                    if (next_child_node_id != first_element_path_node_id) {
                        finished_paths.emplace_back(path);
                        finished_paths.back().push_back(next_child);
                    }
            
                } else {
                    //TODO: Look for the cycle sooner
            
                    if (cycle) { // Case where we find a loop
                        return false;
                    }
                    paths.emplace_back(path);
                    paths.back().push_back(next_child);
                }
                return true;
            };
        }
    }
    
    if (break_snarl) {
        finished_paths.clear();
    }
    
    // Transform the paths into the format we want
    return fill_pretty_paths(distance_index, graph, finished_paths);
    
}

std::tuple<std::vector<stoat::Path_traversal_t>, std::vector<std::string>, std::vector<std::string>> SnarlDataCollection::get_walks_and_sequences_from_partitions(
        const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index, const net_handle_t& snarl,
        const std::vector<std::set<sample_hap_t>>& sample_partitions, bool sequence_requested) const {
    std::vector<stoat::Path_traversal_t> snarl_walks;
    std::vector<std::string> walk_lengths;
    std::vector<std::string> snarl_sequences;

    // Get the walks and sequences by tracing the haplotype paths through the snarl


    ///////////////////////////////////////////// Get one step_handle_t's for the each partition

    handlegraph::handle_t start_in = distance_index.get_handle(distance_index.get_bound(snarl, false, true), &graph);
    handlegraph::handle_t end_in = distance_index.get_handle(distance_index.get_bound(snarl, true, true), &graph);

    std::vector<handlegraph::PathSense> senses = {handlegraph::PathSense::GENERIC,
                                                  handlegraph::PathSense::REFERENCE,
                                                  handlegraph::PathSense::HAPLOTYPE};

    std::vector<handlegraph::step_handle_t> first_steps;
    for (const std::set<sample_hap_t>& partition : sample_partitions) {
        bool found_step = false;
        for (const auto& sense : senses) {
            graph.for_each_step_of_sense(start_in, sense, [&](const handlegraph::step_handle_t& step) {
                if (partition.count(stoat::get_sample_and_haplotype(graph, graph.get_path_handle_of_step(step))) ) {
                    first_steps.emplace_back(step);
                    found_step = true;
                    // Return false to stop looping through steps
                    return false;
                } else{
                    return true;
                }
            }
            if (found_step) {
                // Break out of the inner loop and continue to the next partition
                break;
            }
        }
    }
    #ifdef DEBUG_SNARL_DATA_COLLECTION
    assert(first_steps.size() == sample_partitions.size());
    #endif

    //////////////////////////////////// Follow each of the paths through the snarl and create a walk and sequence

    for (const handlegraph::step_handle_t& first_step : first_steps) {
        //Accumulate the min and max lengths of this walk, based on the lengths of snarls
        size_t min_length = 0;
        size_t max_length = 0;

        // Get the step of this path on both the start and end nodes.
        std::vector<handlegraph::step_handle_t> boundary_steps;
        for (const auto& sense : senses) {
            graph.for_each_step_of_sense(start_in, sense, [&](const handlegraph::step_handle_t& this_step) {
                if (graph.get_path_handle_of_step(first_step) == graph.get_path_handle_of_step(this_step)) {
                    boundary_steps.empalce_back(this_step);
                 }
                return true;
            }); 
            graph.for_each_step_of_sense(end_in, sense, [&](const handlegraph::step_handle_t& this_step) {
                if (graph.get_path_handle_of_step(first_step) == graph.get_path_handle_of_step(this_step)) {
                    boundary_steps.empalce_back(this_step);
                 }
                return true;
            });
        }

        // Now sort the list of steps on the two boundary nodes
        boundary_steps.sort([&](const handlegraph::step_handle_t& a, const handlegraph::step_handle_t& b) {
            return graph.get_position_of_step(a) < graph.get_position_of_step(b);
        });

        // Now walk through the path until we hit all the steps in boundary_steps
        handlegraph::step_handle_t step = graph.get_next_step(start_step);


    }// end for each first step (per partition)

    return std::make_tuple(snarl_walks, walk_lengths, snarl_sequences);
}



bool SnarlDataCollection::snarl_is_eligible(const handlegraph::net_handle_t& snarl) const {
    bool pass = true;

    // If we have distances in the index, make sure that the snarl's maximum lenght is small enough
    if (check_distances) {
        pass &= allele_size_limit <= distance_index.maximum_length(snarl);
    }
    //TODO: Once the libbdsg branch is merged we can use this instead of going through all the children to count them
    //pass &= snarl_child_limit <= distance_index.get_snarl_child_count(snarl);

    // Count children
    size_t children = 0;
    distance_index.for_each_child(snarl, [&](const handlegraph::net_handle_t& snarl) {
        children++;
        return true;
    });
    pass &= snarl_child_limit <= children;

}




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

Each snarl will also contain a vector of Path_traversal_t's, (optionally) a vector of sets of haplotypes that take each Path_traversal_t, 
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
The 8th item is all of the walks, comma separated.
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

        // Next, always include the walks as a single comma-separated string
        stoat::vectorPathToString(snarl_to_walks[snarl_data.start_node]);

        // Next add the partitions, if there are any
        // Each partition is saved as the number of items in it then the items
        if (!snarl_to_partitions.empty()) {
            #ifdef DEBUG_SNARL_DATA_COLLECTION
            assert(snarl_to_walks[snarl_data.start_node].size() == snarl_to_partitions[snarl_data.start_node].size());
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
            assert(snarl_to_walks[snarl_data.start_node].size() == snarl_to_sequences[snarl_data.start_node].size());
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
        snarl_to_walks[all_snarl_data.back().start_node] = stoat::stringToVectorPath(part);
        size_t walk_count = snarl_to_walks[all_snarl_data.back().start_node].size();

        // This may be the end, or there may be partitions and/or sequences
        bool keep_going = std::getline(linestream, part, '\t');
        if (!keep_going) {
            continue;
        }

        if (!(part == "," || part[0] == "A" || part[0] == "C" || part[0] == "G" || part[0] == "T")) {
            // If the next thing isn't a sequence then we need to get the partitions next
            std::vector<std::set<size_t>> partitions;

            for (size_t i = 0 ; i < walk_count ; i++) {

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
            assert(partitions.size() == walk_count);
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
            assert(sequences.size() == walk_count);
            #endif

            snarl_to_sequences[all_snarl_data.back().start_node] = std::move(sequences);
        }
    }


    instream.close();
    return;
}
}

