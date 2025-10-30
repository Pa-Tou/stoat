#include "partitioner.hpp"
#include "log.hpp"
#include <fstream>

using namespace std;
using namespace stoat;
namespace stoat_graph {

std::vector<std::set<sample_hap_t>> SnarlTraverserAndPathPartitioner::partition_samples_in_snarl(const handlegraph::PathPositionHandleGraph& graph, 
                                                                               const handlegraph::net_handle_t& snarl) const {
    stoat::LOG_TRACE("Get sample partitions of " + distance_index->net_handle_as_string(snarl) + " by its paths" );

    //Get the partition of paths. If the snarl is regular, then only check the edges leaving the start node
    std::vector<std::set<stoat::sample_hap_t>> sample_sets = get_walk_sets(graph, snarl, distance_index->is_regular_snarl(snarl, true, &graph));

    stoat::LOG_TRACE((string) "Found sets of paths using " + ( distance_index->is_regular_snarl(snarl, true, &graph) ? "edges from the start node" : "walk sets"));
    for (const std::set<stoat::sample_hap_t>& sample_set : sample_sets) {
        stoat::LOG_TRACE("SET ");
        for (const stoat::sample_hap_t& sample : sample_set) {
            stoat::LOG_TRACE("\t" + sample.sample );
        }
    }

    return sample_sets;
}

// This is supposed to partition the paths in the snarl by the walks they take through the netgraph.
// Instead of explicitly enumerating the paths, it actually finds the sets of edges that each path takes.
// But since paths may loop, it actually finds the order and number of outgoing edges from each node.
// I think this is equivalent to partitioning by the actual sets of unique walks.
std::vector<std::set<stoat::sample_hap_t>> SnarlTraverserAndPathPartitioner::get_walk_sets(const handlegraph::PathPositionHandleGraph& graph, 
                                                                   const handlegraph::net_handle_t& snarl,
                                                                   bool only_bound) const {

    stoat::LOG_TRACE((string) "Get walk sets of " + distance_index->net_handle_as_string(snarl));

    // Make a vector of the paths 
    std::vector<stoat::sample_hap_t> all_samples(all_sample_haplotypes.begin(), all_sample_haplotypes.end());

    // Map each sample(plus haplotype) to its index in all_samples
    std::map<stoat::sample_hap_t, std::size_t> sample_to_index;
    for (size_t i = 0 ; i < all_samples.size() ; i++) {
        sample_to_index[all_samples[i]] = i;
    }

    // For each path (along all_samples), the index of the set it is currently in
    // Everything starts out in the same set, 0, which represents not being in the snarl
    std::vector<size_t> old_sets (all_samples.size(), 0);
    size_t old_set_count = 1;

    // TODO: I think this is unnecessary
    //// First, check for paths that don't go through this snarl
    //std::vector<handlegraph::PathSense> senses = {handlegraph::PathSense::GENERIC,
    //                                              handlegraph::PathSense::REFERENCE,
    //                                              handlegraph::PathSense::HAPLOTYPE};
    //
    ////TODO: This could also use steps_of_handle() but it doesn't seem to work 
    //for (const auto& sense : senses) {
    //    graph.for_each_step_of_sense(distance_index->get_handle(distance_index->get_node_from_sentinel(distance_index->get_bound(snarl, false, true)), &graph),
    //        sense, [&](const handlegraph::step_handle_t& step) {
    //        old_sets[sample_to_index[stoat::get_sample_and_haplotype(graph, graph.get_path_handle_of_step(step))]] = 1;
    //        old_set_count = 2;
    //    });
    //}


    // A struct representing a path's edge going out from child
    // Represents the next node and orientation and the offset along the path of the edge
    // Use (0,0,0,false) as a default value to indicate that the path didn't traverse this child
    // This becomes a linked list when a path traverses the child multiple times (this doesn't happen often) 
    struct path_edge_t {
        size_t offset;           // The offset along the path
        size_t additional_edge;  // If this path traverses the child multiple times, index into additional_steps for others. max() for end of linked list 
        handlegraph::nid_t id;   // The id of the next node
        bool rev;                // Is the edge going to the right side of the next node?

        path_edge_t() : offset(0), additional_edge(std::numeric_limits<size_t>::max()), id(0), rev(false){};
        path_edge_t(size_t offset, size_t edge, handlegraph::nid_t id, bool rev) : 
            offset(offset), additional_edge(edge), id(id), rev(rev){};
    };


    // Helper function to check the edges leaving child
    // This gets called in for_each_child, and also from the start node going in
    auto check_outgoing_edges = [&] (const handlegraph::net_handle_t& child, bool go_left) {

        stoat::LOG_TRACE( (string) "At snarl child " + distance_index->net_handle_as_string(child) + " going " + (go_left ? "left" : "right"));
        
        std::vector<path_edge_t> next_steps (all_samples.size());
        // For paths that go through multiple steps 
        std::vector<path_edge_t> additional_steps;
        
        // Get a handle to the end of the child
        handlegraph::handle_t handle;
        if ( distance_index->is_trivial_chain(child)) {
            handle =  (go_left ? graph.flip(distance_index->get_handle(child, &graph)) 
                               : distance_index->get_handle(child, &graph));
        } else if (distance_index->is_sentinel(child)) {
            handle = distance_index->get_handle(child, &graph);
        } else {
            handle = distance_index->get_handle(distance_index->get_bound(child, go_left, false), &graph);
        }
        
        std::vector<handlegraph::PathSense> senses = {handlegraph::PathSense::GENERIC,
                                                      handlegraph::PathSense::REFERENCE,
                                                      handlegraph::PathSense::HAPLOTYPE};

        stoat::LOG_TRACE( (std::stringstream) "" << "\tgraph handle " << graph.get_id(handle) << " going " << (graph.get_is_reverse(handle) ? "left" : "right"));

        for (const auto& sense : senses) {
            graph.for_each_step_of_sense(handle, sense, [&](const handlegraph::step_handle_t& step) {
                // For each step on the node handle, keep track of which paths take different steps

                stoat::LOG_TRACE( "\ton path " + graph.get_path_name(graph.get_path_handle_of_step(step)));

                sample_hap_t step_sample_haplotype = stoat::get_sample_and_haplotype(graph, graph.get_path_handle_of_step(step));

                // If this is not a sample of interest, skip it
                if (!all_sample_haplotypes.count(step_sample_haplotype)) {
                    return true;
                }
        
                //Do we go forwards in the path? We need to check the direction of the handle in the path
                bool go_forwards = graph.get_is_reverse(handle) == graph.get_is_reverse(graph.get_handle_of_step(step));
        
                //In the case where a path doesn't go all the way through the snarl, stop when the path stops
                if ((go_forwards && !graph.has_next_step(step)) || (!go_forwards && !graph.has_previous_step(step))){
                    return true;
                }
        
                //Get the next step and make an edge
                handlegraph::step_handle_t next_step = go_forwards ? graph.get_next_step(step) : graph.get_previous_step(step);
                handlegraph::handle_t next_handle = graph.get_handle_of_step(next_step);
        
                stoat::LOG_TRACE((std::stringstream) "" << "\t\tgoing to " << graph.get_id(next_handle));
        
                path_edge_t edge (graph.get_position_of_step(step), 
                                  std::numeric_limits<size_t>::max(),
                                  graph.get_id(next_handle), 
                                  graph.get_is_reverse(next_handle));
        
                size_t sample_num = sample_to_index[step_sample_haplotype];
        
                if (next_steps[sample_num].id == 0) {
                    // If this path hasn't been seen before
                    next_steps[sample_num] = std::move(edge);
                } else if (next_steps[sample_num].offset > edge.offset) {
                    // If the new edge comes before the edge stored in the vector, replace it
                    edge.additional_edge = additional_steps.size();
                    additional_steps.emplace_back(std::move(next_steps[sample_num]));
                    next_steps[sample_num] = std::move(edge);
                } else {
                    // If the new edge comes after something in additional_steps, walk through the linked list to find its place

                    // The index into the linked list. Need a bool to know if it's an index into next_steps or additional_steps
                    bool old_edge_first = true;
                    size_t old_edge_index = sample_num;
                    size_t old_next_edge = next_steps[old_edge_index].additional_edge;
                    while (old_next_edge != std::numeric_limits<size_t>::max() &&
                           additional_steps[old_next_edge].offset < edge.offset) {
                        // Step through the linked list
                        old_edge_first = false;
                        old_edge_index = old_next_edge;
                        old_next_edge = additional_steps[old_edge_index].additional_edge;
                    }
                    //old_edge_index now points to the item just smaller than edge
                    if (old_edge_first) {
                        next_steps[old_edge_index].additional_edge = additional_steps.size();
                    } else {
                        additional_steps[old_edge_index].additional_edge = additional_steps.size();
                    }
                    edge.additional_edge = old_next_edge;
                    additional_steps.emplace_back(std::move(edge));
                }
        
                return true;
            });
        }
        
        // We now have the edges for all paths going out of the node in one direction
        // Now we want get an "intermediate set" for each path based on the edge(s) it took from this node/direction
        // Equality in this case is that the edges (as node id and orientation) are in the same order along the path
        // This is later used to split the paths into sets for each pair of old_set and intermediate_set
        
        //This maps each edge list for this node/direction to the intermediate set index
        std::map<std::vector<std::pair<handlegraph::nid_t, bool>>, size_t> edge_to_intermediate_set;
        //Everything starts in the same set, representing not going through this node
        vector<size_t> intermediate_sets (old_sets.size(), 0);
        size_t intermediate_set_count = 1;
        for (size_t path_i = 0 ; path_i < next_steps.size() ; path_i++) {
            const path_edge_t& edge = next_steps[path_i];
            if (edge.id != 0) {
                //If this is a real edge
                std::vector<std::pair<handlegraph::nid_t, bool>> edge_list;
                edge_list.emplace_back(edge.id, edge.rev);
                size_t next_i = edge.additional_edge;
                while (next_i != std::numeric_limits<size_t>::max()) {
                    edge_list.emplace_back(additional_steps[next_i].id, additional_steps[next_i].rev);
                    auto& new_edge = additional_steps[next_i]; 
                    next_i = new_edge.additional_edge;
                }
                if (edge_to_intermediate_set.count(edge_list) == 0) {
                    edge_to_intermediate_set[edge_list] = intermediate_set_count;
                    intermediate_set_count++;
                }
                intermediate_sets[path_i] = edge_to_intermediate_set[edge_list];
            }
        }
        stoat::LOG_TRACE( "Intermediate sets: ");
        for (size_t i = 0 ; i < all_samples.size() ; i++) {
            stoat::LOG_TRACE((std::stringstream) "" << "\t" << all_samples[i] << ": " << intermediate_sets[i]);
        } 
        
        // We now have an old set and an intermediate set for each path
        // Assign the path to a new set. Everything gets a new set
        vector<size_t> new_sets (intermediate_sets.size(), std::numeric_limits<size_t>::max());

        // Map pairs of <old_set, intermediate_set> to new set number
        std::map<std::pair<size_t, size_t>, size_t>  old_to_new_set;

        // Reserve 0 for paths that didn't go through this snarl
        old_to_new_set[std::make_pair(0,0)] = 0;
        size_t new_set_count = 1;

        for (size_t path_i = 0 ; path_i < new_sets.size() ; path_i++) {
            std::pair<size_t, size_t> old_set (old_sets[path_i], intermediate_sets[path_i]);
            if (old_to_new_set.count(old_set) == 0) {
                old_to_new_set[old_set] = new_set_count;
                new_set_count++;
            }
            new_sets[path_i] = old_to_new_set[old_set];
        }
        
        old_sets = std::move(new_sets);
        old_set_count = new_set_count;
        
        stoat::LOG_TRACE("New sets: ");
        for (size_t i = 0 ; i < all_samples.size() ; i++) {
            stoat::LOG_TRACE((std::stringstream) "" << "\t" << all_samples[i] << ": " << old_sets[i]);
        } 
    };

    // Now do the work of going through the edges for the start bound and each child in both directions

    check_outgoing_edges(distance_index->get_bound(snarl, false, true), false);

    if (!only_bound) {
        // Go through each child of the snarl and check the paths on outgoing edges.
        // Split up sets if the paths have different edges leaving this child
        // TODO: This is doubling the work because each edges is looked at twice
        distance_index->for_each_child(snarl, [&] (const handlegraph::net_handle_t& child) {
            for (bool go_left : {true, false}) {
                check_outgoing_edges(child, go_left);
            }
            return true;
        });// end for_each_child of the snarl
    }

    // We have now partitioned the paths into equivalence sets based on the edges they take in this netgraph,
    // stored in old_sets.
    // Return actual sets of paths.
    // Don't return anything with a set number of 0, which means that the path didn't go through this snarl

    std::vector<std::set<stoat::sample_hap_t>> sample_sets (old_set_count-1);
    for (size_t i = 0 ; i < all_samples.size() ; i++) {
        if (old_sets[i] != 0) {
            sample_sets[old_sets[i]-1].emplace(all_samples[i]);
        }
    }
    stoat::LOG_TRACE("Found walk sets ");
    for (const auto& s : sample_sets) {
        stoat::LOG_TRACE( "Set");
        for (const auto& x : s) {
            stoat::LOG_TRACE((std::stringstream) "" << "\t" << x );
        }
    }
    return sample_sets;
}



// Run iteratee on all snarls, either from the distance index or in snarl_partitions
void SnarlTraverserAndPartitioner::for_each_snarl_partition(const handlegraph::PathPositionHandleGraph& graph,
                                                  const std::function<void(const snarl_partition_t& snarl_info)>& iteratee) {
    if (!use_loaded_partitions && distance_index != nullptr) {
        // If the distance index is given, then use that

        // Get a list of all chains in root
        std::vector<handlegraph::net_handle_t> chains;
        chains.reserve(50);
        handlegraph::net_handle_t root = distance_index->get_root();
        distance_index->for_each_child(root, [&] (handlegraph::net_handle_t chain) {
            chains.emplace_back(chain);
            return true;
        });
        // Go through the contents of chains in parallel
        #pragma omp parallel shared(chains, all_references, snarl_partitions)
        {
            // The actual while loop is run on a single thread
            #pragma omp single
            {
                while (!chains.empty()) {
                    handlegraph::net_handle_t chain = chains.back();
                    chains.pop_back();

                    // Everything in here is parallelized
                    #pragma omp task
                    {
                        
                        distance_index->for_each_child(chain, [&] (handlegraph::net_handle_t snarl) {
                        
                            //TODO: Actually use is_eligible
                            //TODO: For now it's fine to check is_eligible here because it's only checking size and we don't want to look at small chains anyway
                            if (distance_index->is_snarl(snarl) && snarl_is_eligible(snarl)) {

                                stoat::LOG_TRACE( "Test snarl " + distance_index->net_handle_as_string(snarl));

                                // Make a snarl_partition_t and call iteratee on it
                                snarl_partition_t snarl_info(snarl, graph, *distance_index);
                                if (check_distances) {
                                    snarl_info.min_length = distance_index->minimum_length(snarl);
                                    snarl_info.max_length = distance_index->maximum_length(snarl);
                                }

                                // Get the offsets of the start and end nodes along the reference
                                std::vector<stoat::path_range_t> ranges = stoat::get_coordinates_of_snarl(graph, *distance_index, snarl, true, reference_sample, false);
                                if (ranges.size() != 0) {
                                    std::tie(snarl_info.ref_path, snarl_info.start_positions, snarl_info.end_positions) = get_name_and_offsets_of_snarl_path_range(graph, ranges.front());
                                } else {
                                    snarl_info.ref_path = "NA";
                                    snarl_info.start_positions = 0;
                                    snarl_info.end_positions = 0;
                                }

                                // Now get the partitions
                                snarl_info.partitions = partition_samples_in_snarl(graph, snarl);

                                // And call iteratee
                                iteratee(snarl_info);

                                #pragma omp critical 
                                {
                                    if (save_partitions) {
                                        // If we are going to serialize the snarls, then save the snarl to snarl_partitions
                                        // TODO :or just write it directly
                                        all_references.emplace(snarl_info.ref_path);
                                        snarl_partitions.emplace_back(std::move(snarl_info));
                                    }

                                    // Add the child chains to the stack
                                    distance_index->for_each_child(snarl, [&] (handlegraph::net_handle_t child) {

                                        chains.emplace_back(child);
                                        return true;
                                    });
                                }

                            }
                            return true;
                        });
                    }
                }
            }
        }

    } else {
        // If the distance index is not given, then go through snarl_partitions
        for (const snarl_partition_t& snarl_info : snarl_partitions) {
            iteratee(snarl_info);
        }
    }
}
bool SnarlTraverserAndPartitioner::snarl_is_eligible(const handlegraph::net_handle_t& snarl) const {
    if (!check_distances) {
        // If the distance index doesn't let us check distances, just return true
        return true;
    } else {
        return distance_index->maximum_length(snarl) >= allele_size_limit;
    }
}

/////////////////////////////////////////// Serializing and de-serializing the snarl partitions
/*
This needs to hold snarl_partition_t's which contain:
net_handle_t snarl; 
std::pair<size_t, size_t> snarl_ids;
std::string ref_path (as an index);
size_t min_length;
size_t max_length;
size_t start_positions;
size_t end_positions;
size_t depth;
std::vector<std::set<sample_hap_t>> partitions (as indices);
std::vector<std::string> type_variants; (can get from partitions I think) 

Start with the header (SnarlTraverserAndPartitioner::file_header)

Next the limits for testing snarls, in case we skipped small snarls for example:
allele_size_limit:N

Next all sample/haplotypes. The order will be the index for sample/haplotypes used when storing partitions
This section starts with #SAMPLES

Then a second list of reference path names. The order will be the index for ref_path
This section starts with #REFS

Each snarl_partition_t is then stored, one per line, as a vector of integers.
the first 9 items are fixed length.
The 10th item is the number of sample_hap_t's in the first partition, followed by that number of entries for the sample_hap_t. And so on
This section starts with #SNARLS
//TODO: Everything is an int so store the bytes

*/
void SnarlTraverserAndPartitioner::serialize(const std::string& filename) {
    ofstream outstream;
    outstream.open(filename);
    if (!outstream.good()) {
        throw std::runtime_error("stoat: could not open file " + filename + " for writing");
    }
    // Write the header
    outstream << file_header << endl;

    outstream << "allele_size_limit:" << allele_size_limit << endl;
    
    // Next will be a list of sample/haplotypes. Also keep a dict mapping sample/haplotype to the index that is used to store it
    outstream << "#SAMPLES" << endl;
    std::unordered_map<stoat::sample_hap_t, size_t> sample_to_index;
    size_t i = 0;
    for (const auto& sample : all_sample_haplotypes) {
        outstream << sample.sample << "\t" << sample.haplotype << endl;
        sample_to_index[sample] = i;
        i++;
    }
    // Next will be a list of reference path names. Also keep a dict mapping reference path names the index that is used to store it
    outstream << "#REFS" << endl;
    std::unordered_map<string, size_t> ref_to_index;
    size_t ref_index = 0;
    i = 0;
    for (const auto& ref : all_references) {
        outstream << ref << endl;
        ref_to_index[ref] = i;
        i++;
    }
    //Finally the snarls
    outstream << "#SNARLS" << endl;
    for (const snarl_partition_t& snarl_partition : snarl_partitions) {
        outstream << handlegraph::as_integer(snarl_partition.snarl) << "\t"
                  << handlegraph::as_integer(snarl_partition.start_handle) << "\t"
                  << handlegraph::as_integer(snarl_partition.end_handle) << "\t"
                  << ref_to_index[snarl_partition.ref_path] << "\t"
                  << snarl_partition.min_length << "\t"
                  << snarl_partition.max_length << "\t"
                  << snarl_partition.start_positions << "\t"
                  << snarl_partition.end_positions << "\t"
                  << snarl_partition.depth << "\t";
        for (const std::set<sample_hap_t> partition : snarl_partition.partitions) {
            outstream << partition.size() << "\t";
            for (const sample_hap_t& sample : partition) {
                outstream << sample_to_index[sample] << "\t"; 
            }
        }
        outstream << endl;
    }


    outstream.close();
    return;
}

void SnarlTraverserAndPartitioner::deserialize(const std::string& filename, const handlegraph::PathPositionHandleGraph& graph) {
    ifstream instream;
    instream.open(filename);

    string line;

    // Read the first line, which must match the header
    std::getline(instream, line);
    if (line != file_header) {
        throw std::runtime_error("stoat: Snarl partitions file " +filename+ " contains the wrong header: " + line);
    }

    // Next should be the allele size limit
    std::getline(instream, line);
    std::stringstream linestream(line);
    string allele_limit_str;
    std::getline(linestream, allele_limit_str, ':');
    std::getline(linestream, allele_limit_str, ':');
    if (allele_size_limit < std::stoull(allele_limit_str)) {
        cerr << "warning [stoat]: The allele_size_limit of the saved snarls file is larger than the given allele_size_limit. Some snarls may be missed" << endl;;
    }


    std::getline(instream, line);
    if (line != "#SAMPLES") {
        throw std::runtime_error("stoat: Snarl partitions file " +filename+ " is not formatted correctly");
    }

    // Get the sample/haplotype from the index
    std::vector<stoat::sample_hap_t> samples;
    std::getline(instream, line);
    while (line != "#REFS") {
        std::stringstream linestream(line);
        string sample_name;
        std::getline(linestream, sample_name, '\t');
        string hap_str;
        std::getline(linestream, hap_str, '\t');
        samples.emplace_back(sample_name, std::stoull(hap_str));
        std::getline(instream, line);
    }

    // Get the reference path names from the index
    std::vector<std::string> refs;
    std::getline(instream, line);
    while (line != "#SNARLS") {
        refs.emplace_back(std::move(line));
        std::getline(instream, line);
    }

    // Get the snarls
    while (std::getline(instream,line)) {

        snarl_partitions.emplace_back();

        std::stringstream linestream(line);
        string part;
        //Snarl net handle
        std::getline(linestream, part, '\t');
        snarl_partitions.back().snarl = handlegraph::as_net_handle(std::stoll(part));
        // snarl start id
        std::getline(linestream, part, '\t');
        snarl_partitions.back().start_handle = handlegraph::as_handle(std::stoull(part));
        //snarl end id
        std::getline(linestream, part, '\t');
        snarl_partitions.back().end_handle = handlegraph::as_handle(std::stoull(part));
        // Put together the ids from the handles
        snarl_partitions.back().snarl_ids = std::make_pair<size_t, size_t>(graph.get_id(snarl_partitions.back().start_handle),
                                                                           graph.get_id(snarl_partitions.back().end_handle));
        //reference path
        std::getline(linestream, part, '\t');
        snarl_partitions.back().ref_path = refs[std::stoull(part)];
        // min lenght
        std::getline(linestream, part, '\t');
        snarl_partitions.back().min_length = std::stoull(part);
        // max length
        std::getline(linestream, part, '\t');
        snarl_partitions.back().max_length = std::stoull(part);
        // start pos
        std::getline(linestream, part, '\t');
        snarl_partitions.back().start_positions = std::stoull(part);
        // end pos 
        std::getline(linestream, part, '\t');
        snarl_partitions.back().end_positions = std::stoull(part);
        // depth
        std::getline(linestream, part, '\t');
        snarl_partitions.back().depth = std::stoull(part);

        //Next is the list of partitions
        std::vector<std::set<sample_hap_t>>& partitions = snarl_partitions.back().partitions;
        while (std::getline(linestream, part, '\t')) {
            partitions.emplace_back();
            size_t sample_count = std::stoull(part);
            for (size_t i = 0 ; i < sample_count ; i++) {
                std::getline(linestream, part, '\t');
                partitions.back().emplace(samples[std::stoull(part)]);
            }
        }
    }


    instream.close();
    return;
}
}
