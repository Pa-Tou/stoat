#include "snarl_data_collection.hpp"
#include <fstream>

//#define DEBUG_SNARL_DATA_COLLECTION

using namespace std;
namespace stoat {


// Constructor
SnarlDataCollection::SnarlDataCollection(size_t allele_size_limit, size_t snarl_child_limit, size_t walk_steps_limit) :
                    allele_size_limit(allele_size_limit),
                    snarl_child_limit(snarl_child_limit),
                    walk_steps_limit(walk_steps_limit) {}


// This goes through all the snarls and fills in the data
void SnarlDataCollection::fill_in_snarl_info(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                                             bool find_sample_sets_first,
                                             bool walks_requested,
                                             const std::function<void(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                                                                      const net_handle_t& snarl, const snarl_info_t& snarl_data,
                                                                      std::vector<std::vector<handlegraph::net_handle_t>>& walks)>& find_walks,
                                             bool sample_set_requested,
                                             const std::function<void(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,

                                                                       const net_handle_t& snarl, const snarl_info_t& snarl_data,
                                                                       std::vector<std::set<sample_hap_t>>& sample_sets_by_allele)>& find_sample_sets,
                                             bool sequence_requested,
                                             const string& reference_sample, bool check_distances) {
    
    // Get a list of all chains in root
    std::vector<handlegraph::net_handle_t> chains;
    chains.reserve(50);
    handlegraph::net_handle_t root = distance_index.get_root();
    distance_index.for_each_child(root, [&] (handlegraph::net_handle_t chain) {
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
    // Keep track of which samples we've seen and their index in sample_haplotypes
    std::unordered_map<stoat::sample_hap_t, size_t> sample_haplotype_to_index;
    
    // Go through the contents of chains in parallel
    // Everything touching chains needs to be in an omp critical block so they don't collide. 
    #pragma omp parallel shared(chains, keep_going, chains_added, chains_processed, all_snarl_data, snarl_to_walks, snarl_to_sample_sets, snarl_to_sequences, reference_names, sample_haplotypes)
    {
        // The actual while loop is run on a single thread
        #pragma omp single
        {
            while (keep_going) {
                handlegraph::net_handle_t chain;
                #pragma omp critical(snarl_collection)
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
                    
                    distance_index.for_each_child(chain, [&] (handlegraph::net_handle_t snarl) {
                    
                        //TODO: Actually use is_eligible
                        //TODO: For now it's fine to check is_eligible here because it's only checking size and we don't want to look at small chains anyway
                        if (distance_index.is_snarl(snarl)) {
    
                            #ifdef DEBUG_SNARL_DATA_COLLECTION
                            cerr << "At snarl " << distance_index.net_handle_as_string(snarl) << endl;
                            #endif
                            if (snarl_is_eligible(distance_index, snarl, check_distances)) {

                                // Make the snarl_info_internal_t to fill in. Since it's multithreaded its better to move() it instead of adding it here
                                snarl_info_internal_t snarl_data;

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
                                std::vector<stoat::path_range_t> ranges = stoat::get_coordinates_of_snarl(graph, distance_index, snarl, true, reference_sample, false);
                                if (ranges.size() != 0) {
                                    // Check if we have already seen the reference path and if not add it
                                    size_t ref_index;

                                    auto reference_range = get_name_and_offsets_of_snarl_path_range(graph, ranges.front());
                                    snarl_data.start_position = std::get<1>(reference_range);
                                    snarl_data.end_position = std::get<2>(reference_range);


                                    #pragma omp critical(snarl_collection)
                                    {
                                        if (reference_name_to_index.count(std::get<0>(reference_range)) == 0) {
                                            ref_index = reference_names.size();
                                            reference_name_to_index[std::get<0>(reference_range)] = ref_index;
                                            reference_names.emplace_back(std::move(std::get<0>(reference_range)));
                                        } else {
                                            ref_index = reference_name_to_index[std::get<0>(reference_range)];
                                        }
                                    }
                                    snarl_data.reference_index = ref_index;

                                } else {
                                    snarl_data.start_position = 0;
                                    snarl_data.end_position = 0;
                                    snarl_data.reference_index = std::numeric_limits<size_t>::max();
                                }

                                // Optionally fill in the walks, sample sets, and sequences.
                                // All three of these vectors must have the same number of entries because they correspond to the same alleles
                                std::vector<std::vector<handlegraph::net_handle_t>> walks_by_allele; 
                                std::vector<stoat::Path_traversal_t> walks_by_allele_as_paths; 
                                std::vector<std::set<sample_hap_t>> sample_sets_by_allele; 
                                std::vector<std::string> snarl_sequences;

                                // TODO decide how to deal with snarls without walks
                                bool save_snarl = true;

                                // Make the snarl_info_t passed to the sample set/walk finders. They don't need to have all the information yet
                                // the snarl_info is const in the finders so it won't change the walks/sample sets/sequences
                                // except when references to them are passed as the thing we're filling in
                                snarl_info_t new_snarl_info (snarl_data.start_node, 
                                                             snarl_data.end_node, 
                                                             snarl_data.reference_index == std::numeric_limits<size_t>::max() ? "NA" 
                                                                                        : reference_names.at(snarl_data.reference_index),
                                                             snarl_data.start_position,
                                                             snarl_data.end_position,
                                                             snarl_data.depth,
                                                             "", //No variant type yet
                                                             walks_by_allele_as_paths,
                                                             sample_sets_by_allele,
                                                             snarl_sequences);

                                if (find_sample_sets_first) {
                                    // Find sample_sets_by_allele then walks
                                    if (sample_set_requested) {
                                       find_sample_sets(graph, distance_index, snarl, new_snarl_info, sample_sets_by_allele); 
                                    }
                                    if (walks_requested) {
                                        find_walks(graph, distance_index, snarl, new_snarl_info, walks_by_allele);
                                    }
                                } else {
                                    // Find walks then sample_sets_by_allele
                                    if (walks_requested) {
                                        find_walks(graph, distance_index, snarl, new_snarl_info, walks_by_allele);
                                    }
                                    if (sample_set_requested) {
                                       find_sample_sets(graph, distance_index, snarl, new_snarl_info, sample_sets_by_allele); 
                                    }
                                }
                                if (walks_requested) {
                                    std::vector<std::string> walk_lengths;
                                    std::tie(walks_by_allele_as_paths, walk_lengths) = stoat::fill_pretty_paths(distance_index, graph, walks_by_allele);
                                    snarl_data.variant_type = stoat::vectorToString(walk_lengths);
                                    if (sequence_requested) {
                                        // Find the sequences
                                        snarl_sequences = get_sequences_from_walks(graph, distance_index, walks_by_allele_as_paths);
                                    }
                                }
                                if (sequence_requested && !walks_requested) {
                                    throw std::runtime_error("stoat: Snarl data collection requested sequences without walks");
                                }

                                if (save_snarl) { 
   
                                    // Add the snarl to the collection
                                    if (walks_requested) {
                                        #pragma omp critical(snarl_collection)
                                        {
                                        snarl_to_walks.emplace(snarl_data.start_node, std::move(walks_by_allele_as_paths));
                                        }
                                    }
                                    if (sample_set_requested) {
                                        // This has its own omp guards
                                        add_sample_sets_to_collection(sample_haplotype_to_index, sample_sets_by_allele, snarl_data.start_node); 
                                    }
                                    if (sequence_requested) {
                                        #pragma omp critical(snarl_collection)
                                        {
                                        snarl_to_sequences.emplace(snarl_data.start_node, std::move(snarl_sequences));
                                        }
                                    }
                                    #pragma omp critical(snarl_collection)
                                    {
                                    all_snarl_data.emplace_back(std::move(snarl_data));
                                    }



                                }
                            }// end if snarl_is_eligible

    
                            #pragma omp critical(snarl_collection)
                            {
   
                                // Add the child chains to the stack
                                distance_index.for_each_child(snarl, [&] (handlegraph::net_handle_t child) {
    
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
                #pragma omp critical(snarl_collection)
                {
                    keep_going = !chains.empty();
                }
    
                if (!keep_going) {
                    // Wait for tasks to complete
                    #pragma omp taskwait
    
                    #pragma omp critical(snarl_collection)
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
void SnarlDataCollection::add_snarl_sample_sets(std::unordered_map<stoat::sample_hap_t, size_t>& sample_haplotype_to_index, const std::function<std::vector<std::set<sample_hap_t>>(const snarl_info_t& snarl_data)>& find_sample_sets,
                          std::string chr){

    #pragma omp parallel for schedule(dynamic)
    for (const snarl_info_internal_t& snarl_info : all_snarl_data) {

        // If we are limiting to a chromosome (by reference path), then skip anything not on this chromosome
        if (chr != "" && (snarl_info.reference_index == std::numeric_limits<size_t>::max() || chr != reference_names.at(snarl_info.reference_index))) {
            continue;
        }

        std::vector<std::set<sample_hap_t>> empty_sample_sets;

        // Make the snarl_info_t from the information we have
        std::vector<Path_traversal_t> empty_walks (0); 
        std::vector<std::string> empty_sequences (0);
        snarl_info_t new_snarl_info (snarl_info.start_node, 
                                snarl_info.end_node, 
                                snarl_info.reference_index == std::numeric_limits<size_t>::max() ? "NA" : reference_names.at(snarl_info.reference_index),
                                snarl_info.start_position,
                                snarl_info.end_position,
                                snarl_info.depth,
                                snarl_info.variant_type,
                                snarl_to_walks.count(snarl_info.start_node) ? snarl_to_walks.at(snarl_info.start_node) : empty_walks,
                                empty_sample_sets,
                                snarl_to_sequences.count(snarl_info.start_node) ? snarl_to_sequences.at(snarl_info.start_node) : empty_sequences);

        std::vector<std::set<sample_hap_t>> new_sample_sets = find_sample_sets(new_snarl_info);

        // Add it to the collection
        add_sample_sets_to_collection(sample_haplotype_to_index, new_sample_sets, snarl_info.start_node);
    }

}

void SnarlDataCollection::add_sample_sets_to_collection(std::unordered_map<stoat::sample_hap_t, size_t>& sample_haplotype_to_index, 
                                                        const std::vector<std::set<sample_hap_t>>& sample_sets, const Node_traversal_t& snarl_start) {

    // Remake sample_sets by index
    std::vector<std::set<size_t>> sample_sets_by_index;
    for (size_t i = 0 ; i < sample_sets.size() ; i++) {
        sample_sets_by_index.emplace_back();
        for (const sample_hap_t& sample : sample_sets[i]) {
    
            // Get the index into sample_haplotypes
            size_t sample_index;
            // Needs to be in an omp block because it accesses sample_haplotypes
            #pragma omp critical(snarl_collection)
            {
                if (sample_haplotype_to_index.count(sample)) {
                    sample_index = sample_haplotype_to_index[sample];
                } else {
                    sample_index = sample_haplotypes.size();
                    sample_haplotype_to_index[sample] = sample_index;
                    sample_haplotypes.emplace_back(sample);
                } 
            }
    
            sample_sets_by_index.back().emplace(sample_index);
        }
    }
    #pragma omp critical(snarl_collection)
    {
        snarl_to_sample_sets.emplace(snarl_start, std::move(sample_sets_by_index));
    }
}



// Call interatee for all snarls
void SnarlDataCollection::for_each_snarl(const std::function<void(const snarl_info_t& snarl_info)>& iteratee) const {
    for (const snarl_info_internal_t& snarl_info : all_snarl_data) {

        std::vector<std::set<sample_hap_t>> sample_sets;

        if (snarl_to_sample_sets.count(snarl_info.start_node)) {
            const std::vector<std::set<size_t>>& set_ints = snarl_to_sample_sets.at(snarl_info.start_node);
            for (const std::set<size_t>& sample_set : set_ints) {
                sample_sets.emplace_back();
                for (const size_t& i : sample_set) {
                    sample_sets.back().emplace(sample_haplotypes[i]);
                }
            }
        }


        std::vector<Path_traversal_t> empty_walks (0); 
        std::vector<std::string> empty_sequences (0);
        snarl_info_t new_snarl_info (snarl_info.start_node, 
                              snarl_info.end_node, 
                              snarl_info.reference_index == std::numeric_limits<size_t>::max() ? "NA" : reference_names.at(snarl_info.reference_index),
                              snarl_info.start_position,
                              snarl_info.end_position,
                              snarl_info.depth,
                              snarl_info.variant_type,
                              snarl_to_walks.count(snarl_info.start_node) ? snarl_to_walks.at(snarl_info.start_node) : empty_walks,
                              sample_sets,
                              snarl_to_sequences.count(snarl_info.start_node) ? snarl_to_sequences.at(snarl_info.start_node) : empty_sequences);
        iteratee(new_snarl_info);
    }
}


void SnarlDataCollection::get_all_walks_through_snarl(
        const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
        const net_handle_t& snarl, const snarl_info_t& snarl_data, std::vector<std::vector<handlegraph::net_handle_t>>& walks ) {


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
    
        // TODO: Put this back
        //for (const auto& net : path) {
        //    if (++dict_path_occ[net] > walk_cycle_limit + 1) {
        //        cycle = true;
        //        break;
        //    }
        //}
    
        // TODO: Add out_fail
        // TODO: Get the child count properly
        //if (steps_taken++ > walk_steps_limit) {
        //    #pragma omp critical(out_fail)
        //    out_fail << distance_index.node_id(distance_index.get_bound(snarl, false, true)) << "_" << distance_index.node_id(distance_index.get_bound(snarl, true, true)) 
        //             << "\titeration_calculation_out = " << endl;// << children << " children\n";
        //    break_snarl = true;
        //    break;
        //}

        // Follow edges from the last element in path
        if (!path.empty()) {
            distance_index.follow_net_edges(path.back(), &graph, false, [&](const handlegraph::net_handle_t& next_child) {
                // If this is the bound of the snarl then we're done
                if (distance_index.is_sentinel(next_child)) {

                    size_t next_child_node_id = distance_index.node_id(distance_index.get_node_from_sentinel(next_child));
                    size_t first_element_path_node_id = distance_index.node_id(distance_index.get_node_from_sentinel(path[0]));

                    // Only keep the walk if it entered and exited the snarl at opposite sides
                    if (next_child_node_id != first_element_path_node_id) {
                        walks.emplace_back(path);
                        walks.back().push_back(next_child);
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
            });
        }
    }
    
    if (break_snarl) {
        walks.clear();
    }
    
    return ;
    
}

void SnarlDataCollection::get_walks_from_sample_sets(
        const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
        const net_handle_t& snarl, const snarl_info_t& snarl_data, std::vector<std::vector<handlegraph::net_handle_t>>& walks ) {

    // Get the walks by tracing the haplotype paths through the snarl


    ///////////////////////////////////////////// Get one step_handle_t's for the each sample set

    handlegraph::handle_t start_in = distance_index.get_handle(distance_index.get_bound(snarl, false, true), &graph);
    handlegraph::handle_t end_in = distance_index.get_handle(distance_index.get_bound(snarl, true, true), &graph);

    std::vector<handlegraph::PathSense> senses = {handlegraph::PathSense::GENERIC,
                                                  handlegraph::PathSense::REFERENCE,
                                                  handlegraph::PathSense::HAPLOTYPE};

    std::vector<handlegraph::step_handle_t> first_steps;
    for (const std::set<sample_hap_t>& sample_set : snarl_data.sample_sets_by_allele) {
        bool found_step = false;
        for (const auto& sense : senses) {
            graph.for_each_step_of_sense(start_in, sense, [&](const handlegraph::step_handle_t& step) {
                if (sample_set.count(stoat::get_sample_and_haplotype(graph, graph.get_path_handle_of_step(step))) ) {
                    first_steps.emplace_back(step);
                    found_step = true;
                    // Return false to stop looping through steps
                    return false;
                } else{
                    return true;
                }
            });
            if (found_step) {
                // Break out of the inner loop and continue to the next sample_set
                break;
            }
        }
    }
    #ifdef DEBUG_SNARL_DATA_COLLECTION
    assert(first_steps.size() == snarl_data.sample_sets_by_allele.size());
    #endif

    //////////////////////////////////// Follow each of the paths through the snarl and create a walk

    for (const handlegraph::step_handle_t& first_step : first_steps) {
        //Accumulate the min and max lengths of this walk, based on the lengths of snarls
        size_t min_length = 0;
        size_t max_length = 0;


        // Get the step of this path on both the start and end nodes.
        std::vector<handlegraph::step_handle_t> boundary_steps;
        for (const auto& sense : senses) {
            graph.for_each_step_of_sense(start_in, sense, [&](const handlegraph::step_handle_t& this_step) {
                if (graph.get_path_handle_of_step(first_step) == graph.get_path_handle_of_step(this_step)) {
                    boundary_steps.emplace_back(this_step);
                 }
                return true;
            }); 
            graph.for_each_step_of_sense(end_in, sense, [&](const handlegraph::step_handle_t& this_step) {
                if (graph.get_path_handle_of_step(first_step) == graph.get_path_handle_of_step(this_step)) {
                    boundary_steps.emplace_back(this_step);
                 }
                return true;
            });
        }

        // Now sort the list of steps on the two boundary nodes
        sort(boundary_steps.begin(), boundary_steps.end(), [&](const handlegraph::step_handle_t& a, const handlegraph::step_handle_t& b) {
            return graph.get_position_of_step(a) < graph.get_position_of_step(b);
        });
        if (boundary_steps.size() > 1) {
            // Go through the boundary nodes up to the next-to-last one and find the walk to the next one
            //TODO: This assumes that the path goes into the snarl then exits, it will fail if the path started 
            // inside the snarl and the first traversal of a boundary node is leaving it
            for (size_t boundary_i = 0 ; boundary_i < boundary_steps.size()-1 ; boundary_i+= 2) {
                // Start a new walk
                //TODO: It would be more efficient to calculate the Path_traversal_t and length counts directly here but I don't want to copy code too much
                walks.emplace_back();
                std::vector<bdsg::net_handle_t>& current_walk = walks.back();
                current_walk.emplace_back(distance_index.get_net(graph.get_handle_of_step(boundary_steps[boundary_i]), &graph));



                handlegraph::step_handle_t step = graph.get_next_step(boundary_steps[boundary_i]);
                while (step != boundary_steps[boundary_i+1]) {
                    // Step is a node inside the snarl, and it may be nested inside children of the snarl
                    // If the node is a child of the snarl, then its parent is a trivial chain and its grandparent is the snarl
                    // If it is the first node in a chain that is the child of the snarl, then it is one of its parent's boundary nodes pointing
                    //   into the chain and its grandparent is the snarl.
                    // In any other case we want to ignore it
                    // TODO: I think this will be faster than skipping to the end of the chain and looking for the right path
                    
                    handlegraph::handle_t node = graph.get_handle_of_step(step);
                    handlegraph::net_handle_t node_net = distance_index.get_net(node, &graph);
                    handlegraph::net_handle_t parent = distance_index.get_parent(node_net);

                    // For the path, add an empty node for each time we leave and re-enter the snarl
                    if (boundary_i != 0) {
                        current_walk.emplace_back(distance_index.get_root());
                    }


                    // Add to the path, depending on if this is a node or chain
                    if (distance_index.get_parent(parent) == snarl) {
                        if (distance_index.is_trivial_chain(parent)) {
                            // This node is really a node in the snarl, then add it to the path 

                            current_walk.emplace_back(parent);

                        } else if (node_net == distance_index.get_bound(parent, false, true) ||  node_net == distance_index.get_bound(parent, true, true)) {
                            // This node is going into the child chain going backward

                            current_walk.emplace_back(parent);
                        }
                    }

                    step = graph.get_next_step(step);

                }//end while loop going through a traversal of the snarl
                //Add the end boundary node
                current_walk.emplace_back(distance_index.get_net(graph.get_handle_of_step(boundary_steps[boundary_i+1]), &graph));

            }// end for loop through each traversal of the snarl in one sample_set
        }// end if there are enough boundary nodes


    }// end for each first step (per sample_set)

    return ;
}

std::vector<std::string> SnarlDataCollection::get_sequences_from_walks(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                 const std::vector<stoat::Path_traversal_t>& paths) const {
    std::vector<std::string> sequences;
    for (const stoat::Path_traversal_t& path : paths) {
        sequences.emplace_back();
        for (const stoat::Node_traversal_t& node : path.get_paths()) {
            if (node.get_node_id() != 0) {
                //TODO: Does this take into account the reverse complement?
                sequences.back() += graph.get_sequence(graph.get_handle(node.get_node_id(), node.get_is_reverse()));
            }
        }
    }
    return sequences;
}



bool SnarlDataCollection::snarl_is_eligible( const bdsg::SnarlDistanceIndex& distance_index, const handlegraph::net_handle_t& snarl, bool check_distances) const {
    bool pass = true;

    // If we have distances in the index, make sure that the snarl's maximum lenght is small enough
    if (check_distances) {
        pass =  pass && (allele_size_limit <= distance_index.maximum_length(snarl));
    }
    //TODO: Once the libbdsg branch is merged we can use this instead of going through all the children to count them
    //pass &= snarl_child_limit <= distance_index.get_snarl_child_count(snarl);

    // Count children
    size_t children = 0;
    distance_index.for_each_child(snarl, [&](const handlegraph::net_handle_t& snarl) {
        children++;
        return true;
    });
    pass = pass && (snarl_child_limit >= children);

    return pass;
}




/////////////////////////////////////////// Writing/reading the snarl data
/*
This needs to hold snarl_info_internal_t's which contain:

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

Next all sample/haplotypes. These get stored in sample_haplotypes and the order will be the index for sample/haplotypes used when storing sample_sets
This section starts with #SAMPLES


Each snarl_info_t is then stored, one per line, as a tab-separated vector of integers/strings.
the first 7 items are the contents of the snarl_info_internal_t.
The 8th item is all of the walks, comma separated.
The sample sets are stored next (if present). Each sample set starts with the number of sample/haplotypes for that sample set, followed with
that number of entries for the sample_hap_t. There should be N sample sets.
The next N items (if present) are the sequences of the Path_traversal_t's. These can be distinguished from the haplotype sample sets because they can only contain A/C/G/T/N, or may be empty
This section starts with #SNARLS

*/
void SnarlDataCollection::write_snarl_data_collection(std::ostream& outstream) const {

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
        outstream << samp.sample << "\t" << samp.haplotype << endl;
    }

    //Finally the snarls
    outstream << "#SNARLS" << endl;
    for (const snarl_info_internal_t& snarl_data : all_snarl_data) {

        // Start with just the contents of the snarl_info_internal_t
        outstream << snarl_data.start_node.to_string() << "\t"
                  << snarl_data.end_node.to_string() << "\t"
                  << snarl_data.reference_index << "\t"
                  << snarl_data.start_position << "\t"
                  << snarl_data.end_position << "\t"
                  << snarl_data.depth << "\t"
                  << snarl_data.variant_type << "\t";

        // Next, always include the walks as a single comma-separated string
        outstream << stoat::vectorPathToString(snarl_to_walks.at(snarl_data.start_node)) << "\t";

        // Next add the sample sets, if there are any
        // Each sample set is saved as the number of items in it then the items
        if (!snarl_to_sample_sets.empty()) {
            #ifdef DEBUG_SNARL_DATA_COLLECTION
            assert(snarl_to_walks.at(snarl_data.start_node).size() == snarl_to_sample_sets.at(snarl_data.start_node).size());
            #endif
            for (const std::set<size_t>& sample_set : snarl_to_sample_sets.at(snarl_data.start_node)){
                outstream << sample_set.size() << "\t";
                for (const size_t& x : sample_set) {
                    outstream << x << "\t";
                }
            }
        }

        // Add the sequences, if any, as comma separated strings
        if (!snarl_to_sequences.empty()) {
            #ifdef DEBUG_SNARL_DATA_COLLECTION
            assert(snarl_to_walks.at(snarl_data.start_node).size() == snarl_to_sequences.at(snarl_data.start_node).size());
            #endif
            for (const std::string& seq : snarl_to_sequences.at(snarl_data.start_node)){
                outstream << seq << ",";
            }
        }
        
        outstream << endl;
    }


    return;
}

void SnarlDataCollection::load_snarl_data_collection(std::istream& instream) {

    string line;

    // Read the first line, which must match the header
    std::getline(instream, line);
    if (line != file_header) {
        throw std::runtime_error("stoat: Snarl data file contains the wrong header: " + line);
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
        throw std::runtime_error("stoat: Snarl file is not formatted correctly");
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
        all_snarl_data.back().reference_index = std::stoull(part);

        // Start offset along reference
        std::getline(linestream, part, '\t');
        all_snarl_data.back().start_position = std::stoull(part);

        // End offset along reference
        std::getline(linestream, part, '\t');
        all_snarl_data.back().end_position = std::stoull(part);

        // Snarl depth
        std::getline(linestream, part, '\t');
        all_snarl_data.back().depth = std::stoull(part);

        // variant type
        std::getline(linestream, part, '\t');
        all_snarl_data.back().variant_type = part;

        // Next are the paths, all as one comma-separated string
        std::getline(linestream, part, '\t');
        snarl_to_walks[all_snarl_data.back().start_node] = stoat::stringToVectorPath(part);
        size_t walk_count = snarl_to_walks[all_snarl_data.back().start_node].size();

        // This may be the end, or there may be sample_sets and/or sequences
        if (!std::getline(linestream, part, '\t')) {
            continue;
        }
        bool finished_line = false;

        if (!(part.at(0) == ',' || part.at(0) == 'A' || part.at(0) == 'C' || part.at(0) == 'G' || part.at(0) == 'T')) {
            // If the next thing isn't a sequence then we need to get the sample_sets next
            std::vector<std::set<size_t>> sample_sets;

            for (size_t i = 0 ; i < walk_count ; i++) {

                sample_sets.emplace_back();
                
                // Since we had to check for the end of the line, part is currently the count of samples
                // and we update it at the end of the loop 
                size_t sample_count = std::stoull(part);
                for (size_t i = 0 ; i < sample_count ; i++) {
                    std::getline(linestream, part, '\t');
                    sample_sets.back().emplace(std::stoull(part));
                }

                // Check if this is the end of the line
                if (!std::getline(linestream, part, '\t')) {
                    finished_line = true;
                }
            }

            #ifdef DEBUG_SNARL_DATA_COLLECTION
            assert(sample_sets.size() == walk_count);
            #endif

            snarl_to_sample_sets[all_snarl_data.back().start_node] = std::move(sample_sets);
        }

        // The last thing we got is now either the end of the line or the first sequence
        if (!finished_line) {
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


    return;
}
}

