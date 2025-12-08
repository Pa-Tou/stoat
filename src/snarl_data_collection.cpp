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
                                             const std::vector<stoat::sample_hap_t>& sample_haplotypes,
                                             bool find_alleles_first,
                                             bool walks_requested,
                                             const std::function<void(const net_handle_t& snarl, const snarl_info_t& snarl_data,
                                                                      std::vector<stoat::PathTraversal>& walks)>& find_walks,
                                             bool alleles_requested,
                                             const std::function<std::vector<size_t>(const net_handle_t& snarl, const snarl_info_t& snarl_data,
                                                                                     const std::vector<stoat::sample_hap_t>& all_sample_haplotypes)>& find_alleles_by_sample,
                                             bool sequence_requested,
                                             const string& reference_sample, bool check_distances) {

    // Make a copy of the samples. Not a reference since it has to be able to be loaded from a file and stored
    all_sample_haplotypes = sample_haplotypes;

    size_t sample_index = 0;
    // Fill in sample_to_index
    for (const sample_hap_t& sample_hap : all_sample_haplotypes) {
        if (!sample_to_index.count(sample_hap.sample)) {
            sample_to_index.emplace(sample_hap.sample, sample_index++);
        }
    }
    
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
    
    // Go through the contents of chains in parallel
    // Everything touching chains needs to be in an omp critical block so they don't collide. 
    #pragma omp parallel shared(chains, keep_going, chains_added, chains_processed, all_snarl_data, snarl_to_walks, snarl_to_alleles_by_sample, snarl_to_sequences, reference_names, all_sample_haplotypes)
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

                                // Optionally fill in the walks, alleles, and sequences.
                                // walks_by_allele and snarl_sequences must have the same number of entries because they correspond to the same alleles
                                // The entries in alleles_by_sample correspond to these alleles so its max value (that is not inf) must be the length of the others
                                std::vector<stoat::PathTraversal> walks_by_allele; 
                                std::vector<size_t> alleles_by_sample; 
                                std::vector<std::string> snarl_sequences;

                                //This might cause problems because it is a reference but it doesn't get used so I think its fine
                                // I don't want to use the actual samples_to_index because then the empty genotype table with allocate memory for the vector
                                GenotypeTable empty_genotypes(std::unordered_map<std::string, size_t>(), 0);

                                // TODO decide how to deal with snarls without walks
                                bool save_snarl = true;

                                // Make the snarl_info_t passed to the sample set/walk finders. They don't need to have all the information yet
                                // the snarl_info is const in the finders so it won't change the walks/alleles/sequences
                                // except when references to them are passed as the thing we're filling in
                                snarl_info_t new_snarl_info (snarl_data.start_node, 
                                                             snarl_data.end_node, 
                                                             snarl_data.reference_index == std::numeric_limits<size_t>::max() ? "NA" 
                                                                                        : reference_names.at(snarl_data.reference_index),
                                                             snarl_data.start_position,
                                                             snarl_data.end_position,
                                                             snarl_data.depth,
                                                             empty_genotypes,
                                                             walks_by_allele,
                                                             snarl_sequences);

                                if (find_alleles_first) {
                                    // Find alleles_by_sample then walks
                                    if (alleles_requested) {
                                       alleles_by_sample = find_alleles_by_sample(snarl, new_snarl_info, all_sample_haplotypes); 
                                    }
                                    if (walks_requested) {
                                        find_walks(snarl, new_snarl_info, walks_by_allele);
                                    }
                                } else {
                                    // Find walks then alleles
                                    if (walks_requested) {
                                        find_walks(snarl, new_snarl_info, walks_by_allele);
                                    }
                                    if (alleles_requested) {
                                       alleles_by_sample = find_alleles_by_sample(snarl, new_snarl_info, all_sample_haplotypes); 
                                    }
                                }
                                if (walks_requested) {
                                    //
                                    //walks_by_allele_as_paths = stoat::convert_path_traversals(distance_index, graph, walks_by_allele);
                                    #ifdef DEBUG_SNARL_DATA_COLLECTION
                                    cerr << "got path from walk" << endl;
                                    for (const auto& path : walks_by_allele) {
                                        cerr << path.to_string() << endl;
                                    }
                                    #endif

                                    if (sequence_requested) {
                                        // Find the sequences
                                        snarl_sequences = get_sequences_from_walks(graph, distance_index, walks_by_allele);
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
                                        snarl_to_walks.emplace(snarl_data.start_node, std::move(walks_by_allele));
                                        }
                                    }
                                    if (alleles_requested) {
                                        // TODO: This could be done earlier but I don't think it's a big deal
                                        size_t max_allele = 0;
                                        for (size_t x : alleles_by_sample) {
                                            if (x != std::numeric_limits<size_t>::max()) {
                                                max_allele = std::max(max_allele, x);
                                            }
                                        }
                                        #pragma omp critical(snarl_collection)
                                        {
                                        snarl_to_alleles_by_sample.emplace(snarl_data.start_node, allele_by_sample_t(max_allele+1, std::move(alleles_by_sample)));
                                        }
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
void SnarlDataCollection::add_alleles_by_sample(const std::function<std::vector<size_t>(const snarl_info_t& snarl_data, 
                                                                                        const std::vector<stoat::sample_hap_t>& all_sample_haplotypes)>& find_alleles_by_sample,
                                                std::string chr){

    #pragma omp parallel for schedule(dynamic)
    for (const snarl_info_internal_t& snarl_info : all_snarl_data) {

        // If we are limiting to a chromosome (by reference path), then skip anything not on this chromosome
        if (chr != "" && (snarl_info.reference_index == std::numeric_limits<size_t>::max() || chr != reference_names.at(snarl_info.reference_index))) {
            continue;
        }

        std::vector<size_t> empty_alleles_by_sample;
        
        //This might cause problems because it is a reference but it doesn't get used so I think its fine
        // I don't want to use the actual samples_to_index because then the empty genotype table with allocate memory for the vector
        GenotypeTable empty_genotypes(std::unordered_map<std::string, size_t>(), 0);

        // Make the snarl_info_t from the information we have
        std::vector<PathTraversal> empty_walks (0); 
        std::vector<std::string> empty_sequences (0);
        snarl_info_t new_snarl_info (snarl_info.start_node, 
                                snarl_info.end_node, 
                                snarl_info.reference_index == std::numeric_limits<size_t>::max() ? "NA" : reference_names.at(snarl_info.reference_index),
                                snarl_info.start_position,
                                snarl_info.end_position,
                                snarl_info.depth,
                                empty_genotypes,
                                snarl_to_walks.count(snarl_info.start_node) ? snarl_to_walks.at(snarl_info.start_node) : empty_walks,
                                snarl_to_sequences.count(snarl_info.start_node) ? snarl_to_sequences.at(snarl_info.start_node) : empty_sequences);

        // Get the alleles from the snarl_info_t
        std::vector<size_t> new_alleles_by_sample = find_alleles_by_sample(new_snarl_info, all_sample_haplotypes);
        size_t max_allele = 0;
        for (size_t x : new_alleles_by_sample) {
            if (x != std::numeric_limits<size_t>::max()) {
                max_allele = std::max(x, max_allele);
            }
        }

        // Save it
        #pragma omp critical(snarl_collection)
        {
            snarl_to_alleles_by_sample.emplace(snarl_info.start_node, allele_by_sample_t(max_allele+1, std::move(new_alleles_by_sample)));
        }
    }

}


// Call interatee for all snarls
void SnarlDataCollection::for_each_snarl(const std::function<void(const snarl_info_t& snarl_info)>& iteratee) const {
    for (const snarl_info_internal_t& snarl_info : all_snarl_data) {

        // GenotypeTable constructor takes a map from sample to index, and the number of alleles
        GenotypeTable genotypes(sample_to_index, snarl_to_alleles_by_sample.count(snarl_info.start_node) ? snarl_to_alleles_by_sample.at(snarl_info.start_node).allele_count
                                                                                                         : 0);

        // Go through the alleles_by_sample vector for this snarl and add the counts to the genotype table
        // alleles_by_sample is a vector with the allele for each sample in all_sample_haplotypes
        if (snarl_to_alleles_by_sample.count(snarl_info.start_node)) {
            const allele_by_sample_t alleles_by_sample = snarl_to_alleles_by_sample.at(snarl_info.start_node);
            for (size_t sample_hap_i = 0 ; sample_hap_i < alleles_by_sample.alleles.size() ; sample_hap_i++) {
                genotypes.increment_count(all_sample_haplotypes[sample_hap_i].sample, alleles_by_sample.alleles[sample_hap_i]);
            }
        }


        std::vector<PathTraversal> empty_walks (0); 
        std::vector<std::string> empty_sequences (0);
        snarl_info_t new_snarl_info (snarl_info.start_node, 
                              snarl_info.end_node, 
                              snarl_info.reference_index == std::numeric_limits<size_t>::max() ? "NA" : reference_names.at(snarl_info.reference_index),
                              snarl_info.start_position,
                              snarl_info.end_position,
                              snarl_info.depth,
                              genotypes,
                              snarl_to_walks.count(snarl_info.start_node) ? snarl_to_walks.at(snarl_info.start_node) : empty_walks,
                              snarl_to_sequences.count(snarl_info.start_node) ? snarl_to_sequences.at(snarl_info.start_node) : empty_sequences);
        iteratee(new_snarl_info);
    }
}


void SnarlDataCollection::get_all_walks_through_snarl(
        const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
        const net_handle_t& snarl, const snarl_info_t& snarl_data, std::vector<stoat::PathTraversal>& walks, 
        size_t walk_cycle_limit ) {


    // Path exploration
    std::vector<std::vector<handlegraph::net_handle_t>> paths = {
        {distance_index.get_bound(snarl, false, true)}
    };
    
    std::vector<std::vector<handlegraph::net_handle_t>> walks_as_net_handles;
    
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
        for (const auto& net : path) {
            if (++dict_path_occ[net] > walk_cycle_limit + 1) {
                cycle = true;
                break;
            }
        }
    
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
                        walks_as_net_handles.emplace_back(path);
                        walks_as_net_handles.back().push_back(next_child);
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
        walks_as_net_handles.clear();
    }
    walks = stoat::convert_path_traversals(distance_index, graph, walks_as_net_handles);   
    return ;
}

void SnarlDataCollection::get_walks_from_alleles(
        const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
        const net_handle_t& snarl, const snarl_info_t& snarl_data, std::vector<stoat::PathTraversal>& walks) {

    // Get the walks by tracing the haplotype paths through the snarl


    ///////////////////////////////////////////// Get one step_handle_t for the each allele

    handlegraph::handle_t start_in = distance_index.get_handle(distance_index.get_bound(snarl, false, true), &graph);
    handlegraph::handle_t end_in = distance_index.get_handle(distance_index.get_bound(snarl, true, true), &graph);

    std::vector<handlegraph::PathSense> senses = {handlegraph::PathSense::GENERIC,
                                                  handlegraph::PathSense::REFERENCE,
                                                  handlegraph::PathSense::HAPLOTYPE};
    const allele_by_sample_t& alleles_by_sample = snarl_to_alleles_by_sample.at(snarl_data.start_node);

    // For each allele, did we find a step for it yet?
    std::vector<bool> got_step(alleles_by_sample.allele_count, false);

    // One step per allele
    std::vector<handlegraph::step_handle_t> first_steps(alleles_by_sample.allele_count);

    for (size_t sample_i = 0 ; sample_i < alleles_by_sample.alleles.size() ; sample_i++) {
        size_t allele_num = alleles_by_sample.alleles[sample_i];
        if (got_step.at(allele_num)) {
            // If we already found a step for this allele, skip it
            continue;
        }

        // Look for a step on this sample haplotype
        const sample_hap_t& target_sample = all_sample_haplotypes[sample_i];
        bool found_step = false;
        for (const auto& sense : senses) {
            graph.for_each_step_of_sense(start_in, sense, [&](const handlegraph::step_handle_t& step) {

                if (target_sample == stoat::sample_hap_t(graph, graph.get_path_handle_of_step(step))) {
                    // If this step is on the haplotype we want, then remember it and remember that we got it
                    first_steps[allele_num] = step;
                    got_step[allele_num] = true;
                    found_step = true;
                    // Return false to stop looping through steps
                    return false;
                } else {
                    // Continue to other haplotypes
                    return true;
                }
            });
            if (found_step) {
                // Break out of the inner loop and continue to the next sample and allele
                break;
            }
        }
        #ifdef DEBUG_SNARL_DATA_COLLECTION
        // Right now, the path partitioner only assigns a path to a partition if it traverses both bounds of the snarl.
        assert(found_step);
        #endif
    }
    #ifdef DEBUG_SNARL_DATA_COLLECTION
    // Check that we got a step for each allele
    for (bool check : got_step) {
        assert(check);
    }
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
        if (boundary_steps.size() > 1) {
            // Start a new walk
            walks.emplace_back();
            stoat::PathTraversal& current_walk = walks.back();

            // Now sort the list of steps on the two boundary nodes
            sort(boundary_steps.begin(), boundary_steps.end(), [&](const handlegraph::step_handle_t& a, const handlegraph::step_handle_t& b) {
                return graph.get_position_of_step(a) < graph.get_position_of_step(b);
            });

            // Go through the boundary nodes up to the next-to-last one and find the walk to the next one
            //TODO: This assumes that the path goes into the snarl then exits, it will fail if the path started 
            // inside the snarl and the first traversal of a boundary node is leaving it

            for (size_t boundary_i = 0 ; boundary_i < boundary_steps.size()-1 ; boundary_i+= 2) {
                //TODO: It would be more efficient to calculate the PathTraversal and length counts directly here but I don't want to copy code too much

                // For the path, add an empty node for each time we leave and re-enter the snarl
                if (boundary_i != 0) {
                    stoat::Node_traversal_t traversal (0, true);
                    current_walk.add_node_traversal_t(traversal);
                }

                // Add the step of the boundary node going into the snarl for each time it re-enters the snarl
                // TODO: Make the PathTraversal add from a handle_t or net_handle_t
                //current_walk.emplace_back(distance_index.get_net(graph.get_handle_of_step(boundary_steps[boundary_i]), &graph));
                handlegraph::handle_t start_handle = graph.get_handle_of_step(boundary_steps[boundary_i]);
                stoat::Node_traversal_t start_traversal (graph.get_id(start_handle), graph.get_is_reverse(start_handle));
                current_walk.add_node_traversal_t(start_traversal);


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



                    // Add to the path, depending on if this is a node or chain
                    if (distance_index.get_parent(parent) == snarl) {
                        if (distance_index.is_trivial_chain(parent)) {
                            // This node is really a node in the snarl, then add it to the path 

                            handlegraph::handle_t node_handle = distance_index.get_handle(parent, &graph);
                            stoat::Node_traversal_t node_traversal (graph.get_id(node_handle), graph.get_is_reverse(node_handle));
                            current_walk.add_node_traversal_t(node_traversal);

                            min_length += distance_index.minimum_length(parent);
                            max_length += distance_index.maximum_length(parent);

                        } else if (node_net == distance_index.get_bound(parent, false, true)) {
                            // This node is going into the child chain going forward

                            // Add the start bound going in
                            handlegraph::handle_t chain_start_handle = distance_index.get_handle(distance_index.get_bound(parent, false, true), &graph);
                            stoat::Node_traversal_t chain_start_traversal (graph.get_id(chain_start_handle), graph.get_is_reverse(chain_start_handle));
                            current_walk.add_node_traversal_t(chain_start_traversal);

                            // Add the interior of the chain as a fake node
                            stoat::Node_traversal_t chain_traversal (0, false);
                            current_walk.add_node_traversal_t(chain_traversal);

                            // Add the end bound going out
                            handlegraph::handle_t chain_end_handle = distance_index.get_handle(distance_index.get_bound(parent, true, false), &graph);
                            stoat::Node_traversal_t chain_end_traversal (graph.get_id(chain_end_handle), graph.get_is_reverse(chain_end_handle));
                            current_walk.add_node_traversal_t(chain_end_traversal);

                            min_length += distance_index.minimum_length(parent);
                            max_length += distance_index.maximum_length(parent);

                        } else if (node_net == distance_index.get_bound(parent, true, true)) {
                            // This node is going into the child chain going backward

                            // Add the end bound going in
                            handlegraph::handle_t chain_start_handle = distance_index.get_handle(distance_index.get_bound(parent, true, true), &graph);
                            stoat::Node_traversal_t chain_start_traversal (graph.get_id(chain_start_handle), graph.get_is_reverse(chain_start_handle));
                            current_walk.add_node_traversal_t(chain_start_traversal);

                            // Add the interior of the chain as a fake node
                            stoat::Node_traversal_t chain_traversal (0, false);
                            current_walk.add_node_traversal_t(chain_traversal);

                            // Add the start bound going out
                            handlegraph::handle_t chain_end_handle = distance_index.get_handle(distance_index.get_bound(parent, false, false), &graph);
                            stoat::Node_traversal_t chain_end_traversal (graph.get_id(chain_end_handle), graph.get_is_reverse(chain_end_handle));
                            current_walk.add_node_traversal_t(chain_end_traversal);

                            min_length += distance_index.minimum_length(parent);
                            max_length += distance_index.maximum_length(parent);
                        }
                    }

                    step = graph.get_next_step(step);

                }//end while loop going through a traversal of the snarl

                // Add the bound
                handlegraph::handle_t end_handle = graph.get_handle_of_step(step);
                stoat::Node_traversal_t end_traversal (graph.get_id(end_handle), graph.get_is_reverse(end_handle));
                current_walk.add_node_traversal_t(end_traversal);

            }// end for loop for a pair of boundary nodes for one traversal of the snarl
            assert(current_walk.get_path().size() >= 2);
            current_walk.add_min_allele_len(min_length);
            current_walk.add_max_allele_len(max_length);

        }// end if there are enough boundary nodes
    }// end for each first step (per allele)

    #ifdef DEBUG_SNARL_DATA_COLLECTION
    cerr << "Found walks from sets" << endl;
    for (const stoat::PathTraversal& walk : walks) {
        std::cerr << walk.to_string();
        std::cerr << std::endl;
    }
    #endif

    return ;
}

std::vector<std::string> SnarlDataCollection::get_sequences_from_walks(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
             const std::vector<stoat::PathTraversal>& paths) const {

    std::vector<std::string> sequences;
    for (const stoat::PathTraversal& path : paths) {
        sequences.emplace_back();
        const std::vector<stoat::Node_traversal_t>& nodes = path.get_path(); 
        if (nodes.size() > 0) {
            handlegraph::nid_t start_id = nodes.front().get_node_id();
            handlegraph::nid_t end_id = nodes.back().get_node_id(); 
            for (size_t i = 0 ; i < nodes.size() ; i++) {
                const stoat::Node_traversal_t& node = nodes[i];
                if (node.get_node_id() != start_id && node.get_node_id() != end_id) {
                    if (node.get_node_id() == 0) {
                        sequences.back() += "N";
                    } else {
                        //TODO: Does this take into account the reverse complement?
                        sequences.back() += graph.get_sequence(graph.get_handle(node.get_node_id(), node.get_is_reverse()));
                    }
                }
            }
        }
    }
    return sequences;
}



bool SnarlDataCollection::snarl_is_eligible( const bdsg::SnarlDistanceIndex& distance_index, const handlegraph::net_handle_t& snarl, bool check_distances) const {

    // If we have distances in the index, make sure that the snarl's maximum length is big enough
    if (check_distances && (allele_size_limit > distance_index.maximum_length(snarl))) {

        return false;
    }
    //TODO: Once the libbdsg branch is merged we can use this instead of going through all the children to count them
    //pass &= snarl_child_limit <= distance_index.get_snarl_child_count(snarl);

    // Count children
    size_t children = 0;
    distance_index.for_each_child(snarl, [&](const handlegraph::net_handle_t& snarl) {
        children++;
        return true;
    });
    if (snarl_child_limit < children) {
        cerr << "Snarl had too many children" << endl;

        return false;
    }

    return true;
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

Each snarl will also contain a vector of PathTraversal's, (optionally) a vector of sets of haplotypes that take each PathTraversal, 
and (optionally) the sequences of each PathTraversal


Start with the header (SnarlDataCollection::file_header)

Next the limits for testing snarls, in case we skipped small snarls for example:
allele_size_limit:N
snarl_child_limit:N

Next a list of reference path names. These get stored in reference_names and the order will be the index for reference_index
This section starts with #REFS

Each snarl_info_t is then stored, one per line, as a tab-separated vector of integers/strings.
The next line is the header defining the 9 items that are always there, and the sample names for the remaining columns, as sample(|haplotype num). The entry for each column is the allele number
the first 7 items are the contents of the snarl_info_internal_t.
The 8th item is all of the walks, comma separated. "." if not present 
The 9th item is all of the sequences, comma separated. "." if not present
The remaining items are the allele number for each sample. "." if the sample is not present in the snarl

*/
void SnarlDataCollection::write_snarl_data_collection(std::ostream& outstream) const {

    // Write the header
    outstream << file_header << endl;

    outstream << "#allele_size_limit:" << allele_size_limit << endl;
    outstream << "#snarl_child_limit:" << snarl_child_limit << endl;
    outstream << "#walk_steps_limit: " << walk_steps_limit << endl;

    // Next will be a list of reference path names.
    outstream << "#REFS" << endl;
    for (const auto& ref : reference_names) {
        outstream << "#" << ref << endl;
    }
    

    //Finally the snarls
    // Start with a header that will contain the names of all samples
    outstream << "#SNARLS" << endl;
    outstream << "#START_NODE\tEND_NODE\tREF\tSTART_OFFSET\tEND_OFFSET\tDEPTH\tALLELE_LENGTHS\tWALKS\tSEQUENCES";

    // The header also includes a list of sample/haplotypes
    size_t sample_count = 0;
    for (const auto& samp : all_sample_haplotypes) {
        outstream << "\t" << samp.sample << "#" + samp.haplotype;
    }
    outstream << endl;

    // Now write the snarls, one per line
    for (const snarl_info_internal_t& snarl_data : all_snarl_data) {

        // Start with just the contents of the snarl_info_internal_t
        outstream << snarl_data.start_node.to_string() << "\t"
                  << snarl_data.end_node.to_string() << "\t"
                  << snarl_data.reference_index << "\t"
                  << snarl_data.start_position << "\t"
                  << snarl_data.end_position << "\t"
                  << snarl_data.depth << "\t";

        // Next, optionally include the walks as a single comma-separated string
        if (snarl_to_walks.count(snarl_data.start_node) == 0) {
            outstream << ".\t.\t";
        } else {
            outstream << stoat::vectorPathToString(snarl_to_walks.at(snarl_data.start_node), true) << "\t";
            outstream << stoat::vectorPathToString(snarl_to_walks.at(snarl_data.start_node)) << "\t";
        }


        // Add the sequences, if any, as comma separated strings
        if (snarl_to_sequences.empty() || snarl_to_sequences.count(snarl_data.start_node) == 0) {

            outstream << ".";
        } else {

            #ifdef DEBUG_SNARL_DATA_COLLECTION
            assert(snarl_to_walks.at(snarl_data.start_node).size() == snarl_to_sequences.at(snarl_data.start_node).size());
            #endif
            const std::vector<std::string>& sequences = snarl_to_sequences.at(snarl_data.start_node);
            for (size_t i = 0 ; i < sequences.size() ; i++){
                const std::string& seq = sequences[i];
                outstream << seq ;
                if (i != sequences.size()-1) {
                    outstream << ",";
                }
            }
        }

        // Next add the allele assignments, if there are any
        if (snarl_to_alleles_by_sample.empty()) {
            for (size_t i = 0 ; i < all_sample_haplotypes.size() ; i++) {
                outstream << "\t.";
            }
        } else {

            const allele_by_sample_t& alleles_by_sample = snarl_to_alleles_by_sample.at(snarl_data.start_node);

            #ifdef DEBUG_SNARL_DATA_COLLECTION
            assert(snarl_to_walks.at(snarl_data.start_node).size() == alleles_by_sample.allele_count);
            assert(alleles_by_sample.alleles.size() == all_sample_haplotypes.size());
            #endif
            for (size_t allele_num : alleles_by_sample.alleles) {
                if (allele_num == std::numeric_limits<size_t>::max()) {
                    // . for a haplotype that doesn't traverse the snarl
                    outstream << "\t.";
                } else {
                    outstream << "\t" << allele_num;
                }
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
    {
        std::stringstream linestream(line);
        string limit_str;
        std::getline(linestream, limit_str, ':');
        #ifdef DEBUG_SNARL_DATA_COLLECTION
        assert(limit_str == "#allele_size_limit");
        #endif
        std::getline(linestream, limit_str, ':');
        if (allele_size_limit < std::stoull(limit_str)) {
            cerr << "warning [stoat]: The allele_size_limit of the saved snarls file is larger than the given allele_size_limit. Some snarls may be missed" << endl;;
        }
    }


    // And the snarl child limit
    std::getline(instream, line);
    {
        std::stringstream linestream(line);
        string limit_str;
        std::getline(linestream, limit_str, ':');
        #ifdef DEBUG_SNARL_DATA_COLLECTION
        assert(limit_str == "#snarl_child_limit");
        #endif
        std::getline(linestream, limit_str, ':');
        if (snarl_child_limit < std::stoull(limit_str)) {
            cerr << "warning [stoat]: The snarl_child_limit of the saved snarls file is larger than the given snarl_child_limit. Some snarls may be missed" << endl;;
        }
    }

    // And the walk steps limit
    std::getline(instream, line);
    {
        std::stringstream linestream(line);
        string limit_str;
        std::getline(linestream, limit_str, ':');
        #ifdef DEBUG_SNARL_DATA_COLLECTION
        assert(limit_str == "#walk_steps_limit");
        #endif
        std::getline(linestream, limit_str, ':');
        if (walk_steps_limit > std::stoull(limit_str)) {
            cerr << "warning [stoat]: The walk_steps_limit of the saved snarls file is smaller than the given walk_steps_limit. Some snarls may be missed" << endl;;
        }
    }

    // Next are the references
    std::getline(instream, line);
    if (line != "#REFS") {
        throw std::runtime_error("stoat: Snarl file is not formatted correctly");
    }

    std::getline(instream, line);
    
    // Get the reference path names from the file
    while (line != "#SNARLS") {
        std::string ref (line.begin()+1, line.end());
        reference_names.emplace_back(std::move(ref));
        std::getline(instream, line);
    }

    // The next header is  "#START_NODE\tEND_NODE\tREF\tSTART_OFFSET\tEND_OFFSET\tDEPTH\tALLELE_LENGTHS\tWALKS\tSEQUENCES", plus all of the sample/haplotypes
    std::getline(instream, line);
    {
        std::stringstream linestream(line);
        string sample_name;
        // First go through the first 9 things and ignore them
        for (size_t i = 0 ; i < 9 ; i++) {
            std::getline(linestream, sample_name, '\t');
        }
        while (std::getline(linestream, sample_name, '\t')) {
            // The sample_hap_t constructor will take care of finding the proper sample name and haplotype 
            all_sample_haplotypes.emplace_back(sample_name);
        }
    }

    // Fill in sample_to_index
    size_t sample_index = 0;
    for (const sample_hap_t& sample_hap : all_sample_haplotypes) {
        if (!sample_to_index.count(sample_hap.sample)) {
            sample_to_index.emplace(sample_hap.sample, sample_index++);
        }
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

        // Next are the paths, first the lengths then the paths themselves
        std::string path_lengths;
        std::string paths;
        std::getline(linestream, path_lengths, '\t');
        std::getline(linestream, paths, '\t');

        snarl_to_walks[all_snarl_data.back().start_node] = stoat::string_to_path_traversals(paths, path_lengths);

        size_t walk_count = snarl_to_walks[all_snarl_data.back().start_node].size();

        // Next are the sequences, or "." if there are no sequences
        std::getline(linestream, part, '\t');
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
            #ifdef DEBUG_SNARL_DATA_COLLECTION
            assert(sequences.size() == walk_count);
            #endif

            snarl_to_sequences[all_snarl_data.back().start_node] = std::move(sequences);
        }


        // The rest of the line will be the allele assignment of each sample
        // "." means that that this sample didn't have an allele in this snarl. If all samples have "." then we just didn't store the samples

        // Fill this in with the contents of the line
        bool has_samples = false;
        std::vector<size_t> allele_assignments;
        allele_assignments.reserve(all_sample_haplotypes.size());
        size_t max_allele = 0;
        while (std::getline(linestream, part, '\t')) {
            if (part == ".") {
                allele_assignments.emplace_back(std::numeric_limits<size_t>::max());
            } else {
                allele_assignments.emplace_back(std::stoull(part));
                max_allele = std::max(max_allele, allele_assignments.back());
                has_samples = true;
            }
        };

        #ifdef DEBUG_SNARL_DATA_COLLECTION
        assert(max_allele == walk_count);
        assert(allele_assignments.size() == all_sample_haplotypes.size()); 
        #endif
        if (has_samples) {
            snarl_to_alleles_by_sample[all_snarl_data.back().start_node] = allele_by_sample_t(max_allele+1, allele_assignments);
        }
    }

    return;
}
}

