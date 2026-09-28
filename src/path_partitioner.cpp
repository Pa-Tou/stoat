#include "path_partitioner.hpp"
#include "log.hpp"
#include <fstream>
#include "snarl_traversals.hpp"

//#define DEBUG_PATH_PARTITIONER

using namespace stoat;
namespace stoat_graph {


// This is supposed to partition the paths in the snarl by the walks they take through the netgraph.
// Instead of explicitly enumerating the paths, it actually finds the sets of edges that each path takes.
// But since paths may loop, it actually finds the order and number of outgoing edges from each node.
// I think this is equivalent to partitioning by the actual sets of unique walks.
std::vector<size_t> partition_embedded_paths_in_snarl(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                          const net_handle_t& snarl,
                          const std::vector<stoat::sample_hap_t>& all_sample_haplotypes) {

    #ifdef DEBUG_PATH_PARTITIONER
    std::cerr <<  "Get walk sets of " << distance_index.net_handle_as_string(snarl) << std::endl;;
    #endif

    // Map each sample(plus haplotype) to its index in all_sample_haplotypes
    std::unordered_map<stoat::sample_hap_t, std::size_t> sample_to_index;
    for (size_t i = 0 ; i < all_sample_haplotypes.size() ; i++) {
        sample_to_index[all_sample_haplotypes[i]] = i;
    }

    // For each path (along all_sample_haplotypes), the index of the set it is currently in
    // Everything starts out in the same set, 0, which represents not being in the snarl
    std::vector<size_t> old_sets (all_sample_haplotypes.size(), 0);
    size_t old_set_count = 1;

    // TODO: I think this is unnecessary
    //// First, check for paths that don't go through this snarl
    //std::vector<handlegraph::PathSense> senses = {handlegraph::PathSense::GENERIC,
    //                                              handlegraph::PathSense::REFERENCE,
    //                                              handlegraph::PathSense::HAPLOTYPE};
    //
    ////TODO: This could also use steps_of_handle() but it doesn't seem to work 
    //for (const auto& sense : senses) {
    //    graph.for_each_step_of_sense(distance_index.get_handle(distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, true)), &graph),
    //        sense, [&](const handlegraph::step_handle_t& step) {
    //        old_sets[sample_to_index[stoat::get_sample_and_haplotype(graph, graph.get_path_handle_of_step(step))]] = 1;
    //        old_set_count = 2;
    //    });
    //}


    // A struct representing a path's edge going out from child
    // Represents the next node and orientation and the offset along the path of the edge
    // Use (0,0,0,false) as a default value to indicate that the path didn't traverse this child
    // This becomes a linked list when a path traverses the child multiple times (this doesn't happen often) 
    // The linked list is sorted by order in the traversal
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
    // Goes through the edges leaving a child and checks with paths take which edges.
    // old_sets and old_set_count are updated to split up sets that take different edges at this node 
    // If require_node is true (used for the bounds going in), then any path that doesn't go through this node
    // is considered to not be in the snarl and its set number is reset to 0.
    // Second required node is true if this is the boundary node checked, in which case the path should be skipped
    // if it didn't pass through the first boundary node 
    // For this to work, check_outgoing_edges must be called on the two bounds last
    auto check_outgoing_edges = [&] (const handlegraph::net_handle_t& child, bool go_left, bool require_node, bool second_required_node) {

        #ifdef DEBUG_PATH_PARTITIONER
        std::cerr << "At snarl child " << distance_index.net_handle_as_string(child) << " going " << (go_left ? "left" : "right") << std::endl;
        #endif
        
        std::vector<path_edge_t> next_steps (all_sample_haplotypes.size());
        // For paths that go through multiple steps 
        std::vector<path_edge_t> additional_steps;
        
        // Get a handle to the end of the child
        handlegraph::handle_t handle;
        if ( distance_index.is_trivial_chain(child)) {
            handle =  (go_left ? graph.flip(distance_index.get_handle(child, &graph)) 
                               : distance_index.get_handle(child, &graph));
        } else if (distance_index.is_sentinel(child)) {
            handle = distance_index.get_handle(child, &graph);
        } else {
            handle = distance_index.get_handle(distance_index.get_bound(child, go_left, false), &graph);
        }
        
        std::vector<handlegraph::PathSense> senses = {handlegraph::PathSense::GENERIC,
                                                      handlegraph::PathSense::REFERENCE,
                                                      handlegraph::PathSense::HAPLOTYPE};

        #ifdef DEBUG_PATH_PARTITIONER
        std::cerr << "\tgraph handle " << graph.get_id(handle) << " going " << (graph.get_is_reverse(handle) ? "left" : "right") << std::endl;
        #endif

        for (const auto& sense : senses) {
            graph.for_each_step_of_sense(handle, sense, [&](const handlegraph::step_handle_t& step) {
                // For each step on the node handle, keep track of which paths take different steps

                #ifdef DEBUG_PATH_PARTITIONER
                std::cerr << "\ton path " << graph.get_path_name(graph.get_path_handle_of_step(step)) << std::endl;
                #endif

                sample_hap_t step_sample_haplotype (graph, graph.get_path_handle_of_step(step));

                // If this is not a sample of interest, skip it
                if (!sample_to_index.count(step_sample_haplotype)) {
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
        
                #ifdef DEBUG_PATH_PARTITIONER
                std::cerr << "" << "\t\tgoing to " << graph.get_id(next_handle) << std::endl;
                #endif
        
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

        #ifdef DEBUG_PATH_PARTITIONER
        // To check that we properly used all additional edges, count how many we follow
        size_t additional_edge_count = 0;
        #endif
        
        //This maps each edge list for this node/direction to the intermediate set index
        std::map<std::vector<std::pair<handlegraph::nid_t, bool>>, size_t> edge_to_intermediate_set;

        //Everything starts in the same set, representing not going through this node
        std::vector<size_t> intermediate_sets (old_sets.size(), 0);
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
                    #ifdef DEBUG_PATH_PARTITIONER
                    additional_edge_count++;
                    #endif
                }
                if (edge_to_intermediate_set.count(edge_list) == 0) {
                    edge_to_intermediate_set[edge_list] = intermediate_set_count;
                    intermediate_set_count++;
                }
                intermediate_sets[path_i] = edge_to_intermediate_set[edge_list];
            }
        }
        #ifdef DEBUG_PATH_PARTITIONER
        std::cerr << "Intermediate sets: " << std::endl;
        for (size_t i = 0 ; i < all_sample_haplotypes.size() ; i++) {
            std::cerr << "" << "\t" << all_sample_haplotypes[i] << ": " << intermediate_sets[i] << std::endl;
        } 
        assert(additional_edge_count == additional_steps.size());
        #endif
        
        // We now have an old set and an intermediate set for each path
        // Assign the path to a new set. Everything gets a new set
        std::vector<size_t> new_sets (intermediate_sets.size(), std::numeric_limits<size_t>::max());

        // Map pairs of <old_set, intermediate_set> to new set number
        std::map<std::pair<size_t, size_t>, size_t>  old_to_new_set;

        // Reserve 0 for paths that didn't go through this snarl
        old_to_new_set[std::make_pair(0,0)] = 0;
        size_t new_set_count = 1;

        for (size_t path_i = 0 ; path_i < new_sets.size() ; path_i++) {
            if (require_node && (intermediate_sets[path_i] == 0 || (second_required_node && old_sets[path_i] == 0 ))) {
                // If this is a boundary node and the path never traverses it, then add it to the set
                // that doesn't go through the snarl
                new_sets[path_i] = 0;
            } else {
                std::pair<size_t, size_t> old_set (old_sets[path_i], intermediate_sets[path_i]);
                if (old_to_new_set.count(old_set) == 0) {
                    old_to_new_set[old_set] = new_set_count;
                    new_set_count++;
                }
                new_sets[path_i] = old_to_new_set[old_set];
            }
        }
        
        old_sets = std::move(new_sets);
        old_set_count = new_set_count;
        
        #ifdef DEBUG_PATH_PARTITIONER
        std::cerr << "New sets: " << std::endl;
        for (size_t i = 0 ; i < all_sample_haplotypes.size() ; i++) {
            std::cerr << "" << "\t" << all_sample_haplotypes[i] << ": " << old_sets[i] << std::endl;
        } 
        #endif
    };

    // Now do the work of going through the edges for the start bound and each child in both directions


    if (!distance_index.is_regular_snarl(snarl, true, &graph)) {
        // Go through each child of the snarl and check the paths on outgoing edges.
        // Split up sets if the paths have different edges leaving this child
        // TODO: This is doubling the work because each edges is looked at twice
        distance_index.for_each_child(snarl, [&] (const handlegraph::net_handle_t& child) {
            for (bool go_left : {true, false}) {
                check_outgoing_edges(child, go_left, false, false);
            }
            return true;
        });// end for_each_child of the snarl
    }

    check_outgoing_edges(distance_index.get_bound(snarl, false, true), false, true, false);
    check_outgoing_edges(distance_index.get_bound(snarl, true, true), false, true, true);

    // Now go through the sets and change 0 to inf, and decrement all others
    for (size_t i = 0 ; i < old_sets.size() ; i++) {
        if (old_sets[i] == 0) {
            old_sets[i] = std::numeric_limits<size_t>::max();
        } else {
            --old_sets[i];
        }
    }



    #ifdef DEBUG_PATH_PARTITIONER
    std::vector<std::set<sample_hap_t>> sample_sets_by_allele(old_set_count-1);
    std::cerr << "Check for " << all_sample_haplotypes.size() << " sample" << std::endl;
    for (size_t i = 0 ; i < all_sample_haplotypes.size() ; i++) {
        if (old_sets[i] != std::numeric_limits<size_t>::max()) {
            sample_sets_by_allele[old_sets[i]].emplace(all_sample_haplotypes[i]);
        }
    }
    std::cerr << "Found walk sets " << std::endl;
    for (const auto& s : sample_sets_by_allele) {
        std::cerr <<  "Set" << std::endl;
        for (const auto& x : s) {
            std::cerr << "" << "\t" << x << std::endl;
        }
    }
    #endif

    // old_sets is now a vector that assigns each sample/haplotype in all_sample_haplotypes to a allele number
    // Return this
    return old_sets;
}




std::vector<size_t> partition_embedded_paths_in_snarl_with_gbwt(const handlegraph::PathPositionHandleGraph& graph, 
                          const gbwt::GBWT& gbwt,
                          const gbwt::FastLocate& r_index,
                          const bdsg::SnarlDistanceIndex& distance_index,
                          const net_handle_t& snarl,
                          const std::vector<stoat::sample_hap_t>& all_sample_haplotypes,
                          std::vector<PathTraversal>& paths_per_allele) {



    // The final boundary-to-boundary paths
    // The size_t is an identifier for the path, since some paths may be identical
    std::vector<gbwt_path_t> finished_paths;

    // Get all start-end or start-start traversals and put them in  finished_paths 
    // Path count will be an upper bound of the number of distinct paths
    // TODO: WHy is it an upper bound and not the exact value?
    size_t path_count = get_gbwt_traversals(graph, gbwt, r_index, distance_index, snarl, finished_paths);

    // At this point, finished paths holds a path for each distinct, haplotype-supported walk through th esnarl netgraph
    // Since a sample may go through the snarl multiple times, we now want to partition the samples based on how many
    // of each walk they take.

    // Map each sample(plus haplotype) to its index in all_sample_haplotypes
    std::map<stoat::sample_hap_t, std::size_t> sample_to_index;
    for (size_t i = 0 ; i < all_sample_haplotypes.size() ; i++) {
        sample_to_index[all_sample_haplotypes[i]] = i;
    }

    // For each sample, remember which paths through the snarl it takes. Since it may take multiple, save this as a linked list.
    // Since it doesn't happen very often, there will be one vector of all the heads of the linked list, and a separate vector of
    // the subsequent nodes/repeats.

    // A struct representing the allele/path taken by a sample.
    // This becomes a linked list when a sample traverses the snarl multiple times (this doesn't happen often) 
    // A value of {max(), max()} indicates that no path was taken
    struct path_link_t {
        size_t path_id;          // The path
        size_t offset;           // The offset of this path (traversal of the snarl) along the haplotype, so we can order in additional_paths 
        size_t additional_path; // Index into a list of additional paths

        path_link_t() : path_id(std::numeric_limits<size_t>::max()), offset(std::numeric_limits<size_t>::max()), additional_path(std::numeric_limits<size_t>::max()){};
        path_link_t(size_t path_id, size_t offset) : 
            path_id(path_id), offset(offset), additional_path(std::numeric_limits<size_t>::max()){};
    };
    std::vector<path_link_t> path_by_sample (all_sample_haplotypes.size());
    std::vector<path_link_t> additional_paths;

    auto add_path_to_sample = [&](size_t sample_id, size_t path_id, size_t offset) {
        path_link_t this_path_link (path_id, offset);
        if (path_by_sample.at(sample_id).path_id == std::numeric_limits<size_t>::max()) {
            // If this is the first time we've found a path for this sample
            path_by_sample.at(sample_id) = std::move(this_path_link);
        } else {
            // If this sample already has a path, place this path in the correct place in the linked list

            if (offset < path_by_sample.at(sample_id).offset) {
                // If this path is the first for this sample, just shift the linked list
                this_path_link.additional_path = additional_paths.size();
                additional_paths.emplace_back(std::move(path_by_sample.at(sample_id)));
                path_by_sample.at(sample_id) = std::move(this_path_link);
            } else {

                // The current and next indices in the linked list
                // since current_index may be an index into the original path_by_sample, remember this 
                size_t current_index = sample_id; 
                bool current_index_is_first = true;
                size_t next_index = path_by_sample.at(sample_id).additional_path;

                while (next_index != std::numeric_limits<size_t>::max() &&
                       additional_paths[next_index].offset < offset) {
                    current_index = next_index;
                    next_index = additional_paths.at(current_index).additional_path;
                    current_index_is_first = false;
                }
                //current_index now points to the item just smaller than the current path
                // and next_index is the item just larger than the current path
                this_path_link.additional_path = next_index;
                if (current_index_is_first) {
                    path_by_sample.at(current_index).additional_path = additional_paths.size();
                } else {
                    additional_paths.at(current_index).additional_path = additional_paths.size();
                }
                additional_paths.emplace_back(std::move(this_path_link)); 
            }
        }
                   
    };

    // Remember the paths themselves, once per path_id
    std::vector<std::vector<handlegraph::net_handle_t>> paths_per_path_id (path_count);


    // Sort the finished paths so by path id so that identical paths are consecutive in the vector.
    // Since actually sorting the vector would be slow, make a vector of indices and sort that
    std::vector<size_t> sort_order(finished_paths.size(), 0);
    for (size_t i = 0 ; i < sort_order.size() ; i++) {
        sort_order[i] = i;
    }
    std::sort(sort_order.begin(), sort_order.end(), [&](size_t a, size_t b) { return finished_paths.at(a).identifier < finished_paths.at(b).identifier; });

    auto gbwt_reference_samples = gbwtgraph::parse_reference_samples_tag(gbwt);


    // Go through all finished paths by path id
    for (size_t sorted_i = 0 ; sorted_i < sort_order.size() ; sorted_i++) {
        const gbwt_path_t& current_path = finished_paths.at(sort_order.at(sorted_i));
        #ifdef DEBUG_PATH_PARTITIONER
            std::cerr << "At path id " << current_path.identifier << ":\t";
            for (const auto& net : current_path.path) {
                std::cerr << distance_index.net_handle_as_string(net) << ",";
            }
            std::cerr << std::endl;
        #endif

        // Now use the r-index to find the occurrences of this path - the range in the search state
        // We use the suffix array offset to find this information for each haplotype in the range
        gbwt::size_type current_sa_index = current_path.sa_index;
        for (size_t i = 0 ; i < current_path.search_state.size() ; i++) {
            gbwt::size_type gbwt_path_id = gbwt::Path::id(r_index.seqId(current_sa_index));
            gbwt::size_type path_offset = r_index.seqOffset(current_sa_index);

            gbwt::size_type next_sa_index = r_index.locateNext(current_sa_index);

            handlegraph::PathSense sense = gbwtgraph::get_path_sense(gbwt, gbwt_path_id, gbwt_reference_samples);

            std::string path_name = handlegraph::PathMetadata::create_path_name(sense,
                                        gbwtgraph::get_path_sample_name(gbwt, gbwt_path_id, sense),
                                        gbwtgraph::get_path_locus_name(gbwt, gbwt_path_id, sense),
                                        gbwtgraph::get_path_haplotype(gbwt, gbwt_path_id, sense),
                                        gbwtgraph::get_path_phase_block(gbwt, gbwt_path_id, sense),
                                        gbwtgraph::get_path_subrange(gbwt, gbwt_path_id, sense)); 
            #ifdef DEBUG_PATH_PARTITIONER
                std::cerr << "\t\tpath " << path_name << " takes this snarl traversal: " << gbwt_path_id << std::endl;
            #endif


            add_path_to_sample(sample_to_index.at(sample_hap_t(path_name)), current_path.identifier, path_offset);

            current_sa_index = next_sa_index;
        }

        if (sorted_i == sort_order.size()-1 || current_path.identifier != finished_paths.at(sort_order.at(sorted_i+1)).identifier) {
            // If we finished going through the last set of paths for this path id, remember the path itself
            paths_per_path_id.at(current_path.identifier) = std::move(current_path.path);
        }
    }

    // we have now filled in path_by_sample to assign each sample to a path (or sequence of paths in additional_paths)

    // Now that we've gone though all paths and assigned paths to samples, do the partitioning by assigning an allele id to each path/linked list of paths
    std::vector<size_t> allele_assignments (all_sample_haplotypes.size(), std::numeric_limits<size_t>::max());
    size_t current_allele_num = 0;
    // For each distinct path (snarl traversal, which may include multiple traversals of the snarl), map it to the new allele number
    std::map<std::vector<size_t>, size_t> path_list_to_allele_num;
    // For each allele, get one sample with this allele. Used to get the paths
    std::vector<size_t> sample_per_allele;
    paths_per_allele.clear();
    for (size_t sample_num = 0 ; sample_num < path_by_sample.size() ; sample_num++) {
        size_t first_path_id = path_by_sample[sample_num].path_id;
        if (first_path_id == std::numeric_limits<size_t>::max()) {
            // If this sample didn't have a path through this snarl
            continue;
        }

        // First, get the list of path ids for this sample
        std::vector<size_t> current_path ({first_path_id});
        size_t next_link = path_by_sample[sample_num].additional_path;
        while (next_link != std::numeric_limits<size_t>::max()) {
            current_path.push_back(additional_paths[next_link].path_id);
            next_link = additional_paths[next_link].additional_path;
        }

        // If we've seen this before, then we have the allele number so assign it
        const auto& allele_assignment = path_list_to_allele_num.find(current_path);
        if (allele_assignment != path_list_to_allele_num.end()) {
            allele_assignments[sample_num] = allele_assignment->second;
        } else {
            allele_assignments[sample_num] = current_allele_num;
            path_list_to_allele_num[current_path] = current_allele_num;

            // Record the path for this allele
            paths_per_allele.emplace_back();
            #ifdef DEBUG_PATH_PARTITIONER
            assert(current_allele_num == sample_per_allele.size());
            #endif
            sample_per_allele.push_back(sample_num);        
            // Get the path from current_path, which has the path id for each sample
            for (size_t path_num = 0 ; path_num < current_path.size() ; path_num++) {
                size_t path_id =  current_path[path_num];
                if (path_num > 0) {
                    // Add an out-of-snarl marker
                    paths_per_allele.back().add_out_of_snarl_walk();
                }
                const std::vector<handlegraph::net_handle_t>& path_as_net_handles = paths_per_path_id.at(path_id);
                
                // Add the path through the snarl
                for (size_t i = 0 ; i < path_as_net_handles.size() ; i++) {
                    const handlegraph::net_handle_t& net = path_as_net_handles.at(i);
                    paths_per_allele.back().add_net_handle(net, distance_index, i == 0 || i == path_as_net_handles.size()-1);
                }
            }

            ++current_allele_num;
        }
    }

    return allele_assignments;
}


}
