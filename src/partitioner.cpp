#include "partitioner.hpp"

//#define DEBUG_PATH_PARTITIONER

using namespace std;
namespace stoat_graph {

std::vector<std::set<std::string>> PathPartitioner::partition_samples_in_snarl(const handlegraph::PathPositionHandleGraph& graph, 
                                                                               const bdsg::SnarlDistanceIndex& distance_index,
                                                                               const handlegraph::net_handle_t& snarl) const {
    #ifdef DEBUG_PATH_PARTITIONER
        cerr << "Get sample partitions of " << distance_index.net_handle_as_string(snarl) << " by its paths" << endl;
    #endif

    //Get the partition of paths, depending on if the snarl is simple or not
    std::vector<std::set<stoat::sample_hap_t>> sample_sets = distance_index.is_regular_snarl(snarl) 
                                                                ? get_start_edge_sets(graph, distance_index, snarl)
                                                                : get_walk_sets(graph, distance_index, snarl);

    #ifdef DEBUG_PATH_PARTITIONER
        cerr << "Found sets of paths using " << ( distance_index.is_regular_snarl(snarl) ? "edges from the start node" : "walk sets") << endl;
        for (const std::set<stoat::sample_hap_t>& sample_set : sample_sets) {
            cerr << "SET "<< endl;
            for (const stoat::sample_hap_t& sample : sample_set) {
                cerr << "\t" << sample.sample << endl;
            }
        }
        cerr << "TRUTH" << endl;
        for ( const std::string& x : samples_of_interest) {
            cerr << "\t" << x << endl;
        }
    #endif

    std::vector<std::set<std::string>> sample_name_sets (sample_sets.size());
    for (size_t i = 0 ; i < sample_sets.size() ; i++) {
        for (const stoat::sample_hap_t& sample : sample_sets[i]) {
            sample_name_sets[i].emplace(sample.sample);
        }
    }
    return sample_name_sets;
}

// This is supposed to partition the paths in the snarl by the walks they take through the netgraph.
// Instead of explicitly enumerating the paths, it actually finds the sets of edges that each path takes.
// But since paths may loop, it also takes into account the order and number of outgoing edges from each node.
// I think this is equivalent to partitioning by the actual sets of unique walks.
std::vector<std::set<stoat::sample_hap_t>> PathPartitioner::get_walk_sets(const handlegraph::PathPositionHandleGraph& graph, 
                                                                   const bdsg::SnarlDistanceIndex& distance_index,
                                                                   const handlegraph::net_handle_t& snarl) const {
    #ifdef DEBUG_PATH_PARTITIONER
    cerr << "Get walk sets of " << distance_index.net_handle_as_string(snarl) << endl;
    #endif

    // Make a vector of the paths 
    std::vector<stoat::sample_hap_t> all_samples(all_sample_haplotypes.begin(), all_sample_haplotypes.end());

    // Map each sample(plus haplotype) to its index in all_samples
    std::map<stoat::sample_hap_t, std::size_t> sample_to_index;
    for (size_t i = 0 ; i < all_samples.size() ; i++) {
        sample_to_index[all_samples[i]] = i;
    }

    // For each path (along all_samples), the index of the set it is currently in
    // Everything starts out in the same set
    std::vector<size_t> old_sets (all_samples.size(), 0);
    size_t old_set_count = 1;

    // First, check for paths that don't go through this snarl
    std::vector<handlegraph::PathSense> senses = {handlegraph::PathSense::GENERIC,
                                                  handlegraph::PathSense::REFERENCE,
                                                  handlegraph::PathSense::HAPLOTYPE};
    
    //TODO: This could also use steps_of_handle() but it doesn't seem to work 
    for (const auto& sense : senses) {
        graph.for_each_step_of_sense(distance_index.get_handle(distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, true)), &graph),
            sense, [&](const handlegraph::step_handle_t& step) {
            old_sets[sample_to_index[stoat::get_sample_and_haplotype(graph, graph.get_path_handle_of_step(step))]] = 1;
            old_set_count = 2;
        });
    }

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


    // Go through each child of the snarl and check the paths on outgoing edges.
    // Split up sets if the paths have different edges leaving this child
    // TODO: This is doubling the work because each edges is looked at twice
    distance_index.for_each_child(snarl, [&] (const handlegraph::net_handle_t& child) {
        for (bool go_left : {true, false}) {
            #ifdef DEBUG_PATH_PARTITIONER
            cerr << "At snarl child " << distance_index.net_handle_as_string(child) << " going " << (go_left ? "left" : "right") << endl;
            #endif

            std::vector<path_edge_t> next_steps (all_samples.size());
            // For paths that 
            std::vector<path_edge_t> additional_steps;

            // Get a handle to the end of the child
            handlegraph::handle_t handle = distance_index.is_trivial_chain(child) ? (go_left ? graph.flip(distance_index.get_handle(child, &graph)) 
                                                                       : distance_index.get_handle(child, &graph))
                                                            : distance_index.get_handle(distance_index.get_bound(child, go_left, false), &graph);

            std::vector<handlegraph::PathSense> senses = {handlegraph::PathSense::GENERIC,
                                                          handlegraph::PathSense::REFERENCE,
                                                          handlegraph::PathSense::HAPLOTYPE};
            for (const auto& sense : senses) {
                graph.for_each_step_of_sense(handle, sense, [&](const handlegraph::step_handle_t& step) {
                    // For each step on the node handle, keep track of which paths take different steps
                    #ifdef DEBUG_PATH_PARTITIONER
                    cerr << "\ton path " << graph.get_path_name(graph.get_path_handle_of_step(step)) << endl;
                    #endif

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
                    cerr << "\t\tgoing to " << graph.get_id(next_handle) << endl;
                    #endif

                    path_edge_t edge (graph.get_position_of_step(step), 
                                      std::numeric_limits<size_t>::max(),
                                      graph.get_id(next_handle), 
                                      graph.get_is_reverse(next_handle));

                    size_t sample_num = sample_to_index[stoat::get_sample_and_haplotype(graph, graph.get_path_handle_of_step(step))];

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
                        path_edge_t& old_edge = next_steps[sample_num];
                        while (old_edge.additional_edge != std::numeric_limits<size_t>::max() &&
                               additional_steps[old_edge.additional_edge].offset < edge.offset) {
                            size_t next = old_edge.additional_edge; 
                            old_edge = additional_steps[next];
                        }
                        //Old_edge_i now points to the item just smaller than edge
                        size_t old_additional_edge = old_edge.additional_edge;
                        old_edge.additional_edge = additional_steps.size();
                        edge.additional_edge = old_additional_edge;
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
            #ifdef DEBUG_PATH_PARTITIONER
            cerr << "Intermediate sets: " << endl;
            for (size_t i = 0 ; i < all_samples.size() ; i++) {
                cerr << "\t" << all_samples[i] << ": " << intermediate_sets[i] << endl;
            } 
            #endif

            // We now have an old set and an intermediate set for each path
            // Assign the path to a new set. Everything gets a new set
            vector<size_t> new_sets (intermediate_sets.size(), std::numeric_limits<size_t>::max());
            size_t new_set_count = 0;
            // Map pairs of <old_set, intermediate_set> to new set number
            std::map<std::pair<size_t, size_t>, size_t>  old_to_new_set;
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

            #ifdef DEBUG_PATH_PARTITIONER
            cerr << "New sets: " << endl;
            for (size_t i = 0 ; i < all_samples.size() ; i++) {
                cerr << "\t" << all_samples[i] << ": " << old_sets[i] << endl;
            } 
            #endif

        } //end for each direction going out of the node
        return true;
    });// end for_each_child of the snarl

    // We have now partitioned the paths into equivalence sets based on the edges they take in this netgraph,
    // stored in old_sets.
    // Return actual sets of paths.

    std::vector<std::set<stoat::sample_hap_t>> sample_sets (old_set_count);
    for (size_t i = 0 ; i < all_samples.size() ; i++) {
        sample_sets[old_sets[i]].emplace(all_samples[i]);
    }
    #ifdef DEBUG_PATH_PARTITIONER
    cerr << "Found walk sets " << endl;
    for (const auto& s : sample_sets) {
        cerr << "Set" << endl;
        for (const auto& x : s) {
            cerr << "\t" << x << endl;;
        }
    }
    #endif
    return sample_sets;
}

std::vector<std::set<stoat::sample_hap_t>> PathPartitioner::get_start_edge_sets(const handlegraph::PathPositionHandleGraph& graph, 
                                                                         const bdsg::SnarlDistanceIndex& distance_index,
                                                                         const bdsg::net_handle_t& snarl) const {

    #ifdef DEBUG_PATH_PARTITIONER
    cerr << "Get start edge sets of " << distance_index.net_handle_as_string(snarl) << endl;
    #endif

    // Map an edge (as the handle reached from the start of the snarl) to a set of paths that took that edge
    std::map<handlegraph::handle_t, std::set<stoat::sample_hap_t>> edge_to_sample_set;

    // The start node going into the snarl
    handlegraph::handle_t start_node = distance_index.get_handle(distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, true)), &graph);


    std::vector<handlegraph::PathSense> senses = {handlegraph::PathSense::GENERIC,
                                                  handlegraph::PathSense::REFERENCE,
                                                  handlegraph::PathSense::HAPLOTYPE};
    for (const auto& sense : senses) {

        graph.for_each_step_of_sense(start_node, sense, [&](const handlegraph::step_handle_t& step) {

            handlegraph::path_handle_t path = graph.get_path_handle_of_step(step);

            // is the step handle is going in the same direction as the original handle?
            bool go_forward = graph.get_is_reverse(graph.get_handle_of_step(step)) == graph.get_is_reverse(start_node);

            if ((go_forward && graph.has_next_step(step)) ||  (!go_forward && graph.has_previous_step(step))) {
                handlegraph::handle_t next_node = go_forward ? graph.get_handle_of_step(graph.get_next_step(step))
                                                             : graph.get_handle_of_step(graph.get_previous_step(step));
                if (edge_to_sample_set.count(next_node) == 0) {
                    edge_to_sample_set[next_node] = std::set<stoat::sample_hap_t>();
                }
                edge_to_sample_set[next_node].emplace(stoat::get_sample_and_haplotype(graph, path));
            }
        });
    }


    std::vector<std::set<stoat::sample_hap_t>> edge_sets;
    for (auto& edge_to_set : edge_to_sample_set) {
        //TODO: idk if this will mess up the map by moving things inside it
        edge_sets.emplace_back(std::move(edge_to_set.second)); 
    }
    return edge_sets;
}

}
