#include "path_partitioner.hpp"
#include "log.hpp"
#include <fstream>

//#define DEBUG_PATH_PARTITIONER

using namespace std;
using namespace stoat;
namespace stoat_graph {


// This is supposed to partition the paths in the snarl by the walks they take through the netgraph.
// Instead of explicitly enumerating the paths, it actually finds the sets of edges that each path takes.
// But since paths may loop, it actually finds the order and number of outgoing edges from each node.
// I think this is equivalent to partitioning by the actual sets of unique walks.
std::vector<size_t> partition_embedded_paths_in_snarl(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                          const net_handle_t& snarl,
                          const std::vector<stoat::sample_hap_t>& all_sample_haplotypes){

    #ifdef DEBUG_PATH_PARTITIONER
    cerr <<  "Get walk sets of " << distance_index.net_handle_as_string(snarl) << endl;;
    #endif

    // Map each sample(plus haplotype) to its index in all_sample_haplotypes
    std::map<stoat::sample_hap_t, std::size_t> sample_to_index;
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
        cerr << "At snarl child " << distance_index.net_handle_as_string(child) << " going " << (go_left ? "left" : "right") << endl;
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
        cerr << "\tgraph handle " << graph.get_id(handle) << " going " << (graph.get_is_reverse(handle) ? "left" : "right") << endl;
        #endif

        for (const auto& sense : senses) {
            graph.for_each_step_of_sense(handle, sense, [&](const handlegraph::step_handle_t& step) {
                // For each step on the node handle, keep track of which paths take different steps

                #ifdef DEBUG_PATH_PARTITIONER
                cerr << "\ton path " << graph.get_path_name(graph.get_path_handle_of_step(step)) << endl;
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
                cerr << "" << "\t\tgoing to " << graph.get_id(next_handle) << endl;
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
        cerr << "Intermediate sets: " << endl;
        for (size_t i = 0 ; i < all_sample_haplotypes.size() ; i++) {
            cerr << "" << "\t" << all_sample_haplotypes[i] << ": " << intermediate_sets[i] << endl;
        } 
        assert(additional_edge_count == additional_steps.size());
        #endif
        
        // We now have an old set and an intermediate set for each path
        // Assign the path to a new set. Everything gets a new set
        vector<size_t> new_sets (intermediate_sets.size(), std::numeric_limits<size_t>::max());

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
        cerr << "New sets: " << endl;
        for (size_t i = 0 ; i < all_sample_haplotypes.size() ; i++) {
            cerr << "" << "\t" << all_sample_haplotypes[i] << ": " << old_sets[i] << endl;
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



    #ifdef DEBUG_PATH_PARTITIONER,
    std::vector<std::set<sample_hap_t>> sample_sets_by_allele(old_set_count-1);
    for (size_t i = 0 ; i < all_sample_haplotypes.size() ; i++) {
        if (old_sets[i] != 0) {
            sample_sets_by_allele[old_sets[i]-1].emplace(all_sample_haplotypes[i]);
        }
    }
    cerr << "Found walk sets " << endl;
    for (const auto& s : sample_sets_by_allele) {
        cerr <<  "Set" << endl;
        for (const auto& x : s) {
            cerr << "" << "\t" << x << endl;
        }
    }
    #endif

    // old_sets is now a vector that assigns each sample/haplotype in all_sample_haplotypes to a allele number
    // Return this
    return old_sets;
}

}
