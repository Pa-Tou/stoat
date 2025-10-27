import bdsg
import re

#Given a handle_t to a node, return a set of the names of paths that traverse that node
def get_paths_on_node(graph, node, ignore_reference, clean_path_name = lambda x : x):
    #Helper function that returns the name of the path of a step
    path_names = set()
    def get_path_name(step):
        path_handle = graph.get_path_handle_of_step(step)
        if not (ignore_reference and (graph.get_sense(path_handle) == bdsg.handlegraph.PathSense.REFERENCE)):
            path_names.add(clean_path_name(graph.get_path_name(path_handle)))
        return True
    graph.for_each_step_on_handle(node, get_path_name)
    return path_names

#Given a handle_t to a node, return a set of the names of paths that traverse that node
def get_paths_and_offsets_on_node(graph, node, ignore_reference, clean_path_name = lambda x : x):
    #Helper function that returns the name of the path of a step
    position_overlay=bdsg.bdsg.PositionOverlay(graph)
    path_names = set()
    def get_path_name(step):
        path_handle = graph.get_path_handle_of_step(step)
        if not (ignore_reference and (graph.get_sense(path_handle) == bdsg.handlegraph.PathSense.REFERENCE)):
            path_names.add((clean_path_name(graph.get_path_name(path_handle)), position_overlay.get_position_of_step(step)))
        return True
    graph.for_each_step_on_handle(node, get_path_name)
    return path_names



def get_sample_from_path_handle(graph, path_handle):
    if graph.get_sense(path_handle) == bdsg.handlegraph.PathSense.GENERIC:
        # Generic paths only have a locus, so return whatever that is
        return graph.get_locus_name(path_handle)
    else:
        return graph.get_sample_name(path_handle)


#Given a handle_t to a node, return a set of the samples that traverse that node
def get_samples_on_node(graph, node, ignore_reference):
    #Helper function that returns the name of the path of a step
    sample_names = set()
    def get_sample_name(step):
        path_handle = graph.get_path_handle_of_step(step)
        if not (ignore_reference and (graph.get_sense(path_handle) == bdsg.handlegraph.PathSense.REFERENCE)):
            sample_names.add(get_sample_from_path_handle(graph, path_handle))
        return True
    graph.for_each_step_on_handle(node, get_sample_name)
    return sample_names

def get_reference_range_of_snarl(graph, position_graph, distance_index, snarl):
    start = distance_index.get_handle(distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, False, True)), graph)
    end = distance_index.get_handle(distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, True, True)), graph)

    path_name = ""
    start_offset = 0
    end_offset = 0

    def get_path_offset_start(step):
        nonlocal path_name
        nonlocal start_offset
        path = graph.get_path_handle_of_step(step)
        if graph.get_sense(path) == bdsg.handlegraph.PathSense.REFERENCE:
            start_offset = position_graph.get_position_of_step(step)
            path_name = graph.get_path_name(path)
        return True

    def get_path_offset_end(step):
        nonlocal end_offset
        path = graph.get_path_handle_of_step(step)
        if graph.get_sense(path) == bdsg.handlegraph.PathSense.REFERENCE:
            end_offset = position_graph.get_position_of_step(step)
            path_name = graph.get_path_name(path)
        return True
        

    graph.for_each_step_on_handle(start, get_path_offset_start)
    graph.for_each_step_on_handle(end, get_path_offset_end)

    return (path_name, min(start_offset, end_offset), max(start_offset, end_offset))

#Given a handle_t to a node, return a dictionary of adjacent nodes, following the handle in its forward orientation
# mapping node id to a dictionary of paths that took the edge, mapped to a count of times the path traversed the edge
# TODO: This should probably be mapping the handle but it doesn't work in python
def get_next_nodes_and_paths(graph, node, ignore_reference, clean_path_name = lambda x : x):
    next_nodes = dict()

    def add_next_step(step):
        path_handle = graph.get_path_handle_of_step(step)
        if not (ignore_reference and (graph.get_sense(path_handle) == bdsg.handlegraph.PathSense.REFERENCE)):
            if graph.get_is_reverse(graph.get_handle_of_step(step)) == graph.get_is_reverse(node): 
                # if the step handle is going in the same direction as the original handle
                if graph.has_next_step(step):
                    next_node = graph.get_handle_of_step(graph.get_next_step(step))
                    if graph.get_id(next_node) not in next_nodes:
                        next_nodes[graph.get_id(next_node)] = set()
                    path_name = clean_path_name(graph.get_path_name(path_handle))
                    next_nodes[graph.get_id(next_node)].add(path_name)
            else: 
                # Otherwise, the step handle got flipped so we go backwards
                if graph.has_previous_step(step):
                    next_node = graph.get_handle_of_step(graph.get_previous_step(step))
                    if graph.get_id(next_node) not in next_nodes:
                        next_nodes[graph.get_id(next_node)] = set()
                    path_name = clean_path_name(graph.get_path_name(path_handle))
                    next_nodes[graph.get_id(next_node)].add(path_name) 


        return True

    graph.for_each_step_on_handle(node, add_next_step)

    return next_nodes

def get_next_nodes_and_samples(graph, node, ignore_reference):
    next_nodes = dict()

    def add_next_step(step):
        path_handle = graph.get_path_handle_of_step(step)
        if not (ignore_reference and (graph.get_sense(path_handle) == bdsg.handlegraph.PathSense.REFERENCE)):
            if graph.get_is_reverse(graph.get_handle_of_step(step)) == graph.get_is_reverse(node): 
                # if the step handle is going in the same direction as the original handle
                if graph.has_next_step(step):
                    next_node = graph.get_handle_of_step(graph.get_next_step(step))
                    if graph.get_id(next_node) not in next_nodes:
                        next_nodes[graph.get_id(next_node)] = set()
                    next_nodes[graph.get_id(next_node)].add(get_sample_from_path_handle(graph, path_handle))
            else: 
                # Otherwise, the step handle got flipped so we go backwards
                if graph.has_previous_step(step):
                    next_node = graph.get_handle_of_step(graph.get_previous_step(step))
                    if graph.get_id(next_node) not in next_nodes:
                        next_nodes[graph.get_id(next_node)] = set()
                    next_nodes[graph.get_id(next_node)].add(get_sample_from_path_handle(graph, path_handle)) 


        return True

    graph.for_each_step_on_handle(node, add_next_step)

    return next_nodes


# Given a child in a chain, get the next child in the forward direction of the net (unless go_left is true)
def get_next_net_in_chain(graph, distance_index, net, go_left=False):
    next_net = distance_index.get_root()
    def set_next(next_net_handle):
        nonlocal next_net
        next_net = next_net_handle
        return True
    distance_index.follow_net_edges(net, graph, go_left, set_next) 
    return next_net

# Given a child in a snarl, get the next children in the forward direction of the net (unless go_left is true)
def get_next_nets_in_snarl(graph, distance_index, net, go_left=False):
    next_chains = []
    def set_next(next_net_handle):
        nonlocal next_chains
        next_chains.apppend(next_net_handle)
        return True
    distance_index.follow_net_edges(net, graph, go_left, set_next)
    return next_chains

# Return a list of children of a net handle
def get_children_of_net(distance_index, net):
    children = []
    def get_child(child):
        nonlocal children
        children.append(child)
        return True
    distance_index.for_each_child(net, get_child)
    return children

# Count the children of a net
def get_child_count(distance_index, net):    
    if distance_index.is_node(net):
        return 0

    count = 0

    def increment_count(n):
        nonlocal count
        count += 1
        return True

    distance_index.for_each_child(net, increment_count)
    return count

def print_ancestors(distance_index, net):
    while not distance_index.is_root(net):
        print(distance_index.net_handle_as_string(net))
    print(distance_index.net_handle_as_string(net))

def print_ancestors_and_child_count(distance_index, net):
    while not distance_index.is_root(net):
        print(distance_index.net_handle_as_string(net), get_child_count(distance_index, net))
        net = distance_index.get_parent(net)
    print(distance_index.net_handle_as_string(net), get_child_count(distance_index, net))

def print_ancestors_and_lengths(distance_index, net):
    while not distance_index.is_root(net):
        print(distance_index.net_handle_as_string(net), distance_index.minimum_length(net))
        net = distance_index.get_parent(net)

def print_node_descendants(distance_index, net):
    to_print = [net]
    while len(to_print) != 0:
        next_net = to_print.pop()
        if distance_index.is_node(next_net):
            print(distance_index.node_id(next_net))
            continue
        def add_net(n):
            nonlocal to_print
            to_print.append(n)
            return True
        distance_index.for_each_child(next_net, add_net)
    if distance_index.is_snarl(net):
        print(distance_index.node_id(distance_index.get_node_from_sentinel(distance_index.get_bound(net, False, True))))
        print(distance_index.node_id(distance_index.get_node_from_sentinel(distance_index.get_bound(net, True, True))))

def print_children_and_lengths(distance_index, net):
    def print_child(n):
        print(distance_index.net_handle_as_string(n), distance_index.minimum_length(n))
        return True
    distance_index.for_each_child(net, print_child)

def print_snarl_edges(graph, distance_index, snarl):
    def print_edge_to(n):
        print("\t->"+distance_index.net_handle_as_string(n))
        return True
    def print_outgoing_edges(n):
        print(distance_index.net_handle_as_string(n))
        distance_index.follow_net_edges(n, graph, False, print_edge_to)

        print(distance_index.net_handle_as_string(distance_index.flip(n)))
        distance_index.follow_net_edges(n, graph, True, print_edge_to)
        return True

    distance_index.for_each_child(snarl, print_outgoing_edges)
