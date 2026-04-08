import bdsg
import random
import argparse


# For distance indexing, it helps to have nodes that are on the top-level chain upweighted to prevent weird loopy chain backbones
# This returns a list of nodes that all haplotypes in the graph pass through
# TODO: idk how this will deal with fragmented paths
# TODO: I think this only works with hg
def find_nodes_with_all_haplotypes(graph, node_count):


    # First find all the paths in the graph
    path_names = set()
    def add_path_name(path_handle):
        path_names.add(graph.get_path_name(path_handle))
        return True
    graph.for_each_path_handle(add_path_name)
    path_count = len(path_names)

    nodes = []

    def check_handle(node_handle):

        return True
            
 
    # Since it would be super slow to check all nodes, we're going to pick random nodes from the graph
    # and check them until we've got 20, or until we've checked 20% of the nodes in the graph
    min_node_id = graph.min_node_id()
    max_node_id = graph.max_node_id()        
    i = 0
    while i < ((max_node_id-min_node_id)*0.2) and len(nodes) < node_count:
        node_id = random.randrange(min_node_id, max_node_id)

        if graph.has_node(node_id): 
            i+=1
            node_handle = graph.get_handle(node_id)

            # For each node handle, get all paths on the handle
            node_path_names = set()
            def add_node_path_names(step_handle):
                # For each step on the node handle, get the path
                node_path_names.add(graph.get_path_name(graph.get_path_handle_of_step(step_handle)))
                return True

            graph.for_each_step_on_handle(node_handle, add_node_path_names)

            # Now that we have all paths on this node, check if we have all paths we wanted
            if len(node_path_names) == path_count:
                nodes.append(graph.get_id(node_handle))
    return nodes

# using the nodes from find_nodes_with_all_haplotypes, print the -w options for the vg index -j command
def print_distance_index_node_command(graph, node_count):
    nodes = find_nodes_with_all_haplotypes(graph, node_count)
    string_list = ["-w " + str(x) for x in nodes] 
    print (" ".join(string_list))


def main():
    parser = argparse.ArgumentParser(description="For snarl finding, it can sometimes help to give the decomposition finder hints about which nodes should be on the top-level chain backbone. This randomly finds and prints nodes through which all paths in the graph pass. \n \
    This can be used with vg index: \n \
    vg index -j [graph.dist] $(python3 get_snarl_backbone_nodes.py [graph.(pg|hg)]) [graph.(pg|hg)]")
    parser.add_argument('graph_name', help="the graph. The file extension must be .hg or .pg for a HashGraph or PackedGraph respectively.")
    parser.add_argument('--node_count', '-n', help="How many node ids do we want? ", default=100, type=int) 
    args = parser.parse_args()

    graph = bdsg.bdsg.HashGraph()
    if args.graph_name.endswith(".pg"):
        graph = bdsg.bdsg.PackedGraph()
    elif not args.graph_name.endswith(".hg"):
        raise OSError("Invalid graph file: "+ args.graph_name)
        print()

    graph.deserialize(args.graph_name)

    print_distance_index_node_command(graph, args.node_count)

if __name__ == "__main__":
    main()
