#include <iostream>
#include <unordered_set>
#include <string>

#include <handlegraph/path_position_handle_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include <bdsg/overlays/overlay_helper.hpp>
#include <arpa/inet.h>
#include <vg/io/registry.hpp>
#include <vg/io/vpkg.hpp>

#include "bdsg/packed_graph.hpp"

using namespace std;
using namespace handlegraph;
using namespace bdsg;


using namespace vg::io;
using namespace handlegraph;

void register_loader_saver_packed_graph() {

    // Convert the bdsg::PackedGraph SerializableHandleGraph magic number to a std::string
    bdsg::PackedGraph empty;
    // Make sure it is in network byte order
    uint32_t new_magic_number = htonl(empty.get_magic_number());
    // Load all 4 characters of it into a std::string
    std::string new_magic((char*)&new_magic_number, 4);
    
    Registry::register_bare_loader_saver_with_magic<bdsg::PackedGraph, MutablePathDeletableHandleGraph, MutablePathMutableHandleGraph, MutableHandleGraph, PathHandleGraph, HandleGraph>("PackedGraph", new_magic, [](istream& input) -> void* {
        // Allocate a bdsg::PackedGraph
         bdsg::PackedGraph* packed_graph = new bdsg::PackedGraph();
        
        // Load it
        packed_graph->deserialize(input);
        
        // Return it so the caller owns it.
        return (void*) packed_graph;
    }, [](const void* packed_graph_void, ostream& output) {
        // Cast to bdsg::PackedGraph and serialize to the stream.
        assert(packed_graph_void != nullptr);
        ((const bdsg::PackedGraph*) packed_graph_void)->serialize(output);
    });
}

int main(int argc, char* argv[]) {

    // Load the graph and make it a PathPositionHandleGraph
    std::string graph_path = argv[1];
    std::string distance_index_path = argv[2];

    register_loader_saver_packed_graph();
    std::unique_ptr<handlegraph::PathHandleGraph> graph = std::move(vg::io::VPKG::load_one<handlegraph::PathHandleGraph>(graph_path));
    bdsg::PathPositionOverlayHelper overlay_helper;
    bdsg::PathPositionHandleGraph* path_position_graph;
    path_position_graph = overlay_helper.apply(graph.get());

    std::unique_ptr<bdsg::SnarlDistanceIndex> distance_index = std::make_unique<bdsg::SnarlDistanceIndex>();
    distance_index->deserialize(distance_index_path);

    cout << "node_id" << "\t"
         << "pos" << "\t"
         << "ref_name" << "\t"
         << "depth" << "\n";

    path_position_graph->for_each_handle([&](const handle_t& handle) {

        nid_t node_id = path_position_graph->get_id(handle);

        // Compute depth
        // net_handle_t net = path_position_graph->get_net(handle);
        net_handle_t net = distance_index->get_net(handle, path_position_graph);
        size_t depth = distance_index->get_depth(net) - 1;

        // Iterate over each occurrence of the node
        path_position_graph->for_each_step_on_handle(handle, [&](const step_handle_t& step) {

            path_handle_t path = path_position_graph->get_path_handle_of_step(step);
            string ref_name = path_position_graph->get_path_name(path);
            size_t pos = path_position_graph->get_position_of_step(step);

            cout << node_id << "\t"
                 << pos << "\t"
                 << ref_name << "\t"
                 << depth << "\n";

            return false;
        });

    });

    return 0;
}

// g++ -std=c++17 node_parsing.cpp -o node_parsing \
  -I../../deps/libvgio/include \
  -L../../build/deps/libvgio \
  -Xpreprocessor -fopenmp \
  -I/usr/local/opt/libomp/include \
  -L/usr/local/opt/libomp/lib -lomp \
  -lvgio -lbdsg -lhandlegraph -lsdsl -lz -lpthread

// ./node_parsing pg.full.pg pg.full.dist > node_parsing_output.txt
