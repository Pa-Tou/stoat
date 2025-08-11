/**
 * \file register_loader_saver_packed_graph.cpp
 * Defines IO for a bdsg::PackedGraph from stream files.
 * copied from vg
 */
#include <arpa/inet.h>
#include <vg/io/registry.hpp>
#include "register_loader_saver_packed_graph.hpp"

#include "bdsg/packed_graph.hpp"

namespace stoat {

namespace io {

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

}

}

