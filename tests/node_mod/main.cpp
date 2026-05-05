#include <bdsg/packed_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include <bdsg/overlays/packed_path_position_overlay.hpp>

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>

#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <vector>

using namespace std;
using namespace bdsg;
using namespace handlegraph;

int main() {

    // Load graph
    auto graph = std::make_unique<PackedGraph>();

    ifstream graph_in("pg.full.pg");
    graph->deserialize(graph_in);

    // Load snarl distance index
    auto sdi = std::make_unique<SnarlDistanceIndex>();

    ifstream sdi_in("pg.full.dist");
    sdi->deserialize(sdi_in);

    // Overlay to compute path positions efficiently
    PackedPositionOverlay position_overlay(graph.get());

    set<string> set_reference_name = {"ref"};
    vector<string> vector_reference;

    // Iterate nodes
    graph->for_each_handle([&](const handle_t& h) {

        nid_t nid_t_id = graph->get_id(h);
        size_t node_id = nid_t_id;
        bool is_reverse = graph->get_is_reverse(h);
        string reverse = is_reverse ? "<" : ">";
        string node_str = reverse + to_string(node_id);
        std::string node_sequence = graph->get_sequence(h);
        // size_t depth = sdi->get_depth(h);
        size_t pos;
        string reference_name;

        cout << "Node: " << node_str << ", node_sequence: " << node_sequence;

        graph->for_each_step_on_handle(h, [&](const step_handle_t& step) {
            reference_name = graph->get_path_name(graph->get_path_handle_of_step(step));
            pos = position_overlay.get_position_of_step(step);
            if (set_reference_name.find(reference_name) != set_reference_name.end()) {
                return false;
            } else {
                reference_name = "NA";
                return true;
            }
        });

        cout << ", reference_name: " << reference_name << ", position: " << pos << endl;
    });
}

// g++ -std=c++17 -O2 main.cpp -o graph_analyse -lbdsg -lhandlegraph -lsdsl -lz -lpthread



