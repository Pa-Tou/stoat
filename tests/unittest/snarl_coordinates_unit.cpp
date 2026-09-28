#include <catch.hpp>
#include <bdsg/overlays/overlay_helper.hpp>
#include <bdsg/hash_graph.hpp>

#include "../../src/snarl_coordinates.hpp"

using namespace stoat; 


TEST_CASE("Snarl coordinates simple nested chain", "[snarl_coords]") {

    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/test_data/test_graphs/simple_nested_chain.dist");
    
    bdsg::HashGraph hash_graph;
    hash_graph.deserialize("../tests/test_data/test_graphs/simple_nested_chain.hg");
    bdsg::PathPositionOverlayHelper overlay_helper;
    auto graph = overlay_helper.apply(&hash_graph);



    handlegraph::net_handle_t snarl1 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(2)));
    handlegraph::net_handle_t snarl2 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(5)));
    handlegraph::net_handle_t snarl3 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(6)));
    handlegraph::net_handle_t snarl4 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(9)));

    SECTION("reference path") {
        SnarlCoordinates snarl_coordinate_finder;
        std::tuple<std::string, size_t, size_t> snarl1_coords = snarl_coordinate_finder.get_reference_coordinates_as_string(*graph, distance_index, snarl1);
        REQUIRE(std::get<0>(snarl1_coords) == "path0#0#path0");
        REQUIRE(std::get<1>(snarl1_coords) == 1);
        REQUIRE(std::get<2>(snarl1_coords) == 2);

        std::tuple<std::string, size_t, size_t> snarl3_coords = snarl_coordinate_finder.get_reference_coordinates_as_string(*graph, distance_index,  snarl3);
        REQUIRE(std::get<0>(snarl3_coords) == "path0#0#path0");
        REQUIRE(std::get<1>(snarl3_coords) == 4);
        REQUIRE(std::get<2>(snarl3_coords) == 5);
    }
    SECTION("named path") {
        SnarlCoordinates snarl_coordinate_finder;
        snarl_coordinate_finder.add_reference_path("path2");

        std::tuple<std::string, size_t, size_t> snarl1_coords = snarl_coordinate_finder.get_reference_coordinates_as_string(*graph, distance_index, snarl1);
        REQUIRE(std::get<0>(snarl1_coords) == "path2");
        REQUIRE(std::get<1>(snarl1_coords) == 1);
        REQUIRE(std::get<2>(snarl1_coords) == 2);

        // Nested snarl off the reference
        std::tuple<std::string, size_t, size_t> snarl3_coords = snarl_coordinate_finder.get_reference_coordinates_as_string(*graph, distance_index, snarl3);
        REQUIRE(std::get<0>(snarl3_coords) == "path2");
        REQUIRE(std::get<1>(snarl3_coords) == 3);
        REQUIRE(std::get<2>(snarl3_coords) == 3);
    }

}
TEST_CASE("Snarl coordinates deeply nested snarls with deletion in reference", "[snarl_coords]") {
    /*
                      3
                    /   \
                   2 ----4  
                 /         \
               1  ----------5
              /              \
            0 ----------------6

   */

    bdsg::HashGraph hash_graph;

    std::vector<std::string> sequences = { "C", "C", "C", "A", "T", "C", "A"};

    std::vector<handlegraph::handle_t> nodes;
    for (auto& seq : sequences) {
        nodes.emplace_back(hash_graph.create_handle(seq));
    }

    hash_graph.create_edge(nodes[0], nodes[1]);
    hash_graph.create_edge(nodes[0], nodes[6]);
    hash_graph.create_edge(nodes[1], nodes[2]);
    hash_graph.create_edge(nodes[1], nodes[5]);
    hash_graph.create_edge(nodes[2], nodes[3]);
    hash_graph.create_edge(nodes[2], nodes[4]);
    hash_graph.create_edge(nodes[3], nodes[4]);
    hash_graph.create_edge(nodes[4], nodes[5]);
    hash_graph.create_edge(nodes[5], nodes[6]);

    std::vector<std::vector<std::size_t>> paths_seqs = { {0, 6}, {0, 1, 5, 6}, {0, 1, 2, 3, 4, 5, 6}};
    std::vector<handlegraph::path_handle_t> paths;

    for (int path_i = 0 ; path_i < paths_seqs.size() ; path_i++) {
        if (path_i == 0) {
            // Set first path with deletion as reference
            paths.emplace_back(hash_graph.create_path_handle("path"+std::to_string(path_i)+"#0#0"));
        } else {
            paths.emplace_back(hash_graph.create_path_handle("path"+std::to_string(path_i)+"#0#0#0"));
        }
        for (size_t node_i : paths_seqs[path_i]) {
            hash_graph.append_step(paths.back(), nodes[node_i]);
        }
    }

    //// vg isn't included so the distance index can only be built from the command line
    //hash_graph.serialize("../tests/test_data/test_graphs/deeply_nested_snarl.hg");
    //int built = system("vg index -j ../tests/test_data/test_graphs/deeply_nested_snarl.dist ../tests/test_data/test_graphs/deeply_nested_snarl.hg"); 

    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/test_data/test_graphs/deeply_nested_snarl.dist");
    
    //bdsg::HashGraph hash_graph;
    //hash_graph.deserialize("../tests/test_data/test_graphs/deeply_nested_snarl.hg");

    bdsg::PathPositionOverlayHelper overlay_helper;
    auto graph = overlay_helper.apply(&hash_graph);

    handlegraph::net_handle_t snarl1 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(6)));
    handlegraph::net_handle_t snarl2 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(5)));
    handlegraph::net_handle_t snarl3 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(4)));

    SECTION("reference with deletion") {
        SnarlCoordinates snarl_coordinate_finder;
        snarl_coordinate_finder.add_reference_path("path0#0#0");

        std::tuple<std::string, size_t, size_t> snarl1_coords = snarl_coordinate_finder.get_reference_coordinates_as_string(*graph, distance_index, snarl1);
        REQUIRE(std::get<0>(snarl1_coords) == "path0#0#0");
        REQUIRE(std::get<1>(snarl1_coords) == 1);
        REQUIRE(std::get<2>(snarl1_coords) == 1);

        // Nested snarl off the reference
        std::tuple<std::string, size_t, size_t> snarl3_coords = snarl_coordinate_finder.get_reference_coordinates_as_string(*graph, distance_index, snarl3);
        REQUIRE(std::get<0>(snarl3_coords) == "path0#0#0");
        REQUIRE(std::get<1>(snarl3_coords) == 1);
        REQUIRE(std::get<2>(snarl3_coords) == 1);
    }
}

TEST_CASE("Snarl coordinates deeply nested snarls with no reference sense", "[snarl_coords]") {
    /*
                      3
                    /   \
                   2 ----4  
                 /         \
               1  ----------5
              /              \
            0 ----------------6

   */

    bdsg::HashGraph hash_graph;

    std::vector<std::string> sequences = { "C", "C", "C", "A", "T", "C", "A"};

    std::vector<handlegraph::handle_t> nodes;
    for (auto& seq : sequences) {
        nodes.emplace_back(hash_graph.create_handle(seq));
    }

    hash_graph.create_edge(nodes[0], nodes[1]);
    hash_graph.create_edge(nodes[0], nodes[6]);
    hash_graph.create_edge(nodes[1], nodes[2]);
    hash_graph.create_edge(nodes[1], nodes[5]);
    hash_graph.create_edge(nodes[2], nodes[3]);
    hash_graph.create_edge(nodes[2], nodes[4]);
    hash_graph.create_edge(nodes[3], nodes[4]);
    hash_graph.create_edge(nodes[4], nodes[5]);
    hash_graph.create_edge(nodes[5], nodes[6]);

    std::vector<std::vector<std::size_t>> paths_seqs = { {0, 6}, {0, 1, 5, 6}, {0, 1, 2, 3, 4, 5, 6}};
    std::vector<handlegraph::path_handle_t> paths;

    for (int path_i = 0 ; path_i < paths_seqs.size() ; path_i++) {
        paths.emplace_back(hash_graph.create_path_handle("path"+std::to_string(path_i)+"#0#0#0"));
        for (size_t node_i : paths_seqs[path_i]) {
            hash_graph.append_step(paths.back(), nodes[node_i]);
        }
    }

    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/test_data/test_graphs/deeply_nested_snarl.dist");


    bdsg::PathPositionOverlayHelper overlay_helper;
    auto graph = overlay_helper.apply(&hash_graph);

    handlegraph::net_handle_t snarl1 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(6)));
    handlegraph::net_handle_t snarl2 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(5)));
    handlegraph::net_handle_t snarl3 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(4)));

    SECTION("reference with deletion") {

        SnarlCoordinates snarl_coordinate_finder;
        snarl_coordinate_finder.add_reference_path("path0#0#0#0");

        std::tuple<std::string, size_t, size_t> snarl1_coords = snarl_coordinate_finder.get_reference_coordinates_as_string(*graph, distance_index, snarl1);
        REQUIRE(std::get<0>(snarl1_coords) == "path0#0#0#0");
        REQUIRE(std::get<1>(snarl1_coords) == 1);
        REQUIRE(std::get<2>(snarl1_coords) == 1);

        // Nested snarl off the reference
        std::tuple<std::string, size_t, size_t> snarl3_coords = snarl_coordinate_finder.get_reference_coordinates_as_string(*graph, distance_index, snarl3);
        REQUIRE(std::get<0>(snarl3_coords) == "path0#0#0#0");
        REQUIRE(std::get<1>(snarl3_coords) == 1);
        REQUIRE(std::get<2>(snarl3_coords) == 1);
    }
}

TEST_CASE("Snarl coordinates on given path instead of reference", "[snarl_coords]") {
    /*
                      3
                    /   \
                   2 ----4  
                 /         \
               1  ----------5
              /              \
            0 ----------------6

   */

    bdsg::HashGraph hash_graph;

    std::vector<std::string> sequences = { "C", "C", "C", "A", "T", "C", "A"};

    std::vector<handlegraph::handle_t> nodes;
    for (auto& seq : sequences) {
        nodes.emplace_back(hash_graph.create_handle(seq));
    }

    hash_graph.create_edge(nodes[0], nodes[1]);
    hash_graph.create_edge(nodes[0], nodes[6]);
    hash_graph.create_edge(nodes[1], nodes[2]);
    hash_graph.create_edge(nodes[1], nodes[5]);
    hash_graph.create_edge(nodes[2], nodes[3]);
    hash_graph.create_edge(nodes[2], nodes[4]);
    hash_graph.create_edge(nodes[3], nodes[4]);
    hash_graph.create_edge(nodes[4], nodes[5]);
    hash_graph.create_edge(nodes[5], nodes[6]);

    std::vector<std::vector<std::size_t>> paths_seqs = {{0, 1, 5, 6}, {0, 6}, {0, 1, 2, 3, 4, 5, 6}, {}};
    std::vector<handlegraph::path_handle_t> paths;

    for (int path_i = 0 ; path_i < paths_seqs.size() ; path_i++) {
        if (path_i == 0) {
            // Set first path with deletion as reference
            paths.emplace_back(hash_graph.create_path_handle("path"+std::to_string(path_i)+"#0#0"));
        } else {
            paths.emplace_back(hash_graph.create_path_handle("path"+std::to_string(path_i)+"#0#0#0"));
        }
        for (size_t node_i : paths_seqs[path_i]) {
            hash_graph.append_step(paths.back(), nodes[node_i]);
        }
    }

    //// vg isn't included so the distance index can only be built from the command line
    //hash_graph.serialize("../tests/test_data/test_graphs/deeply_nested_snarl.hg");
    //int built = system("vg index -j ../tests/test_data/test_graphs/deeply_nested_snarl.dist ../tests/test_data/test_graphs/deeply_nested_snarl.hg"); 

    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/test_data/test_graphs/deeply_nested_snarl.dist");
    
    //bdsg::HashGraph hash_graph;
    //hash_graph.deserialize("../tests/test_data/test_graphs/deeply_nested_snarl.hg");

    bdsg::PathPositionOverlayHelper overlay_helper;
    auto graph = overlay_helper.apply(&hash_graph);

    handlegraph::net_handle_t snarl1 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(6)));
    handlegraph::net_handle_t snarl2 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(5)));
    handlegraph::net_handle_t snarl3 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(4)));

    SECTION("given sample with deletion") {
        SnarlCoordinates snarl_coordinate_finder;
        snarl_coordinate_finder.add_reference_path("path1#0#0#0");

        std::tuple<std::string, size_t, size_t> snarl1_coords = snarl_coordinate_finder.get_reference_coordinates_as_string(*graph, distance_index, snarl1);
        REQUIRE(std::get<0>(snarl1_coords) == "path1#0#0#0");
        REQUIRE(std::get<1>(snarl1_coords) == 1);
        REQUIRE(std::get<2>(snarl1_coords) == 1);

        // Nested snarl off the reference
        std::tuple<std::string, size_t, size_t> snarl3_coords = snarl_coordinate_finder.get_reference_coordinates_as_string(*graph, distance_index, snarl3);
        REQUIRE(std::get<0>(snarl3_coords) == "path1#0#0#0");
        REQUIRE(std::get<1>(snarl3_coords) == 1);
        REQUIRE(std::get<2>(snarl3_coords) == 1);
    }
    SECTION("reference") {
        SnarlCoordinates snarl_coordinate_finder;

        std::tuple<std::string, size_t, size_t> snarl1_coords = snarl_coordinate_finder.get_reference_coordinates_as_string(*graph, distance_index, snarl1);
        REQUIRE(std::get<0>(snarl1_coords) == "path0#0#0");
        REQUIRE(std::get<1>(snarl1_coords) == 1);
        REQUIRE(std::get<2>(snarl1_coords) == 3);

        // Nested snarl off the reference
        std::tuple<std::string, size_t, size_t> snarl3_coords = snarl_coordinate_finder.get_reference_coordinates_as_string(*graph, distance_index, snarl3);
        REQUIRE(std::get<0>(snarl3_coords) == "path0#0#0");
        REQUIRE(std::get<1>(snarl3_coords) == 2);
        REQUIRE(std::get<2>(snarl3_coords) == 2);
    }
    SECTION("given sample not on snarl") {
        SnarlCoordinates snarl_coordinate_finder;
        snarl_coordinate_finder.add_reference_path("path3#0#0#0");

        std::tuple<std::string, size_t, size_t> snarl1_coords = snarl_coordinate_finder.get_reference_coordinates_as_string(*graph, distance_index, snarl1);
        REQUIRE(std::get<0>(snarl1_coords) == "path0#0#0");
        REQUIRE(std::get<1>(snarl1_coords) == 1);
        REQUIRE(std::get<2>(snarl1_coords) == 3);

        // Nested snarl off the reference
        std::tuple<std::string, size_t, size_t> snarl3_coords = snarl_coordinate_finder.get_reference_coordinates_as_string(*graph, distance_index, snarl3);
        REQUIRE(std::get<0>(snarl3_coords) == "path0#0#0");
        REQUIRE(std::get<1>(snarl3_coords) == 2);
        REQUIRE(std::get<2>(snarl3_coords) == 2);
    }
}
