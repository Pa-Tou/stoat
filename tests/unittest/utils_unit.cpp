#include <catch.hpp>
#include <bdsg/overlays/overlay_helper.hpp>
#include <bdsg/hash_graph.hpp>

#include "../../src/types_and_structs.hpp"
#include "../../src/utils.hpp"

using namespace stoat; 

using boost::multiprecision::cpp_dec_float_50;

TEST_CASE("set_precision handles small and normal values") {
    REQUIRE(set_precision(0.00001234) == "1.2340e-05");
    REQUIRE(set_precision(0.123456) == "0.1235");
}

TEST_CASE("set_precision handles small and normal values 2") {
    REQUIRE(set_precision(0.00001234567890123456789) == "1.2346e-05");
    REQUIRE(set_precision(0.34567890123456789) == "0.3457");
}

TEST_CASE("set_precision handles small and normal values 3") {
    REQUIRE(set_precision(0.333333333) == "0.3333");
    REQUIRE(set_precision(1.0000) == "1");
    REQUIRE(set_precision(1.000000000) == "1");
}

TEST_CASE("string_to_pvalue converts valid p-values or returns 1.0 for NA") {
    REQUIRE(stoat::string_to_pvalue("0.01") == 0.01);
    REQUIRE(stoat::string_to_pvalue("NA") == 1.0);
    REQUIRE(stoat::string_to_pvalue("") == 1.0);
}

TEST_CASE("remove_prefix cases") {

    REQUIRE(remove_prefix("", "foo") == "");
    REQUIRE(remove_prefix("foo", "") == "foo");
    REQUIRE(remove_prefix("$$#$value", "$$") == "#$value");

    std::string prefix = "prefix_";
    std::string long_str = prefix + std::string(1000, 'a');
    REQUIRE(remove_prefix(long_str, prefix) == std::string(1000, 'a'));
}

TEST_CASE("adjusted_hochberg basic example") {
    std::vector<double> raw = {0.02, 0.15, 0.03, 0.001, 0.25, 0.05};
    std::pair<double, size_t> expected = {0.006, 3};

    auto [pvalue, index] = adjusted_hochberg(raw);
    REQUIRE(pvalue == expected.first);
    REQUIRE(index == expected.second);
}

TEST_CASE("adjusted_hochberg with descending input") {
    std::vector<double> raw = {0.5, 0.04, 0.03, 0.02, 0.001};
    std::pair<double, size_t> expected = {0.005, 4};

    auto [pvalue, index] = adjusted_hochberg(raw);
    REQUIRE(pvalue == expected.first);
    REQUIRE(index == expected.second);
}

TEST_CASE("adjusted_hochberg with one low p-value and others near 0.1") {
    std::vector<double> raw = {0.000001, 0.09, 0.08, 0.07};
    std::pair<double, size_t> expected = {0.000004, 0};

    auto [pvalue, index] = adjusted_hochberg(raw);
    REQUIRE(pvalue == expected.first);
    REQUIRE(index == expected.second);
}

TEST_CASE("adjusted_hochberg same adjustement values") {
    std::vector<double> raw = {0.745386, 0.425089};
    std::pair<double, size_t> expected = {0.745386, 0.745386};

    auto [pvalue, index] = adjusted_hochberg(raw);
    REQUIRE(pvalue == expected.first);
    REQUIRE(index == expected.second);
}

TEST_CASE("stringToVector parses comma-separated values to vector") {
    std::string s = "4,578,6";
    std::vector<size_t> v = stringToVector<size_t>(s);
    REQUIRE(v.size() == 3);
    REQUIRE(v[0] == 4);
    REQUIRE(v[1] == 578);
    REQUIRE(v[2] == 6);
}

TEST_CASE("stringToVector throws on invalid input") {
    std::string s = "4,abc,6";
    REQUIRE_THROWS_AS(stringToVector<size_t>(s), std::runtime_error);
}
TEST_CASE("Snarl coordinates simple nested chain", "[snarl_info]") {

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
        std::vector<stoat::path_range_t> snarl1_paths = get_coordinates_of_snarl(*graph, distance_index, snarl1, true, std::unordered_set<std::string>(), false);
        REQUIRE(snarl1_paths.size() == 1);
        std::tuple<std::string, size_t, size_t> snarl1_coords = get_name_and_offsets_of_snarl_path_range(*graph, snarl1_paths.front());
        REQUIRE(std::get<0>(snarl1_coords) == "path0#0#path0");
        REQUIRE(std::get<1>(snarl1_coords) == 1);
        REQUIRE(std::get<2>(snarl1_coords) == 2);

        std::vector<stoat::path_range_t> snarl3_paths = get_coordinates_of_snarl(*graph, distance_index, snarl3, true, std::unordered_set<std::string>(), false);
        REQUIRE(snarl3_paths.size() == 1);
        std::tuple<std::string, size_t, size_t> snarl3_coords = get_name_and_offsets_of_snarl_path_range(*graph, snarl3_paths.front());
        REQUIRE(std::get<0>(snarl3_coords) == "path0#0#path0");
        REQUIRE(std::get<1>(snarl3_coords) == 4);
        REQUIRE(std::get<2>(snarl3_coords) == 5);
    }
    SECTION("named path") {
        std::vector<stoat::path_range_t> snarl1_paths = get_coordinates_of_snarl(*graph, distance_index, snarl1, false, std::unordered_set<std::string>({"path2"}), false);
        REQUIRE(snarl1_paths.size() == 1);
        std::tuple<std::string, size_t, size_t> snarl1_coords = get_name_and_offsets_of_snarl_path_range(*graph, snarl1_paths.front());
        REQUIRE(std::get<0>(snarl1_coords) == "path2");
        REQUIRE(std::get<1>(snarl1_coords) == 1);
        REQUIRE(std::get<2>(snarl1_coords) == 2);

        // Nested snarl off the reference
        std::vector<stoat::path_range_t> snarl3_paths = get_coordinates_of_snarl(*graph, distance_index, snarl3, false, std::unordered_set<std::string>({"path2"}), false);
        REQUIRE(snarl3_paths.size() == 1);
        std::tuple<std::string, size_t, size_t> snarl3_coords = get_name_and_offsets_of_snarl_path_range(*graph, snarl3_paths.front());
        REQUIRE(std::get<0>(snarl3_coords) == "path2");
        REQUIRE(std::get<1>(snarl3_coords) == 3);
        REQUIRE(std::get<2>(snarl3_coords) == 3);
    }

}
TEST_CASE("Snarl coordinates deeply nested snarls with deletion in reference", "[snarl_info]") {
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
        std::vector<stoat::path_range_t> snarl1_paths = get_coordinates_of_snarl(*graph, distance_index, snarl1, true, std::unordered_set<std::string>({"path0#0#0"}), false);
        REQUIRE(snarl1_paths.size() == 1);
        std::tuple<std::string, size_t, size_t> snarl1_coords = get_name_and_offsets_of_snarl_path_range(*graph, snarl1_paths.front());
        REQUIRE(std::get<0>(snarl1_coords) == "path0#0#0");
        REQUIRE(std::get<1>(snarl1_coords) == 1);
        REQUIRE(std::get<2>(snarl1_coords) == 1);

        // Nested snarl off the reference
        std::vector<stoat::path_range_t> snarl3_paths = get_coordinates_of_snarl(*graph, distance_index, snarl3, true, std::unordered_set<std::string>({"path0#0#0"}), false);
        REQUIRE(snarl3_paths.size() == 1);
        std::tuple<std::string, size_t, size_t> snarl3_coords = get_name_and_offsets_of_snarl_path_range(*graph, snarl3_paths.front());
        REQUIRE(std::get<0>(snarl3_coords) == "path0#0#0");
        REQUIRE(std::get<1>(snarl3_coords) == 1);
        REQUIRE(std::get<2>(snarl3_coords) == 1);
    }
}

TEST_CASE("Snarl coordinates deeply nested snarls with no reference sense", "[snarl_info]") {
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
        std::vector<stoat::path_range_t> snarl1_paths = get_coordinates_of_snarl(*graph, distance_index, snarl1, true, std::unordered_set<std::string>({"path0#0#0#0"}), false);
        REQUIRE(snarl1_paths.size() == 1);
        std::tuple<std::string, size_t, size_t> snarl1_coords = get_name_and_offsets_of_snarl_path_range(*graph, snarl1_paths.front());
        REQUIRE(std::get<0>(snarl1_coords) == "path0#0#0#0");
        REQUIRE(std::get<1>(snarl1_coords) == 1);
        REQUIRE(std::get<2>(snarl1_coords) == 1);

        // Nested snarl off the reference
        std::vector<stoat::path_range_t> snarl3_paths = get_coordinates_of_snarl(*graph, distance_index, snarl3, true, std::unordered_set<std::string>({"path0#0#0#0"}), false);
        REQUIRE(snarl3_paths.size() == 1);
        std::tuple<std::string, size_t, size_t> snarl3_coords = get_name_and_offsets_of_snarl_path_range(*graph, snarl3_paths.front());
        REQUIRE(std::get<0>(snarl3_coords) == "path0#0#0#0");
        REQUIRE(std::get<1>(snarl3_coords) == 1);
        REQUIRE(std::get<2>(snarl3_coords) == 1);
    }
}

TEST_CASE("Snarl coordinates on given path instead of reference", "[snarl_info]") {
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
        std::vector<stoat::path_range_t> snarl1_paths = get_coordinates_of_snarl(*graph, distance_index, snarl1, true, std::unordered_set<std::string>({"path1#0#0#0"}), false);
        REQUIRE(snarl1_paths.size() == 1);
        std::tuple<std::string, size_t, size_t> snarl1_coords = get_name_and_offsets_of_snarl_path_range(*graph, snarl1_paths.front());
        REQUIRE(std::get<0>(snarl1_coords) == "path1#0#0#0");
        REQUIRE(std::get<1>(snarl1_coords) == 1);
        REQUIRE(std::get<2>(snarl1_coords) == 1);

        // Nested snarl off the reference
        std::vector<stoat::path_range_t> snarl3_paths = get_coordinates_of_snarl(*graph, distance_index, snarl3, true, std::unordered_set<std::string>({"path1#0#0#0"}), false);
        REQUIRE(snarl3_paths.size() == 1);
        std::tuple<std::string, size_t, size_t> snarl3_coords = get_name_and_offsets_of_snarl_path_range(*graph, snarl3_paths.front());
        REQUIRE(std::get<0>(snarl3_coords) == "path1#0#0#0");
        REQUIRE(std::get<1>(snarl3_coords) == 1);
        REQUIRE(std::get<2>(snarl3_coords) == 1);
    }
    SECTION("reference") {
        std::vector<stoat::path_range_t> snarl1_paths = get_coordinates_of_snarl(*graph, distance_index, snarl1, true, std::unordered_set<std::string>(), false);
        REQUIRE(snarl1_paths.size() == 1);
        std::tuple<std::string, size_t, size_t> snarl1_coords = get_name_and_offsets_of_snarl_path_range(*graph, snarl1_paths.front());
        REQUIRE(std::get<0>(snarl1_coords) == "path0#0#0");
        REQUIRE(std::get<1>(snarl1_coords) == 1);
        REQUIRE(std::get<2>(snarl1_coords) == 3);

        // Nested snarl off the reference
        std::vector<stoat::path_range_t> snarl3_paths = get_coordinates_of_snarl(*graph, distance_index, snarl3, true, std::unordered_set<std::string>(), false);
        REQUIRE(snarl3_paths.size() == 1);
        std::tuple<std::string, size_t, size_t> snarl3_coords = get_name_and_offsets_of_snarl_path_range(*graph, snarl3_paths.front());
        REQUIRE(std::get<0>(snarl3_coords) == "path0#0#0");
        REQUIRE(std::get<1>(snarl3_coords) == 2);
        REQUIRE(std::get<2>(snarl3_coords) == 2);
    }
    SECTION("given sample not on snarl") {
        std::vector<stoat::path_range_t> snarl1_paths = get_coordinates_of_snarl(*graph, distance_index, snarl1, true, std::unordered_set<std::string>({"path3#0#0#0"}), false);
        REQUIRE(snarl1_paths.size() == 1);
        std::tuple<std::string, size_t, size_t> snarl1_coords = get_name_and_offsets_of_snarl_path_range(*graph, snarl1_paths.front());
        REQUIRE(std::get<0>(snarl1_coords) == "path0#0#0");
        REQUIRE(std::get<1>(snarl1_coords) == 1);
        REQUIRE(std::get<2>(snarl1_coords) == 3);

        // Nested snarl off the reference
        std::vector<stoat::path_range_t> snarl3_paths = get_coordinates_of_snarl(*graph, distance_index, snarl3, true, std::unordered_set<std::string>({"path3#0#0#0"}), false);
        REQUIRE(snarl3_paths.size() == 1);
        std::tuple<std::string, size_t, size_t> snarl3_coords = get_name_and_offsets_of_snarl_path_range(*graph, snarl3_paths.front());
        REQUIRE(std::get<0>(snarl3_coords) == "path0#0#0");
        REQUIRE(std::get<1>(snarl3_coords) == 2);
        REQUIRE(std::get<2>(snarl3_coords) == 2);
    }
}
