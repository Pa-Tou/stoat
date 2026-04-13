#include <catch.hpp>
#include "../../src/types_and_structs.hpp"
#include "../../src/snarl_analyzer.hpp"
#include "../../src/utils.hpp"
#include "../../src/matrix.hpp"

#include <limits>

using namespace stoat;
using namespace stoat_vcf;

TEST_CASE("stoat::node_traversal_t Basic Functionality") {
    stoat::node_traversal_t node(42, true);
    REQUIRE(node.get_node_id() == 42);
    REQUIRE(node.get_is_reverse() == true);
    REQUIRE(node.to_string() == "<42");
}

TEST_CASE("edge_t Functionality") {
    stoat::node_traversal_t a(1, false);
    stoat::node_traversal_t b(2, true);
    stoat::edge_t edge(a, b);

    auto pair = edge.get_node_pair();
    REQUIRE(pair.first == 1);
    REQUIRE(pair.second == 2);

    REQUIRE(edge.to_string() == ">1<2");
}

TEST_CASE("PathTraversal Add and Get") {
    PathTraversal path;
    path.add_node_traversal_t({1, false});
    path.add_node_traversal_t({2, true});
    const auto& paths = path.get_path();
    REQUIRE(paths.size() == 2);
    REQUIRE(paths[0].get_node_id() == 1);
    REQUIRE(paths[1].get_is_reverse() == true);
}

TEST_CASE("decompose_path_str_to_edge handles basic and complex strings") {
    SECTION("Basic case") {
        auto edges = decompose_path_str_to_edge(">1<2>3");
        REQUIRE(edges.size() == 2);
        REQUIRE(edges[0].to_string() == ">1<2");
        REQUIRE(edges[1].to_string() == "<2>3");
    }

    SECTION("Special case with zeros (complex path)") {
        auto edges = decompose_path_str_to_edge(">1<324<323<0<213>214<0<213");
        REQUIRE(edges.size() == 7);

        REQUIRE(edges[0].to_string() == ">1<324");
        REQUIRE(edges[1].to_string() == "<324<323");
        REQUIRE(edges[2].to_string() == "<323<0");
        REQUIRE(edges[3].to_string() == "<0<213");
        REQUIRE(edges[4].to_string() == "<213>214");
        REQUIRE(edges[5].to_string() == ">214<0");
        REQUIRE(edges[6].to_string() == "<0<213");
    }
}

TEST_CASE("identify_path with AlleleBySampleMatrix") {
    stoat::node_traversal_t a(1, false), b(2, false), c(3, false);
    stoat::edge_t edge1(a, b);
    stoat::edge_t edge2(b, c);

    stoat::PathTraversal path;
    path.add_node_traversal_t(a);
    path.add_node_traversal_t(b);
    path.add_node_traversal_t(c);
    
    std::vector<std::string> samples = {"sample1", "sample2", "sample3"};
    AlleleBySampleMatrix<stoat::edge_t> matrix(samples, 2);

    matrix.add_sample_key(edge1, 0); // Set true at [row for edge1][0]
    matrix.add_sample_key(edge2, 0); // Set true at [row for edge2][0]

    matrix.add_sample_key(edge1, 2); // Also at [][2]
    matrix.add_sample_key(edge2, 2);

    auto result = matrix.get_samples_on_path(path);
    REQUIRE(result == std::vector<size_t>({0, 2}));
}
