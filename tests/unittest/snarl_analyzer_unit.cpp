#include <catch.hpp>
#include "../../src/snarl_data_t.hpp"
#include "../../src/snarl_analyzer.hpp"
#include "../../src/utils.hpp"
#include "../../src/matrix.hpp"

#include <limits>

using namespace stoat;
using namespace stoat_vcf;

TEST_CASE("stoat::Node_traversal_t Basic Functionality") {
    stoat::Node_traversal_t node(42, true);
    REQUIRE(node.get_node_id() == 42);
    REQUIRE(node.get_is_reverse() == true);
    REQUIRE(node.to_string() == "<42");
}

TEST_CASE("Edge_t Functionality") {
    stoat::Node_traversal_t a(1, false);
    stoat::Node_traversal_t b(2, true);
    stoat::Edge_t edge(a, b);

    auto pair = edge.print_pair_edge();
    REQUIRE(pair.first == 1);
    REQUIRE(pair.second == 2);

    REQUIRE(edge.print_string_edge() == ">1<2");
}

TEST_CASE("Path_traversal_t Add and Get") {
    Path_traversal_t path;
    path.add_node_traversal_t({1, false});
    path.add_node_traversal_t({2, true});
    const auto& paths = path.get_paths();
    REQUIRE(paths.size() == 2);
    REQUIRE(paths[0].get_node_id() == 1);
    REQUIRE(paths[1].get_is_reverse() == true);
}

TEST_CASE("decompose_path_to_edges") {
    Path_traversal_t path;
    path.add_node_traversal_t({1, false});
    path.add_node_traversal_t({2, false});
    path.add_node_traversal_t({3, false});

    auto edges = decompose_path_to_edges(path);
    REQUIRE(edges.size() == 2);
    REQUIRE(edges[0].print_pair_edge() == std::make_pair(1ul, 2ul));
    REQUIRE(edges[1].print_pair_edge() == std::make_pair(2ul, 3ul));
}

TEST_CASE("decompose_path_str_to_edge handles basic and complex strings") {
    SECTION("Basic case") {
        auto edges = decompose_path_str_to_edge(">1<2>3");
        REQUIRE(edges.size() == 2);
        REQUIRE(edges[0].print_string_edge() == ">1<2");
        REQUIRE(edges[1].print_string_edge() == "<2>3");
    }

    SECTION("Special case with zeros (complex path)") {
        auto edges = decompose_path_str_to_edge(">1<324<323<0<213>214<0<213");
        REQUIRE(edges.size() == 7);

        REQUIRE(edges[0].print_string_edge() == ">1<324");
        REQUIRE(edges[1].print_string_edge() == "<324<323");
        REQUIRE(edges[2].print_string_edge() == "<323<0");
        REQUIRE(edges[3].print_string_edge() == "<0<213");
        REQUIRE(edges[4].print_string_edge() == "<213>214");
        REQUIRE(edges[5].print_string_edge() == ">214<0");
        REQUIRE(edges[6].print_string_edge() == "<0<213");
    }
}

TEST_CASE("decompose_path_list_str supports multiple strings including ones with zero nodes") {
    std::vector<std::string> input = {
        ">1>2",
        "<3<4",
        ">1<324<323<0<213>214<0<213" // complex path
    };

    auto decomposed = decompose_path_list_str(input);
    
    REQUIRE(decomposed.size() == 3);

    // First path: >1>2
    REQUIRE(decomposed[0].size() == 1);
    REQUIRE(decomposed[0][0].print_string_edge() == ">1>2");

    // Second path: <3<4
    REQUIRE(decomposed[1].size() == 1);
    REQUIRE(decomposed[1][0].print_string_edge() == "<3<4");

    // Third path: complex path with 0s
    REQUIRE(decomposed[2].size() == 7);
    REQUIRE(decomposed[2][0].print_string_edge() == ">1<324");
    REQUIRE(decomposed[2][6].print_string_edge() == "<0<213");
}

TEST_CASE("identify_path with EdgeBySampleMatrix") {
    stoat::Node_traversal_t a(1, false), b(2, false), c(3, false);
   stoat::Edge_t edge1(a, b);
   stoat::Edge_t edge2(b, c);

    std::vector<stoat::Edge_t> path = {edge1, edge2};

    std::vector<std::string> samples = {"sample1", "sample2", "sample3"};
    EdgeBySampleMatrix matrix(samples, 2, 3);

    matrix.push_matrix(edge1, 0); // Set true at [row for edge1][0]
    matrix.push_matrix(edge2, 0); // Set true at [row for edge2][0]

    matrix.push_matrix(edge1, 2); // Also at [][2]
    matrix.push_matrix(edge2, 2);

    auto result = identify_path(path, matrix, 3);
    REQUIRE(result == std::vector<size_t>({0, 2}));
}

TEST_CASE("filtration_quantitative_table basic filtering") {
    std::vector<std::vector<double>> df = {
        {1.0, 1.0},
        {0.5, 0.5}
    };

    SECTION("Valid data, should NOT be filtered") {
        bool result = filtration_quantitative_table(df, 2, 2, 0.1);
        REQUIRE(result == false);
    }

    SECTION("Not enough individuals, should be filtered") {
        bool result = filtration_quantitative_table({{1.0, 1.0}}, 2, 2, 0.1);
        REQUIRE(result == true);
    }

    SECTION("Not enough haplotypes, should be filtered") {
        std::vector<std::vector<double>> low_sum = {
            {0.1, 0.1},
            {0.1, 0.1}
        };
        bool result = filtration_quantitative_table(low_sum, 2, 2, 0.1);
        REQUIRE(result == true);
    }

    SECTION("Low MAFs, should be filtered") {
        std::vector<std::vector<double>> low_maf = {
            {1.9, 0.1},
            {1.9, 0.1}
        };
        bool result = filtration_quantitative_table(low_maf, 2, 2, 0.4);
        REQUIRE(result == true);
    }
}

TEST_CASE("remove_empty_columns_quantitative_table removes zero") {
    std::vector<std::vector<double>> df = {
        {1.0, 0.0, 0.0, 0.0},
        {0.5, 0.0, 0.5, 0.0},
    };

    remove_empty_columns_quantitative_table(df);
    
    REQUIRE(df[0].size() == 2);
    REQUIRE(df[1].size() == 2);
    REQUIRE(df[0][0] == 1.0);
    REQUIRE(df[1][1] == 0.5);
}

TEST_CASE("remove_last_columns_quantitative_table works correctly") {
    std::vector<std::vector<double>> df = {
        {0.5, 0.5},
        {0, 0.5}
    };

    remove_last_columns_quantitative_table(df);

    REQUIRE(df[0].size() == 1);
    REQUIRE(df[1].size() == 1);
    REQUIRE(df[0][0] == 0.5);
    REQUIRE(df[1][0] == 0);
}

TEST_CASE("remove_empty_columns_binary_table filters empty columns") {
    std::vector<size_t> g0 = {1, 0, 2, 0};
    std::vector<size_t> g1 = {0, 0, 1, 0};

    remove_empty_columns_binary_table(g0, g1);

    REQUIRE(g0 == std::vector<size_t>({1, 2}));
    REQUIRE(g1 == std::vector<size_t>({0, 1}));
}

TEST_CASE("filtration_binary_table logic correctness") {
    SECTION("Valid case, should NOT be filtered") {
        std::vector<size_t> g0 = {5, 1};
        std::vector<size_t> g1 = {1, 5};
        bool result = filtration_binary_table(g0, g1, 12, 2, 4, 0.1);
        REQUIRE(result == false);
    }

    SECTION("Too few individuals, should be filtered") {
        std::vector<size_t> g0 = {5, 1};
        std::vector<size_t> g1 = {1, 5};
        bool result = filtration_binary_table(g0, g1, 2, 2, 4, 0.1);
        REQUIRE(result == true);
    }

    SECTION("Too few haplotypes, should be filtered") {
        std::vector<size_t> g0 = {1, 1};
        std::vector<size_t> g1 = {1, 1};
        bool result = filtration_binary_table(g0, g1, 4, 2, 10, 0.1);
        REQUIRE(result == true);
    }

    SECTION("Low MAFs, should be filtered") {
        std::vector<size_t> g0 = {10, 1};
        std::vector<size_t> g1 = {0, 10};
        bool result = filtration_binary_table(g0, g1, 21, 2, 4, 0.45);
        REQUIRE(result == true);
    }
}
