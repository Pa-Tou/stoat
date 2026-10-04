#include <catch.hpp>
#include <bdsg/overlays/overlay_helper.hpp>
#include <bdsg/hash_graph.hpp>

#include "../../src/types_and_structs.hpp"
#include "../../src/utils.hpp"

using namespace stoat; 

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

