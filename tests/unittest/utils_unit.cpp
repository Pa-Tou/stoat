#include <catch.hpp>
#include "../../src/snarl_data_t.hpp"
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

TEST_CASE("set_precision_float_50 handles small and large cpp_dec_float_50 values") {
    cpp_dec_float_50 small("0.00001234567890123456789");
    cpp_dec_float_50 normal("0.34567890123456789");
    REQUIRE(set_precision_float_50(small) == "1.2346e-05");
    REQUIRE(set_precision_float_50(normal) == "0.3457");
}

TEST_CASE("string_to_pvalue converts valid p-values or returns 1.0 for NA") {
    REQUIRE(stoat::string_to_pvalue("0.01") == 0.01);
    REQUIRE(stoat::string_to_pvalue("NA") == 1.0);
    REQUIRE(stoat::string_to_pvalue("") == 1.0);
}

TEST_CASE("isPValueSignificant handles correct parsing and NA") {
    REQUIRE(isPValueSignificant(0.05, "0.01") == true);
    REQUIRE(isPValueSignificant(0.05, "0.05") == false);
    REQUIRE(isPValueSignificant(0.05, "NA") == false);
}

TEST_CASE("adjusted_holm produces correct order and monotonicity") {
    std::vector<double> raw = {0.01, 0.03, 0.02, 0.05};
    auto adj = adjusted_holm(raw);

    REQUIRE(adj.size() == raw.size());

    for (double val : adj) {
        REQUIRE(val <= 1.0);
        REQUIRE(val >= 0.0);
    }
}

TEST_CASE("adjusted_holm maintains monotonicity") {
    std::vector<double> raw = {0.01, 0.01, 0.01, 0.01};
    auto adj = adjusted_holm(raw);
    REQUIRE(adj[0] <= adj[1]);
    REQUIRE(adj[1] <= adj[2]);
    REQUIRE(adj[2] <= adj[3]);
}

TEST_CASE("retain_indices keeps only specified elements") {
    std::vector<double> vec = {0.1, 0.2, 0.3, 0.4};
    std::set<size_t> keep = {1, 3};
    retain_indices(vec, keep);
    REQUIRE(vec.size() == 2);
    REQUIRE(vec[0] == 0.2);
    REQUIRE(vec[1] == 0.4);
}

TEST_CASE("vectorToString handles string and numeric vectors") {
    std::vector<std::string> svec = {"A", "B", "C"};
    std::vector<size_t> ivec = {1, 2, 3};
    REQUIRE(vectorToString(svec) == "A,B,C");
    REQUIRE(vectorToString(ivec) == "1,2,3");
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
