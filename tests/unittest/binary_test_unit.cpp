#include <catch.hpp>
#include "../../src/binary_table.hpp"
#include "../../src/stats_test.hpp"

// using namespace stoat_vcf; 
// TEST_CASE("Chi-square & Fisher test function", "[chi2_2xN]") {
//     SECTION("Valid chi-square test & valid Fisher test calculation") {
//         std::vector<size_t> g0 = {10, 20};
//         std::vector<size_t> g1 = {20, 10};
//         size_t a = g0[0];
//         size_t b = g0[1];
//         size_t c = g1[0];
//         size_t d = g1[1];
//         REQUIRE(chi2_2x2(a, b, c, d)  == "0.0098");
//         REQUIRE(fastFishersExactTest(a, b, c, d)  == "0.0194");
//     }

//     SECTION("Chi-square & Fisher test (significatif)") {
//         std::vector<size_t> g0 = {30, 5};
//         std::vector<size_t> g1 = {2, 25};
//         size_t a = g0[0];
//         size_t b = g0[1];
//         size_t c = g1[0];
//         size_t d = g1[1];
//         REQUIRE(chi2_2x2(a, b, c, d)  == "9.5037e-10");
//         REQUIRE(fastFishersExactTest(a, b, c, d)  == "3.5379e-10");
//     }

//     SECTION("Chi-square fail (N row)") {
//         std::vector<size_t> g0 = {10, 15, 5};
//         std::vector<size_t> g1 = {20, 10, 10};
//         REQUIRE(chi2_2xN(g0, g1) == "0.0970");
//     }

//     SECTION("Chi-square N row significatif") {
//         std::vector<size_t> g0 = {5, 10, 15, 20};
//         std::vector<size_t> g1 = {20, 15, 10, 5};
//         REQUIRE(chi2_2xN(g0, g1) == "0.0002");
//     }

//     SECTION("Chi-square fail (N row 1.0000)") {
//         std::vector<size_t> g0 = {10, 10, 10, 10, 10};
//         std::vector<size_t> g1 = {10, 10, 10, 10, 10};
//         REQUIRE(chi2_2xN(g0, g1) == "1.0000");
//     }

//     SECTION("Chi-square fail & Fisher test fail (full zero row)") {
//         std::vector<size_t> g0 = {0, 0};
//         std::vector<size_t> g1 = {0, 0};
//         size_t a = g0[0];
//         size_t b = g0[1];
//         size_t c = g1[0];
//         size_t d = g1[1];
//         REQUIRE(chi2_2x2(a, b, c, d)  == "NA");
//         REQUIRE(fastFishersExactTest(a, b, c, d)  == "NA");
//     }

//     SECTION("Chi-square fail (zero row)") {
//         std::vector<size_t> g0 = {0, 0, 0};
//         std::vector<size_t> g1 = {10, 20, 30};
//         REQUIRE(chi2_2xN(g0, g1) == "NA");
//     }

//     SECTION("Chi-square fail & (zero column)") {
//         std::vector<size_t> g0 = {0, 10, 5};
//         std::vector<size_t> g1 = {0, 20, 15};
//         REQUIRE(chi2_2xN(g0, g1) == "NA");
//     }

//     SECTION("Chi-square fail & Fisher test valid (zero row + column)") {
//         std::vector<size_t> g0 = {0, 0};
//         std::vector<size_t> g1 = {0, 1};
//         size_t a = g0[0];
//         size_t b = g0[1];
//         size_t c = g1[0];
//         size_t d = g1[1];
//         REQUIRE(chi2_2x2(a, b, c, d) == "NA");
//         REQUIRE(fastFishersExactTest(a, b, c, d) == "NA");
//     }

//     SECTION("Chi-square & Fisher test (1/0 0/1)") {
//         std::vector<size_t> g0 = {1, 0};
//         std::vector<size_t> g1 = {0, 1};
//         size_t a = g0[0];
//         size_t b = g0[1];
//         size_t c = g1[0];
//         size_t d = g1[1];
//         REQUIRE(chi2_2x2(a, b, c, d)  == "0.1573");
//         REQUIRE(fastFishersExactTest(a, b, c, d)  == "1.0000");
//     }

//     SECTION("Chi-square & Fisher test (strange but correct)") {
//         std::vector<size_t> g0 = {79, 18};
//         std::vector<size_t> g1 = {96, 23};
//         size_t a = g0[0];
//         size_t b = g0[1];
//         size_t c = g1[0];
//         size_t d = g1[1];
//         REQUIRE(chi2_2x2(a, b, c, d)  == "0.8857");
//         REQUIRE(fastFishersExactTest(a, b, c, d)  == "1.0000");
//     }

//     SECTION("Chi-square & Fisher test (very significative)") {
//         std::vector<size_t> g0 = {122, 78};
//         std::vector<size_t> g1 = {27, 173};
//         size_t a = g0[0];
//         size_t b = g0[1];
//         size_t c = g1[0];
//         size_t d = g1[1];
//         REQUIRE(chi2_2x2(a, b, c, d)  == "8.8051e-23");
//         REQUIRE(fastFishersExactTest(a, b, c, d)  == "1.4799e-23");
//     }
// }
