#include <catch.hpp>
#include "../../src/stats_test.hpp"

using namespace stoat; 

TEST_CASE("Chi-square & Fisher test function", "[fchi.chi2_2xN]") {
    FisherChi2 fchi;

    SECTION("Valid chi-square test & valid Fisher test calculation") {
        std::vector<size_t> g0 = {10, 20};
        std::vector<size_t> g1 = {20, 10};
        size_t a = g0[0];
        size_t b = g0[1];
        size_t c = g1[0];
        size_t d = g1[1];
        REQUIRE(stoat::set_precision(fchi.chi2_2x2(a, b, c, d)) == "9.8233e-03");
        REQUIRE(stoat::set_precision(fchi.fastFishersExactTest(a, b, c, d)) == "1.9383e-02");
    }

    SECTION("Chi-square & Fisher test (significatif)") {
        std::vector<size_t> g0 = {30, 5};
        std::vector<size_t> g1 = {2, 25};
        size_t a = g0[0];
        size_t b = g0[1];
        size_t c = g1[0];
        size_t d = g1[1];
        REQUIRE(stoat::set_precision(fchi.chi2_2x2(a, b, c, d)) == "9.5037e-10");
        REQUIRE(stoat::set_precision(fchi.fastFishersExactTest(a, b, c, d)) == "3.5379e-10");
    }

    SECTION("Chi-square fail (N row)") {
        std::vector<size_t> g0 = {10, 15, 5};
        std::vector<size_t> g1 = {20, 10, 10};
        REQUIRE(stoat::set_precision(fchi.chi2_2xN(g0, g1)) == "9.6972e-02");
    }

    SECTION("Chi-square N row significatif") {
        std::vector<size_t> g0 = {5, 10, 15, 20};
        std::vector<size_t> g1 = {20, 15, 10, 5};
        REQUIRE(stoat::set_precision(fchi.chi2_2xN(g0, g1)) == "1.6974e-04");
    }

    SECTION("Chi-square fail (N row 1)") {
        std::vector<size_t> g0 = {10, 10, 10, 10, 10};
        std::vector<size_t> g1 = {10, 10, 10, 10, 10};
        REQUIRE(stoat::set_precision(fchi.chi2_2xN(g0, g1)) == "1");
    }

    SECTION("Chi-square fail & Fisher test fail (full zero row)") {
        std::vector<size_t> g0 = {0, 0};
        std::vector<size_t> g1 = {0, 0};
        size_t a = g0[0];
        size_t b = g0[1];
        size_t c = g1[0];
        size_t d = g1[1];
        REQUIRE(std::isnan(fchi.chi2_2x2(a, b, c, d)));
        REQUIRE(std::isnan(fchi.fastFishersExactTest(a, b, c, d)));
    }

    SECTION("Chi-square fail (zero row)") {
        std::vector<size_t> g0 = {0, 0, 0};
        std::vector<size_t> g1 = {10, 20, 30};
        REQUIRE(std::isnan(fchi.chi2_2xN(g0, g1)));
    }

    SECTION("Chi-square fail & Fisher test valid (zero row + column)") {
        std::vector<size_t> g0 = {0, 0};
        std::vector<size_t> g1 = {0, 1};
        size_t a = g0[0];
        size_t b = g0[1];
        size_t c = g1[0];
        size_t d = g1[1];
        REQUIRE(std::isnan(fchi.chi2_2x2(a, b, c, d)));
        REQUIRE(std::isnan(fchi.fastFishersExactTest(a, b, c, d)));
    }

    SECTION("Chi-square & Fisher test (1/0 0/1)") {
        std::vector<size_t> g0 = {1, 0};
        std::vector<size_t> g1 = {0, 1};
        size_t a = g0[0];
        size_t b = g0[1];
        size_t c = g1[0];
        size_t d = g1[1];
        REQUIRE(stoat::set_precision(fchi.chi2_2x2(a, b, c, d)) == "0.1573");
        REQUIRE(stoat::set_precision(fchi.fastFishersExactTest(a, b, c, d)) == "1");
    }

    SECTION("Chi-square & Fisher test (strange but correct)") {
        std::vector<size_t> g0 = {79, 18};
        std::vector<size_t> g1 = {96, 23};
        size_t a = g0[0];
        size_t b = g0[1];
        size_t c = g1[0];
        size_t d = g1[1];
        REQUIRE(stoat::set_precision(fchi.chi2_2x2(a, b, c, d)) == "0.8857");
        REQUIRE(stoat::set_precision(fchi.fastFishersExactTest(a, b, c, d)) == "1");
    }

    SECTION("Chi-square & Fisher test (very significative)") {
        std::vector<size_t> g0 = {122, 78};
        std::vector<size_t> g1 = {27, 173};
        size_t a = g0[0];
        size_t b = g0[1];
        size_t c = g1[0];
        size_t d = g1[1];
        REQUIRE(stoat::set_precision(fchi.chi2_2x2(a, b, c, d)) == "8.8051e-23");
        REQUIRE(stoat::set_precision(fchi.fastFishersExactTest(a, b, c, d)) == "1.4799e-23");
    }
    SECTION("Chi-square & Fisher test (very significative) with empty counts") {
        std::vector<size_t> g0 = {122, 0, 78};
        std::vector<size_t> g1 = {27, 0, 173};
        size_t a = g0[0];
        size_t b = g0[1];
        size_t c = g1[0];
        size_t d = g1[1];
        REQUIRE(stoat::set_precision(fchi.chi2_2xN(g0, g1)) == "8.8051e-23");
    }

    SECTION("four alleles and some very low counts") {
        //     A0 A1 A2 A3
        // S1 {1, 0, 1, 0}
        // S2 {1, 0, 1, 0}
        // S3 {1, 1, 0, 0}
        // S4 {0, 1, 1, 0}
        // S5 {0, 0, 1, 1}

        std::vector<size_t> g0 = {2, 2, 2, 0};
        std::vector<size_t> g1 = {1, 0, 2, 1};
        
        std::pair<double, double> pvs = fchi.fisher_chi2(g0, g1);

        INFO("p_value fisher = " << pvs.first);
        INFO("p_value chi2 = " << pvs.second);
        
        REQUIRE(std::isnan(pvs.first));
        REQUIRE(std::abs(pvs.second - 0.3831) < 0.01);
    }


}
