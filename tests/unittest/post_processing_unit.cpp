#include <catch.hpp>
#include "../../src/post_processing.hpp"

#include <vector>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <limits>

using namespace stoat;

TEST_CASE("adjust_pvalues_with_BH basic behavior") {
    std::vector<std::tuple<double, double, size_t>> data = {
        {0.01, 0.0, 0},
        {0.04, 0.0, 1},
        {0.03, 0.0, 2},
        {0.002, 0.0, 3},
        {0.05, 0.0, 4}
    };

    adjust_pvalues_with_BH(data);

    // Extract adjusted values in original order
    std::vector<double> adjusted;
    for (const auto& tup : data) {
        adjusted.push_back(std::get<1>(tup));
    }

    REQUIRE(adjusted.size() == 5);

    // Expected monotonic BH adjustments after reordering:
    // sorted p: 0.002, 0.01, 0.03, 0.04, 0.05
    // BH:       0.002*5/1 = 0.01
    //           0.01*5/2 = 0.025
    //           0.03*5/3 ≈ 0.05
    //           0.04*5/4 = 0.05
    //           0.05*5/5 = 0.05
    // Monotonicity enforced: 0.01, 0.025, 0.05, 0.05, 0.05

    // Ordered back to indices: 0, 1, 2, 3, 4
    REQUIRE(adjusted[0] == Approx(0.025));
    REQUIRE(adjusted[1] == Approx(0.05));
    REQUIRE(adjusted[2] == Approx(0.05));
    REQUIRE(adjusted[3] == Approx(0.01));
    REQUIRE(adjusted[4] == Approx(0.05));
}

TEST_CASE("adjust_pvalues_with_BH handles empty input") {
    std::vector<std::tuple<double, double, size_t>> data;
    adjust_pvalues_with_BH(data);
    REQUIRE(data.empty());
}

TEST_CASE("adjust_pvalues_with_BH clamps above 1.0") {
    std::vector<std::tuple<double, double, size_t>> data = {
        {0.9, 0.0, 0},
        {0.95, 0.0, 1},
        {0.99, 0.0, 2}
    };

    adjust_pvalues_with_BH(data);

    for (const auto& tup : data) {
        REQUIRE(std::get<1>(tup) <= 1.0);
    }
}

TEST_CASE("adjust_pvalues_with_BH preserves index ordering") {
    std::vector<std::tuple<double, double, size_t>> data = {
        {0.05, 0.0, 2},
        {0.01, 0.0, 0},
        {0.03, 0.0, 1}
    };

    adjust_pvalues_with_BH(data);

    REQUIRE(std::get<2>(data[0]) == 0);
    REQUIRE(std::get<2>(data[1]) == 1);
    REQUIRE(std::get<2>(data[2]) == 2);
}
