#include <catch.hpp>
#include "../../src/matrix.hpp"

using namespace stoat;

class TestEdgeBySampleMatrix : stoat_vcf::EdgeBySampleMatrix {
    public:
    TestEdgeBySampleMatrix(const std::vector<std::string>& sampleNames, size_t rows) : EdgeBySampleMatrix(sampleNames, rows) {}

    using stoat_vcf::EdgeBySampleMatrix::matrix_1D;
    using stoat_vcf::EdgeBySampleMatrix::operator();
    using stoat_vcf::EdgeBySampleMatrix::max_edges;
    using stoat_vcf::EdgeBySampleMatrix::expand_matrix;
    using stoat_vcf::EdgeBySampleMatrix::shrink;
    using stoat_vcf::EdgeBySampleMatrix::set_edge;
};

TEST_CASE("stoat_vcf::EdgeBySampleMatrix Constructor and Basic Properties", "[stoat_vcf::EdgeBySampleMatrix]") {
    SECTION("stoat_vcf::EdgeBySampleMatrix initializes correctly") {
        std::vector<string> sample_names = {"1", "2"};
        TestEdgeBySampleMatrix mat(sample_names, 4);
        REQUIRE(mat.matrix_1D.size() > 0);  // Ensure matrix is allocated
        REQUIRE_FALSE(mat(0, 0));  // Initially, all elements should be false
    }
}

TEST_CASE("stoat_vcf::EdgeBySampleMatrix Expansion", "[stoat_vcf::EdgeBySampleMatrix]") {
    SECTION("stoat_vcf::EdgeBySampleMatrix expands properly") {
        std::vector<string> sample_names = {"1", "2"};
        TestEdgeBySampleMatrix mat(sample_names, 4);
        size_t original_size = mat.matrix_1D.size();
        mat.expand_matrix();
        REQUIRE(mat.matrix_1D.size() > original_size);
    }
}
TEST_CASE("stoat_vcf::EdgeBySampleMatrix Set and Access Elements", "[stoat_vcf::EdgeBySampleMatrix]") {
    SECTION("stoat_vcf::EdgeBySampleMatrix correctly sets and retrieves values") {
        std::vector<string> sample_names = {"1", "2"};
        TestEdgeBySampleMatrix mat(sample_names, 4);
        REQUIRE_FALSE(mat(1, 3));  // Initially, should be false
        mat.set_edge(1, 3);
        REQUIRE(mat(1, 3));  // Should now be true
    }
}
// TODO: Make this shrink to a specific size
TEST_CASE("stoat_vcf::EdgeBySampleMatrix Shrink", "[stoat_vcf::EdgeBySampleMatrix]") {
    SECTION("stoat_vcf::EdgeBySampleMatrix correctly shrinks") {
        std::vector<string> sample_names = {"1", "2"};
        TestEdgeBySampleMatrix mat(sample_names, 10);
        size_t original_size = mat.matrix_1D.size();
        mat.shrink();  // Reduce row count
        REQUIRE(mat.matrix_1D.size() < original_size);
    }
}

TEST_CASE("stoat_vcf::EdgeBySampleMatrix Maximum Element", "[stoat_vcf::EdgeBySampleMatrix]") {
    SECTION("stoat_vcf::EdgeBySampleMatrix tracks maximum element correctly") {
        std::vector<string> sample_names = {"1", "2"};
        TestEdgeBySampleMatrix mat(sample_names, 4);
        REQUIRE(mat.max_edges == 4);  // Initially zero
        mat.expand_matrix();
        REQUIRE(mat.max_edges == 8);  // Should be updated
    }
}
