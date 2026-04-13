#include <catch.hpp>
#include "../../src/matrix.hpp"

using namespace stoat;

class TestAlleleBySampleMatrix : stoat_vcf::AlleleBySampleMatrix<stoat::edge_t> {
    public:
    TestAlleleBySampleMatrix(const std::vector<std::string>& sampleNames, size_t rows) : AlleleBySampleMatrix<stoat::edge_t>(sampleNames, rows) {}

    using stoat_vcf::AlleleBySampleMatrix<stoat::edge_t>::matrix_1D;
    using stoat_vcf::AlleleBySampleMatrix<stoat::edge_t>::operator();
    using stoat_vcf::AlleleBySampleMatrix<stoat::edge_t>::max_keys;
    using stoat_vcf::AlleleBySampleMatrix<stoat::edge_t>::expand_matrix;
    using stoat_vcf::AlleleBySampleMatrix<stoat::edge_t>::shrink;
    using stoat_vcf::AlleleBySampleMatrix<stoat::edge_t>::set_key;
};

TEST_CASE("stoat_vcf::AlleleBySampleMatrix Constructor and Basic Properties", "[stoat_vcf::AlleleBySampleMatrix]") {
    SECTION("stoat_vcf::AlleleBySampleMatrix initializes correctly") {
        std::vector<string> sample_names = {"1", "2"};
        TestAlleleBySampleMatrix mat(sample_names, 4);
        REQUIRE(mat.matrix_1D.size() > 0);  // Ensure matrix is allocated
        REQUIRE_FALSE(mat(0, 0));  // Initially, all elements should be false
    }
}

TEST_CASE("stoat_vcf::AlleleBySampleMatrix Expansion", "[stoat_vcf::AlleleBySampleMatrix]") {
    SECTION("stoat_vcf::AlleleBySampleMatrix expands properly") {
        std::vector<string> sample_names = {"1", "2"};
        TestAlleleBySampleMatrix mat(sample_names, 4);
        size_t original_size = mat.matrix_1D.size();
        mat.expand_matrix();
        REQUIRE(mat.matrix_1D.size() > original_size);
    }
}
TEST_CASE("stoat_vcf::AlleleBySampleMatrix Set and Access Elements", "[stoat_vcf::AlleleBySampleMatrix]") {
    SECTION("stoat_vcf::AlleleBySampleMatrix correctly sets and retrieves values") {
        std::vector<string> sample_names = {"1", "2"};
        TestAlleleBySampleMatrix mat(sample_names, 4);
        REQUIRE_FALSE(mat(1, 3));  // Initially, should be false
        mat.set_key(1, 3);
        REQUIRE(mat(1, 3));  // Should now be true
    }
}
// TODO: Make this shrink to a specific size
TEST_CASE("stoat_vcf::AlleleBySampleMatrix Shrink", "[stoat_vcf::AlleleBySampleMatrix]") {
    SECTION("stoat_vcf::AlleleBySampleMatrix correctly shrinks") {
        std::vector<string> sample_names = {"1", "2"};
        TestAlleleBySampleMatrix mat(sample_names, 10);
        size_t original_size = mat.matrix_1D.size();
        mat.shrink();  // Reduce row count
        REQUIRE(mat.matrix_1D.size() < original_size);
    }
}

TEST_CASE("stoat_vcf::AlleleBySampleMatrix Maximum Element", "[stoat_vcf::AlleleBySampleMatrix]") {
    SECTION("stoat_vcf::AlleleBySampleMatrix tracks maximum element correctly") {
        std::vector<string> sample_names = {"1", "2"};
        TestAlleleBySampleMatrix mat(sample_names, 4);
        REQUIRE(mat.max_keys == 4);  // Initially zero
        mat.expand_matrix();
        REQUIRE(mat.max_keys == 8);  // Should be updated
    }
}
