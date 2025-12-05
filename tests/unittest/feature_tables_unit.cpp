#include <catch.hpp>
#include "../../src/feature_tables.hpp"

using namespace stoat;

TEST_CASE( "BinaryPhenotypeTable with three samples", "[table]" ) {

    // Get sample_to_index 
    std::unordered_map<std::string, size_t> sample_to_index;
    for (size_t i = 0 ; i < 3 ; i++) {
        sample_to_index.emplace("sample" + std::to_string(i+1), i);
    }

    // Make the table
    BinaryPhenotypeTable table(sample_to_index);

    SECTION("Set and get values") {
        table.set_value_for_sample("sample1", false);
        table.set_value_for_sample("sample2", false);
        table.set_value_for_sample("sample3", true);
        table.set_value_for_sample("sample1", true);

        REQUIRE(table.get_value_for_sample("sample1") == true);
        REQUIRE(table.get_value_for_sample("sample2") == false);
        REQUIRE(table.get_value_for_sample("sample3") == true);
        
    }
}

TEST_CASE( "QuantitativePhenotypeTable with five samples", "[table]" ) {

    // Get sample_to_index 
    std::unordered_map<std::string, size_t> sample_to_index;
    for (size_t i = 0 ; i < 5 ; i++) {
        sample_to_index.emplace("sample" + std::to_string(i+1), i);
    }

    // Make the table
    QuantitativePhenotypeTable table(sample_to_index);

    SECTION("Set and get values") {
        table.set_value_for_sample("sample1", 1.0);
        table.set_value_for_sample("sample2", std::numeric_limits<double>::max());
        table.set_value_for_sample("sample3", 0.0);
        table.set_value_for_sample("sample4", 4.43243);
        table.set_value_for_sample("sample5", 4.43243);

        REQUIRE(table.get_value_for_sample("sample1") == 1.0);
        REQUIRE(table.get_value_for_sample("sample2") == std::numeric_limits<double>::max());
        REQUIRE(table.get_value_for_sample("sample3") == 0.0);
        REQUIRE(table.get_value_for_sample("sample4") == 4.43243);
        REQUIRE(table.get_value_for_sample("sample5") == 4.43243);
        
    }
}


TEST_CASE( "GenotypeTable with three samples and three alleles", "[table]" ) {

    // Get sample_to_index 
    std::unordered_map<std::string, size_t> sample_to_index;
    for (size_t i = 0 ; i < 3 ; i++) {
        sample_to_index.emplace("sample" + std::to_string(i+1), i);
    }

    // Get feature_to_index
    std::unordered_map<std::string, size_t> feature_to_index;
    for (size_t i = 0 ; i < 3 ; i++) {
        feature_to_index.emplace("allele" + std::to_string(i+1), i);
    }

    // Make the table
    GenotypeTable table(sample_to_index, feature_to_index);

    SECTION("Set and get values") {
        table.set_value_for_sample_and_feature("sample1", "allele1", 0);
        table.set_value_for_sample_and_feature("sample1", "allele2", 3);
        table.set_value_for_sample_and_feature("sample1", "allele3", 1);
        table.set_value_for_sample_and_feature("sample2", "allele1", std::numeric_limits<size_t>::max());
        table.set_value_for_sample_and_feature("sample2", "allele2", 0);
        table.set_value_for_sample_and_feature("sample2", "allele3", 1);
        table.set_value_for_sample_and_feature("sample3", "allele1", 0);
        table.set_value_for_sample_and_feature("sample3", "allele2", 10);
        table.set_value_for_sample_and_feature("sample3", "allele3", 8);

        REQUIRE(table.get_value_for_sample_and_feature("sample1", "allele1") == 0);
        REQUIRE(table.get_value_for_sample_and_feature("sample1", "allele2") == 3);
        REQUIRE(table.get_value_for_sample_and_feature("sample1", "allele3") == 1);
        REQUIRE(table.get_value_for_sample_and_feature("sample2", "allele1") == std::numeric_limits<size_t>::max());
        REQUIRE(table.get_value_for_sample_and_feature("sample2", "allele2") == 0);
        REQUIRE(table.get_value_for_sample_and_feature("sample2", "allele3") == 1);
        REQUIRE(table.get_value_for_sample_and_feature("sample3", "allele1") == 0);
        REQUIRE(table.get_value_for_sample_and_feature("sample3", "allele2") == 10);
        REQUIRE(table.get_value_for_sample_and_feature("sample3", "allele3") == 8);
        
    }
}

TEST_CASE( "GeneExpressionTable with three samples and three genes", "[table]" ) {

    // Get sample_to_index 
    std::unordered_map<std::string, size_t> sample_to_index;
    for (size_t i = 0 ; i < 3 ; i++) {
        sample_to_index.emplace("sample" + std::to_string(i+1), i);
    }

    // Get feature_to_index
    std::unordered_map<std::string, size_t> feature_to_index;
    for (size_t i = 0 ; i < 3 ; i++) {
        feature_to_index.emplace("gene" + std::to_string(i+1), i);
    }

    // Make the table
    GeneExpressionTable table(sample_to_index, feature_to_index);

    SECTION("Set and get values") {
        table.set_value_for_sample_and_feature("sample1", "gene1", 0.1);
        table.set_value_for_sample_and_feature("sample1", "gene2", 3.0);
        table.set_value_for_sample_and_feature("sample1", "gene3", 1.34);
        table.set_value_for_sample_and_feature("sample2", "gene1", std::numeric_limits<double>::max());
        table.set_value_for_sample_and_feature("sample2", "gene2", 0.0);
        table.set_value_for_sample_and_feature("sample2", "gene3", 1.0932);
        table.set_value_for_sample_and_feature("sample3", "gene1", 0.1234);
        table.set_value_for_sample_and_feature("sample3", "gene2", 10.509345);
        table.set_value_for_sample_and_feature("sample3", "gene3", 8.234234);

        REQUIRE(table.get_value_for_sample_and_feature("sample1", "gene1") == 0.1);
        REQUIRE(table.get_value_for_sample_and_feature("sample1", "gene2") == 3.0);
        REQUIRE(table.get_value_for_sample_and_feature("sample1", "gene3") == 1.34);
        REQUIRE(table.get_value_for_sample_and_feature("sample2", "gene1") == std::numeric_limits<double>::max());
        REQUIRE(table.get_value_for_sample_and_feature("sample2", "gene2") == 0.0);
        REQUIRE(table.get_value_for_sample_and_feature("sample2", "gene3") == 1.0932);
        REQUIRE(table.get_value_for_sample_and_feature("sample3", "gene1") == 0.1234);
        REQUIRE(table.get_value_for_sample_and_feature("sample3", "gene2") == 10.509345);
        REQUIRE(table.get_value_for_sample_and_feature("sample3", "gene3") == 8.234234);
        
    }
}
TEST_CASE( "CovariateTable with three samples and three covariates", "[table]" ) {

    // Get sample_to_index 
    std::unordered_map<std::string, size_t> sample_to_index;
    for (size_t i = 0 ; i < 3 ; i++) {
        sample_to_index.emplace("sample" + std::to_string(i+1), i);
    }

    // Get feature_to_index
    std::unordered_map<std::string, size_t> feature_to_index;
    for (size_t i = 0 ; i < 3 ; i++) {
        feature_to_index.emplace("covar" + std::to_string(i+1), i);
    }

    // Make the table
    CovariateTable table(sample_to_index, feature_to_index);

    SECTION("Set and get values") {
        table.set_value_for_sample_and_feature("sample1", "covar1", 0.1);
        table.set_value_for_sample_and_feature("sample1", "covar2", 3.0);
        table.set_value_for_sample_and_feature("sample1", "covar3", 1.34);
        table.set_value_for_sample_and_feature("sample2", "covar1", std::numeric_limits<double>::max());
        table.set_value_for_sample_and_feature("sample2", "covar2", 0.0);
        table.set_value_for_sample_and_feature("sample2", "covar3", 1.0932);
        table.set_value_for_sample_and_feature("sample3", "covar1", 0.1234);
        table.set_value_for_sample_and_feature("sample3", "covar2", 10.509345);
        table.set_value_for_sample_and_feature("sample3", "covar3", 8.234234);

        REQUIRE(table.get_value_for_sample_and_feature("sample1", "covar1") == 0.1);
        REQUIRE(table.get_value_for_sample_and_feature("sample1", "covar2") == 3.0);
        REQUIRE(table.get_value_for_sample_and_feature("sample1", "covar3") == 1.34);
        REQUIRE(table.get_value_for_sample_and_feature("sample2", "covar1") == std::numeric_limits<double>::max());
        REQUIRE(table.get_value_for_sample_and_feature("sample2", "covar2") == 0.0);
        REQUIRE(table.get_value_for_sample_and_feature("sample2", "covar3") == 1.0932);
        REQUIRE(table.get_value_for_sample_and_feature("sample3", "covar1") == 0.1234);
        REQUIRE(table.get_value_for_sample_and_feature("sample3", "covar2") == 10.509345);
        REQUIRE(table.get_value_for_sample_and_feature("sample3", "covar3") == 8.234234);
        
    }
}

