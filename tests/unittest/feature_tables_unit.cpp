#include <catch.hpp>
#include "../../src/feature_tables.hpp"

using namespace stoat;

class TestBinaryPhenotypeTable : BinaryPhenotypeTable { 
    public:
    TestBinaryPhenotypeTable(const std::unordered_map<std::string, size_t>& sample_to_index) : BinaryPhenotypeTable(sample_to_index){}
    using BinaryPhenotypeTable::set_value_for_sample;
    using BinaryPhenotypeTable::get_value_for_sample;
    using BinaryPhenotypeTable::values_per_sample;
};
class TestQuantitativePhenotypeTable : QuantitativePhenotypeTable { 
    public:
    TestQuantitativePhenotypeTable(const std::unordered_map<std::string, size_t>& sample_to_index) : QuantitativePhenotypeTable(sample_to_index){}
    using QuantitativePhenotypeTable::set_value_for_sample;
    using QuantitativePhenotypeTable::get_value_for_sample;
    using QuantitativePhenotypeTable::values_per_sample;
};
class TestGenotypeTable : GenotypeTable { 
    public:
    TestGenotypeTable(const std::unordered_map<std::string, size_t>& sample_to_index, size_t allele_count) : GenotypeTable(sample_to_index, allele_count){}
    using GenotypeTable::get_count_for_sample_and_allele;
    using GenotypeTable::increment_count;
    using GenotypeTable::values_per_sample;
};
class TestGeneExpressionTable : GeneExpressionTable { 
    public:
    TestGeneExpressionTable(const std::unordered_map<std::string, size_t>& sample_to_index, 
                      const std::unordered_map<std::string, size_t>& feature_to_index) : GeneExpressionTable(sample_to_index, feature_to_index){}
    using GeneExpressionTable::set_value_for_sample_and_feature;
    using GeneExpressionTable::get_value_for_sample_and_feature;
    using GeneExpressionTable::values_per_sample;
};
class TestCovariateTable : CovariateTable { 
    public:
    TestCovariateTable(const std::unordered_map<std::string, size_t>& sample_to_index, 
                      const std::unordered_map<std::string, size_t>& feature_to_index) : CovariateTable(sample_to_index, feature_to_index){}
    using CovariateTable::set_value_for_sample_and_feature;
    using CovariateTable::get_value_for_sample_and_feature;
    using CovariateTable::values_per_sample;
};

TEST_CASE( "BinaryPhenotypeTable with three samples", "[table]" ) {

    // Get sample_to_index 
    std::unordered_map<std::string, size_t> sample_to_index;
    for (size_t i = 0 ; i < 3 ; i++) {
        sample_to_index.emplace("sample" + std::to_string(i+1), i);
    }

    // Make the table
    TestBinaryPhenotypeTable table(sample_to_index);

    SECTION("Table is the right size") {
        REQUIRE(table.values_per_sample.size() == 3);
    }

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
    TestQuantitativePhenotypeTable table(sample_to_index);

    SECTION("Table is the right size") {
        REQUIRE(table.values_per_sample.size() == 5);
    }

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

    // Make the table
    TestGenotypeTable table(sample_to_index, 3);


    SECTION("Table is the right size") {
        REQUIRE(table.values_per_sample.size() == 3);
        for (size_t i = 3 ; i < table.values_per_sample.size() ; i++) {
            REQUIRE(table.values_per_sample.at(i).size() == 3);
        }
    }

    SECTION("Set and get values") {

        table.increment_count("sample1", 1);
        table.increment_count("sample1", 1);
        table.increment_count("sample1", 1);

        table.increment_count("sample1", 2);

        table.increment_count("sample2", 0);
        table.increment_count("sample2", 0);

        table.increment_count("sample2", 1);

        table.increment_count("sample2", 2);
        table.increment_count("sample2", 2);

        table.increment_count("sample3", 2);


        REQUIRE(table.get_count_for_sample_and_allele("sample1", 0) == 0);
        REQUIRE(table.get_count_for_sample_and_allele("sample1", 1) == 3);
        REQUIRE(table.get_count_for_sample_and_allele("sample1", 2) == 1);
        REQUIRE(table.get_count_for_sample_and_allele("sample2", 0) == 2);
        REQUIRE(table.get_count_for_sample_and_allele("sample2", 1) == 1);
        REQUIRE(table.get_count_for_sample_and_allele("sample2", 2) == 2);
        REQUIRE(table.get_count_for_sample_and_allele("sample3", 0) == 0);
        REQUIRE(table.get_count_for_sample_and_allele("sample3", 1) == 0);
        REQUIRE(table.get_count_for_sample_and_allele("sample3", 2) == 1);
        
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
    TestGeneExpressionTable table(sample_to_index, feature_to_index);

    SECTION("Table is the right size") {
        REQUIRE(table.values_per_sample.size() == 3);
        for (size_t i = 3 ; i < table.values_per_sample.size() ; i++) {
            REQUIRE(table.values_per_sample.at(i).size() == 3);
        }
    }
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
    TestCovariateTable table(sample_to_index, feature_to_index);

    SECTION("Table is the right size") {
        REQUIRE(table.values_per_sample.size() == 3);
        for (size_t i = 3 ; i < table.values_per_sample.size() ; i++) {
            REQUIRE(table.values_per_sample.at(i).size() == 3);
        }
    }
    SECTION("Set and get values") {
        table.set_value_for_sample_and_feature("sample1", "covar1", 0.1);
        table.set_value_for_sample_and_feature("sample1", "covar2", 3.0);
        table.set_value_for_sample_and_feature("sample1", "covar3", 1.34);
        table.set_value_for_sample_and_feature("sample2", "covar1", std::numeric_limits<double>::max());
        table.set_value_for_sample_and_feature("sample2", "covar3", 1.0932);
        table.set_value_for_sample_and_feature("sample3", "covar1", 0.1234);
        table.set_value_for_sample_and_feature("sample3", "covar2", 10.509345);
        table.set_value_for_sample_and_feature("sample3", "covar3", 8.234234);

        REQUIRE(table.get_value_for_sample_and_feature("sample1", "covar1") == 0.1);
        REQUIRE(table.get_value_for_sample_and_feature("sample1", "covar2") == 3.0);
        REQUIRE(table.get_value_for_sample_and_feature("sample1", "covar3") == 1.34);
        REQUIRE(table.get_value_for_sample_and_feature("sample2", "covar1") == std::numeric_limits<double>::max());
        REQUIRE(table.get_value_for_sample_and_feature("sample2", "covar2") == std::numeric_limits<double>::max());
        REQUIRE(table.get_value_for_sample_and_feature("sample2", "covar3") == 1.0932);
        REQUIRE(table.get_value_for_sample_and_feature("sample3", "covar1") == 0.1234);
        REQUIRE(table.get_value_for_sample_and_feature("sample3", "covar2") == 10.509345);
        REQUIRE(table.get_value_for_sample_and_feature("sample3", "covar3") == 8.234234);
        
    }
}

TEST_CASE( "Combined GenotypeTable", "[table]" ) {


    //     A0 A1 A2 A3
    // S1 {1, 0, 1, 1}
    // S2 {1, 0, 1, 1}
    // S3 {1, 0, 0, 0}
    // S4 {1, 0, 0, 0}
    // S5 {0, 0, 1, 1}
    // S6 {0, 0, 0, 0}

    std::unordered_map<std::string, size_t> sample_to_index = {{"S1", 0}, {"S2", 1}, {"S3", 2},
                                                               {"S4", 3}, {"S5", 4}, {"S6", 5}};
    stoat::GenoTable table(sample_to_index, 4);
    table.increment_count(0, 0);
    table.increment_count(0, 2);
    table.increment_count(0, 3);
    table.increment_count(1, 0);
    table.increment_count(1, 2);
    table.increment_count(1, 3);
    table.increment_count(2, 0);
    table.increment_count(3, 0);
    table.increment_count(4, 2);
    table.increment_count(4, 3);


    SECTION("Table is the right size") {
        REQUIRE(table.get_n_active_alleles() == 4);
        REQUIRE(table.get_n_active_samples() == 6);
    }

    // add covariates
    stoat::BinaryPhenotypeTable pheno(sample_to_index);
    pheno.set_value_for_sample("S1", 1);
    pheno.set_value_for_sample("S2", 0);
    pheno.set_value_for_sample("S3", 0);
    pheno.set_value_for_sample("S4", 0);
    pheno.set_value_for_sample("S5", 1);
    pheno.set_value_for_sample("S6", 1);
    
    std::unordered_map<std::string, size_t> covar_to_index = {{"covar1", 0}};
    stoat::CovariateTable covar(sample_to_index, covar_to_index);
    covar.set_value_for_sample_and_feature("S1", "covar1", 0.2);
    covar.set_value_for_sample_and_feature("S2", "covar1", 10);
    covar.set_value_for_sample_and_feature("S3", "covar1", 2);
    covar.set_value_for_sample_and_feature("S4", "covar1", 5);
    covar.set_value_for_sample_and_feature("S5", "covar1", 11);
    covar.set_value_for_sample_and_feature("S6", "covar1", 1);

    // link to them in the GenoTable
    table.link_to_binary_phenotype(pheno);
    table.link_to_covariates(covar);

    SECTION("Table is the right size with covariates") {
        REQUIRE(table.get_n_active_alleles() == 4);
        REQUIRE(table.get_n_active_columns() == 5);
        REQUIRE(table.get_n_active_samples() == 6);
    }

    SECTION("Accessor works on alleles and covariates") {
        REQUIRE(table.get_value(0, 2) == 1);
        REQUIRE(table.get_value(1, 4) == 10);
    }

    SECTION("Remove noncovered samples") {
        table.remove_noncovered_samples();
        REQUIRE(table.get_n_active_samples() == 5);
    }

    SECTION("Remove constant allele") {
        table.remove_constant_predictors();
        REQUIRE(table.get_n_active_alleles() == 3);
        REQUIRE(table.get_n_active_columns() == 4);
    }

    SECTION("Remove first allele") {
        table.remove_one_allele();
        REQUIRE(table.get_n_active_alleles() == 3);
        REQUIRE(table.get_n_active_columns() == 4);
    }

    SECTION("Remove duplicated predictor") {
        table.remove_duplicated_predictors();
        REQUIRE(table.get_n_active_alleles() == 3);
        REQUIRE(table.get_n_active_columns() == 4);
    }

    SECTION("Add new total allele count column") {
        Eigen::MatrixXd mat = table.make_matrixXd_features();
        REQUIRE(mat.cols() == 6);
        table.add_total_allele_count_covariable();
        mat = table.make_matrixXd_features();
        REQUIRE(mat.cols() == 7);
    }

    SECTION("Make the correct feature matrix, no filtering but with total allele counts") {
        table.add_total_allele_count_covariable();
        Eigen::MatrixXd mat = table.make_matrixXd_features();
        // intercept
        REQUIRE(mat(2,0) == 1);
        // alleles
        REQUIRE(mat(4,3) == 1);
        REQUIRE(mat(4,4) == 1);
        // covariates
        REQUIRE(mat(4,5) == 11);
        REQUIRE(mat(2,5) == 2);
        // total allele count
        REQUIRE(mat(0,6) == 3);
        REQUIRE(mat(2,6) == 1);
    }

    SECTION("Make the correct feature matrix, after filtering") {
        table.remove_noncovered_samples();
        table.remove_constant_predictors();
        table.remove_duplicated_predictors();
        table.remove_one_allele();
        Eigen::MatrixXd mat = table.make_matrixXd_features();
        // alleles becomes
        //    A2
        // S1 1
        // S2 1
        // S3 0
        // S4 0
        // S5 1

        // intercept
        REQUIRE(mat(2,0) == 1);
        // alleles
        REQUIRE(mat(2,1) == 0);
        REQUIRE(mat(4,1) == 1);
    }

    SECTION("Samples were masked also in the phenotype") {
        Eigen::VectorXd Y = table.make_vectorxd_phenotype();
        REQUIRE(Y.rows() == 6);
        table.remove_noncovered_samples();
        Y = table.make_vectorxd_phenotype();
        REQUIRE(Y.rows() == 5);
        REQUIRE(Y(0) == 1);
        REQUIRE(Y(1) == 0);
    }

    SECTION("Check that the phenotype is correctly matched") {
        //     A0 A1
        // S1 {1, 0}
        // S2 {1, 0}
        // S3 {1, 0}
        // S4 {0, 1}
        // S5 {0, 1}

        std::unordered_map<std::string, size_t> sample_to_index = {{"S1", 0}, {"S2", 1}, {"S3", 2},
                                                                   {"S4", 3}, {"S5", 4}};
        stoat::GenoTable geno(sample_to_index, 3);
        geno.increment_count(0, 0);
        geno.increment_count(1, 0);
        geno.increment_count(2, 0);
        geno.increment_count(3, 1);
        geno.increment_count(4, 1);

        stoat::QuantitativePhenotypeTable pheno(sample_to_index);
        pheno.set_value_for_sample("S1", 11);
        pheno.set_value_for_sample("S2", 10.1);
        pheno.set_value_for_sample("S3", 5.2);
        pheno.set_value_for_sample("S4", -0.3);
        pheno.set_value_for_sample("S5", 2);

        geno.link_to_quantitative_phenotype(pheno);

        Eigen::VectorXd Y = geno.make_vectorxd_phenotype();
        REQUIRE(Y(0) == 11);
        REQUIRE(Y(3) == -0.3);        
    }

    SECTION("Check that filters work") {
        //     A0 A1
        // S1 {1, 0}
        // S2 {1, 0}
        // S3 {1, 0}
        // S4 {0, 1}
        // S5 {0, 1}

        std::unordered_map<std::string, size_t> sample_to_index = {{"S1", 0}, {"S2", 1}, {"S3", 2},
                                                                   {"S4", 3}, {"S5", 4}};
        stoat::GenoTable geno(sample_to_index, 3);
        geno.increment_count(0, 0);
        geno.increment_count(1, 0);
        geno.increment_count(2, 0);
        geno.increment_count(3, 1);
        geno.increment_count(4, 1);

        stoat::QuantitativePhenotypeTable pheno(sample_to_index);
        pheno.set_value_for_sample("S1", 11);
        pheno.set_value_for_sample("S2", 10.1);
        pheno.set_value_for_sample("S3", 5.2);
        pheno.set_value_for_sample("S4", -0.3);
        pheno.set_value_for_sample("S5", 2);

        geno.link_to_quantitative_phenotype(pheno);

        // remove non-variable allele, e.g. absent in both groups
        geno.remove_noncovered_samples();    
        geno.remove_constant_predictors();
        
        REQUIRE(geno.passes_filters(0, 0));
        REQUIRE(geno.passes_filters(.2, 5));
        REQUIRE(!geno.passes_filters(0.5, 0));
        REQUIRE(!geno.passes_filters(0, 6));
    }

}

