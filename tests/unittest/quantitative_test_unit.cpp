#include <catch.hpp>

#include "../../src/quantitative_table.hpp"
#include "../../src/stats_test.hpp"
#include "../../src/arg_parser.hpp"

using namespace stoat_vcf;

TEST_CASE("Quantitative table modification") {
    SECTION("Remove empty columns quantitative table") {

        //     A0 A1 A2
        // S1 {1, 0, 1}
        // S2 {1, 0, 1}
        // S3 {1, 0, 0}
        // S4 {1, 0, 0}
        // S5 {0, 0, 1}

        std::unordered_map<std::string, size_t> sample_to_index = {{"S1", 0}, {"S2", 1}, {"S3", 2},
                                                                   {"S4", 3}, {"S5", 4}};
        stoat::GenoTable geno(sample_to_index, 3);
        geno.increment_count(0, 0);
        geno.increment_count(0, 2);
        geno.increment_count(1, 0);
        geno.increment_count(1, 2);
        geno.increment_count(2, 0);
        geno.increment_count(3, 0);
        geno.increment_count(4, 2);
        
        // only 3 alleles originally
        REQUIRE(geno.get_n_active_alleles() == 3);
        // second allele should have an allele count of 0
        Eigen::MatrixXd mat = geno.make_matrixXd_features();
        REQUIRE(mat(0, 2) == 0);
        // note: second allele is at index 2 because of the added intercept column in the matrix

        geno.remove_constant_predictors();

        // only 2 alleles now
        REQUIRE(geno.get_n_active_alleles() == 2);
        // second allele should have an allele count of 1 now
        mat = geno.make_matrixXd_features();
        REQUIRE(mat(0, 2) == 1);
    }

    SECTION("Combine identical columns quantitative table") {
        
        //     A0 A1 A2
        // S1 {0, 1, 1}
        // S2 {1, 0, 0}
        // S3 {1, 0, 0}
        // S4 {1, 0, 0}
        // S5 {0, 1, 1}

        std::unordered_map<std::string, size_t> sample_to_index = {{"S1", 0}, {"S2", 1}, {"S3", 2},
                                                                   {"S4", 3}, {"S5", 4}};
        stoat::GenoTable geno(sample_to_index, 3);
        geno.increment_count(0, 1);
        geno.increment_count(0, 2);
        geno.increment_count(1, 0);
        geno.increment_count(2, 0);
        geno.increment_count(3, 0);
        geno.increment_count(4, 1);
        geno.increment_count(4, 2);
        
        // 3 alleles originally
        REQUIRE(geno.get_n_active_alleles() == 3);
        
        geno.remove_duplicated_predictors();

        // only 2 alleles now
        REQUIRE(geno.get_n_active_alleles() == 2);
    }

    SECTION("Filters monoallelic variants") {
        //     A0
        // S1 {1}
        // S2 {1}
        // S3 {1}
        // S4 {1}
        // S5 {1}

        std::unordered_map<std::string, size_t> sample_to_index = {{"S1", 0}, {"S2", 1}, {"S3", 2},
                                                                   {"S4", 3}, {"S5", 4}};
        stoat::GenoTable geno(sample_to_index, 1);
        geno.increment_count(0, 0);
        geno.increment_count(1, 0);
        geno.increment_count(2, 0);
        geno.increment_count(3, 0);
        geno.increment_count(4, 0);
        
        REQUIRE(!geno.passes_filters(0, 0));
    }

    SECTION("Remove last columns quantitative table") {
        //     A0 A1 A2
        // S1 {0, 1, 1}
        // S2 {1, 0, 0}
        // S3 {1, 0, 0}
        // S4 {1, 0, 0}
        // S5 {0, 1, 1}

        std::unordered_map<std::string, size_t> sample_to_index = {{"S1", 0}, {"S2", 1}, {"S3", 2},
                                                                   {"S4", 3}, {"S5", 4}};
        stoat::GenoTable geno(sample_to_index, 3);
        geno.increment_count(0, 1);
        geno.increment_count(0, 2);
        geno.increment_count(1, 0);
        geno.increment_count(2, 0);
        geno.increment_count(3, 0);
        geno.increment_count(4, 1);
        geno.increment_count(4, 2);
        
        // 3 alleles originally
        REQUIRE(geno.get_n_active_alleles() == 3);
        
        geno.remove_one_allele();

        // 2 alleles now
        REQUIRE(geno.get_n_active_alleles() == 2);
    }
}

TEST_CASE("Logistic Regression") {
     SECTION("Logistic Regression 2 paths - Perfect Relationship") {
        //     A0 A1
        // S1  {0, 1}
        // S2  {1, 0}
        // S3  {1, 0}
        // S4  {1, 0}
        // S5  {0, 1}
        // S11 {0, 1}
        // S12 {1, 0}
        // S13 {1, 0}
        // S14 {1, 0}
        // S15 {0, 1}

        Eigen::MatrixXd X(10, 2);
        X << 1, 1,
            1, 0, 
            1, 0, 
            1, 0, 
            1, 1, 
            1, 1, 
            1, 0, 
            1, 0, 
            1, 0, 
            1, 1;

        Eigen::VectorXd Y(10);
        Y << 1, 0, 0, 0, 1, 1, 0, 0, 0, 1;
        
        stoat::LogisticRegression lr;
        std::string pv = lr.logistic_regression(X, Y, 1);

        INFO("p_value = " << pv);

        REQUIRE(std::abs(std::stod(pv) - 0.002842742) < 0.01);
    }

    SECTION("Logistic Regression 2 paths - Moderate") {
        //     A0 A1
        // S1 {0, 1}
        // S2 {1, 0}
        // S3 {1, 0}
        // S4 {1, 0}
        // S5 {0, 1}

        Eigen::MatrixXd X(5, 2);
        X << 1, 1,
            1, 0, 
            1, 0, 
            1, 0, 
            1, 1;
        
        Eigen::VectorXd Y(5);
        Y << 0, 1, 0, 0, 1;
        
        stoat::LogisticRegression lr;
        std::string pv = lr.logistic_regression(X, Y, 1);

        INFO("p_value = " << pv);

        REQUIRE(std::stod(pv) == 0.7098);
    }

    SECTION("Logistic Regression 3 paths - Moderate") {

        //     A0 A1 A2
        // S1 {1, 0, 1}
        // S2 {1, 0, 1}
        // S3 {1, 1, 0}
        // S4 {0, 1, 1}
        // S5 {0, 1, 1}

        Eigen::MatrixXd X(5, 3);
        X << 1, 0, 1, 
            1, 0, 1, 
            1, 1, 0, 
            1, 1, 1, 
            1, 1, 1; 
        
        Eigen::VectorXd Y(5);
        Y << 1, 1, 1, 0, 0;
        
        stoat::LogisticRegression lr;
        std::string pv = lr.logistic_regression(X, Y, 2);

        INFO("p_value = " << pv);

        REQUIRE(std::abs(std::stod(pv) - 0.8732353) < .01);
    }
}

TEST_CASE("Linear Regression Test without cov", "[linear_regression]") {
    SECTION("Linear Regression 2 paths - good linear relationship") {
        //     A0 A1
        // S1 {1, 0}
        // S2 {1, 0}
        // S3 {1, 0}
        // S4 {0, 1}
        // S5 {0, 1}

        Eigen::MatrixXd X(5, 2);
        X << 1, 0,
            1, 0, 
            1, 0, 
            1, 1, 
            1, 1;

        Eigen::VectorXd Y(5);
        Y << 11, 10.1, 5.2, -0.3, 2;

        stoat::LinearRegression lr;
        std::string pv = lr.linear_regression(X, Y, 1);

        INFO("p_value = " << pv);

        REQUIRE(std::stod(pv) == 0.0496);
    }

    SECTION("Linear Regression 3 paths - Moderate") {
        //     A0 A1 A2
        // S1 {1, 0, 1}
        // S2 {1, 0, 1}
        // S3 {1, 1, 0}
        // S4 {0, 1, 1}
        // S5 {0, 1, 1}

        Eigen::MatrixXd X(5, 3);
        X << 1, 0, 1,
            1, 0, 1, 
            1, 1, 0, 
            1, 1, 1, 
            1, 1, 1; 

        Eigen::VectorXd Y(5);
        Y << 11, 10.1, 5.2, -0.3, 2;

        stoat::LinearRegression lr;
        std::string pv = lr.linear_regression(X, Y, 2);

        INFO("p_value = " << pv);

        REQUIRE(std::abs(std::stod(pv) - 0.03133) < 0.01);
    }

    SECTION("Linear Regression pseudo inversion 3 paths - Moderate") {
        //     A0 A1 A2
        // S1 {1, 0, 0}
        // S2 {0, 0, 1}
        // S3 {0, 0, 1}
        // S4 {0, 1, 0}
        // S5 {0, 1, 0}

        // we need to add the extra column with the total allele counts JEAN
        Eigen::MatrixXd X(5, 3);
        X << 1, 0, 0,
            1, 0, 1, 
            1, 0, 1, 
            1, 1, 0, 
            1, 1, 0; 

        Eigen::VectorXd Y(5);
        Y << 11, 10.1, 5.2, -0.3, 2;

        stoat::LinearRegression lr;
        std::string pv = lr.linear_regression(X, Y, 2);

        INFO("p_value = " << pv);

        REQUIRE(std::stod(pv) == 0.1505);
    }
}

TEST_CASE("Linear Regression Test with covariates", "[linear_regression]") {
    SECTION("Linear Regression 2 paths - Perfect Linear Relationship with Covariate") {
        //     A0 A1
        // S1 {1, 0}
        // S2 {1, 0}
        // S3 {1, 0}
        // S4 {0, 1}
        // S5 {0, 1}

        // we add the covariate column
        Eigen::MatrixXd X(5, 3);
        X << 1, 0, 0.2,
            1, 0, 10, 
            1, 0, 2, 
            1, 1, 5, 
            1, 1, 11; 

        Eigen::VectorXd Y(5);
        Y << 11, 10.1, 5.2, -0.3, 2;

        stoat::LinearRegression lr;
        std::string pv = lr.linear_regression(X, Y, 1);

        INFO("p_value = " << pv);

        REQUIRE(std::abs(std::stod(pv) - 0.1140625) < .01);
    }

    SECTION("Linear Regression pseudo inversion 3 paths - Moderate with Covariate") {
        //     A0 A1 A2
        // S1 {1, 0, 0}
        // S2 {0, 0, 1}
        // S3 {0, 0, 1}
        // S4 {0, 1, 0}
        // S5 {0, 1, 0}

        Eigen::MatrixXd X(5, 4);
        X << 1, 0, 0, 0.2,
            1, 0, 1, 10,  
            1, 0, 1, 2,   
            1, 1, 0, 5,   
            1, 1, 0, 11;  

        Eigen::VectorXd Y(5);
        Y << 11, 10.1, 5.2, -0.3, 2;

        stoat::LinearRegression lr;
        std::string pv = lr.linear_regression(X, Y, 2);

        INFO("p_value = " << pv);

        REQUIRE(std::abs(std::stod(pv) - 0.08149071) < .0001);
    }
}
