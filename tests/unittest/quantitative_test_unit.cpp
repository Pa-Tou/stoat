#include <catch.hpp>

#include "../../src/quantitative_table.hpp"
#include "../../src/stats_test.hpp"
#include "../../src/arg_parser.hpp"  // Qtl_data

using namespace stoat_vcf;

TEST_CASE("Quantitative table creation") {
    SECTION("Creation 1") {

    }
}

TEST_CASE("Quantitative table filtration") {
    SECTION("Remove empty columns quantitative table") {
        
        std::vector<std::vector<double>> df = {
            {0.5, 0, 0.5},
            {0.5, 0, 0.5},
            {1, 0, 0},
            {1, 0, 0},
            {0, 0, 1}
        };
        
        stoat_vcf::remove_empty_columns_quantitative_table(df);

        std::vector<std::vector<double>> df_expected = {
            {0.5, 0.5},
            {0.5, 0.5},
            {1, 0},
            {1, 0},
            {0, 1}
        };
        REQUIRE(df == df_expected);
    }

    SECTION("MAF filtration") {

        std::vector<std::vector<double>> df = {
            {0.5, 0},
            {0.5, 0},
            {1, 0},
            {1, 0},
            {1, 0},
            {1, 0},
            {1, 0},
            {1, 0},
            {1, 0},
            {0, 1},
        };

        double maf_threshold = 0.01;
        bool filtration = stoat_vcf::filtration_quantitative_table(df, 3, 5, maf_threshold);
        REQUIRE(filtration == false);

        double maf_threshold2 = 0.2;
        bool filtration2 = stoat_vcf::filtration_quantitative_table(df, 3, 5, maf_threshold2);
        REQUIRE(filtration2 == true);
    }

    SECTION("Remove last columns quantitative table") {
        
        std::vector<std::vector<double>> df = {
            {0.5, 0.5},
            {0.5, 0.5},
            {1, 0},
            {1, 0},
            {0, 1}
        };
        
        stoat_vcf::remove_last_columns_quantitative_table(df);

        std::vector<std::vector<double>> df_expected = {
            {0.5},
            {0.5},
            {1},
            {1},
            {0}
        };
        REQUIRE(df == df_expected);
    }
}

// TEST_CASE("Linear Regression Test without cov", "[linear_regression]") {
//     SECTION("Linear Regression 1 - Perfect Linear Relationship") {

//         std::vector<std::vector<double>> df = {
//             {0},
//             {1},
//             {0}
//         };

//         std::vector<double> quantitative_phenotype = {2.0, 4.0, 6.0};
//         std::vector<std::vector<double>> covar;  // No covariates

//         LinearRegression lr;
//         auto [r2, beta, se, p_value] = lr.linear_regression(df, quantitative_phenotype, covar);

//         INFO("se = " << se);
//         INFO("beta = " << beta);
//         INFO("p_value = " << p_value);
//         INFO("r2 = " << r2);

//         REQUIRE(se == "3.464");
//         REQUIRE(beta == "1.332e-15");
//         REQUIRE(p_value == "1.0000");
//         REQUIRE(r2 == "0.0000");
//     }

//     SECTION("Linear Regression 2 - Moderate") {

//         std::vector<std::vector<double>> df = {
//             {0.5, 0},
//             {0, 0.5},
//             {1, 0},
//             {0, 1},
//             {0, 0.5}
//         };

//         std::vector<double> quantitative_phenotype = {10.5, 13.0, 15.8, 19.7, 21.5};
//         std::vector<std::vector<double>> covar;
        
//         LinearRegression lr;
//         auto [r2, beta, se, p_value] = lr.linear_regression(df, quantitative_phenotype, covar);

//         INFO("se = " << se);
//         INFO("beta = " << beta);
//         INFO("p_value = " << p_value);
//         INFO("r2 = " << r2);

//         REQUIRE(std::stod(se) == 9.762);
//         REQUIRE(std::stod(beta) == 9.7000);
//         REQUIRE(std::stod(p_value) == 0.425);
//         REQUIRE(std::stod(r2) == 0.427);
//     }

//     SECTION("Linear Regression 3 - Weaker Correlation") {

//         std::vector<std::vector<double>> df = {
//             {1, 0},
//             {1, 0},
//             {1, 0},
//             {1, 0},
//             {1, 0},
//             {1, 0},
//             {1, 0},
//             {0, 1},
//             {0, 0},

//         };

//         std::vector<double> quantitative_phenotype = {4.5, 7.0, 9.2, 10.9, 13.0, 14.0, 11.0, 15.0, 16.0};
//         std::vector<std::vector<double>> covar;  // No covariates
//         LinearRegression lr;
//         auto [r2, beta, se, p_value] = lr.linear_regression(df, quantitative_phenotype, covar);
        
//         INFO("se = " << se);
//         INFO("beta = " << beta);
//         INFO("p_value = " << p_value);
//         INFO("r2 = " << r2);

//         REQUIRE(std::stod(se) == 3.564);
//         REQUIRE(std::stod(beta) == -6.0571);
//         REQUIRE(std::stod(p_value) == 0.140);
//         REQUIRE(std::stod(r2) == 0.421);
//     }
// }

// TEST_CASE("Linear Regression Test with covariates", "[linear_regression]") {
//     SECTION("Linear Regression 1 - Perfect Linear Relationship with Covariate") {

//         std::vector<std::vector<double>> df = {
//             {0},
//             {1},
//             {0}
//         };

//         std::vector<double> quantitative_phenotype = {2.0, 4.0, 6.0};

//         std::vector<std::vector<double>> covar = {
//             {1.0}, {1.0}, {1.0}  // Intercept-only
//         };

//         LinearRegression lr;
//         auto [r2, beta, se, p_value] = lr.linear_regression(df, quantitative_phenotype, covar);

//         INFO("se = " << se);
//         INFO("beta = " << beta);
//         INFO("p_value = " << p_value);
//         INFO("r2 = " << r2);

//         REQUIRE(se == "NA");
//         REQUIRE(beta == "-8.0000");
//         REQUIRE(p_value == "NA");
//         REQUIRE(r2 == "1.0000");
//     }

//     SECTION("Linear Regression 2 - Moderate with Covariate") {

//         std::vector<std::vector<double>> df = {
//             {0.5, 0},
//             {0, 0.5},
//             {1, 0},
//             {0, 1},
//             {0, 0.5}
//         };

//         std::vector<double> quantitative_phenotype = {10.5, 13.0, 15.8, 19.7, 21.5};

//         std::vector<std::vector<double>> covar = {
//             {1.0}, {1.0}, {1.0}, {1.0}, {1.0}
//         };

//         LinearRegression lr;
//         auto [r2, beta, se, p_value] = lr.linear_regression(df, quantitative_phenotype, covar);

//         INFO("se = " << se);
//         INFO("beta = " << beta);
//         INFO("p_value = " << p_value);
//         INFO("r2 = " << r2);

//         REQUIRE(se == "0.880");
//         REQUIRE(beta == "-17.4400");
//         REQUIRE(p_value == "0.0320");
//         REQUIRE(r2 == "0.999");
//     }

//     SECTION("Linear Regression 3 - Weaker Correlation with Covariate") {

//         std::vector<std::vector<double>> df = {
//             {1, 0},
//             {1, 0},
//             {1, 0},
//             {1, 0},
//             {1, 0},
//             {1, 0},
//             {1, 0},
//             {0, 1},
//             {0, 0},
//         };

//         std::vector<double> quantitative_phenotype = {
//             4.5, 7.0, 9.2, 10.9, 13.0, 14.0, 11.0, 15.0, 16.0
//         };

//         std::vector<std::vector<double>> covar = {
//             {1.0}, {1.0}, {1.0}, {1.0}, {1.0},
//             {1.0}, {1.0}, {1.0}, {1.0}
//         };

//         LinearRegression lr;
//         auto [r2, beta, se, p_value] = lr.linear_regression(df, quantitative_phenotype, covar);

//         INFO("se = " << se);
//         INFO("beta = " << beta);
//         INFO("p_value = " << p_value);
//         INFO("r2 = " << r2);

//         REQUIRE(se == "3.033");
//         REQUIRE(beta == "6.5878");
//         REQUIRE(p_value == "0.0730");
//         REQUIRE(r2 == "0.4210");
//     }
// }
