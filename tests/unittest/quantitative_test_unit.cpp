#include <catch.hpp>

#include "../../src/quantitative_table.hpp"
#include "../../src/stats_test.hpp"
#include "../../src/arg_parser.hpp"  // Qtl_data

using namespace stoat_vcf;

TEST_CASE("Quantitative table creation") {
    SECTION("Creation 1") {

    }
}

TEST_CASE("Quantitative table modification") {
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

    SECTION("Combine identical columns quantitative table") {
        
        std::vector<std::vector<double>> df = {
            {0.5, 0.5, 0.5},
            {0.5, 0.5, 0.5},
            {1, 0, 1},
            {1, 0, 1},
            {0, 1, 0}
        };
        
        stoat_vcf::combine_identical_columns_quantitative_table(df);

        std::vector<std::vector<double>> df_expected = {
            {1.0, 0.5},
            {1.0, 0.5},
            {2.0, 0},
            {2.0, 0},
            {0.0, 1.0}
        };
        REQUIRE(df == df_expected);
    }

    SECTION("Check last columns quantitative table") {
        std::vector<std::vector<double>> df = {
            {0.5},
            {0.5},
            {0.5},
            {0.5},
            {0.5}
        };
        
        bool to_filter = stoat_vcf::check_last_columns_quantitative_table(df);
        REQUIRE(to_filter == true);

        std::vector<std::vector<double>> df2 = {
            {0.5},
            {0.5},
            {1.0},
            {0.5},
            {0.5}
        };
        
        bool to_filter2 = stoat_vcf::check_last_columns_quantitative_table(df2);
        REQUIRE(to_filter2 == false);
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

TEST_CASE("Linear Regression Test without cov", "[linear_regression]") {
    SECTION("Linear Regression 1 path - Perfect Linear Relationship") {

        std::vector<std::vector<double>> df = {
            {0.5},
            {1},
            {0.5}
        };

        std::vector<double> quantitative_phenotype = {39.83, 49.92, 34.56};
        std::vector<std::vector<double>> covar;  // No covariates

        stoat::LinearRegression lr;
        auto [p_value, beta, se, r2] = lr.linear_regression(df, quantitative_phenotype, covar);

        INFO("se = " << se);
        INFO("beta = " << beta);
        INFO("p_value = " << p_value);
        INFO("r2 = " << r2);

        REQUIRE(se == "9.128");
        REQUIRE(beta == "25.45");
        REQUIRE(p_value == "0.2192");
        REQUIRE(r2 == "0.886");
    }

    SECTION("Linear Regression 2 paths - Moderate") {

        std::vector<std::vector<double>> df = {
            {0.5, 0},
            {0, 0.5},
            {1, 0},
            {0, 1},
            {0, 0.5}
        };

        std::vector<double> quantitative_phenotype = {10.5, 13.0, 15.8, 19.7, 21.5};
        std::vector<std::vector<double>> covar;
        
        stoat::LinearRegression lr;
        auto [p_value, beta, se, r2] = lr.linear_regression(df, quantitative_phenotype, covar);

        INFO("se = " << se);
        INFO("beta = " << beta);
        INFO("p_value = " << p_value);
        INFO("r2 = " << r2);

        REQUIRE(std::stod(se) == 9.131);
        REQUIRE(std::stod(beta) == 3.4);
        REQUIRE(std::stod(p_value) == 0.7454);
        REQUIRE(std::stod(r2) == 0.4272);
    }

    SECTION("Linear Regression speudo inversion 3 paths - Moderate") {

        std::vector<std::vector<double>> df = {
            {0, 0.5, 0.5},
            {0.5, 0, 0},
            {0.25, 0.5, 0},
            {0.5, 0, 0},
            {0.333333, 0.333333, 0}
        };

        std::vector<double> quantitative_phenotype = {39.83, 30.24, 49, 72, 54.36};
        std::vector<std::vector<double>> covar;
        
        stoat::LinearRegression lr;
        auto [p_value, beta, se, r2] = lr.linear_regression(df, quantitative_phenotype, covar);

        INFO("se = " << se);
        INFO("beta = " << beta);
        INFO("p_value = " << p_value);
        INFO("r2 = " << r2);

        REQUIRE(std::stod(se) == 24.56);
        REQUIRE(std::stod(beta) == 21.5);
        REQUIRE(std::stod(p_value) == 0.5422);
        REQUIRE(std::stod(r2) == 0.108);
    }
}

TEST_CASE("Linear Regression Test with covariates", "[linear_regression]") {
    SECTION("Linear Regression 1 path - Perfect Linear Relationship with Covariate") {

        std::vector<std::vector<double>> df = {
            {0.5},
            {1},
            {0.5},
            {1}
        };

        std::vector<double> quantitative_phenotype = {39.83, 49.92, 34.56, 33.46};

        std::vector<std::vector<double>> covar = {
            {1.2}, {0.8}, {1.0}, {0.9}
        };

        stoat::LinearRegression lr;
        auto [p_value, beta, se, r2] = lr.linear_regression(df, quantitative_phenotype, covar);

        INFO("se = " << se);
        INFO("beta = " << beta);
        INFO("p_value = " << p_value);
        INFO("r2 = " << r2);

        REQUIRE(se == "45.19");
        REQUIRE(beta == "3.07");
        REQUIRE(p_value == "0.9568");
        REQUIRE(r2 == "0.1398");
    }

    SECTION("Linear Regression 2 paths - Moderate with Covariate") {

        std::vector<std::vector<double>> df = {
            {0.5, 0},
            {0, 0.5},
            {1, 0},
            {0, 1},
            {0, 0.5}
        };

        std::vector<double> quantitative_phenotype = {10.5, 13.0, 15.8, 19.7, 21.5};

        std::vector<std::vector<double>> covar = {
            {1.1}, {1.2}, {0.8}, {0.9}, {1.0}
        };

        stoat::LinearRegression lr;
        auto [p_value, beta, se, r2] = lr.linear_regression(df, quantitative_phenotype, covar);

        INFO("se = " << se);
        INFO("beta = " << beta);
        INFO("p_value = " << p_value);
        INFO("r2 = " << r2);

        REQUIRE(se == "7.447");
        REQUIRE(beta == "-20.56");
        REQUIRE(p_value == "0.2212");
        REQUIRE(r2 == "0.9595");

        // se = 7.447
        // beta = -20.56
        // p_value = 0.2212
        // r2 = 0.9595

    }

    SECTION("Linear Regression speudo inversion 3 paths - Moderate with Covariate") {

        std::vector<std::vector<double>> df = {
            {0, 0.5, 0.5},
            {0.5, 0, 0},
            {0.25, 0.5, 0},
            {0.5, 0, 0},
            {0.333333, 0.333333, 0}
        };

        std::vector<double> quantitative_phenotype = {39.83, 30.24, 49, 72, 54.36};
        std::vector<std::vector<double>> covar = {
            {1.1}, {1.2}, {0.8}, {0.9}, {1.0}
        };
        
        stoat::LinearRegression lr;
        auto [p_value, beta, se, r2] = lr.linear_regression(df, quantitative_phenotype, covar);

        INFO("se = " << se);
        INFO("beta = " << beta);
        INFO("p_value = " << p_value);
        INFO("r2 = " << r2);

        REQUIRE(std::stod(se) == 22.85);
        REQUIRE(std::stod(beta) == 52.34);
        REQUIRE(std::stod(p_value) == 0.2621);
        REQUIRE(std::stod(r2) == 0.7563);

        // se = 22.85
        // beta = 52.34
        // p_value = 0.2621
        // r2 = 0.7563
    }
}
