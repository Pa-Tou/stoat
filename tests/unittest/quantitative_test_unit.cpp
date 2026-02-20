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
        // S1 {0, 1}
        // S2 {1, 0}
        // S3 {1, 0}
        // S4 {1, 0}
        // S5 {0, 1}

        std::unordered_map<std::string, size_t> sample_to_index = {{"S1", 0}, {"S2", 1}, {"S3", 2},
                                                                   {"S4", 3}, {"S5", 4},
                                                                   {"S11", 5}, {"S12", 6}, {"S13", 7},
                                                                   {"S14", 8}, {"S15", 9}};
        stoat::GenoTable geno(sample_to_index, 2);
        geno.increment_count(0, 1);
        geno.increment_count(1, 0);
        geno.increment_count(2, 0);
        geno.increment_count(3, 0);
        geno.increment_count(4, 1);
        geno.increment_count(5, 1);
        geno.increment_count(6, 0);
        geno.increment_count(7, 0);
        geno.increment_count(8, 0);
        geno.increment_count(9, 1);

        
        stoat::BinaryPhenotypeTable pheno(sample_to_index);
        pheno.set_value_for_sample("S1", 1);
        pheno.set_value_for_sample("S2", 0);
        pheno.set_value_for_sample("S3", 0);
        pheno.set_value_for_sample("S4", 0);
        pheno.set_value_for_sample("S5", 1);
        pheno.set_value_for_sample("S11", 1);
        pheno.set_value_for_sample("S12", 0);
        pheno.set_value_for_sample("S13", 0);
        pheno.set_value_for_sample("S14", 0);
        pheno.set_value_for_sample("S15", 1);

        std::unordered_map<std::string, size_t> covar_to_index;
        stoat::CovariateTable covar(sample_to_index, covar_to_index);

        stoat::LogisticRegression lr;
        stoat::test_result_t tres = lr.logistic_regression(pheno, geno, covar, 0, 0);

        INFO("p_value = " << tres.pv);

        REQUIRE(std::abs(std::stod(tres.pv) - 0.002842742) < 0.01);
    }

    SECTION("Logistic Regression 2 paths - Moderate") {
        //     A0 A1
        // S1 {0, 1}
        // S2 {1, 0}
        // S3 {1, 0}
        // S4 {1, 0}
        // S5 {0, 1}

        std::unordered_map<std::string, size_t> sample_to_index = {{"S1", 0}, {"S2", 1}, {"S3", 2},
                                                                   {"S4", 3}, {"S5", 4}};
        stoat::GenoTable geno(sample_to_index, 2);
        geno.increment_count(0, 1);
        geno.increment_count(1, 0);
        geno.increment_count(2, 0);
        geno.increment_count(3, 0);
        geno.increment_count(4, 1);

        stoat::BinaryPhenotypeTable pheno(sample_to_index);
        pheno.set_value_for_sample("S1", 0);
        pheno.set_value_for_sample("S2", 1);
        pheno.set_value_for_sample("S3", 0);
        pheno.set_value_for_sample("S4", 0);
        pheno.set_value_for_sample("S5", 1);

        std::unordered_map<std::string, size_t> covar_to_index;
        stoat::CovariateTable covar(sample_to_index, covar_to_index);

        stoat::LogisticRegression lr;
        stoat::test_result_t tres = lr.logistic_regression(pheno, geno, covar, 0, 0);

        INFO("p_value = " << tres.pv);

        REQUIRE(std::stod(tres.pv) == 0.7098);
    }

    SECTION("Logistic Regression 3 paths - Moderate") {

        //     A0 A1 A2
        // S1 {1, 0, 1}
        // S2 {1, 0, 1}
        // S3 {1, 1, 0}
        // S4 {0, 1, 1}
        // S5 {0, 1, 1}

        std::unordered_map<std::string, size_t> sample_to_index = {{"S1", 0}, {"S2", 1}, {"S3", 2},
                                                                   {"S4", 3}, {"S5", 4}};
        stoat::GenoTable geno(sample_to_index, 3);
        geno.increment_count(0, 0);
        geno.increment_count(0, 2);
        geno.increment_count(1, 0);
        geno.increment_count(1, 2);
        geno.increment_count(2, 0);
        geno.increment_count(2, 1);
        geno.increment_count(3, 1);
        geno.increment_count(3, 2);
        geno.increment_count(4, 1);
        geno.increment_count(4, 2);

        stoat::BinaryPhenotypeTable pheno(sample_to_index);
        pheno.set_value_for_sample("S1", 1);
        pheno.set_value_for_sample("S2", 1);
        pheno.set_value_for_sample("S3", 1);
        pheno.set_value_for_sample("S4", 0);
        pheno.set_value_for_sample("S5", 0);

        std::unordered_map<std::string, size_t> covar_to_index;
        stoat::CovariateTable covar(sample_to_index, covar_to_index);
        
        stoat::LogisticRegression lr;
        stoat::test_result_t tres = lr.logistic_regression(pheno, geno, covar, 0, 0);

        INFO("p_value = " << tres.pv);

        REQUIRE(std::abs(std::stod(tres.pv) - 0.8732353) < .01);
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

        std::unordered_map<std::string, size_t> covar_to_index;
        stoat::CovariateTable covar(sample_to_index, covar_to_index);
        
        stoat::LinearRegression lr;
        stoat::test_result_t tres = lr.linear_regression(pheno, geno, covar, 0, 0);

        INFO("p_value = " << tres.pv);

        REQUIRE(std::stod(tres.pv) == 0.0496);
    }

    SECTION("Linear Regression 3 paths - Moderate") {
        //     A0 A1 A2
        // S1 {1, 0, 1}
        // S2 {1, 0, 1}
        // S3 {1, 1, 0}
        // S4 {0, 1, 1}
        // S5 {0, 1, 1}

        std::unordered_map<std::string, size_t> sample_to_index = {{"S1", 0}, {"S2", 1}, {"S3", 2},
                                                                   {"S4", 3}, {"S5", 4}};
        stoat::GenoTable geno(sample_to_index, 3);
        geno.increment_count(0, 0);
        geno.increment_count(0, 2);
        geno.increment_count(1, 0);
        geno.increment_count(1, 2);
        geno.increment_count(2, 0);
        geno.increment_count(2, 1);
        geno.increment_count(3, 1);
        geno.increment_count(3, 2);
        geno.increment_count(4, 1);
        geno.increment_count(4, 2);

        stoat::QuantitativePhenotypeTable pheno(sample_to_index);
        pheno.set_value_for_sample("S1", 11);
        pheno.set_value_for_sample("S2", 10.1);
        pheno.set_value_for_sample("S3", 5.2);
        pheno.set_value_for_sample("S4", -0.3);
        pheno.set_value_for_sample("S5", 2);

        std::unordered_map<std::string, size_t> covar_to_index;
        stoat::CovariateTable covar(sample_to_index, covar_to_index);
        
        stoat::LinearRegression lr;
        stoat::test_result_t tres = lr.linear_regression(pheno, geno, covar, 0, 0);

        INFO("p_value = " << tres.pv);

        REQUIRE(std::abs(std::stod(tres.pv) - 0.03133) < 0.01);
    }

    SECTION("Linear Regression pseudo inversion 3 paths - Moderate") {
        //     A0 A1 A2
        // S1 {1, 0, 0}
        // S2 {0, 0, 1}
        // S3 {0, 0, 1}
        // S4 {0, 1, 0}
        // S5 {0, 1, 0}
        std::unordered_map<std::string, size_t> sample_to_index = {{"S1", 0}, {"S2", 1}, {"S3", 2},
                                                                   {"S4", 3}, {"S5", 4}};
        stoat::GenoTable geno(sample_to_index, 3);
        geno.increment_count(0, 0);
        geno.increment_count(1, 2);
        geno.increment_count(2, 2);
        geno.increment_count(3, 1);
        geno.increment_count(4, 1);

        stoat::QuantitativePhenotypeTable pheno(sample_to_index);
        pheno.set_value_for_sample("S1", 11);
        pheno.set_value_for_sample("S2", 10.1);
        pheno.set_value_for_sample("S3", 5.2);
        pheno.set_value_for_sample("S4", -0.3);
        pheno.set_value_for_sample("S5", 2);

        std::unordered_map<std::string, size_t> covar_to_index;
        stoat::CovariateTable covar(sample_to_index, covar_to_index);
        
        stoat::LinearRegression lr;
        stoat::test_result_t tres = lr.linear_regression(pheno, geno, covar, 0, 0);

        INFO("p_value = " << tres.pv);

        REQUIRE(std::stod(tres.pv) == 0.1505);
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

        std::unordered_map<std::string, size_t> covar_to_index = {{"covar1", 0}};
        stoat::CovariateTable covar(sample_to_index, covar_to_index);
        covar.set_value_for_sample_and_feature("S1", "covar1", 0.2);
        covar.set_value_for_sample_and_feature("S2", "covar1", 10);
        covar.set_value_for_sample_and_feature("S3", "covar1", 2);
        covar.set_value_for_sample_and_feature("S4", "covar1", 5);
        covar.set_value_for_sample_and_feature("S5", "covar1", 11);
        
        stoat::LinearRegression lr;
        stoat::test_result_t tres = lr.linear_regression(pheno, geno, covar, 0, 0);

        INFO("p_value = " << tres.pv);

        REQUIRE(std::abs(std::stod(tres.pv) - 0.1140625) < .01);
    }

    SECTION("Linear Regression pseudo inversion 3 paths - Moderate with Covariate") {
        //     A0 A1 A2
        // S1 {1, 0, 0}
        // S2 {0, 0, 1}
        // S3 {0, 0, 1}
        // S4 {0, 1, 0}
        // S5 {0, 1, 0}

        std::unordered_map<std::string, size_t> sample_to_index = {{"S1", 0}, {"S2", 1}, {"S3", 2},
                                                                   {"S4", 3}, {"S5", 4}};
        stoat::GenoTable geno(sample_to_index, 3);
        geno.increment_count(0, 0);
        geno.increment_count(1, 2);
        geno.increment_count(2, 2);
        geno.increment_count(3, 1);
        geno.increment_count(4, 1);

        stoat::QuantitativePhenotypeTable pheno(sample_to_index);
        pheno.set_value_for_sample("S1", 11);
        pheno.set_value_for_sample("S2", 10.1);
        pheno.set_value_for_sample("S3", 5.2);
        pheno.set_value_for_sample("S4", -0.3);
        pheno.set_value_for_sample("S5", 2);

        std::unordered_map<std::string, size_t> covar_to_index = {{"covar1", 0}};
        stoat::CovariateTable covar(sample_to_index, covar_to_index);
        covar.set_value_for_sample_and_feature("S1", "covar1", 0.2);
        covar.set_value_for_sample_and_feature("S2", "covar1", 10);
        covar.set_value_for_sample_and_feature("S3", "covar1", 2);
        covar.set_value_for_sample_and_feature("S4", "covar1", 5);
        covar.set_value_for_sample_and_feature("S5", "covar1", 11);
        
        stoat::LinearRegression lr;
        std::cerr << "Start test\n";
        stoat::test_result_t tres = lr.linear_regression(pheno, geno, covar, 0, 0);

        INFO("p_value = " << tres.pv);

        REQUIRE(std::abs(std::stod(tres.pv) - 0.08149071) < .0001);
    }
}

TEST_CASE("Quantitative phenotype filters") {
    SECTION("minor allele frequency filters") {
        //     A0 A1
        // S1 {2, 0}
        // S2 {2, 0}
        // S3 {2, 0}
        // S4 {2, 0}
        // S5 {0, 1}

        std::unordered_map<std::string, size_t> sample_to_index = {{"S1", 0}, {"S2", 1}, {"S3", 2},
                                                                   {"S4", 3}, {"S5", 4}};
        stoat::GenoTable geno(sample_to_index, 2);
        geno.increment_count(0, 0);
        geno.increment_count(0, 0);
        geno.increment_count(1, 0);
        geno.increment_count(1, 0);
        geno.increment_count(2, 0);
        geno.increment_count(2, 0);
        geno.increment_count(3, 0);
        geno.increment_count(3, 0);
        geno.increment_count(4, 1);

        stoat::QuantitativePhenotypeTable pheno(sample_to_index);
        pheno.set_value_for_sample("S1", 11);
        pheno.set_value_for_sample("S2", 10.1);
        pheno.set_value_for_sample("S3", 5.2);
        pheno.set_value_for_sample("S4", -0.3);
        pheno.set_value_for_sample("S5", 2);

        std::unordered_map<std::string, size_t> covar_to_index;
        stoat::CovariateTable covar(sample_to_index, covar_to_index);
        
        stoat::LinearRegression lr;
        stoat::test_result_t tres = lr.linear_regression(pheno, geno, covar, 0.5, 0);

        INFO("p_value = " << tres.pv);

        REQUIRE(tres.pv == "NA");

        tres = lr.linear_regression(pheno, geno, covar, 0.001, 0);

        INFO("p_value = " << tres.pv);
        
        REQUIRE(tres.pv != "NA");
    }

    SECTION("minimum individual filters") {
        //     A0 A1
        // S1 {1, 0}
        // S2 {1, 0}
        // S3 {1, 0}
        // S4 {1, 0}
        // S5 {0, 1}

        std::unordered_map<std::string, size_t> sample_to_index = {{"S1", 0}, {"S2", 1}, {"S3", 2},
                                                                   {"S4", 3}, {"S5", 4}};
        stoat::GenoTable geno(sample_to_index, 2);
        geno.increment_count(0, 0);
        geno.increment_count(1, 0);
        geno.increment_count(2, 0);
        geno.increment_count(3, 0);
        geno.increment_count(4, 1);

        stoat::QuantitativePhenotypeTable pheno(sample_to_index);
        pheno.set_value_for_sample("S1", 11);
        pheno.set_value_for_sample("S2", 10.1);
        pheno.set_value_for_sample("S3", 5.2);
        pheno.set_value_for_sample("S4", -0.3);
        pheno.set_value_for_sample("S5", 2);

        std::unordered_map<std::string, size_t> covar_to_index;
        stoat::CovariateTable covar(sample_to_index, covar_to_index);
        
        stoat::LinearRegression lr;
        stoat::test_result_t tres = lr.linear_regression(pheno, geno, covar, 0, 10);

        INFO("p_value = " << tres.pv);

        REQUIRE(tres.pv == "NA");

        tres = lr.linear_regression(pheno, geno, covar, 0, 3);
        
        INFO("p_value = " << tres.pv);
        
        REQUIRE(tres.pv != "NA");
    }

}
