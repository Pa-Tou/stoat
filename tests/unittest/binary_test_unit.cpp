#include <catch.hpp>
#include "../../src/quantitative_table.hpp"
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
        REQUIRE(fchi.chi2_2x2(a, b, c, d)  == "9.8233e-03");
        REQUIRE(fchi.fastFishersExactTest(a, b, c, d)  == "1.9383e-02");
    }

    SECTION("Chi-square & Fisher test (significatif)") {
        std::vector<size_t> g0 = {30, 5};
        std::vector<size_t> g1 = {2, 25};
        size_t a = g0[0];
        size_t b = g0[1];
        size_t c = g1[0];
        size_t d = g1[1];
        REQUIRE(fchi.chi2_2x2(a, b, c, d)  == "9.5037e-10");
        REQUIRE(fchi.fastFishersExactTest(a, b, c, d)  == "3.5379e-10");
    }

    SECTION("Chi-square fail (N row)") {
        std::vector<size_t> g0 = {10, 15, 5};
        std::vector<size_t> g1 = {20, 10, 10};
        REQUIRE(fchi.chi2_2xN(g0, g1) == "9.6972e-02");
    }

    SECTION("Chi-square N row significatif") {
        std::vector<size_t> g0 = {5, 10, 15, 20};
        std::vector<size_t> g1 = {20, 15, 10, 5};
        REQUIRE(fchi.chi2_2xN(g0, g1) == "1.6974e-04");
    }

    SECTION("Chi-square fail (N row 1)") {
        std::vector<size_t> g0 = {10, 10, 10, 10, 10};
        std::vector<size_t> g1 = {10, 10, 10, 10, 10};
        REQUIRE(fchi.chi2_2xN(g0, g1) == "1");
    }

    SECTION("Chi-square fail & Fisher test fail (full zero row)") {
        std::vector<size_t> g0 = {0, 0};
        std::vector<size_t> g1 = {0, 0};
        size_t a = g0[0];
        size_t b = g0[1];
        size_t c = g1[0];
        size_t d = g1[1];
        REQUIRE(fchi.chi2_2x2(a, b, c, d)  == "NA");
        REQUIRE(fchi.fastFishersExactTest(a, b, c, d)  == "NA");
    }

    SECTION("Chi-square fail (zero row)") {
        std::vector<size_t> g0 = {0, 0, 0};
        std::vector<size_t> g1 = {10, 20, 30};
        REQUIRE(fchi.chi2_2xN(g0, g1) == "NA");
    }

    SECTION("Chi-square fail & Fisher test valid (zero row + column)") {
        std::vector<size_t> g0 = {0, 0};
        std::vector<size_t> g1 = {0, 1};
        size_t a = g0[0];
        size_t b = g0[1];
        size_t c = g1[0];
        size_t d = g1[1];
        REQUIRE(fchi.chi2_2x2(a, b, c, d) == "NA");
        REQUIRE(fchi.fastFishersExactTest(a, b, c, d) == "NA");
    }

    SECTION("Chi-square & Fisher test (1/0 0/1)") {
        std::vector<size_t> g0 = {1, 0};
        std::vector<size_t> g1 = {0, 1};
        size_t a = g0[0];
        size_t b = g0[1];
        size_t c = g1[0];
        size_t d = g1[1];
        REQUIRE(fchi.chi2_2x2(a, b, c, d)  == "0.1573");
        REQUIRE(fchi.fastFishersExactTest(a, b, c, d)  == "1");
    }

    SECTION("Chi-square & Fisher test (strange but correct)") {
        std::vector<size_t> g0 = {79, 18};
        std::vector<size_t> g1 = {96, 23};
        size_t a = g0[0];
        size_t b = g0[1];
        size_t c = g1[0];
        size_t d = g1[1];
        REQUIRE(fchi.chi2_2x2(a, b, c, d)  == "0.8857");
        REQUIRE(fchi.fastFishersExactTest(a, b, c, d)  == "1");
    }

    SECTION("Chi-square & Fisher test (very significative)") {
        std::vector<size_t> g0 = {122, 78};
        std::vector<size_t> g1 = {27, 173};
        size_t a = g0[0];
        size_t b = g0[1];
        size_t c = g1[0];
        size_t d = g1[1];
        REQUIRE(fchi.chi2_2x2(a, b, c, d)  == "8.8051e-23");
        REQUIRE(fchi.fastFishersExactTest(a, b, c, d)  == "1.4799e-23");
    }
    SECTION("Chi-square & Fisher test (very significative) with empty counts") {
        std::vector<size_t> g0 = {122, 0, 78};
        std::vector<size_t> g1 = {27, 0, 173};
        size_t a = g0[0];
        size_t b = g0[1];
        size_t c = g1[0];
        size_t d = g1[1];
        REQUIRE(fchi.chi2_2xN(g0, g1)  == "8.8051e-23");
    }
}

TEST_CASE("Binary phenotype filters") {
    FisherChi2 fchi;

    SECTION("minor allele frequency filters") {
        //     A0 A1
        // S1 {2, 0}
        // S2 {2, 0}
        // S3 {2, 0}
        // S4 {2, 0}
        // S5 {0, 1}

        std::unordered_map<std::string, size_t> sample_to_index = {{"S1", 0}, {"S2", 1}, {"S3", 2},
                                                                   {"S4", 3}, {"S5", 4}};
        stoat::GenotypeTable geno(sample_to_index, 3);
        geno.increment_count("S1", 0);
        geno.increment_count("S1", 0);
        geno.increment_count("S2", 0);
        geno.increment_count("S2", 0);
        geno.increment_count("S3", 0);
        geno.increment_count("S3", 0);
        geno.increment_count("S4", 0);
        geno.increment_count("S4", 0);
        geno.increment_count("S5", 1);

        stoat::BinaryPhenotypeTable pheno(sample_to_index);
        pheno.set_value_for_sample("S1", 1);
        pheno.set_value_for_sample("S2", 0);
        pheno.set_value_for_sample("S3", 0);
        pheno.set_value_for_sample("S4", 0);
        pheno.set_value_for_sample("S5", 1);
        
        stoat::test_result_t tres = fchi.fisher_chi2(pheno, geno, 0.5, 0);

        INFO("p_value fisher = " << tres.pv);
        INFO("p_value chi2 = " << tres.second_pv);

        REQUIRE(tres.pv == "NA");
        REQUIRE(tres.second_pv == "NA");

        tres = fchi.fisher_chi2(pheno, geno, 0.001, 0);

        INFO("p_value fisher = " << tres.pv);
        INFO("p_value chi2 = " << tres.second_pv);
        
        REQUIRE(tres.pv != "NA");
        REQUIRE(tres.second_pv != "NA");
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
        stoat::GenotypeTable geno(sample_to_index, 3);
        geno.increment_count("S1", 0);
        geno.increment_count("S2", 0);
        geno.increment_count("S3", 0);
        geno.increment_count("S4", 0);
        geno.increment_count("S5", 1);

        stoat::BinaryPhenotypeTable pheno(sample_to_index);
        pheno.set_value_for_sample("S1", 1);
        pheno.set_value_for_sample("S2", 0);
        pheno.set_value_for_sample("S3", 0);
        pheno.set_value_for_sample("S4", 0);
        pheno.set_value_for_sample("S5", 1);
        
        stoat::test_result_t tres = fchi.fisher_chi2(pheno, geno, 0, 10);

        INFO("p_value fisher = " << tres.pv);
        INFO("p_value chi2 = " << tres.second_pv);

        REQUIRE(tres.pv == "NA");
        REQUIRE(tres.second_pv == "NA");

        tres = fchi.fisher_chi2(pheno, geno, 0, 3);
        
        INFO("p_value fisher = " << tres.pv);
        INFO("p_value chi2 = " << tres.second_pv);

        REQUIRE(tres.pv != "NA");
        REQUIRE(tres.second_pv != "NA");
    }

}
