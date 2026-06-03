#include <catch.hpp>
#include <filesystem>
#include <fstream>
#include <unordered_map>
#include <string>
#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <tuple>
#include <vector>
#include "../../src/arg_parser.hpp"  // adjust path as needed
#include "../../src/log.hpp"  // adjust path as needed

// Helper to write minimal VCF content to file
void write_vcf_file(const std::string& path, const std::string& content) {
    std::filesystem::create_directories("test_data");
    std::ofstream file(path);
    file << content;
    file.close();
}

// Helper to create a temporary test file
std::string create_test_pheno_file(const std::string& content, const std::string& filename = "test_pheno.txt") {
    std::ofstream out(filename);
    out << content;
    out.close();
    return filename;
}

// Helper function to create a test file
std::string create_test_covar_file(const std::string& content, const std::string& filename = "test_covariate.txt") {
    std::ofstream out(filename);
    out << content;
    out.close();
    return filename;
}

TEST_CASE("Binary phenotype parsing", "[stoat_vcf::parse_binary_pheno_table]") {
    SECTION("Valid phenotype file with mixed cases and controls") {
        
        std::string file_content =
            "SAMPLE PHENO\n"
            "I1 0\n"
            "I2 1\n"
            "I3 0\n"
            "I4 1\n";
        std::string file_path = create_test_pheno_file(file_content);

        std::unordered_map<std::string, size_t> sample_to_index;
        std::unique_ptr<stoat::BinaryPhenotypeTable> result = stoat_vcf::parse_binary_pheno_table(file_path, sample_to_index);

        REQUIRE(result->get_sample_names().size() == 4);
        REQUIRE(result->get_value_for_sample("I1") == false);
        REQUIRE(result->get_value_for_sample("I2") == true);
        REQUIRE(result->get_value_for_sample("I3") == false);
        REQUIRE(result->get_value_for_sample("I4") == true);
    }

    SECTION("Invalid header") {
        std::string file_content =
            "ID PHENO\n"
            "I1 1\n";
        std::string file_path = create_test_pheno_file(file_content);
        std::unordered_map<std::string, size_t> sample_to_index;
       
        REQUIRE_THROWS_AS(stoat_vcf::parse_binary_pheno_table(file_path, sample_to_index), std::invalid_argument);
    }

    SECTION("Non-binary phenotype value") {
        std::string file_content =
            "SAMPLE PHENO\n"
            "I1 3\n";
        std::string file_path = create_test_pheno_file(file_content);
        std::unordered_map<std::string, size_t> sample_to_index;

        REQUIRE_THROWS_AS(stoat_vcf::parse_binary_pheno_table(file_path, sample_to_index), std::invalid_argument);
    }

    SECTION("Malformed line") {
        std::string file_content =
            "SAMPLE PHENO\n"
            "I1\n";
        std::string file_path = create_test_pheno_file(file_content);
        std::unordered_map<std::string, size_t> sample_to_index;

        REQUIRE_THROWS_AS(stoat_vcf::parse_binary_pheno_table(file_path, sample_to_index), std::invalid_argument);
    }

    SECTION("Non-integer phenotype") {
        std::string file_content =
            "SAMPLE PHENO\n"
            "I1 X\n";
        std::string file_path = create_test_pheno_file(file_content);
        std::unordered_map<std::string, size_t> sample_to_index;

        REQUIRE_THROWS_AS(stoat_vcf::parse_binary_pheno_table(file_path, sample_to_index), std::invalid_argument);
    }
}


TEST_CASE("Quantitative phenotype parsing", "[stoat_vcf::parse_quantitative_pheno]") {
    SECTION("Valid quantitative phenotype file") {
        std::string file_content =
            "SAMPLE PHENO\n"
            "I1 1.5\n"
            "I2 2.0\n"
            "I3 -3.25\n"
            "I4 0.0\n";
        std::string file_path = create_test_pheno_file(file_content);

        std::unordered_map<std::string, size_t> sample_to_index;
        std::unique_ptr<stoat::QuantitativePhenotypeTable> result = stoat_vcf::parse_quantitative_pheno_table(file_path, sample_to_index);

        REQUIRE(result->get_sample_names().size() == 4);
        REQUIRE(result->get_value_for_sample("I1") == 1.5);
        REQUIRE(result->get_value_for_sample("I2") == 2.0);
        REQUIRE(result->get_value_for_sample("I3") == -3.25);
        REQUIRE(result->get_value_for_sample("I4") == 0.0);
    }

    SECTION("Invalid header in quantitative phenotype file") {
        std::string file_content =
            "ID PHENO\n"
            "I1 1.5\n";
        std::string file_path = create_test_pheno_file(file_content);

        std::unordered_map<std::string, size_t> sample_to_index;
        REQUIRE_THROWS_AS(stoat_vcf::parse_quantitative_pheno_table(file_path, sample_to_index), std::invalid_argument);
    }

    SECTION("Non-numeric phenotype value") {
        std::string file_content =
            "SAMPLE PHENO\n"
            "I1 abc\n";
        std::string file_path = create_test_pheno_file(file_content);

        std::unordered_map<std::string, size_t> sample_to_index;
        REQUIRE_THROWS_AS(stoat_vcf::parse_quantitative_pheno_table(file_path, sample_to_index), std::invalid_argument);
    }

    SECTION("Malformed line in quantitative phenotype file") {
        std::string file_content =
            "SAMPLE PHENO\n"
            "I1\n";
        std::string file_path = create_test_pheno_file(file_content);

        std::unordered_map<std::string, size_t> sample_to_index;
        REQUIRE_THROWS_AS(stoat_vcf::parse_quantitative_pheno_table(file_path, sample_to_index), std::invalid_argument);
    }
}

TEST_CASE("VCF Parsing", "[parse_vcf][stoat_vcf::parseHeader]") {
    SECTION("Valid VCF with 3 samples") {
        std::string vcf_content =
            "##fileformat=VCFv4.2\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsampleA\tsampleB\tsampleC\n"
            "1\t12345\trsTest\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\n";

        std::string vcf_path = "test_data/test_valid.vcf";
        write_vcf_file(vcf_path, vcf_content);

        auto [samples, ptr, hdr, rec] = stoat_vcf::parseHeader(vcf_path);

        REQUIRE(samples.size() == 3);
        REQUIRE(samples[0] == "sampleA");
        REQUIRE(samples[1] == "sampleB");
        REQUIRE(samples[2] == "sampleC");

        bcf_destroy(rec);
        bcf_hdr_destroy(hdr);
        bcf_close(ptr);
    }

    SECTION("Malformed VCF should throw error") {
        std::string vcf_content =
            "This is not a VCF\n"
            "Just some junk lines\n";

        std::string vcf_path = "test_data/test_invalid.vcf";
        write_vcf_file(vcf_path, vcf_content);

        REQUIRE_THROWS_AS(stoat_vcf::parseHeader(vcf_path), std::invalid_argument);
    }

    SECTION("Empty VCF should throw error on header read") {
        std::string vcf_path = "test_data/test_empty.vcf";
        write_vcf_file(vcf_path, "");  // Empty file

        REQUIRE_THROWS_AS(stoat_vcf::parseHeader(vcf_path), std::invalid_argument);
    }
}

TEST_CASE("Check sample matching", "[check_match_samples]") {
    SECTION("All keys found, size match (bool)") {
        std::unordered_map<std::string, bool> pheno = {
            {"sample1", true},
            {"sample2", false},
            {"sample3", true}
        };
        std::vector<std::string> vcf_samples = {"sample1", "sample2", "sample3"};

        REQUIRE_NOTHROW(stoat_vcf::check_match_samples<bool>(pheno, vcf_samples));
    }

    SECTION("All keys found, size mismatch (double)") {
        std::unordered_map<std::string, double> pheno = {
            {"sample1", 1.0},
            {"sample2", 2.0},
            {"sample3", 3.0},
            {"extra_sample", 4.0}
        };
        std::vector<std::string> vcf_samples = {"sample1", "sample2", "sample3"};

        // Should not throw, just a warning (which we can't easily assert here)
        REQUIRE_NOTHROW(stoat_vcf::check_match_samples<double>(pheno, vcf_samples));
    }

    SECTION("Missing key throws invalid_argument") {
        std::unordered_map<std::string, bool> pheno = {
            {"sample1", true},
            {"sample2", false}
        };
        std::vector<std::string> vcf_samples = {"sample1", "sample2", "sample3"};

        REQUIRE_THROWS_AS(stoat_vcf::check_match_samples<bool>(pheno, vcf_samples), std::invalid_argument);
    }
}

TEST_CASE("Parse covariates from file", "[stoat_vcf::parse_covariates]") {
    SECTION("Valid covariate file and columns") {
        std::string content =
            "SAMPLE age sex\n"
            "samp0 25 1\n"
            "samp1 30 0\n"
            "samp2 45 1\n";
        std::string path = create_test_covar_file(content);
        std::vector<std::string> covars = {"age", "sex"};

        std::unordered_map<std::string, size_t> sample_to_index;
        std::unordered_map<std::string, size_t> covar_to_index = {{"age", 0}, {"sex", 1}};
        std::unique_ptr<stoat::CovariateTable> result = stoat_vcf::parse_covariate_table(path, sample_to_index, covar_to_index);

        REQUIRE(result->get_sample_names().size() == 3);
        REQUIRE(result->get_value_for_sample_and_feature("samp0", "age") == 25);
        REQUIRE(result->get_value_for_sample_and_feature("samp0", "sex") == 1);
        REQUIRE(result->get_value_for_sample_and_feature("samp1", "age") == 30);
        REQUIRE(result->get_value_for_sample_and_feature("samp1", "sex") == 0);
        REQUIRE(result->get_value_for_sample_and_feature("samp2", "age") == 45);
        REQUIRE(result->get_value_for_sample_and_feature("samp2", "sex") == 1);
    }

    SECTION("Missing SAMPLE column") {
        std::string content =
            "ID age sex\n"
            "samp0 25 1\n";
        std::string path = create_test_covar_file(content);

        std::unordered_map<std::string, size_t> sample_to_index;
        std::unordered_map<std::string, size_t> covar_to_index = {{"age", 0}, {"sex", 1}};
                
        REQUIRE_THROWS_AS(stoat_vcf::parse_covariate_table(path, sample_to_index, covar_to_index), std::invalid_argument);
    }

    SECTION("Missing covariate column") {
        std::string content =
            "SAMPLE age sex\n"
            "A 25 1\n";
        std::string path = create_test_covar_file(content);

        std::unordered_map<std::string, size_t> sample_to_index;
        std::unordered_map<std::string, size_t> covar_to_index = {{"height", 0}};
                
        REQUIRE_THROWS_AS(stoat_vcf::parse_covariate_table(path, sample_to_index, covar_to_index), std::invalid_argument);
    }

    SECTION("Non-numeric value in covariate field") {
        std::string content =
            "SAMPLE age sex\n"
            "samp0 XX 1\n";
        std::string path = create_test_covar_file(content);

        std::unordered_map<std::string, size_t> sample_to_index;
        std::unordered_map<std::string, size_t> covar_to_index = {{"age", 0}, {"sex", 1}};

        REQUIRE_THROWS_AS(stoat_vcf::parse_covariate_table(path, sample_to_index, covar_to_index), std::invalid_argument);
    }
}
