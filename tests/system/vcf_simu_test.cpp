#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "compare_files_utils.hpp"

namespace fs = std::filesystem;

bool run_test_snarl(
    const std::string& stoat_command,
    const std::string& output_dir,
    const std::string& expected_dir,
    const std::string& data_path) {

    std::string cmd = stoat_command + " vcf"
    + " -p " + data_path + "/pg.full.pg"
    + " -d " + data_path + "/pg.full.dist"
    + " -r " + data_path + "/pg.chromosome"
    + " --output " + output_dir;

    std::cout << "Command run : \n" << cmd << std::endl;

    int command_output = std::system(cmd.c_str());
    if (command_output != 0) {
        std::cerr << "Command failed: " << cmd << "\n";
        return false;
    }

    bool result = compare_output_dirs(output_dir, expected_dir);
    clean_output_dir(output_dir);

    return result;
}

bool run_test(
    const std::string& stoat_command,
    const std::string& output_dir,
    const std::string& expected_dir,
    const std::string& data_path,
    const std::string& phenotype,
    bool use_covariate = false) {

    std::string cmd = stoat_command + " vcf"
    + " -s " + data_path + "/snarl_analyse.tsv"
    + " -v " + data_path + "/merged_output.vcf.gz";

    std::string type;

    if (phenotype == "eqtl") {
        cmd += " -e " + data_path + "/qtl.tsv" 
        + " --gene-position " + data_path + "/gene_position.tsv";

    } else if (phenotype == "binary") {
        cmd += " -b " + data_path + "/phenotype.tsv";

    } else if (phenotype == "quantitative") {
        cmd += " -q " + data_path + "/phenotype.tsv";

    } else {
        std::cerr << "Phenotype Error !" << std::endl;
        return false;
    }

    if (use_covariate) {
        cmd += " --covariate " + data_path + "/covariate.tsv"
             + " --covar-name CP1,SEX,CP3";
    }

    cmd += " --output " + output_dir;

    std::cout << "Command run : \n" << cmd << std::endl;

    int command_output = std::system(cmd.c_str());
    if (command_output != 0) {
        std::cerr << "Command failed: " << cmd << "\n";
        return false;
    }

    bool result = compare_output_dirs(output_dir, expected_dir);
    clean_output_dir(output_dir);

    return result;
}

bool run_test_full(
    const std::string& stoat_command,
    const std::string& output_dir,
    const std::string& expected_dir,
    const std::string& data_path,
    const std::string& phenotype,
    bool use_covariate = false) {

    std::string cmd = stoat_command + " vcf"
    + " -p " + data_path + "/pg.full.pg"
    + " -d " + data_path + "/pg.full.dist"
    + " -r " + data_path + "/pg.chromosome"
    + " -v " + data_path + "/merged_output.vcf.gz";

    std::string type;

    if (phenotype == "eqtl") {
        cmd += " -e " + data_path + "/qtl.tsv" 
        + " --gene-position " + data_path + "/gene_position.tsv";

    } else if (phenotype == "binary") {
        cmd += " -b " + data_path + "/phenotype.tsv";

    } else if (phenotype == "quantitative") {
        cmd += " -q " + data_path + "/phenotype.tsv";

    } else {
        std::cerr << "Phenotype Error !" << std::endl;
        return false;
    }

    if (use_covariate) {
        cmd += " --covariate " + data_path + "/covariate.tsv"
             + " --covar-name CP1,SEX,CP3";
    }

    cmd += " --output " + output_dir;

    std::cout << "Command run : \n" << cmd << std::endl;

    int command_output = std::system(cmd.c_str());
    if (command_output != 0) {
        std::cerr << "Command failed: " << cmd << "\n";
        return false;
    }

    bool result = compare_output_dirs(output_dir, expected_dir);
    clean_output_dir(output_dir);

    return result;
}

TEST_CASE("Snarl decomposition", "[snarl]") {
    const std::string stoat_command = "../bin/stoat";
    const std::string output_dir = "../output_snarl";
    const std::string expected_dir = "../tests/expected_output/vcf/snarl";
    const std::string data_path = "../data/binary";

    SECTION("Binary decomposition") {
        REQUIRE(run_test_snarl(stoat_command, output_dir, expected_dir, data_path));
    }
}

TEST_CASE("Binary association tests vcf", "[binary]") {
    const std::string stoat_command = "../bin/stoat";
    const std::string output_dir = "../output_binary";
    const std::string expected_dir = "../tests/expected_output/vcf/binary";
    const std::string expected_dir_covar = "../tests/expected_output/vcf/binary_covar";
    const std::string data_path = "../data/binary";
    const std::string phenotype = "binary";

    SECTION("Without covariate") {
        REQUIRE(run_test_full(stoat_command, output_dir, expected_dir, data_path, phenotype, false));
    }

    SECTION("With covariate") {
        REQUIRE(run_test_full(stoat_command, output_dir, expected_dir_covar, data_path, phenotype, true));
    }
}

TEST_CASE("Quantitative trait tests vcf", "[quantitative]") {
    const std::string stoat_command = "../bin/stoat";
    const std::string output_dir = "../output_quantitative";
    const std::string expected_dir = "../tests/expected_output/vcf/quantitative";
    const std::string expected_dir_covar = "../tests/expected_output/vcf/quantitative_covar";
    const std::string data_path = "../data/quantitative";
    const std::string phenotype = "quantitative";

    SECTION("Without covariate") {
        REQUIRE(run_test_full(stoat_command, output_dir, expected_dir, data_path, phenotype, false));
    }

    SECTION("With covariate") {
        REQUIRE(run_test_full(stoat_command, output_dir, expected_dir_covar, data_path, phenotype, true));
    }
}

TEST_CASE("eQTL tests vcf", "[eqtl]") {
    const std::string stoat_command = "./stoat";
    const std::string output_dir = "../output_eqtl";
    const std::string expected_dir = "../tests/expected_output/vcf/eqtl";
    const std::string expected_dir_covar = "../tests/expected_output/vcf/eqtl_covar";
    const std::string data_path = "../data/eqtl";
    const std::string phenotype = "eqtl";

    SECTION("Without covariate") {
        REQUIRE(run_test(stoat_command, output_dir, expected_dir, data_path, phenotype, false));
    }

    SECTION("With covariate") {
        REQUIRE(run_test(stoat_command, output_dir, expected_dir_covar, data_path, phenotype, true));
    }
}
