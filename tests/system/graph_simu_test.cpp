#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "compare_files_utils.hpp"

namespace fs = std::filesystem;

bool run_test(
    const std::string& binary,
    const std::string& output_dir,
    const std::string& expected_dir,
    const std::string& data_path,
    const std::string phenotype_command,
    const std::string phenotype,
    bool use_covariate = false) {

    clean_output_dir(output_dir);

    std::string cmd = binary + " graph"
        + " -g " + data_path + "/pg.full.pg"
        + " -d " + data_path + "/pg.full.dist"
        + " -S " + data_path + "/samples.g0.tsv"
        + phenotype_command + " -r ref";

    if (use_covariate) {
        cmd += " --covariate " + data_path + "/covariate.tsv"
             + " --covar-name CP1,SEX,CP3";
    }

    cmd += " --output " + output_dir;

    int command_output = std::system(cmd.c_str());
    if (command_output != 0) {
        std::cerr << "Command failed: " << cmd << "\n";
        return false;
    }

    bool result = compare_output_dirs(output_dir, expected_dir);
    clean_output_dir(output_dir);

    return result;
}

TEST_CASE("Binary association tests graph", "[binary]") {
    const std::string binary = "../bin/stoat";
    const std::string output_dir = "../output_binary";
    const std::string expected_dir = "../tests/expected_output/graph/binary";
    const std::string expected_dir_covar = "../tests/expected_output/graph/binary_covar";
    const std::string data_path = "../data/binary";
    const std::string phenotype_command = " -T chi2 ";
    const std::string phenotype = "binary";

    SECTION("Without covariate") {
        REQUIRE(run_test(binary, output_dir, expected_dir, data_path, phenotype_command, phenotype, false));
        clean_output_dir(output_dir);
    }

    // SECTION("With covariate") {
    //     REQUIRE(run_test(binary, output_dir, expected_dir_covar, data_path, phenotype_command, true));
    // }
}

// TEST_CASE("Quantitative trait tests graph", "[quantitative]") {
//     const std::string binary = "./stoat";
//     const std::string output_dir = "../output_quantitative";
//     const std::string expected_dir = "../tests/expected_output/graph/quantitative";
//     const std::string expected_dir_covar = "../tests/expected_output/graph/quantitative_covar";
//     const std::string data_path = "../data/quantitative";
//     const std::string phenotype_command = " -q ";
//     const std::string phenotype = "quantitative";

//     //TODO: I added this so it would compile, idk what it should be
//     std::string sample_of_interest;

//     SECTION("Without covariate") {
//         REQUIRE(run_test(binary, output_dir, expected_dir, data_path, sample_of_interest, false));
//     }

//     SECTION("With covariate") {
//         REQUIRE(run_test(binary, output_dir, expected_dir_covar, data_path, sample_of_interest, true));
//     }
// }
