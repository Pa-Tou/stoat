#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include <iostream>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

namespace fs = std::filesystem;

// === UTILITY FUNCTIONS ===

void clean_output_dir(const std::string& output_dir) {
    if (fs::exists(output_dir))
        fs::remove_all(output_dir);
    fs::create_directory(output_dir);
}

// if you wander why we doing this instead of a simple compare files like 'diff' in bash 
// the ans is : because multi-threading.
bool files_equal(const std::string& file1, const std::string& file2) {
    std::unordered_map<std::string, std::string> map1, map2;

    auto load_file = [](const std::string& path, std::unordered_map<std::string, std::string>& map) {
        std::ifstream infile(path);
        std::string line;
        while (std::getline(infile, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            std::string token;
            std::vector<std::string> columns;

            while (std::getline(ss, token, '\t')) {
                columns.push_back(token);
            }

            if (columns.size() < 4) {
                throw std::runtime_error("Invalid line with fewer than 4 columns: " + line);
            }

            std::string snarl_id = columns[3]; // 4th column = SNARL
            map[snarl_id] = line;
        }
    };

    load_file(file1, map1);
    load_file(file2, map2);

    // Check missing in file2
    for (const auto& [snarl, line1] : map1) {
        auto it = map2.find(snarl);
        if (it == map2.end()) {
            std::cerr << "Missing SNARL in file2: " << snarl << std::endl;
            return false;
        } else if (line1 != it->second) {
            std::cerr << "Mismatch for SNARL " << snarl << ":\n"
                      << "File1: " << line1 << "\n"
                      << "File2: " << it->second << "\n";
            return false;
        }
    }

    // Check missing in file1
    for (const auto& [snarl, _] : map2) {
        if (map1.find(snarl) == map1.end()) {
            std::cerr << "Missing SNARL in file1: " << snarl << std::endl;
            return false;
        }
    }

    return true;
}

bool compare_output_dirs(const std::string& output_dir, const std::string& expected_dir) {
    for (const auto& file : fs::directory_iterator(expected_dir)) {
        auto expected_file = file.path();
        auto output_file = fs::path(output_dir) / expected_file.filename();

        if (!fs::exists(output_file)) {
            std::cerr << "Missing output file: " << output_file << "\n";
            return false;
        }
        if (!files_equal(expected_file, output_file)) {
            std::cerr << "Mismatch in file: " << expected_file.filename() << "\n";
            return false;
        }
    }
    return true;
}

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
