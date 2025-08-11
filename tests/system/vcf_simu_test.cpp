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
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <string>

bool files_equal(const std::string& file1, const std::string& file2) {
    std::unordered_map<std::string, std::string> map1, map2;

    auto load_file = [](const std::string& path, std::unordered_map<std::string, std::string>& map) {
        std::ifstream infile(path);
        std::string line;
        int snarl_column = -1;
        bool header_processed = false;

        while (std::getline(infile, line)) {
            if (line.empty()) continue;

            std::istringstream ss(line);
            std::string token;
            std::vector<std::string> columns;

            while (std::getline(ss, token, '\t')) {
                columns.push_back(token);
            }

            // First non-empty line is the header
            if (!header_processed) {
                for (size_t i = 0; i < columns.size(); ++i) {
                    if (columns[i] == "SNARL") {
                        snarl_column = static_cast<int>(i);
                        break;
                    }
                }

                if (snarl_column == -1) {
                    std::cerr << "SNARL column not found in header of file: " << path << std::endl;
                    exit(1);
                }

                header_processed = true;
                continue; // skip header line
            }

            if (snarl_column >= columns.size()) {
                std::cerr << "Invalid line (too few columns) in file: " << path << "\nLine: " << line << std::endl;
                exit(1);
            }

            std::string snarl_id = columns[snarl_column];
            map[snarl_id] = line;
        }
    };

    load_file(file1, map1);
    load_file(file2, map2);

    // Compare file1 against file2
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

    // Compare file2 against file1
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
        cmd += " -e " + data_path + "/quantitative/qtl.tsv" 
        + " --gene-position " + data_path + "/quantitative/gene_position.tsv";

    } else if (phenotype == "binary") {
        cmd += " -b " + data_path + "/phenotype.tsv";

    } else if (phenotype == "quantitative") {
        cmd += " -q " + data_path + "/phenotype.tsv";
    }

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

TEST_CASE("Binary association tests vcf", "[binary]") {
    const std::string stoat_command = "../bin/stoat";
    const std::string output_dir = "../output_binary";
    const std::string expected_dir = "../tests/expected_output/vcf/binary";
    const std::string expected_dir_covar = "../tests/expected_output/vcf/binary_covar";
    const std::string data_path = "../data/binary";
    const std::string phenotype = "binary";

    SECTION("Without covariate") {
        REQUIRE(run_test(stoat_command, output_dir, expected_dir, data_path, phenotype, false));
    }

    SECTION("With covariate") {
        REQUIRE(run_test(stoat_command, output_dir, expected_dir_covar, data_path, phenotype, true));
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
        REQUIRE(run_test(stoat_command, output_dir, expected_dir, data_path, phenotype, false));
    }

    SECTION("With covariate") {
        REQUIRE(run_test(stoat_command, output_dir, expected_dir_covar, data_path, phenotype, true));
    }
}

// TEST_CASE("eQTL tests vcf", "[eqtl]") {
//     const std::string stoat_command = "./stoat";
//     const std::string output_dir = "../output_eqtl";
//     const std::string expected_dir = "../tests/expected_output/vcf/eqtl";
//     const std::string expected_dir_covar = "../tests/expected_output/vcf/eqtl_covar";
//     const std::string phenotype_command = " -e ";

//     const std::string data_path = "../data/eqtl";

//     SECTION("Without covariate") {
//         REQUIRE(run_test(stoat_command, output_dir, expected_dir, data_path, phenotype_command, false, true));
//     }

//     SECTION("With covariate") {
//         REQUIRE(run_test(stoat_command, output_dir, expected_dir_covar, data_path, phenotype_command, true, true));
//     }
// }
