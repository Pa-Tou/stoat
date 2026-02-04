#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "compare_files_utils.hpp"
#include "../../src/snarl_data_collection.hpp"

namespace fs = std::filesystem;

//Check if two snarl collections are equivalent
bool compare_snarl_collection(std::string test_name, std::string truth_name) {
    std::cerr << "Compare files " << test_name << " and " << truth_name << std::endl;
    SnarlDataCollection test_snarl(0,std::numeric_limits<size_t>::max(),std::numeric_limits<size_t>::max());
    SnarlDataCollection truth_snarl(0,std::numeric_limits<size_t>::max(),std::numeric_limits<size_t>::max());

    test_snarl.load_snarl_data_collection(test_name);
    truth_snarl.load_snarl_data_collection(truth_name);

    return SnarlDataCollection::is_equivalent(test_snarl, truth_snarl); 
}

bool run_test_snarl(
    const std::string& stoat_command,
    const std::string& output_dir,
    const std::string& expected_dir,
    const std::string& data_path) {

    std::string cmd = stoat_command + " vcf"
    + " -g " + data_path + "/pg.full.pg"
    + " -d " + data_path + "/pg.full.dist"
    + " -r " + data_path + "/pg.chromosome"
    + " --output " + output_dir;

    std::cout << "Command run : \n" << cmd << std::endl;

    int command_output = std::system(cmd.c_str());
    if (command_output != 0) {
        std::cerr << "Command failed: " << cmd << "\n";
        return false;
    }

    bool passed = compare_snarl_collection(output_dir + "/snarl_info.tsv", expected_dir + "/snarl_info.tsv");
    clean_output_dir(output_dir);
    return passed;
}


bool run_test(
    const std::string& stoat_command,
    const std::string& output_dir,
    const std::string& expected_dir,
    const std::string& data_path,
    const std::string& phenotype,
    bool use_covariate = false) {

    std::string cmd = stoat_command + " vcf"
    + " -s " + data_path + "/snarl_info.tsv"
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
             + " --covar-name PC1,SEX,PC3";
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
    + " -g " + data_path + "/pg.full.pg"
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
             + " --covar-name PC1,SEX,PC3";
    }

    cmd += " --output " + output_dir;

    std::cout << "Command run : \n" << cmd << std::endl;

    int command_output = std::system(cmd.c_str());
    if (command_output != 0) {
        std::cerr << "Command failed: " << cmd << "\n";
        return false;
    }

    bool result = compare_output_dirs(output_dir, expected_dir);

    result &= compare_snarl_collection(output_dir + "/snarl_info.tsv", expected_dir + "/snarl_info.tsv");

    clean_output_dir(output_dir);

    return result;
}

TEST_CASE("Snarl decomposition", "[snarl]") {
    const std::string stoat_command = "../bin/stoat";
    const std::string output_dir = "../output_snarl";
    const std::string expected_dir = "../tests/test_data/expected_output/vcf/output_snarl";
    const std::string data_path = "../tests/test_data/input_data/binary";

    SECTION("Binary decomposition") {
        REQUIRE(run_test_snarl(stoat_command, output_dir, expected_dir, data_path));
    }
}

TEST_CASE("Bad paths", "[test]") {
    const std::string stoat_command = "../bin/stoat";
    const std::string output_dir = "../output_bad_paths";
    const std::string expected_dir = "../tests/test_data/expected_output/vcf/output_bad_paths";
    const std::string data_path = "../tests/test_data/input_data/bad_paths";
    const std::string phenotype = "binary";

    SECTION("Bad paths fail - snarl paths given") {
        REQUIRE(run_test(stoat_command, output_dir, expected_dir, data_path, phenotype));
    }
}

TEST_CASE("Binary association tests vcf", "[binary]") {
    const std::string stoat_command = "../bin/stoat";
    const std::string output_dir = "../output_binary";
    const std::string expected_dir = "../tests/test_data/expected_output/vcf/output_binary";
    const std::string expected_dir_covar = "../tests/test_data/expected_output/vcf/output_binary_covar";
    const std::string data_path = "../tests/test_data/input_data/binary";
    const std::string phenotype = "binary";

    SECTION("Without covariate") {
        REQUIRE(run_test_full(stoat_command, output_dir, expected_dir, data_path, phenotype, false));
    }

    SECTION("With covariate") {
        REQUIRE(run_test_full(stoat_command, output_dir, expected_dir_covar, data_path, phenotype, true));
    }
}
TEST_CASE("Binary association tests with snarl resolving vcf", "[binary]") {
    const std::string stoat_command = "../bin/stoat";
    const std::string output_dir = "../output_binary";
    const std::string expected_dir = "../tests/test_data/expected_output/vcf/output_binary";
    const std::string expected_dir_covar = "../tests/test_data/expected_output/vcf/output_binary_covar";
    const std::string data_path = "../tests/test_data/input_data/binary";
    const std::string phenotype = "binary";

    SECTION("Without covariate") {

        std::string cmd = stoat_command + " vcf -R"
        + " -g " + data_path + "/pg.full.pg"
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

        } 

        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }

        bool result = compare_output_dirs(output_dir, expected_dir);

        result &= compare_snarl_collection(output_dir + "/snarl_info.tsv", expected_dir + "/snarl_info.tsv");

        clean_output_dir(output_dir);

        REQUIRE( result);

    }

}

TEST_CASE("Quantitative trait tests vcf", "[quantitative]") {
    const std::string stoat_command = "../bin/stoat";
    const std::string output_dir = "../output_quantitative";
    const std::string expected_dir = "../tests/test_data/expected_output/vcf/output_quantitative";
    const std::string expected_dir_covar = "../tests/test_data/expected_output/vcf/output_quantitative_covar";
    const std::string data_path = "../tests/test_data/input_data/quantitative";
    const std::string phenotype = "quantitative";

    SECTION("Without covariate") {
        REQUIRE(run_test_full(stoat_command, output_dir, expected_dir, data_path, phenotype, false));
    }

    SECTION("With covariate") {
        REQUIRE(run_test_full(stoat_command, output_dir, expected_dir_covar, data_path, phenotype, true));
    }
}

TEST_CASE("eQTL tests vcf", "[eqtl]") {
    const std::string stoat_command = "../bin/stoat";
    const std::string output_dir = "../output_eqtl";
    const std::string expected_dir = "../tests/test_data/expected_output/vcf/output_eqtl";
    const std::string expected_dir_covar = "../tests/test_data/expected_output/vcf/output_eqtl_covar";
    const std::string data_path = "../tests/test_data/input_data/eqtl";
    const std::string phenotype = "eqtl";

    SECTION("Without covariate") {
        REQUIRE(run_test(stoat_command, output_dir, expected_dir, data_path, phenotype, false));
    }

    SECTION("With covariate") {
        REQUIRE(run_test(stoat_command, output_dir, expected_dir_covar, data_path, phenotype, true));
    }
}

TEST_CASE("Output simple nested chain", "[detangle]") {
/*
                        6
                      /   \
             2       5 ----7    9
           /   \   /         \ /  \
         1       4  ----------8---10
           \   /
             3

*/


    const std::string output_dir = "../output_binary";
    const std::string graph_base = "../tests/test_data/test_graphs/simple_nested_chain";
    const std::string samples_filename = "./samples.tsv";
    const std::string reference_filename = "./references.tsv";
    std::string vcf_filename = "./test.vcf";

    std::ofstream vcf_out;
    vcf_out.open(vcf_filename);
    vcf_out << "##fileformat=VCFv4.2" << std::endl;
    vcf_out << "##FILTER=<ID=PASS,Description=\"All filters passed\">" << std::endl;
    vcf_out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
    vcf_out << "##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the snarl tree (0=top level)\">" << std::endl;
    vcf_out << "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph\">" << std::endl;
    vcf_out << "##contig=<ID=path0#0#path0,length=100>" << std::endl;
    vcf_out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2" << std::endl;
    vcf_out << "path0#0#path0\t0\t>1>4\tCGT\tCCT,A\t60\t.\tLV=0;AT=>1>2>4,>1>3>4\tGT\t0/0\t0/1" << std::endl;
    // 0 and 1 both take the insertion but 0 takes the inner insertion too
    vcf_out << "path0#0#path0\t3\t>4>8\tCGT\tCCT,A\t60\t.\tLV=0;AT=>4>5>6>7>8,>4>5>7>8,>4>8\tGT\t0/1\t2/0" << std::endl;
    vcf_out << "path0#0#path0\t4\t>5>7\tCGT\tCCT,A\t60\t.\tLV=1;AT=>5>6>7,>5>7\tGT\t1/0\t1/0" << std::endl;
    // This isn't actually on this reference so it will get skipped
    vcf_out << "path0#0#path0\t5\t>8>10\tCGT\tCCT,A\t60\t.\tLV=0;AT=>8>9>10,>8>10\tGT\t0/1\t1/0" << std::endl;
    vcf_out.close();


    std::vector<std::string> samples_of_interest = {"S1"};
    std::vector<std::string> other_samples = {"S2"};

    string write_cmd = "echo \"SAMPLE\tPHENO\" > " + samples_filename;
    int ignore = std::system(write_cmd.c_str());
    for (auto sample : samples_of_interest) {
        write_cmd = "echo \"" + sample + "\t1\" >> " + samples_filename;
        ignore = std::system(write_cmd.c_str());
    }
    for (auto sample : other_samples) {
        write_cmd = "echo \"" + sample + "\t0\" >> " + samples_filename;
        ignore = std::system(write_cmd.c_str());
    }

    // Write the reference file
    write_cmd = "echo path0#0#path0 > " + reference_filename;
    ignore = std::system(write_cmd.c_str());

    SECTION("Test without untangling") {
        // Make the snarl file

        std::string cmd = (std::string)"stoat vcf"
            + " -g " + graph_base + ".hg"
            + " -d " + graph_base + ".dist"
            + " -r " + reference_filename
            + " -v " + vcf_filename
            + " -b " + samples_filename
            + " -o " + output_dir;
        std::cerr << "Run command " << cmd << std::endl;
        int command_output = std::system(cmd.c_str());

        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE(false);
        }

        // Snarls should now be in output_dir/snarl_info.tsv
        // Genotypes should be in output_dir/snarl_genotypes.tsv
        // Final values should be in output_dir/stoat.assoc.pvalues.tsv
        ifstream in_genotypes;
        in_genotypes.open(output_dir + "/snarl_genotypes.tsv");
        std::string line; 
        while (std::getline(in_genotypes, line)) {
            if (line.at(0) == '#') {
                //skip the header
                continue;
            }
            // Get the start node, which is the first thing in the tab separated line
            std::stringstream linestream(line);
            std::string first_node;
            std::getline(linestream, first_node, '\t');
            // Get to the 10th item, which is the first genotype
            std::string item;
            std::getline(linestream, item, '\t'); // got end node
            std::getline(linestream, item, '\t'); // got ref
            std::getline(linestream, item, '\t'); // got start offset
            std::getline(linestream, item, '\t'); // got end offset
            std::getline(linestream, item, '\t'); // got depth
            std::getline(linestream, item, '\t'); // got allele lengths
            std::string walks;
            std::getline(linestream, walks, '\t'); // got walks
            std::getline(linestream, item, '\t'); // got sequences

            if (first_node == ">1" || first_node == "<4") {
                // Get the walks. should be >1>2>4,1>3>4
                std::stringstream walkstream(walks);
                std::string walk;
                std::getline(walkstream, walk, ',');
                bool flipped = false;
                if (walk == ">1>2>4") {
                    std::getline(walkstream, walk, ',');
                    REQUIRE(walk == ">1>3>4");
                } else if (walk == ">1>3>4") {
                    flipped = true;
                    std::getline(walkstream, walk, ',');
                    REQUIRE(walk == ">1>2>4");
                } else {
                    // Bad walk
                    std::cerr << "Walk shouldn't exist " << walk << std::endl;
                    REQUIRE(false);
                }
                // should be only two walks
                REQUIRE(!std::getline(walkstream, walk, ','));

                // Get the genotypes. Should be 0/0 0/1, assuming the same order
                std::string genotype;
                std::getline(linestream, genotype, '\t');
                REQUIRE(genotype == (flipped ? "1" : "0"));
                std::getline(linestream, genotype, '\t');
                REQUIRE(genotype == (flipped ? "1" : "0"));
                std::getline(linestream, genotype, '\t');
                REQUIRE(genotype == (flipped ? "1" : "0"));
                std::getline(linestream, genotype, '\t');
                REQUIRE(genotype == (flipped ? "0" : "1"));
            } else if (first_node == ">4" || first_node == "<8") {
                // Get the walks. should be >4>5>6>7>8,>4>5>7>8,>4>8
                // genotypes should be 0/0 1/0
                std::stringstream walkstream(walks);
                std::string walk;
                std::vector<size_t> walk_indices(2,0);

                for (size_t i = 0 ; i < 3 ; i++) {
                    std::getline(walkstream, walk, ',');
                    if (walk == ">4>5>0>7>8") {
                        walk_indices[0] = i;
                    } else if (walk == ">4>8") {
                        walk_indices[1] = i;
                    } else {
                        // Bad walk
                        std::cerr << "Walk shouldn't exist " << walk << std::endl;
                        REQUIRE(false);
                    }
                }
                // should be only two walks
                REQUIRE(!std::getline(walkstream, walk, ','));

                // Get the genotypes. Should be 0/1 2/0, assuming the same order
                std::string genotype;
                std::getline(linestream, genotype, '\t');
                REQUIRE(genotype == std::to_string(walk_indices[0]));
                std::getline(linestream, genotype, '\t');
                REQUIRE(genotype == std::to_string(walk_indices[0]));
                std::getline(linestream, genotype, '\t');
                REQUIRE(genotype == std::to_string(walk_indices[1]));
                std::getline(linestream, genotype, '\t');
                REQUIRE(genotype == std::to_string(walk_indices[0]));
            } else if (first_node == ">5" || first_node == "<7") {
                // Get the walks. should be >5>6>7, >5>7
                // genotypes should be 0/1 ./0
                std::stringstream walkstream(walks);
                std::string walk;
                std::getline(walkstream, walk, ',');
                bool flipped = false;
                if (walk == ">5>6>7") {
                    std::getline(walkstream, walk, ',');
                    REQUIRE(walk == ">5>7");
                } else if (walk == ">5>7") {
                    flipped = true;
                    std::getline(walkstream, walk, ',');
                    REQUIRE(walk == ">5>6>7");
                } else {
                    // Bad walk
                    std::cerr << "Walk shouldn't exist " << walk << std::endl;
                    REQUIRE(false);
                }
                // should be only two walks
                REQUIRE(!std::getline(walkstream, walk, ','));

                // Get the genotypes. Should be 0/1 ./0, assuming the same order
                std::string genotype;
                std::getline(linestream, genotype, '\t');
                REQUIRE(genotype == (flipped ? "1" : "0"));
                std::getline(linestream, genotype, '\t');
                REQUIRE(genotype == (flipped ? "0" : "1"));
                std::getline(linestream, genotype, '\t');
                REQUIRE(genotype == ".");
                std::getline(linestream, genotype, '\t');
                REQUIRE(genotype == (flipped ? "1" : "0"));
            } else if (first_node == ">8" || first_node == "<10") {
                // This doesn't matter
            } else {
                // Snarl that shouldn't exist
                std::cerr << "Snarl that shouldn't exist starting at node " << first_node << std::endl;
                REQUIRE(false);
            }
        }

        in_genotypes.close();
        clean_output_dir(output_dir);
    }
    SECTION("Test untangling") {
        // Make the snarl file

        std::string cmd = (std::string)"stoat vcf"
            + " -g " + graph_base + ".hg"
            + " -d " + graph_base + ".dist"
            + " -r " + reference_filename
            + " -v " + vcf_filename
            + " -b " + samples_filename
            + " -o " + output_dir
            + " -R";
        std::cerr << "Run command " << cmd << std::endl;
        int command_output = std::system(cmd.c_str());

        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE(false);
        }

        // Snarls should now be in output_dir/snarl_info.tsv
        // Genotypes should be in output_dir/snarl_genotypes.tsv
        // Final values should be in output_dir/stoat.assoc.pvalues.tsv
        ifstream in_genotypes;
        in_genotypes.open(output_dir + "/snarl_genotypes.tsv");
        std::string line; 
        while (std::getline(in_genotypes, line)) {
            if (line.at(0) == '#') {
                //skip the header
                continue;
            }
            // Get the start node, which is the first thing in the tab separated line
            std::stringstream linestream(line);
            std::string first_node;
            std::getline(linestream, first_node, '\t');
            // Get to the 10th item, which is the first genotype
            std::string item;
            std::getline(linestream, item, '\t'); // got end node
            std::getline(linestream, item, '\t'); // got ref
            std::getline(linestream, item, '\t'); // got start offset
            std::getline(linestream, item, '\t'); // got end offset
            std::getline(linestream, item, '\t'); // got depth
            std::getline(linestream, item, '\t'); // got allele lengths
            std::string walks;
            std::getline(linestream, walks, '\t'); // got walks
            std::getline(linestream, item, '\t'); // got sequences

            if (first_node == ">1" || first_node == "<4") {
                // Get the walks. should be >1>2>4,1>3>4
                std::stringstream walkstream(walks);
                std::string walk;
                std::getline(walkstream, walk, ',');
                bool flipped = false;
                if (walk == ">1>2>4") {
                    std::getline(walkstream, walk, ',');
                    REQUIRE(walk == ">1>3>4");
                } else if (walk == ">1>3>4") {
                    flipped = true;
                    std::getline(walkstream, walk, ',');
                    REQUIRE(walk == ">1>2>4");
                } else {
                    // Bad walk
                    std::cerr << "Walk shouldn't exist " << walk << std::endl;
                    REQUIRE(false);
                }
                // should be only two walks
                REQUIRE(!std::getline(walkstream, walk, ','));

                // Get the genotypes. Should be 0/0 0/1, assuming the same order
                std::string genotype;
                std::getline(linestream, genotype, '\t');
                REQUIRE(genotype == (flipped ? "1" : "0"));
                std::getline(linestream, genotype, '\t');
                REQUIRE(genotype == (flipped ? "1" : "0"));
                std::getline(linestream, genotype, '\t');
                REQUIRE(genotype == (flipped ? "1" : "0"));
                std::getline(linestream, genotype, '\t');
                REQUIRE(genotype == (flipped ? "0" : "1"));
            } else if (first_node == ">4" || first_node == "<8") {
                // Get the walks. should be >4>5>0>5>8 and >4>8
                // genotypes should be 0/0 1/0
                std::stringstream walkstream(walks);
                std::string walk;
                std::getline(walkstream, walk, ',');
                bool flipped = false;
                if (walk == ">4>5>0>7>8") {
                    std::getline(walkstream, walk, ',');
                    REQUIRE(walk == ">4>8");
                } else if (walk == ">4>8") {
                    flipped = true;
                    std::getline(walkstream, walk, ',');
                    REQUIRE(walk == ">4>5>0>7>8");
                } else {
                    // Bad walk
                    std::cerr << "Walk shouldn't exist " << walk << std::endl;
                    REQUIRE(false);
                }
                // should be only two walks
                REQUIRE(!std::getline(walkstream, walk, ','));

                // Get the genotypes. Should be 0/0 1/0, assuming the same order
                std::string genotype;
                std::getline(linestream, genotype, '\t');
                REQUIRE(genotype == (flipped ? "1" : "0"));
                std::getline(linestream, genotype, '\t');
                REQUIRE(genotype == (flipped ? "1" : "0"));
                std::getline(linestream, genotype, '\t');
                REQUIRE(genotype == (flipped ? "0" : "1"));
                std::getline(linestream, genotype, '\t');
                REQUIRE(genotype == (flipped ? "1" : "0"));
            } else if (first_node == ">5" || first_node == "<7") {
                // Get the walks. should be >5>6>7, >5>7
                // genotypes should be 0/1 ./0
                std::stringstream walkstream(walks);
                std::string walk;
                std::getline(walkstream, walk, ',');
                bool flipped = false;
                if (walk == ">5>6>7") {
                    std::getline(walkstream, walk, ',');
                    REQUIRE(walk == ">5>7");
                } else if (walk == ">5>7") {
                    flipped = true;
                    std::getline(walkstream, walk, ',');
                    REQUIRE(walk == ">5>6>7");
                } else {
                    // Bad walk
                    std::cerr << "Walk shouldn't exist " << walk << std::endl;
                    REQUIRE(false);
                }
                // should be only two walks
                REQUIRE(!std::getline(walkstream, walk, ','));

                // Get the genotypes. Should be 1/0 ./0, assuming the same order
                std::string genotype;
                std::getline(linestream, genotype, '\t');
                REQUIRE(genotype == (flipped ? "0" : "1"));
                std::getline(linestream, genotype, '\t');
                REQUIRE(genotype == (flipped ? "1" : "0"));
                std::getline(linestream, genotype, '\t');
                REQUIRE(genotype == ".");
                std::getline(linestream, genotype, '\t');
                REQUIRE(genotype == (flipped ? "1" : "0"));
            } else if (first_node == ">8" || first_node == "<10") {
                // This doesn't matter
            } else {
                // Snarl that shouldn't exist
                std::cerr << "Snarl that shouldn't exist starting at node " << first_node << std::endl;
                REQUIRE(false);
            }
        }

        in_genotypes.close();
        clean_output_dir(output_dir);
    }


    clean_output_dir(output_dir);
    //fs::remove(samples_filename);
    //fs::remove(reference_filename);
    //fs::remove(vcf_filename);
}

