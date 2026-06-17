#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "compare_files_utils.hpp"
#include "load_tables.hpp"

namespace fs = std::filesystem;
using namespace std;


TEST_CASE("Output loop with snarl", "[bh]") {
    /*
              --------
             |   2    |
             \ / \    /
         0 ---1---3--4----5

         Node sequences:
         "AAAAAAAAAA", "A", "G", "C", "T",  "AAAAAAAAA"

         Paths:
         {0, 1, 2, 3, 4, 5}, {0, 1, 3, 4, 1, 3, 4, 5}, {0, 1, 2, 3, 4, 1, 3, 4, 5}
     */


    const std::string output_dir = "../output_binary";
    const std::string graph_base = "../tests/test_data/test_graphs/loop_with_indel";

    std::vector<std::string> samples_of_interest = {"path1", "path2"};
    std::vector<std::string> other_samples = {"path0"};

    const std::string samples_file = "./samples.tsv";
    const std::string refs_file = "./refs.tsv";
    
    string write_cmd = "echo \"SAMPLE\tPHENO\" > " + samples_file;
    int ignore = std::system(write_cmd.c_str());
    for (auto sample : samples_of_interest) {
        write_cmd = "echo \"" + sample + "\t1\" >> " + samples_file;
        ignore = std::system(write_cmd.c_str());
    }
    for (auto sample : other_samples) {
        write_cmd = "echo \"" + sample + "\t0\" >> " + samples_file;
        ignore = std::system(write_cmd.c_str());
    }
    write_cmd = "echo path0 > " + refs_file; 
    ignore = std::system(write_cmd.c_str());

    // Run stoat to get 
    std::string cmd = "../bin/stoat graph -u";
    cmd += " -g " + graph_base + ".hg"
        + " -d " + graph_base + ".dist"
        + " -L"
        + " -r " + refs_file
        + " --output " + output_dir;

    std::cout << "Command run : \n" << cmd << std::endl;
    int command_output = std::system(cmd.c_str());

    std::string cmd_test = "../bin/stoat test -u";
    cmd_test += " -g " + output_dir + "/snarl_genotypes.tsv"
        + " -p " + samples_file
        + " -m chi2"
        + " --output " + output_dir;

    std::cout << "Command run : \n" << cmd_test << std::endl;
    command_output = std::system(cmd_test.c_str());
    
    SECTION("Bh correct") {

        std::string rewrite_cmd = "../bin/stoat BHcorrect";
        rewrite_cmd += " -t " + output_dir+"/stoat.assoc.pvalues.tsv"; 

        std::cout << "Command run : \n" << cmd << std::endl;
        int command_output = std::system(rewrite_cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << rewrite_cmd << "\n";
            REQUIRE( false);
        }
    }


    fs::remove(samples_file);
    fs::remove(refs_file);
    clean_output_dir(output_dir);
}
TEST_CASE("Output loop with snarl with bgzipped output", "[bh]") {
    /*
              --------
             |   2    |
             \ / \    /
         0 ---1---3--4----5

         Node sequences:
         "AAAAAAAAAA", "A", "G", "C", "T",  "AAAAAAAAA"

         Paths:
         {0, 1, 2, 3, 4, 5}, {0, 1, 3, 4, 1, 3, 4, 5}, {0, 1, 2, 3, 4, 1, 3, 4, 5}
     */


    const std::string output_dir = "../output_binary";
    const std::string graph_base = "../tests/test_data/test_graphs/loop_with_indel";

    std::vector<std::string> samples_of_interest = {"path1", "path2"};
    std::vector<std::string> other_samples = {"path0"};

    const std::string samples_file = "./samples.tsv";
    const std::string refs_file = "./refs.tsv";
    
    string write_cmd = "echo \"SAMPLE\tPHENO\" > " + samples_file;
    int ignore = std::system(write_cmd.c_str());
    for (auto sample : samples_of_interest) {
        write_cmd = "echo \"" + sample + "\t1\" >> " + samples_file;
        ignore = std::system(write_cmd.c_str());
    }
    for (auto sample : other_samples) {
        write_cmd = "echo \"" + sample + "\t0\" >> " + samples_file;
        ignore = std::system(write_cmd.c_str());
    }
    write_cmd = "echo path0 > " + refs_file; 
    ignore = std::system(write_cmd.c_str());

    // Run stoat to get 
    std::string cmd = "../bin/stoat graph -u";
    cmd += " -g " + graph_base + ".hg"
        + " -d " + graph_base + ".dist"
        + " -L"
        + " -r " + refs_file
        + " --output " + output_dir;

    std::cout << "Command run : \n" << cmd << std::endl;
    int command_output = std::system(cmd.c_str());

    std::string cmd_test = "../bin/stoat test";
    cmd_test += " -g " + output_dir + "/snarl_genotypes.tsv"
        + " -p " + samples_file
        + " -m chi2"
        + " --output " + output_dir;

    std::cout << "Command run : \n" << cmd_test << std::endl;
    command_output = std::system(cmd_test.c_str());
    
    SECTION("BHcorrect") {

        std::string rewrite_cmd = "../bin/stoat BHcorrect";
        rewrite_cmd +=
            + " -t " + output_dir+"/stoat.assoc.pvalues.tsv.gz"; 

        std::cout << "Command run : \n" << cmd << std::endl;
        int command_output = std::system(rewrite_cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << rewrite_cmd << "\n";
            REQUIRE( false);
        }


    }



    fs::remove(samples_file);
    fs::remove(refs_file);
    clean_output_dir(output_dir);
}


