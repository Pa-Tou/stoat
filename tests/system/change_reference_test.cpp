#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "compare_files_utils.hpp"
#include "load_tables.hpp"

namespace fs = std::filesystem;
using namespace std;


TEST_CASE("Output loop with snarl", "[graph]") {
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
    
    SECTION("Change ref to same path") {

        std::string rewrite_cmd = "../bin/stoat change-ref";
        rewrite_cmd += " -g " + graph_base + ".hg"
            + " -d " + graph_base + ".dist"
            + " -T " + output_dir+"/stoat.assoc.pvalues.tsv"
            + " -r " + refs_file 
            + " -o " + output_dir+"/stoat.assoc.pvalues.path0.tsv"; 

        std::cout << "Command run : \n" << cmd << std::endl;
        int command_output = std::system(rewrite_cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << rewrite_cmd << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/stoat.assoc.pvalues.path0.tsv"));

        // Go through each file line by line and check that they are identical

        ifstream in_original;
        ifstream in_new;
        in_original.open(output_dir+"/stoat.assoc.pvalues.tsv");
        in_new.open(output_dir+"/stoat.assoc.pvalues.path0.tsv");

        std::string line_original;
        std::string line_new;
        while (std::getline(in_original, line_original)) {
            REQUIRE(std::getline(in_new, line_new));
            REQUIRE(line_original == line_new);
        }

        in_original.close();
        in_new.close();

    }
    SECTION("Change ref to path1") {

        write_cmd = "echo path1 > " + refs_file; 
        ignore = std::system(write_cmd.c_str());

        std::string rewrite_cmd = "../bin/stoat change-ref";
        rewrite_cmd += " -g " + graph_base + ".hg"
            + " -d " + graph_base + ".dist"
            + " -T " + output_dir+"/stoat.assoc.pvalues.tsv"
            + " -r " + refs_file 
            + " -o " + output_dir+"/stoat.assoc.pvalues.path1.tsv"; 

        std::cout << "Command run : \n" << rewrite_cmd << std::endl;
        int command_output = std::system(rewrite_cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << rewrite_cmd << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/stoat.assoc.pvalues.path1.tsv"));

        // Go through each file line by line and check that they are identical except for the coordinates

        ifstream in_original;
        ifstream in_new;
        in_original.open(output_dir+"/stoat.assoc.pvalues.tsv");
        in_new.open(output_dir+"/stoat.assoc.pvalues.path1.tsv");

        std::string line_original;
        std::string line_new;
        // Check that the header is the same
        std::getline(in_original, line_original);
        std::getline(in_new, line_new);
        REQUIRE(line_original == line_new);
        while (std::getline(in_original, line_original)) {
            REQUIRE(std::getline(in_new, line_new));


            // Get the line as a vector of strings
            std::vector<std::string> original_values;
            std::stringstream original_linestream(line_original);
            std::string col;
            while (std::getline(original_linestream, col, '\t')) {
                original_values.push_back(col);
            }
            std::vector<std::string> new_values;
            std::stringstream new_linestream(line_new);
            while (std::getline(new_linestream, col, '\t')) {
                new_values.push_back(col);
            }
            REQUIRE(original_values.size() == new_values.size());

            REQUIRE(new_values[0] == "path1");
            if (original_values[3] == "1_6" || original_values[3] == "6_1") {
                REQUIRE(new_values[1] == "10");
                REQUIRE(new_values[2] == "16");
            } else if (original_values[3] == "2_4" || original_values[3] == "4_2") {
                REQUIRE(((new_values[1] == "11" && new_values[2] == "11") ||
                        (new_values[1] == "14" && new_values[2] == "14")));
            }


            // Now check that everything else is the same
            for (size_t i = 3 ; i < original_values.size() ; i++) {
                REQUIRE(original_values[i] == new_values[i]);
            }

        }

        in_original.close();
        in_new.close();

    }


    //fs::remove(samples_file);
    //fs::remove(refs_file);
    //clean_output_dir(output_dir);
}


