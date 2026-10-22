#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "compare_files_utils.hpp"
#include "load_tables.hpp"

namespace fs = std::filesystem;
using namespace std;


TEST_CASE( "simple nested chain","[get-tree]" ) {

    /*
                       5
                     /   \
            1       4 ----6    8
          /   \   /         \ / \
        0       3  ----------7---9
          \   /
            2

   */

 
    const std::string graph_base = "../tests/test_data/test_graphs/simple_nested_chain";
    const std::string outfile = "test_snarls.tsv";   

    SECTION("all snarls") {

        // Run stoat to get 
        std::string cmd = "../bin/stoat get-tree";
        cmd += " -g " + graph_base + ".hg"
            + " -d " + graph_base + ".dist"
            + " -p path0#0#path0:0-5 > " 
            + outfile;

        std::cout << "Command run : \n" << cmd << std::endl;
        int command_output = std::system(cmd.c_str());

        // Go through each file line by line and check that they are identical

        std::unordered_map<std::string, std::string> snarl_to_parent;
        ifstream instream;
        instream.open(outfile);

        std::string line;
        while (std::getline(instream, line)) {
            std::stringstream ss (line);
            std::string snarl;
            std::string parent;
            std::getline(ss, snarl, '\t');
            std::getline(ss, parent, '\t');
            snarl_to_parent.emplace(snarl, parent);
        }
        instream.close();

        REQUIRE(snarl_to_parent.size() == 3);

        // Top-level snarl 1-4
        if (snarl_to_parent.count(">1<4") == 1) {
            REQUIRE(snarl_to_parent[">1<4"] == ".");
        } else if (snarl_to_parent.count("<4>1") == 1) {
            REQUIRE(snarl_to_parent["<4>1"] == ".");
        } else {
            REQUIRE(false);
        }
        // Top-level snarl 4-8
        if (snarl_to_parent.count(">4<8") == 1) {
            REQUIRE(snarl_to_parent[">4<8"] == ".");
        } else if (snarl_to_parent.count("<8>4") == 1) {
            REQUIRE(snarl_to_parent["<8>4"] == ".");
        } else {
            REQUIRE(false);
        }

        // Snarl 5-7 child of 4-8
        if (snarl_to_parent .count(">5<7") == 1) {
            REQUIRE((snarl_to_parent[">5<7"] == ">4<8" || snarl_to_parent[">5<7"] == "<8>4"));
        } else if (snarl_to_parent.count("<4>2") == 1) {
            REQUIRE((snarl_to_parent[">7<5"] == ">4<8" || snarl_to_parent[">7<5"] == "<8>4"));
        } else {
            REQUIRE(false);
        }

    }
    SECTION("only nested snarl") {

        // Run stoat to get 
        std::string cmd = "../bin/stoat get-tree";
        cmd += " -g " + graph_base + ".hg"
            + " -d " + graph_base + ".dist"
            + " -p path0#0#path0:3-4 > " 
            + outfile;

        std::cout << "Command run : \n" << cmd << std::endl;
        int command_output = std::system(cmd.c_str());

        // Go through each file line by line and check that they are identical

        std::unordered_map<std::string, std::string> snarl_to_parent;
        ifstream instream;
        instream.open(outfile);

        std::string line;
        while (std::getline(instream, line)) {
            std::stringstream ss (line);
            std::string snarl;
            std::string parent;
            std::getline(ss, snarl, '\t');
            std::getline(ss, parent, '\t');
            snarl_to_parent.emplace(snarl, parent);
        }
        instream.close();

        REQUIRE(snarl_to_parent.size() == 1);


        // Snarl 5-7 child of 4-8 but since we didn't include 4-8 it doesn't have a parent
        if (snarl_to_parent .count(">5<7") == 1) {
            REQUIRE(snarl_to_parent[">5<7"] == ".");
        } else if (snarl_to_parent.count("<4>2") == 1) {
            REQUIRE(snarl_to_parent[">7<5"] == ".");
        } else {
            REQUIRE(false);
        }

    }



    fs::remove(outfile);

}

TEST_CASE("Loop with snarl", "[get-tree]") {
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

 
    const std::string graph_base = "../tests/test_data/test_graphs/loop_with_indel";
    const std::string outfile = "test_snarls.tsv";   

    SECTION("get-tree") {



        // Run stoat to get 
        std::string cmd = "../bin/stoat get-tree";
        cmd += " -g " + graph_base + ".hg"
            + " -d " + graph_base + ".dist"
            + " -p path0:0-20 > " 
            + outfile;

        std::cout << "Command run : \n" << cmd << std::endl;
        int command_output = std::system(cmd.c_str());

        // Go through each file line by line and check that they are identical

        std::unordered_map<std::string, std::string> snarl_to_parent;
        ifstream instream;
        instream.open(outfile);

        std::string line;
        while (std::getline(instream, line)) {
            std::stringstream ss (line);
            std::string snarl;
            std::string parent;
            std::getline(ss, snarl, '\t');
            std::getline(ss, parent, '\t');
            snarl_to_parent.emplace(snarl, parent);
        }
        instream.close();

        // There should only be two snarls, the outer one with no parent and the inner one with the outer one as a parent
        REQUIRE(snarl_to_parent.size() == 2);

        if (snarl_to_parent.count(">1<6") == 1) {
            REQUIRE(snarl_to_parent[">1<6"] == ".");
        } else if (snarl_to_parent.count("<6>1") == 1) {
            REQUIRE(snarl_to_parent["<6>1"] == ".");
        } else {
            REQUIRE(false);
        }

        if (snarl_to_parent .count(">2<4") == 1) {
            REQUIRE((snarl_to_parent[">2<4"] == ">1<6" || snarl_to_parent[">2<4"] == "<6>1"));
        } else if (snarl_to_parent.count("<4>2") == 1) {
            REQUIRE((snarl_to_parent[">4<2"] == ">1<6" || snarl_to_parent[">4<2"] == "<6>1"));
        } else {
            REQUIRE(false);
        }

    }



    fs::remove(outfile);
}
TEST_CASE("Giant unverified graph just to see if it crashes", "[get-tree]") {
    // Just check that this runs and produces some output

    const std::string graph_base = "../tests/test_data/input_data/binary/pg.full";

    SECTION("test") {


        // Run stoat to get 
        std::string cmd = "../bin/stoat get-tree";
        cmd += " -g " + graph_base + ".pg"
            + " -d " + graph_base + ".dist"
            + " -p ref:0-248073 > /dev/null";

        std::cout << "Command run : \n" << cmd << std::endl;
        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE(false);
        }

    }
    SECTION("test subrange") {


        // Run stoat to get 
        std::string cmd = "../bin/stoat get-tree";
        cmd += " -g " + graph_base + ".pg"
            + " -d " + graph_base + ".dist"
            + " -p ref:2073-183728 > /dev/null";

        std::cout << "Command run : \n" << cmd << std::endl;
        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE(false);
        }

    }
}


