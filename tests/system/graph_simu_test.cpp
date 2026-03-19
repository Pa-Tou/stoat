#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "compare_files_utils.hpp"
#include "load_tables.hpp"

namespace fs = std::filesystem;
using namespace std;


TEST_CASE("Giant unverified binary association tests graph plus test", "[test]") {
    // Just check that this runs and produces some output

    const std::string output_dir = "../output_binary";
    const std::string data_path = "../tests/test_data/input_data/binary";
    const std::string graph_base = "pg.full";

    SECTION("Test tsv output saving snarls") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph -u";

        cmd +=" -g " + data_path + "/" + graph_base + ".pg"
            + " -d " + data_path + "/" + graph_base + ".dist"
            + " -L -r ref --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE(false);
        }

        REQUIRE(std::filesystem::exists(output_dir + "/snarl_genotypes.tsv"));
        std::ifstream snarlsfile;
        snarlsfile.open(output_dir + "/snarl_genotypes.tsv");
        REQUIRE(snarlsfile.peek() != std::ifstream::traits_type::eof());

        size_t line_count = 0;
        std::string line;
        while (std::getline(snarlsfile, line)) {
            // Only start counting snarls after proper header
            if (line[0] != '#') {
                line_count++;
            }
        }
        snarlsfile.close();
        REQUIRE(line_count==1524);

        // TODO: Add something that actually checks this
        //bool passed = compare_output_dirs(output_dir, expected_dir);
        //REQUIRE(passed);

    }
    SECTION("Test tsv output saving snarls multithreaded") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph -u";

        cmd +=" -g " + data_path + "/" + graph_base + ".pg"
            + " -d " + data_path + "/" + graph_base + ".dist"
            + " -L -r ref --output " + output_dir 
            + " -t 4";

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE(false);
        }

        REQUIRE(std::filesystem::exists(output_dir + "/snarl_genotypes.tsv"));
        std::ifstream snarlsfile;
        snarlsfile.open(output_dir + "/snarl_genotypes.tsv");
        REQUIRE(snarlsfile.peek() != std::ifstream::traits_type::eof());

        size_t line_count = 0;
        std::string line;
        while (std::getline(snarlsfile, line)) {
            // Only start counting snarls after proper header
            if (line[0] != '#') {
                line_count++;
            }
        }
        snarlsfile.close();
        REQUIRE(line_count==1524);

        // TODO: Add something that actually checks this
        //bool passed = compare_output_dirs(output_dir, expected_dir);
        //REQUIRE(passed);

    }

    clean_output_dir(output_dir);
}

TEST_CASE("Output simple nested chain stats", "[test]") {
    const std::string output_dir = "../output_binary";
    const std::string graph_base = "../tests/test_data/test_graphs/simple_nested_chain";

    SECTION("Test stoat graph with lengths") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph -u";
           cmd += " -g " + graph_base + ".hg"
            + " -d " + graph_base + ".dist"
            + " -L"
            + " -r path0 -V 4"
            + " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/snarl_genotypes.tsv"));
       

        std::vector<std::string> lengths(2, "1");
        std::vector<std::string> paths;
        paths.emplace_back(">1>2>4");
        paths.emplace_back(">1>3>4");
        std::vector<std::string> sequences;
        std::unordered_map<std::string, size_t> allele_assignment;
        allele_assignment["path0#0#path0"] = 0;
        allele_assignment["path1#0#path1#0"] = 0;
        allele_assignment["path2#"] = 1;
        allele_assignment["path3#"] = 1;
        snarl_genotype_values_t truth1 ({">1", "<4",
                                       (std::string)"path0#0#path0", (size_t)1, (size_t)2, 
                                       (size_t) 1,
                                       lengths,
                                       paths, 
                                       sequences,
                                       allele_assignment});

        lengths.clear();
        lengths.emplace_back("0");
        lengths.emplace_back("2/3");
        paths.clear();
        paths.emplace_back(">4>8");
        paths.emplace_back(">4>5>0>7>8");
        allele_assignment.clear();
        allele_assignment["path0#0#path0"] = 1;
        allele_assignment["path1#0#path1#0"] = 1;
        allele_assignment["path2#"] = 0;
        allele_assignment["path3#"] = 1;
        snarl_genotype_values_t truth2 ({">4", "<8",
                                       (std::string)"path0#0#path0", (size_t)3, (size_t)6, 
                                       (size_t) 1,
                                       lengths,
                                       paths, 
                                       sequences,
                                       allele_assignment});

        lengths.clear();
        lengths.emplace_back("0");
        lengths.emplace_back("1");
        paths.clear();
        paths.emplace_back(">5>7");
        paths.emplace_back(">5>6>7");
        allele_assignment.clear();
        allele_assignment["path0#0#path0"] = 1;
        allele_assignment["path1#0#path1#0"] = 0;
        allele_assignment["path2#"] = std::numeric_limits<size_t>::max();
        allele_assignment["path3#"] = 0;
        snarl_genotype_values_t truth3 ({">5", "<7",
                                       (std::string)"path0#0#path0", (size_t)4, (size_t)5, 
                                       (size_t) 2,
                                       lengths,
                                       paths, 
                                       sequences,
                                       allele_assignment});
        
        lengths.clear();
        paths.clear();
        allele_assignment.clear();
        allele_assignment["path0#0#path0"] = std::numeric_limits<size_t>::max();
        allele_assignment["path1#0#path1#0"] = std::numeric_limits<size_t>::max();
        allele_assignment["path2#"] = std::numeric_limits<size_t>::max();
        allele_assignment["path3#"] = std::numeric_limits<size_t>::max();
        snarl_genotype_values_t truth4 ({">8", "<10",
                                       (std::string)".", (size_t)0, (size_t)0, 
                                       (size_t) 1,
                                       lengths,
                                       paths, 
                                       sequences,
                                       allele_assignment});

        std::vector<snarl_genotype_values_t> loaded_snarls= load_genotype_file(output_dir+"/snarl_genotypes.tsv");
        std::vector<bool> found_snarls (4, false);
        for (const snarl_genotype_values_t& test : loaded_snarls) {
            if (test == truth1) {
                found_snarls[0] = true;
            } else if (test == truth2) {
                found_snarls[1] = true;
            }  else if (test == truth3) {
                found_snarls[2] = true;
            }  else if (test == truth4) {
                found_snarls[3] = true;
            } else {
                std::cerr << "Unknown snarl " << std::endl << test << std::endl;
            }
        }

        if (!found_snarls[0]) {
            std::cerr << "Missing snarl " << std::endl << truth1 << std::endl;
            REQUIRE(false);
        }
        if (!found_snarls[1]) {
            std::cerr << "Missing snarl " << std::endl << truth2 << std::endl;
            REQUIRE(false);
        }
        if (!found_snarls[2]) {
            std::cerr << "Missing snarl " << std::endl << truth3 << std::endl;
            REQUIRE(false);
        }
        if (!found_snarls[3]) {
            std::cerr << "Missing snarl " << std::endl << truth4 << std::endl;
            REQUIRE(false);
        }

    }
    SECTION("Test stoat graph without lengths") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph -u";
           cmd += " -g " + graph_base + ".hg"
            + " -d " + graph_base + ".dist"
            + " -r path0 -V 4"
            + " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/snarl_genotypes.tsv"));
       

        std::vector<std::string> lengths;
        std::vector<std::string> paths;
        std::vector<std::string> sequences;
        std::unordered_map<std::string, size_t> allele_assignment;
        allele_assignment["path0#0#path0"] = 0;
        allele_assignment["path1#0#path1#0"] = 0;
        allele_assignment["path2#"] = 1;
        allele_assignment["path3#"] = 1;
        snarl_genotype_values_t truth1 ({">1", "<4",
                                       (std::string)"path0#0#path0", (size_t)1, (size_t)2, 
                                       (size_t) 1,
                                       lengths,
                                       paths, 
                                       sequences,
                                       allele_assignment});

        allele_assignment.clear();
        allele_assignment["path0#0#path0"] = 1;
        allele_assignment["path1#0#path1#0"] = 1;
        allele_assignment["path2#"] = 0;
        allele_assignment["path3#"] = 1;
        snarl_genotype_values_t truth2 ({">4", "<8",
                                       (std::string)"path0#0#path0", (size_t)3, (size_t)6, 
                                       (size_t) 1,
                                       lengths,
                                       paths, 
                                       sequences,
                                       allele_assignment});

        allele_assignment.clear();
        allele_assignment["path0#0#path0"] = 1;
        allele_assignment["path1#0#path1#0"] = 0;
        allele_assignment["path2#"] = std::numeric_limits<size_t>::max();
        allele_assignment["path3#"] = 0;
        snarl_genotype_values_t truth3 ({">5", "<7",
                                       (std::string)"path0#0#path0", (size_t)4, (size_t)5, 
                                       (size_t) 2,
                                       lengths,
                                       paths, 
                                       sequences,
                                       allele_assignment});
        
        allele_assignment.clear();
        allele_assignment["path0#0#path0"] = std::numeric_limits<size_t>::max();
        allele_assignment["path1#0#path1#0"] = std::numeric_limits<size_t>::max();
        allele_assignment["path2#"] = std::numeric_limits<size_t>::max();
        allele_assignment["path3#"] = std::numeric_limits<size_t>::max();
        snarl_genotype_values_t truth4 ({">8", "<10",
                                       (std::string)".", (size_t)0, (size_t)0, 
                                       (size_t) 1,
                                       lengths,
                                       paths, 
                                       sequences,
                                       allele_assignment});

        std::vector<snarl_genotype_values_t> loaded_snarls= load_genotype_file(output_dir+"/snarl_genotypes.tsv");
        std::vector<bool> found_snarls (4, false);
        for (const snarl_genotype_values_t& test : loaded_snarls) {
            if (test == truth1) {
                found_snarls[0] = true;
            } else if (test == truth2) {
                found_snarls[1] = true;
            }  else if (test == truth3) {
                found_snarls[2] = true;
            }  else if (test == truth4) {
                found_snarls[3] = true;
            } else {
                std::cerr << "Unknown snarl " << std::endl << test << std::endl;
            }
        }

        if (!found_snarls[0]) {
            std::cerr << "Missing snarl " << std::endl << truth1 << std::endl;
            REQUIRE(false);
        }
        if (!found_snarls[1]) {
            std::cerr << "Missing snarl " << std::endl << truth2 << std::endl;
            REQUIRE(false);
        }
        if (!found_snarls[2]) {
            std::cerr << "Missing snarl " << std::endl << truth3 << std::endl;
            REQUIRE(false);
        }
        if (!found_snarls[3]) {
            std::cerr << "Missing snarl " << std::endl << truth4 << std::endl;
            REQUIRE(false);
        }

    }

    SECTION("test stoat graph gbz") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph -u";
           cmd += " -g " + graph_base + ".gbz"
            + " -d " + graph_base + ".dist"
            + " -L"
            + " -r path0 -V 4"
            + " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }


        REQUIRE(std::filesystem::exists(output_dir+"/snarl_genotypes.tsv"));
       

        std::vector<std::string> lengths(2, "1");
        std::vector<std::string> paths;
        paths.emplace_back(">1>2>4");
        paths.emplace_back(">1>3>4");
        std::vector<std::string> sequences;
        std::unordered_map<std::string, size_t> allele_assignment;
        allele_assignment["path0#0#path0"] = 0;
        allele_assignment["path1#0#path1#0"] = 0;
        allele_assignment["path2#"] = 1;
        allele_assignment["path3#"] = 1;
        snarl_genotype_values_t truth1 ({">1", "<4",
                                       (std::string)"path0#0#path0", (size_t)1, (size_t)2, 
                                       (size_t) 1,
                                       lengths,
                                       paths, 
                                       sequences,
                                       allele_assignment});

        lengths.clear();
        lengths.emplace_back("0");
        lengths.emplace_back("2/3");
        paths.clear();
        paths.emplace_back(">4>8");
        paths.emplace_back(">4>5>0>7>8");
        allele_assignment.clear();
        allele_assignment["path0#0#path0"] = 1;
        allele_assignment["path1#0#path1#0"] = 1;
        allele_assignment["path2#"] = 0;
        allele_assignment["path3#"] = 1;
        snarl_genotype_values_t truth2 ({">4", "<8",
                                       (std::string)"path0#0#path0", (size_t)3, (size_t)6, 
                                       (size_t) 1,
                                       lengths,
                                       paths, 
                                       sequences,
                                       allele_assignment});

        lengths.clear();
        lengths.emplace_back("0");
        lengths.emplace_back("1");
        paths.clear();
        paths.emplace_back(">5>7");
        paths.emplace_back(">5>6>7");
        allele_assignment.clear();
        allele_assignment["path0#0#path0"] = 1;
        allele_assignment["path1#0#path1#0"] = 0;
        allele_assignment["path2#"] = std::numeric_limits<size_t>::max();
        allele_assignment["path3#"] = 0;
        snarl_genotype_values_t truth3 ({">5", "<7",
                                       (std::string)"path0#0#path0", (size_t)4, (size_t)5, 
                                       (size_t) 2,
                                       lengths,
                                       paths, 
                                       sequences,
                                       allele_assignment});
        
        lengths.clear();
        paths.clear();
        allele_assignment.clear();
        allele_assignment["path0#0#path0"] = std::numeric_limits<size_t>::max();
        allele_assignment["path1#0#path1#0"] = std::numeric_limits<size_t>::max();
        allele_assignment["path2#"] = std::numeric_limits<size_t>::max();
        allele_assignment["path3#"] = std::numeric_limits<size_t>::max();
        snarl_genotype_values_t truth4 ({">8", "<10",
                                       (std::string)".", (size_t)0, (size_t)0, 
                                       (size_t) 1,
                                       lengths,
                                       paths, 
                                       sequences,
                                       allele_assignment});

        std::vector<snarl_genotype_values_t> loaded_snarls= load_genotype_file(output_dir+"/snarl_genotypes.tsv");
        std::vector<bool> found_snarls (4, false);
        for (const snarl_genotype_values_t& test : loaded_snarls) {
            if (test == truth1) {
                found_snarls[0] = true;
            } else if (test == truth2) {
                found_snarls[1] = true;
            }  else if (test == truth3) {
                found_snarls[2] = true;
            }  else if (test == truth4) {
                found_snarls[3] = true;
            } else {
                std::cerr << "Unknown snarl " << std::endl << test << std::endl;
            }
        }

        if (!found_snarls[0]) {
            std::cerr << "Missing snarl " << std::endl << truth1 << std::endl;
            REQUIRE(false);
        }
        if (!found_snarls[1]) {
            std::cerr << "Missing snarl " << std::endl << truth2 << std::endl;
            REQUIRE(false);
        }
        if (!found_snarls[2]) {
            std::cerr << "Missing snarl " << std::endl << truth3 << std::endl;
            REQUIRE(false);
        }
        if (!found_snarls[3]) {
            std::cerr << "Missing snarl " << std::endl << truth4 << std::endl;
            REQUIRE(false);
        }


    }
    clean_output_dir(output_dir);
}

TEST_CASE("Output loop with snarl test", "[test]") {
    const std::string output_dir = "../output_binary";
    const std::string graph_base = "../tests/test_data/test_graphs/loop_with_indel";

    SECTION("stoat graph") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph -u";
           cmd += " -g " + graph_base + ".hg"
            + " -d " + graph_base + ".dist"
            + " -L"
            + " -r path0"
            + " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }
       

        REQUIRE(std::filesystem::exists(output_dir+"/snarl_genotypes.tsv"));
       

        std::vector<std::string> lengths;
        lengths.emplace_back("3/4");
        lengths.emplace_back("6/8");
        std::vector<std::string> paths;
        paths.emplace_back(">1>2>0>5>6");
        paths.emplace_back(">1>2>0>5>2>0>5>6");
        std::vector<std::string> sequences;
        std::unordered_map<std::string, size_t> allele_assignment;
        allele_assignment["path0#"] = 0;
        allele_assignment["path1#"] = 1;
        allele_assignment["path2#"] = 1;
        snarl_genotype_values_t truth1 ({">1", "<6",
                                       (std::string)"path0", (size_t)10, (size_t)14, 
                                       (size_t) 1,
                                       lengths,
                                       paths, 
                                       sequences,
                                       allele_assignment});

        lengths.clear();
        lengths.emplace_back("1");
        lengths.emplace_back("0");
        lengths.emplace_back("1");
        paths.clear();
        paths.emplace_back(">2>3>4");
        paths.emplace_back(">2>4<0>2>4");
        paths.emplace_back(">2>3>4<0>2>4");
        allele_assignment.clear();
        allele_assignment["path0#"] = 0;
        allele_assignment["path1#"] = 1;
        allele_assignment["path2#"] = 2;
        snarl_genotype_values_t truth2 ({">2", "<4",
                                       (std::string)"path0", (size_t)11, (size_t)12, 
                                       (size_t) 2,
                                       lengths,
                                       paths, 
                                       sequences,
                                       allele_assignment});

        std::vector<snarl_genotype_values_t> loaded_snarls= load_genotype_file(output_dir+"/snarl_genotypes.tsv");
        std::vector<bool> found_snarls (2, false);
        for (const snarl_genotype_values_t& test : loaded_snarls) {
            if (test == truth1) {
                found_snarls[0] = true;
            } else if (test == truth2) {
                found_snarls[1] = true;
            } else {
                std::cerr << "Unknown snarl " << std::endl << test << std::endl;
            }
        }

        if (!found_snarls[0]) {
            std::cerr << "Missing snarl " << std::endl << truth1 << std::endl;
            REQUIRE(false);
        }
        if (!found_snarls[1]) {
            std::cerr << "Missing snarl " << std::endl << truth2 << std::endl;
            REQUIRE(false);
        }


    }
    clean_output_dir(output_dir);
}

TEST_CASE("Multiple connected components test", "[test]") {
   /*
    This graph is duplicated twice 
                       5
                     /   \
            1       4 ----6    8
          /   \   /         \ / \
        0       3  ----------7---9
          \   /
            2

   */


    const std::string output_dir = "../output_binary/";
    const std::string data_path = "../tests/test_data/test_graphs/";
    const std::string graph_base = "simple_nested_chains";

    SECTION("Test snarls output") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph -u";
        cmd += " -g " + data_path + graph_base + ".hg"
            + " -d " + data_path + graph_base + ".dist"
            + " -L"
            + " -r ref"
            + " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE(false);
        }

        REQUIRE(std::filesystem::exists(output_dir + "/snarl_genotypes.tsv"));
        std::ifstream snarlsfile;
        snarlsfile.open(output_dir + "/snarl_genotypes.tsv");
        REQUIRE(snarlsfile.peek() != std::ifstream::traits_type::eof());

        size_t line_count = 0;
        std::string line;
        while (std::getline(snarlsfile, line)) {
            // Only start counting snarls after proper header
            if (line[0] != '#') {
                line_count++;
            }
        }
        snarlsfile.close();
        REQUIRE(line_count==24);

        // TODO: Add something that actually checks this
        //bool passed = compare_output_dirs(output_dir, expected_dir);
        //REQUIRE(passed);

    }
    SECTION("Test snarls output multithreaded") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph -u";
        cmd += " -g " + data_path + graph_base + ".hg"
            + " -d " + data_path + graph_base + ".dist"
            + " -L"
            + " -r ref -t 4"
            + " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE(false);
        }


        REQUIRE(std::filesystem::exists(output_dir + "/snarl_genotypes.tsv"));
        std::ifstream snarlsfile;
        snarlsfile.open(output_dir + "/snarl_genotypes.tsv");
        REQUIRE(snarlsfile.peek() != std::ifstream::traits_type::eof());

        size_t line_count = 0;
        std::string line;
        while (std::getline(snarlsfile, line)) {
            // Only count non-header lines
            if (line[0] != '#') {
                line_count++;
            }
        }
        snarlsfile.close();
        REQUIRE(line_count==24);

        // TODO: Add something that actually checks this
        //bool passed = compare_output_dirs(output_dir, expected_dir);
        //REQUIRE(passed);

    }


    clean_output_dir(output_dir);
}

