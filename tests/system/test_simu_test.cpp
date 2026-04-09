#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "compare_files_utils.hpp"
#include "load_tables.hpp"

namespace fs = std::filesystem;
using namespace std;
// TODO: Add vcf tests. This just copies the stats part of the old graph test

TEST_CASE("Giant unverified binary association tests graph", "[graph]") {

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

TEST_CASE("Output simple nested chain", "[graph]") {
    const std::string output_dir = "../output_binary";
    const std::string graph_base = "../tests/test_data/test_graphs/simple_nested_chain";
    const std::string samples_file = "./samples.tsv";
    
    std::vector<std::string> samples_of_interest = {"path1", "path3"};
    std::vector<std::string> other_samples = {"path0", "path2"};
    
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

    SECTION("Test chi2 tsv output") {

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
        
        std::string cmd_test = "../bin/stoat test -u";
        cmd_test += " -g " + output_dir + "/snarl_genotypes.tsv"
            + " -p " + samples_file
            + " -m chi2"
            + " --output " + output_dir;

        std::cout << "Command run : \n" << cmd_test << std::endl;

        command_output = std::system(cmd_test.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd_test << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/stoat.assoc.pvalues.tsv"));

        binary_table_values_t truth1 ({(std::string)"path0#0#path0", 
                                  (size_t)1, 
                                  (size_t)2, 
                                  std::make_pair<std::string, std::string>(">1","<4"), 
                                  std::vector<std::string>({"1","1"}),
                                  (std::string) "1",
                                  (std::string) "1", 
                                  std::vector<std::string>({"1:1","1:1"}), 
                                  (size_t)1});

        binary_table_values_t truth2 ({(std::string)"path0#0#path0", 
                                  (size_t)3, 
                                  (size_t)6, 
                                  std::make_pair<std::string, std::string>(">4","<8"), 
                                  std::vector<std::string>({"0","2/3"}),
                                  (std::string) "1",
                                  (std::string) "0.2482", 
                                  std::vector<std::string>({"1:0","1:2"}), 
                                  (size_t)1});

        binary_table_values_t truth3 ({(std::string)"path0#0#path0", 
                                  (size_t)4, 
                                  (size_t)5, 
                                  std::make_pair<std::string, std::string>(">5","<7"), 
                                  std::vector<std::string>({"0","1"}),
                                  (std::string) "0.3333",
                                  (std::string) "8.3265e-02", 
                                  std::vector<std::string>({"0:2","1:0"}), 
                                  (size_t)2});
        ifstream in_table;
        in_table.open(output_dir+"/stoat.assoc.pvalues.tsv");
        std::string line; 
        std::getline(in_table, line); // Header
        std::vector<bool> found_snarls (3, false);
        while (std::getline(in_table, line)) {
            binary_table_values_t test = load_binary_snarl_line(line);
            if (test == truth1) {
                found_snarls[0] = true;
            } else if (test == truth2) {
                found_snarls[1] = true;
            } else if (test == truth3) {
                found_snarls[2] = true;
            } else if (test.chr != "NA") { // could include the extra snarl but maybe not
                cerr << "Line doesn't match the truth" << endl;
                cerr << line << endl;
                REQUIRE(false);
            }
        }
        in_table.close();
        REQUIRE(found_snarls[0]);
        REQUIRE(found_snarls[1]);
        REQUIRE(found_snarls[2]);

    }

    SECTION("Test chi2 tsv output without lengths") {

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

        std::string cmd_test = "../bin/stoat test -u";
        cmd_test += " -g " + output_dir + "/snarl_genotypes.tsv"
            + " -p " + samples_file
            + " -m chi2"
            + " --output " + output_dir;

        std::cout << "Command run : \n" << cmd_test << std::endl;

        command_output = std::system(cmd_test.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd_test << "\n";
            REQUIRE( false);
        }
        
        REQUIRE(std::filesystem::exists(output_dir+"/stoat.assoc.pvalues.tsv"));

        binary_table_values_t truth1 ({(std::string)"path0#0#path0", 
                                  (size_t)1, 
                                  (size_t)2, 
                                  std::make_pair<std::string, std::string>(">1","<4"), 
                                  std::vector<std::string>({}),
                                  (std::string) "1",
                                  (std::string) "1", 
                                  std::vector<std::string>({"1:1","1:1"}), 
                                  (size_t)1});

        binary_table_values_t truth2 ({(std::string)"path0#0#path0", 
                                  (size_t)3, 
                                  (size_t)6, 
                                  std::make_pair<std::string, std::string>(">4","<8"), 
                                  std::vector<std::string>({}),
                                  (std::string) "1",
                                  (std::string) "0.2482", 
                                  std::vector<std::string>({"1:0","1:2"}), 
                                  (size_t)1});

        binary_table_values_t truth3 ({(std::string)"path0#0#path0", 
                                  (size_t)4, 
                                  (size_t)5, 
                                  std::make_pair<std::string, std::string>(">5","<7"), 
                                  std::vector<std::string>({}),
                                  (std::string) "0.3333",
                                  (std::string) "8.3265e-02", 
                                  std::vector<std::string>({"0:2","1:0"}), 
                                  (size_t)2});
        ifstream in_table;
        in_table.open(output_dir+"/stoat.assoc.pvalues.tsv");
        std::string line; 
        std::getline(in_table, line); // Header
        std::vector<bool> found_snarls (3, false);
        while (std::getline(in_table, line)) {
            binary_table_values_t test = load_binary_snarl_line(line);
            if (test == truth1) {
                found_snarls[0] = true;
            } else if (test == truth2) {
                found_snarls[1] = true;
            } else if (test == truth3) {
                found_snarls[2] = true;
            } else if (test.chr != "NA") { // could include the extra snarl but maybe not
                cerr << "Line doesn't match the truth" << endl;
                cerr << line << endl;
                REQUIRE(false);
            }
        }
        in_table.close();
        REQUIRE(found_snarls[0]);
        REQUIRE(found_snarls[1]);
        REQUIRE(found_snarls[2]);

    }

    SECTION("Test exact tsv output") {

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

        std::string cmd_test = "../bin/stoat test -u";
        cmd_test += " -g " + output_dir + "/snarl_genotypes.tsv"
            + " -p " + samples_file
            + " -m exact"
            + " --output " + output_dir;

        std::cout << "Command run : \n" << cmd_test << std::endl;

        command_output = std::system(cmd_test.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd_test << "\n";
            REQUIRE( false);
        }

        
        REQUIRE(std::filesystem::exists(output_dir+"/stoat.assoc.pvalues.tsv"));

        binary_table_values_t truth ({(std::string)"path0#0#path0", 
                                  (size_t)4, 
                                  (size_t)5, 
                                  std::make_pair<std::string, std::string>(">5","<7"), 
                                  std::vector<std::string>({"0","1"}),
                                  (std::string) "NA",
                                  (std::string) "NA", 
                                  std::vector<std::string>({"NA"}), 
                                  (size_t)2});
        ifstream in_table;
        in_table.open(output_dir+"/stoat.assoc.pvalues.tsv");
        std::string line; 
        std::getline(in_table, line); // Header
        std::getline(in_table, line);
        binary_table_values_t test = load_binary_snarl_line(line);
        REQUIRE(test == truth);
        REQUIRE((!std::getline(in_table, line)));
        in_table.close();

    }

    clean_output_dir(output_dir);
    fs::remove(samples_file);
}

TEST_CASE("Output simple nested chain gbz", "[graph]") {
    const std::string output_dir = "../output_binary";
    const std::string graph_base = "../tests/test_data/test_graphs/simple_nested_chain";
    const std::string samples_file = "./samples.tsv";
    

    std::vector<std::string> samples_of_interest = {"path1", "path3"};
    std::vector<std::string> other_samples = {"path0", "path2"};
    
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
        write_cmd = "cat " + samples_file;
        ignore = std::system(write_cmd.c_str());



    SECTION("Test chi2 tsv output") {

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

        std::string cmd_test = "../bin/stoat test -u";
        cmd_test += " -g " + output_dir + "/snarl_genotypes.tsv"
            + " -p " + samples_file
            + " -m chi2"
            + " --output " + output_dir;

        std::cout << "Command run : \n" << cmd_test << std::endl;

        command_output = std::system(cmd_test.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd_test << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/stoat.assoc.pvalues.tsv"));

        binary_table_values_t truth1 ({(std::string)"path0#0#path0", 
                                  (size_t)1, 
                                  (size_t)2, 
                                  std::make_pair<std::string, std::string>(">1","<4"), 
                                  std::vector<std::string>({"1","1"}),
                                  (std::string) "1",
                                  (std::string) "1", 
                                  std::vector<std::string>({"1:1","1:1"}), 
                                  (size_t)1});

        binary_table_values_t truth2 ({(std::string)"path0#0#path0", 
                                  (size_t)3, 
                                  (size_t)6, 
                                  std::make_pair<std::string, std::string>(">4","<8"), 
                                  std::vector<std::string>({"0","2/3"}),
                                  (std::string) "1",
                                  (std::string) "0.2482", 
                                  std::vector<std::string>({"1:0","1:2"}), 
                                  (size_t)1});

        binary_table_values_t truth3 ({(std::string)"path0#0#path0", 
                                  (size_t)4, 
                                  (size_t)5, 
                                  std::make_pair<std::string, std::string>(">5","<7"), 
                                  std::vector<std::string>({"0","1"}),
                                  (std::string) "0.3333",
                                  (std::string) "8.3265e-02", 
                                  std::vector<std::string>({"0:2","1:0"}), 
                                  (size_t)2});
        ifstream in_table;
        in_table.open(output_dir+"/stoat.assoc.pvalues.tsv");
        std::string line; 
        std::getline(in_table, line); // Header
        std::vector<bool> found_snarls (3, false);
        while (std::getline(in_table, line)) {
            binary_table_values_t test = load_binary_snarl_line(line);
            if (test == truth1) {
                found_snarls[0] = true;
            } else if (test == truth2) {
                found_snarls[1] = true;
            } else if (test == truth3) {
                found_snarls[2] = true;
            } else if (test.chr != "NA") { // could include the extra snarl but maybe not
                cerr << "Line doesn't match the truth" << endl;
                cerr << line << endl;
                REQUIRE(false);
            }
        }
        in_table.close();
        REQUIRE(found_snarls[0]);
        REQUIRE(found_snarls[1]);
        REQUIRE(found_snarls[2]);


    }

    SECTION("Test exact tsv output") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph -u";
           cmd += " -g " + graph_base + ".gbz"
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

        std::string cmd_test = "../bin/stoat test -u";
        cmd_test += " -g " + output_dir + "/snarl_genotypes.tsv"
            + " -p " + samples_file
            + " -m exact"
            + " --output " + output_dir;

        std::cout << "Command run : \n" << cmd_test << std::endl;

        command_output = std::system(cmd_test.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd_test << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/stoat.assoc.pvalues.tsv"));

        binary_table_values_t truth ({(std::string)"path0#0#path0", 
                                  (size_t)4, 
                                  (size_t)5, 
                                  std::make_pair<std::string, std::string>(">5","<7"), 
                                  std::vector<std::string>({"0","1"}),
                                  (std::string) "NA",
                                  (std::string) "NA", 
                                  std::vector<std::string>({"NA"}), 
                                  (size_t)2});
        ifstream in_table;
        in_table.open(output_dir+"/stoat.assoc.pvalues.tsv");
        std::string line; 
        std::getline(in_table, line); // Header
        std::getline(in_table, line);
        binary_table_values_t test = load_binary_snarl_line(line);
        REQUIRE(test == truth);
        REQUIRE((!std::getline(in_table, line)));
        in_table.close();
    }

    SECTION("Test chi2 tsv output with maf excluding one snarl") {

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

        std::string cmd_test = "../bin/stoat test -u";
        cmd_test += " -g " + output_dir + "/snarl_genotypes.tsv"
            + " -p " + samples_file
            + " -m chi2 -M 0.26"
            + " --output " + output_dir;

        std::cout << "Command run : \n" << cmd_test << std::endl;

        command_output = std::system(cmd_test.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd_test << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/stoat.assoc.pvalues.tsv"));

        binary_table_values_t truth1 ({(std::string)"path0#0#path0", 
                                  (size_t)1, 
                                  (size_t)2, 
                                  std::make_pair<std::string, std::string>(">1","<4"), 
                                  std::vector<std::string>({"1","1"}),
                                  (std::string) "1",
                                  (std::string) "1", 
                                  std::vector<std::string>({"1:1","1:1"}), 
                                  (size_t)1});

        binary_table_values_t truth3 ({(std::string)"path0#0#path0", 
                                  (size_t)4, 
                                  (size_t)5, 
                                  std::make_pair<std::string, std::string>(">5","<7"), 
                                  std::vector<std::string>({"0","1"}),
                                  (std::string) "0.3333",
                                  (std::string) "8.3265e-02", 
                                  std::vector<std::string>({"0:2","1:0"}), 
                                  (size_t)2});

        ifstream in_table;
        in_table.open(output_dir+"/stoat.assoc.pvalues.tsv");
        std::string line; 
        std::getline(in_table, line); // Header
        std::vector<bool> found_snarls (3, false);
        while (std::getline(in_table, line)) {
            binary_table_values_t test = load_binary_snarl_line(line);
            if (test == truth1) {
                found_snarls[0] = true;
            } else if (test == truth3) {
                found_snarls[2] = true;
            } else if (test.chr != "NA") { // could include the extra snarl but maybe not
                cerr << "Line doesn't match the truth" << endl;
                cerr << line << endl;
                REQUIRE(false);
            }
        }
        in_table.close();
        REQUIRE(found_snarls[0]);
        REQUIRE(found_snarls[2]);
    }

    SECTION("Test chi2 tsv output with individual count excluding one snarl") {

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

        std::string cmd_test = "../bin/stoat test -u";
        cmd_test += " -g " + output_dir + "/snarl_genotypes.tsv"
            + " -p " + samples_file
            + " -m chi2 -I 4"
            + " --output " + output_dir;

        std::cout << "Command run : \n" << cmd_test << std::endl;

        command_output = std::system(cmd_test.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd_test << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/stoat.assoc.pvalues.tsv"));


        binary_table_values_t truth1 ({(std::string)"path0#0#path0", 
                                  (size_t)1, 
                                  (size_t)2, 
                                  std::make_pair<std::string, std::string>(">1","<4"), 
                                  std::vector<std::string>({"1","1"}),
                                  (std::string) "1",
                                  (std::string) "1", 
                                  std::vector<std::string>({"1:1","1:1"}), 
                                  (size_t)1});

        binary_table_values_t truth2 ({(std::string)"path0#0#path0", 
                                  (size_t)3, 
                                  (size_t)6, 
                                  std::make_pair<std::string, std::string>(">4","<8"), 
                                  std::vector<std::string>({"0","2/3"}),
                                  (std::string) "1",
                                  (std::string) "0.2482", 
                                  std::vector<std::string>({"1:0","1:2"}), 
                                  (size_t)1});
        ifstream in_table;
        in_table.open(output_dir+"/stoat.assoc.pvalues.tsv");
        std::string line; 
        std::getline(in_table, line); // Header
        std::vector<bool> found_snarls (3, false);
        while (std::getline(in_table, line)) {
            binary_table_values_t test = load_binary_snarl_line(line);
            if (test == truth1) {
                found_snarls[0] = true;
            } else if (test == truth2) {
                found_snarls[1] = true;
            } else if (test.chr != "NA") { // could include the extra snarl but maybe not
                cerr << "Line doesn't match the truth" << endl;
                cerr << line << endl;
                REQUIRE(false);
            }
        }
        in_table.close();
        REQUIRE(found_snarls[0]);
        REQUIRE(found_snarls[1]);


    }

    clean_output_dir(output_dir);
    fs::remove(samples_file);
}

TEST_CASE("Output loop with snarl", "[graph]") {
    const std::string output_dir = "../output_binary";
    const std::string graph_base = "../tests/test_data/test_graphs/loop_with_indel";

    std::vector<std::string> samples_of_interest = {"path1", "path2"};
    std::vector<std::string> other_samples = {"path0"};

    const std::string samples_file = "./samples.tsv";
    
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

    SECTION("Test chi2 tsv output") {

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
        
        std::string cmd_test = "../bin/stoat test -u";
        cmd_test += " -g " + output_dir + "/snarl_genotypes.tsv"
            + " -p " + samples_file
            + " -m chi2"
            + " --output " + output_dir;

        std::cout << "Command run : \n" << cmd_test << std::endl;

        command_output = std::system(cmd_test.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd_test << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/stoat.assoc.pvalues.tsv"));

        binary_table_values_t truth1 ({(std::string)"path0", 
                                  (size_t)10, 
                                  (size_t)14, 
                                  std::make_pair<std::string, std::string>(">1","<6"), 
                                  std::vector<std::string>({"3/4","6/8"}),
                                  (std::string) "0.3333",
                                  (std::string) "8.3265e-02", 
                                  std::vector<std::string>({"1:0","0:2"}), 
                                  (size_t)1});

        binary_table_values_t truth2 ({(std::string)"path0", 
                                  (size_t)11, 
                                  (size_t)12, 
                                  std::make_pair<std::string, std::string>(">2","<4"), 
                                  std::vector<std::string>({"1","0","1"}),
                                  (std::string) "NA",
                                  (std::string) "0.2231", 
                                  std::vector<std::string>({"1:0","0:1", "0:1"}), 
                                  (size_t)2});

        ifstream in_table;
        in_table.open(output_dir+"/stoat.assoc.pvalues.tsv");
        std::string line; 
        std::getline(in_table, line); // Header
        std::vector<bool> found_snarls (2, false);
        while (std::getline(in_table, line)) {
            binary_table_values_t test = load_binary_snarl_line(line);
            cerr << "LINE: " << line << endl;
            if (test == truth1) {
                found_snarls[0] = true;
            } else if (test == truth2) {
                found_snarls[1] = true;
            }
        }
        in_table.close();
        REQUIRE(found_snarls[0]);
        REQUIRE(found_snarls[1]);

    }

    SECTION("Test exact tsv output") {

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
        
        std::string cmd_test = "../bin/stoat test -u";

        cmd_test += " -g " + output_dir + "/snarl_genotypes.tsv"
            + " -p " + samples_file
            + " -m exact"
            + " --output " + output_dir;

        std::cout << "Command run : \n" << cmd_test << std::endl;
        
        command_output = std::system(cmd_test.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd_test << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/stoat.assoc.pvalues.tsv"));


        binary_table_values_t truth1 ({(std::string)"path0", 
                                  (size_t)10, 
                                  (size_t)14, 
                                  std::make_pair<std::string, std::string>(">1","<6"), 
                                  std::vector<std::string>({"3/4","6/8"}),
                                  (std::string) "NA",
                                  (std::string) "NA", 
                                  std::vector<std::string>({"NA"}), 
                                  (size_t)1});

        binary_table_values_t truth2 ({(std::string)"path0", 
                                  (size_t)11, 
                                  (size_t)12, 
                                  std::make_pair<std::string, std::string>(">2","<4"), 
                                  std::vector<std::string>({"1","0","1"}),
                                  (std::string) "NA",
                                  (std::string) "NA", 
                                  std::vector<std::string>({"NA"}), 
                                  (size_t)2});

        ifstream in_table;
        in_table.open(output_dir+"/stoat.assoc.pvalues.tsv");
        std::string line; 
        std::getline(in_table, line); // Header
        std::vector<bool> found_snarls (2, false);
        while (std::getline(in_table, line)) {
            binary_table_values_t test = load_binary_snarl_line(line);
            cerr << "LINE: " << line << endl;
            if (test == truth1) {
                found_snarls[0] = true;
            } else if (test == truth2) {
                found_snarls[1] = true;
            }
        }
        in_table.close();
        REQUIRE(found_snarls[0]);
        REQUIRE(found_snarls[1]);
    }

    fs::remove(samples_file);
    clean_output_dir(output_dir);
}

TEST_CASE("Multiple connected components", "[graph]") {
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

