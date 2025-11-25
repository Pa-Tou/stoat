#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "compare_files_utils.hpp"
#include "load_tables.hpp"

namespace fs = std::filesystem;
using namespace std;

TEST_CASE("Giant unverified binary association tests graph", "[graph]") {
    // Just check that this runs and produces some output

    const std::string output_dir = "../output_binary";
    const std::string data_path = "../data/binary";
    const std::string graph_base = "pg.full";


    SECTION("Test tsv output") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph";

        cmd +=" -g " + data_path + "/" + graph_base + ".pg"
            + " -d " + data_path + "/" + graph_base + ".dist"
            + " -b " + data_path + "/phenotype_samples.tsv"
            + " -T chi2 -r ref";

        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE(false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/binary_table_graph.tsv"));
        std::ifstream testfile;
        testfile.open(output_dir+"/binary_table_graph.tsv");
        REQUIRE(testfile.peek() != std::ifstream::traits_type::eof());
        testfile.close();

        // TODO: Add something that actually checks this
        //bool passed = compare_output_dirs(output_dir, expected_dir);
        //REQUIRE(passed);

    }
    SECTION("Test tsv output saving snarls") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph";

        cmd +=" -g " + data_path + "/" + graph_base + ".pg"
            + " -d " + data_path + "/" + graph_base + ".dist"
            + " -b " + data_path + "/phenotype_samples.tsv"
            + " -M 0 -l 0 -s " + output_dir + "/" + graph_base + ".snarls.txt" 
            + " -T chi2 -r ref";

        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE(false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/binary_table_graph.tsv"));
        std::ifstream testfile;
        testfile.open(output_dir+"/binary_table_graph.tsv");
        REQUIRE(testfile.peek() != std::ifstream::traits_type::eof());
        testfile.close();

        REQUIRE(std::filesystem::exists(output_dir+"/"+graph_base+".snarls.txt"));
        std::ifstream snarlsfile;
        snarlsfile.open(output_dir+"/"+graph_base+".snarls.txt");
        REQUIRE(snarlsfile.peek() != std::ifstream::traits_type::eof());

        size_t line_count = 0;
        bool start_counting= false;
        std::string line;
        while (std::getline(snarlsfile, line)) {
            // Only start counting snarls after proper header
            if (start_counting) {
                line_count++;
            }
            if (line == "#SNARLS") {
                start_counting = true;
            }
        }
        snarlsfile.close();
        REQUIRE(line_count==1524);

        // TODO: Add something that actually checks this
        //bool passed = compare_output_dirs(output_dir, expected_dir);
        //REQUIRE(passed);

    }


    SECTION("Test chi2 fasta output") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph";

        cmd +=  " -g " + data_path + "/" + graph_base + ".pg"
            + " -d " + data_path + "/" + graph_base + ".dist"
            + " -b " + data_path + "/phenotype_samples.tsv"
            + " -T chi2 -r ref -O fasta";

        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/binary_output.fasta"));
        std::ifstream testfile;
        testfile.open(output_dir+"/binary_output.fasta");
        REQUIRE(testfile.peek() != std::ifstream::traits_type::eof());
        testfile.close();


        REQUIRE(is_valid_fasta(output_dir+"/binary_output.fasta"));


    }

    SECTION("Test exact fasta output") {


        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph";

        cmd +=  " -g " + data_path + "/" + graph_base + ".pg"
            + " -d " + data_path + "/" + graph_base + ".dist"
            + " -b " + data_path + "/phenotype_samples.tsv"
            + " -T exact -r ref -O fasta";


        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/binary_output.fasta"));
        std::ifstream testfile;
        testfile.open(output_dir+"/binary_output.fasta");
        //REQUIRE(testfile.peek() != std::ifstream::traits_type::eof());
        testfile.close();

        REQUIRE(is_valid_fasta(output_dir+"/binary_output.fasta"));

    }

    SECTION("Test serialization of snarls") {

        clean_output_dir(output_dir);
        string output_dir_loaded = output_dir + "_loaded";
        clean_output_dir(output_dir_loaded);
        std::filesystem::create_directory(output_dir_loaded);

        std::string cmd = "../bin/stoat graph";

        cmd +=" -g " + data_path + "/" + graph_base + ".pg"
            + " -d " + data_path + "/" + graph_base + ".dist"
            + " -b " + data_path + "/phenotype_samples.tsv"
            + " -s " + output_dir_loaded + "/" + graph_base + ".saved_snarls"
            + " -T chi2 -r ref";

        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE(false);
        }

        REQUIRE(std::filesystem::exists(output_dir_loaded+"/" + graph_base + ".saved_snarls"));
        REQUIRE(std::filesystem::exists(output_dir+"/binary_table_graph.tsv"));
        std::ifstream testfile;
        testfile.open(output_dir+"/binary_table_graph.tsv");
        REQUIRE(testfile.peek() != std::ifstream::traits_type::eof());
        testfile.close();

        // Now run it again but using the saved snarls

        cmd = "../bin/stoat graph";

        cmd +=" -g " + data_path + "/" + graph_base + ".pg"
            + " -b " + data_path + "/phenotype_samples.tsv"
            + " -s " + output_dir_loaded + "/" + graph_base + ".saved_snarls"
            + " -T chi2 -r ref";

        cmd += " --output " + output_dir_loaded;

        std::cout << "Command run : \n" << cmd << std::endl;

        command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE(false);
        }

        REQUIRE(std::filesystem::exists(output_dir_loaded+"/binary_table_graph.tsv"));
        testfile;
        testfile.open(output_dir_loaded+"/binary_table_graph.tsv");
        REQUIRE(testfile.peek() != std::ifstream::traits_type::eof());
        testfile.close();

        REQUIRE(files_equal(output_dir+"/binary_table_graph.tsv", output_dir_loaded+"/binary_table_graph.tsv"));

        // TODO: Add something that actually checks this
        //bool passed = compare_output_dirs(output_dir, expected_dir);
        //REQUIRE(passed);

    }

    //clean_output_dir(output_dir);
    //clean_output_dir(output_dir_loaded);
}

TEST_CASE("Output simple nested chain", "[graph]") {
    const std::string output_dir = "../output_binary";
    const std::string graph_base = "../tests/graph_test/simple_nested_chain";
    const std::string samples_file = "./samples.tsv";
    

    std::vector<std::string> samples_of_interest = {"path1", "path3"};
    std::vector<std::string> other_samples = {"path0", "path2"};
    
    string write_cmd = "echo \"FID\tIID\tPHENO\" > " + samples_file;
    int ignore = std::system(write_cmd.c_str());
    for (auto sample : samples_of_interest) {
        write_cmd = "echo \"" + sample + "\t" + sample + "\t2\" >> " + samples_file;
        ignore = std::system(write_cmd.c_str());
    }
    for (auto sample : other_samples) {
        write_cmd = "echo \"" + sample + "\t" + sample + "\t1\" >> " + samples_file;
        ignore = std::system(write_cmd.c_str());
    }


    SECTION("Test chi2 tsv output") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph";

        cmd +=  " -g " + graph_base + ".hg"
            + " -d " + graph_base + ".dist"
            + " -b " + samples_file
            + " -T chi2 -r path0 -V 4";


        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/binary_table_graph.tsv"));

        binary_table_values_t truth1 ({(std::string)"path0#0#path0", 
                                  (size_t)1, 
                                  (size_t)2, 
                                  std::make_pair<size_t, size_t>(1,4), 
                                  std::vector<std::string>({"1","1"}),
                                  (std::string) "1",
                                  (std::string) "1", 
                                  std::vector<std::string>({"1:1","1:1"}), 
                                  (size_t)1});

        binary_table_values_t truth2 ({(std::string)"path0#0#path0", 
                                  (size_t)3, 
                                  (size_t)6, 
                                  std::make_pair<size_t, size_t>(4,8), 
                                  std::vector<std::string>({"0","2/3"}),
                                  (std::string) "1",
                                  (std::string) "0.2482", 
                                  std::vector<std::string>({"0:1","2:1"}), 
                                  (size_t)1});

        binary_table_values_t truth3 ({(std::string)"path0#0#path0", 
                                  (size_t)4, 
                                  (size_t)5, 
                                  std::make_pair<size_t, size_t>(5,7), 
                                  std::vector<std::string>({"0","1"}),
                                  (std::string) "0.3333",
                                  (std::string) "8.3265e-02", 
                                  std::vector<std::string>({"2:0","0:1"}), 
                                  (size_t)2});
        ifstream in_table;
        in_table.open(output_dir+"/binary_table_graph.tsv");
        std::string line; 
        std::getline(in_table, line); // Header
        std::vector<bool> found_snarls (3, false);
        while (std::getline(in_table, line)) {
            binary_table_values_t test = load_binary_snarl_line(line);
            bool match = false;
            if (test == truth1) {
                match = true;
                found_snarls[0] = true;
            } else if (test == truth2) {
                match = true;
                found_snarls[1] = true;
            } else if (test == truth3) {
                match = true;
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

        std::string cmd = "../bin/stoat graph";

        cmd +=  " -g " + graph_base + ".hg"
            + " -d " + graph_base + ".dist"
            + " -b " + samples_file
            + " -T exact -r path0";


        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/binary_table_graph.tsv"));

        binary_table_values_t truth ({(std::string)"path0#0#path0", 
                                  (size_t)4, 
                                  (size_t)5, 
                                  std::make_pair<size_t, size_t>(5,7), 
                                  std::vector<std::string>({"0","1"}),
                                  (std::string) "NA",
                                  (std::string) "NA", 
                                  std::vector<std::string>({"NA"}), 
                                  (size_t)2});
        ifstream in_table;
        in_table.open(output_dir+"/binary_table_graph.tsv");
        std::string line; 
        std::getline(in_table, line); // Header
        std::getline(in_table, line);
        binary_table_values_t test = load_binary_snarl_line(line);
        REQUIRE(test == truth);
        REQUIRE((!std::getline(in_table, line)));
        in_table.close();

    }

    SECTION("Test chi2 fasta output") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph";

        cmd +=  " -g " + graph_base + ".hg"
            + " -d " + graph_base + ".dist"
            + " -b " + samples_file
            + " -T chi2 -r path0 -O fasta";


        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/binary_output.fasta"));

        REQUIRE(is_valid_fasta(output_dir+"/binary_output.fasta"));

        std::vector<std::tuple<size_t, std::string, std::string>> truth_fasta;

        truth_fasta.emplace_back(1, ">snarl:1_4|path0#0#path0:1-2|path0#0#path0:1-2", "C");
        truth_fasta.emplace_back(1, ">snarl:1_4|path0#0#path0:1-2|path1#0#path1#0:1-2", "C");
        truth_fasta.emplace_back(2, ">snarl:1_4|path0#0#path0:1-2|path2:1-2", "C");
        truth_fasta.emplace_back(2, ">snarl:1_4|path0#0#path0:1-2|path3:1-2", "C");

        truth_fasta.emplace_back(3, ">snarl:4_8|path0#0#path0:3-6|path0#0#path0:3-6", "TNA");
        truth_fasta.emplace_back(3, ">snarl:4_8|path0#0#path0:3-6|path1#0#path1#0:3-6", "TNA");
        truth_fasta.emplace_back(3, ">snarl:4_8|path0#0#path0:3-6|path3:3-5", "TNA");
        truth_fasta.emplace_back(4, ">snarl:4_8|path0#0#path0:3-6|path2:3-3", "");

        truth_fasta.emplace_back(5, ">snarl:5_7|path0#0#path0:4-5|path0#0#path0:4-5", "C");
        truth_fasta.emplace_back(6, ">snarl:5_7|path0#0#path0:4-5|path1#0#path1#0:4-4", "");
        truth_fasta.emplace_back(6, ">snarl:5_7|path0#0#path0:4-5|path3:4-4", "");

        REQUIRE(fasta_equal(output_dir+"/binary_output.fasta", truth_fasta));

    }

    SECTION("Test exact fasta output") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph";

        cmd +=  " -g " + graph_base + ".hg"
            + " -d " + graph_base + ".dist"
            + " -b " + samples_file
            + " -T exact -r path0 -O fasta";


        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/binary_output.fasta"));

        REQUIRE(is_valid_fasta(output_dir+"/binary_output.fasta"));

        std::vector<std::tuple<size_t, std::string, std::string>> truth_fasta;

        truth_fasta.emplace_back(2, ">snarl:5_7|path0#0#path0:4-5|path0#0#path0:4-5", "C");

        truth_fasta.emplace_back(1, ">snarl:5_7|path0#0#path0:4-5|path1#0#path1#0:4-4", "");
        truth_fasta.emplace_back(1, ">snarl:5_7|path0#0#path0:4-5|path3:4-4", "");
        REQUIRE(fasta_equal(output_dir+"/binary_output.fasta", truth_fasta));


    }

    clean_output_dir(output_dir);
    fs::remove(samples_file);
}

TEST_CASE("Output simple nested chain gbz", "[graph]") {
    const std::string output_dir = "../output_binary";
    const std::string graph_base = "../tests/graph_test/simple_nested_chain";
    const std::string samples_file = "./samples.tsv";
    

    std::vector<std::string> samples_of_interest = {"path1", "path3"};
    std::vector<std::string> other_samples = {"path0", "path2"};
    
    string write_cmd = "echo \"FID\tIID\tPHENO\" > " + samples_file;
    int ignore = std::system(write_cmd.c_str());
    for (auto sample : samples_of_interest) {
        write_cmd = "echo \"" + sample + "\t" + sample + "\t2\" >> " + samples_file;
        ignore = std::system(write_cmd.c_str());
    }
    for (auto sample : other_samples) {
        write_cmd = "echo \"" + sample + "\t" + sample + "\t1\" >> " + samples_file;
        ignore = std::system(write_cmd.c_str());
    }
        write_cmd = "cat " + samples_file;
        ignore = std::system(write_cmd.c_str());



    SECTION("Test chi2 tsv output") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph";

        cmd +=  " -g " + graph_base + ".gbz"
            + " -d " + graph_base + ".dist"
            + " -b " + samples_file
            + " -T chi2 -r path0 -V 4";


        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/binary_table_graph.tsv"));

        binary_table_values_t truth1 ({(std::string)"path0#0#path0", 
                                  (size_t)1, 
                                  (size_t)2, 
                                  std::make_pair<size_t, size_t>(1,4), 
                                  std::vector<std::string>({"1","1"}),
                                  (std::string) "1",
                                  (std::string) "1", 
                                  std::vector<std::string>({"1:1","1:1"}), 
                                  (size_t)1});

        binary_table_values_t truth2 ({(std::string)"path0#0#path0", 
                                  (size_t)3, 
                                  (size_t)6, 
                                  std::make_pair<size_t, size_t>(4,8), 
                                  std::vector<std::string>({"0","2/3"}),
                                  (std::string) "1",
                                  (std::string) "0.2482", 
                                  std::vector<std::string>({"0:1","2:1"}), 
                                  (size_t)1});

        binary_table_values_t truth3 ({(std::string)"path0#0#path0", 
                                  (size_t)4, 
                                  (size_t)5, 
                                  std::make_pair<size_t, size_t>(5,7), 
                                  std::vector<std::string>({"0","1"}),
                                  (std::string) "0.3333",
                                  (std::string) "8.3265e-02", 
                                  std::vector<std::string>({"2:0","0:1"}), 
                                  (size_t)2});
        ifstream in_table;
        in_table.open(output_dir+"/binary_table_graph.tsv");
        std::string line; 
        std::getline(in_table, line); // Header
        std::vector<bool> found_snarls (3, false);
        while (std::getline(in_table, line)) {
            binary_table_values_t test = load_binary_snarl_line(line);
            bool match = false;
            if (test == truth1) {
                match = true;
                found_snarls[0] = true;
            } else if (test == truth2) {
                match = true;
                found_snarls[1] = true;
            } else if (test == truth3) {
                match = true;
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

        std::string cmd = "../bin/stoat graph";

        cmd +=  " -g " + graph_base + ".gbz"
            + " -d " + graph_base + ".dist"
            + " -b " + samples_file
            + " -T exact -r path0";


        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/binary_table_graph.tsv"));

        binary_table_values_t truth ({(std::string)"path0#0#path0", 
                                  (size_t)4, 
                                  (size_t)5, 
                                  std::make_pair<size_t, size_t>(5,7), 
                                  std::vector<std::string>({"0","1"}),
                                  (std::string) "NA",
                                  (std::string) "NA", 
                                  std::vector<std::string>({"NA"}), 
                                  (size_t)2});
        ifstream in_table;
        in_table.open(output_dir+"/binary_table_graph.tsv");
        std::string line; 
        std::getline(in_table, line); // Header
        std::getline(in_table, line);
        binary_table_values_t test = load_binary_snarl_line(line);
        REQUIRE(test == truth);
        REQUIRE((!std::getline(in_table, line)));
        in_table.close();
    }

    SECTION("Test chi2 fasta output") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph";

        cmd +=  " -g " + graph_base + ".gbz"
            + " -d " + graph_base + ".dist"
            + " -b " + samples_file
            + " -T chi2 -r path0 -O fasta";


        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/binary_output.fasta"));

        REQUIRE(is_valid_fasta(output_dir+"/binary_output.fasta"));

        std::vector<std::tuple<size_t, std::string, std::string>> truth_fasta;

        truth_fasta.emplace_back(1, ">snarl:1_4|path0#0#path0:1-2|path0#0#path0:1-2", "C");
        truth_fasta.emplace_back(1, ">snarl:1_4|path0#0#path0:1-2|path1#0#path1#0:1-2", "C");
        truth_fasta.emplace_back(2, ">snarl:1_4|path0#0#path0:1-2|path2:1-2", "C");
        truth_fasta.emplace_back(2, ">snarl:1_4|path0#0#path0:1-2|path3:1-2", "C");

        truth_fasta.emplace_back(3, ">snarl:4_8|path0#0#path0:3-6|path0#0#path0:3-6", "TNA");
        truth_fasta.emplace_back(3, ">snarl:4_8|path0#0#path0:3-6|path1#0#path1#0:3-6", "TNA");
        truth_fasta.emplace_back(3, ">snarl:4_8|path0#0#path0:3-6|path3:3-5", "TNA");
        truth_fasta.emplace_back(4, ">snarl:4_8|path0#0#path0:3-6|path2:3-3", "");

        truth_fasta.emplace_back(5, ">snarl:5_7|path0#0#path0:4-5|path0#0#path0:4-5", "C");
        truth_fasta.emplace_back(6, ">snarl:5_7|path0#0#path0:4-5|path1#0#path1#0:4-4", "");
        truth_fasta.emplace_back(6, ">snarl:5_7|path0#0#path0:4-5|path3:4-4", "");

        REQUIRE(fasta_equal(output_dir+"/binary_output.fasta", truth_fasta));

    }

    SECTION("Test exact fasta output") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph";

        cmd +=  " -g " + graph_base + ".gbz"
            + " -d " + graph_base + ".dist"
            + " -b " + samples_file
            + " -T exact -r path0 -O fasta";


        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/binary_output.fasta"));

        REQUIRE(is_valid_fasta(output_dir+"/binary_output.fasta"));

        std::vector<std::tuple<size_t, std::string, std::string>> truth_fasta;

        truth_fasta.emplace_back(2, ">snarl:5_7|path0#0#path0:4-5|path0#0#path0:4-5", "C");

        truth_fasta.emplace_back(1, ">snarl:5_7|path0#0#path0:4-5|path1#0#path1#0:4-4", "");
        truth_fasta.emplace_back(1, ">snarl:5_7|path0#0#path0:4-5|path3:4-4", "");
        REQUIRE(fasta_equal(output_dir+"/binary_output.fasta", truth_fasta));


    }

    SECTION("Test chi2 tsv output with maf excluding one snarl") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph";

        cmd +=  " -g " + graph_base + ".gbz"
            + " -d " + graph_base + ".dist"
            + " -b " + samples_file
            + " -M 0.26"
            + " -T chi2 -r path0 -V 4";


        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/binary_table_graph.tsv"));


        binary_table_values_t truth1 ({(std::string)"path0#0#path0", 
                                  (size_t)1, 
                                  (size_t)2, 
                                  std::make_pair<size_t, size_t>(1,4), 
                                  std::vector<std::string>({"1","1"}),
                                  (std::string) "1",
                                  (std::string) "1", 
                                  std::vector<std::string>({"1:1","1:1"}), 
                                  (size_t)1});

        binary_table_values_t truth3 ({(std::string)"path0#0#path0", 
                                  (size_t)4, 
                                  (size_t)5, 
                                  std::make_pair<size_t, size_t>(5,7), 
                                  std::vector<std::string>({"0","1"}),
                                  (std::string) "0.3333",
                                  (std::string) "8.3265e-02", 
                                  std::vector<std::string>({"2:0","0:1"}), 
                                  (size_t)2});

        ifstream in_table;
        in_table.open(output_dir+"/binary_table_graph.tsv");
        std::string line; 
        std::getline(in_table, line); // Header
        std::vector<bool> found_snarls (3, false);
        while (std::getline(in_table, line)) {
            binary_table_values_t test = load_binary_snarl_line(line);
            bool match = false;
            if (test == truth1) {
                match = true;
                found_snarls[0] = true;
            } else if (test == truth3) {
                match = true;
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

        std::string cmd = "../bin/stoat graph";

        cmd +=  " -g " + graph_base + ".gbz"
            + " -d " + graph_base + ".dist"
            + " -b " + samples_file
            + " -I 4"
            + " -T chi2 -r path0 -V 4";


        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/binary_table_graph.tsv"));


        binary_table_values_t truth1 ({(std::string)"path0#0#path0", 
                                  (size_t)1, 
                                  (size_t)2, 
                                  std::make_pair<size_t, size_t>(1,4), 
                                  std::vector<std::string>({"1","1"}),
                                  (std::string) "1",
                                  (std::string) "1", 
                                  std::vector<std::string>({"1:1","1:1"}), 
                                  (size_t)1});

        binary_table_values_t truth2 ({(std::string)"path0#0#path0", 
                                  (size_t)3, 
                                  (size_t)6, 
                                  std::make_pair<size_t, size_t>(4,8), 
                                  std::vector<std::string>({"0","2/3"}),
                                  (std::string) "1",
                                  (std::string) "0.2482", 
                                  std::vector<std::string>({"0:1","2:1"}), 
                                  (size_t)1});
        ifstream in_table;
        in_table.open(output_dir+"/binary_table_graph.tsv");
        std::string line; 
        std::getline(in_table, line); // Header
        std::vector<bool> found_snarls (3, false);
        while (std::getline(in_table, line)) {
            binary_table_values_t test = load_binary_snarl_line(line);
            bool match = false;
            if (test == truth1) {
                match = true;
                found_snarls[0] = true;
            } else if (test == truth2) {
                match = true;
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

TEST_CASE("Output loop with snarl", "[graph][bug]") {
    const std::string output_dir = "../output_binary";
    const std::string graph_base = "../tests/graph_test/loop_with_indel";

    std::vector<std::string> samples_of_interest = {"path1", "path2"};
    std::vector<std::string> other_samples = {"path0"};

    const std::string samples_file = "./samples.tsv";
    
    string write_cmd = "echo \"FID\tIID\tPHENO\" > " + samples_file;
    int ignore = std::system(write_cmd.c_str());
    for (auto sample : samples_of_interest) {
        write_cmd = "echo \"" + sample + "\t" + sample + "\t2\" >> " + samples_file;
        ignore = std::system(write_cmd.c_str());
    }
    for (auto sample : other_samples) {
        write_cmd = "echo \"" + sample + "\t" + sample + "\t1\" >> " + samples_file;
        ignore = std::system(write_cmd.c_str());
    }

    SECTION("Test chi2 tsv output") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph";
        cmd += " -g " + graph_base + ".hg"
            + " -d " + graph_base + ".dist"
            + " -b " + samples_file
            + " -T chi2 -r path0";


        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/binary_table_graph.tsv"));

        binary_table_values_t truth1 ({(std::string)"path0", 
                                  (size_t)10, 
                                  (size_t)14, 
                                  std::make_pair<size_t, size_t>(1,6), 
                                  std::vector<std::string>({"3/4","6/8"}),
                                  (std::string) "0.3333",
                                  (std::string) "8.3265e-02", 
                                  std::vector<std::string>({"0:1","2:0"}), 
                                  (size_t)1});

        binary_table_values_t truth2 ({(std::string)"path0", 
                                  (size_t)11, 
                                  (size_t)12, 
                                  std::make_pair<size_t, size_t>(2,4), 
                                  std::vector<std::string>({"1","0","1"}),
                                  (std::string) "NA",
                                  (std::string) "0.2231", 
                                  std::vector<std::string>({"0:1","1:0", "1:0"}), 
                                  (size_t)2});

        ifstream in_table;
        in_table.open(output_dir+"/binary_table_graph.tsv");
        std::string line; 
        std::getline(in_table, line); // Header
        std::vector<bool> found_snarls (2, false);
        while (std::getline(in_table, line)) {
            binary_table_values_t test = load_binary_snarl_line(line);
            cerr << "LINE: " << line << endl;
            bool match = false;
            if (test == truth1) {
                match = true;
                found_snarls[0] = true;
            } else if (test == truth2) {
                match = true;
                found_snarls[1] = true;
            }
        }
        in_table.close();
        REQUIRE(found_snarls[0]);
        REQUIRE(found_snarls[1]);

    }

    SECTION("Test exact tsv output") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph";
        cmd += " -g " + graph_base + ".hg"
            + " -d " + graph_base + ".dist"
            + " -b " + samples_file
            + " -T exact -r path0";


        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/binary_table_graph.tsv"));


        binary_table_values_t truth1 ({(std::string)"path0", 
                                  (size_t)10, 
                                  (size_t)14, 
                                  std::make_pair<size_t, size_t>(1,6), 
                                  std::vector<std::string>({"3/4","6/8"}),
                                  (std::string) "NA",
                                  (std::string) "NA", 
                                  std::vector<std::string>({"NA"}), 
                                  (size_t)1});

        binary_table_values_t truth2 ({(std::string)"path0", 
                                  (size_t)11, 
                                  (size_t)12, 
                                  std::make_pair<size_t, size_t>(2,4), 
                                  std::vector<std::string>({"1","0","1"}),
                                  (std::string) "NA",
                                  (std::string) "NA", 
                                  std::vector<std::string>({"NA"}), 
                                  (size_t)2});

        ifstream in_table;
        in_table.open(output_dir+"/binary_table_graph.tsv");
        std::string line; 
        std::getline(in_table, line); // Header
        std::vector<bool> found_snarls (2, false);
        while (std::getline(in_table, line)) {
            binary_table_values_t test = load_binary_snarl_line(line);
            cerr << "LINE: " << line << endl;
            bool match = false;
            if (test == truth1) {
                match = true;
                found_snarls[0] = true;
            } else if (test == truth2) {
                match = true;
                found_snarls[1] = true;
            }
        }
        in_table.close();
        REQUIRE(found_snarls[0]);
        REQUIRE(found_snarls[1]);
    }

    SECTION("Test chi2 fasta output") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph";
        cmd += " -g " + graph_base + ".hg"
            + " -d " + graph_base + ".dist"
            + " -b " + samples_file
            + " -T chi2 -r path0 -O fasta";


        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }


        REQUIRE(std::filesystem::exists(output_dir+"/binary_output.fasta"));

        REQUIRE(is_valid_fasta(output_dir+"/binary_output.fasta"));

        std::vector<std::tuple<size_t, std::string, std::string>> truth_fasta;
        truth_fasta.emplace_back(1, ">snarl:6_1|path0:10-14|path1:10-17", "ANTANT");

        truth_fasta.emplace_back(1, ">snarl:6_1|path0:10-14|path2:10-17", "ANTANT");

        truth_fasta.emplace_back(2, ">snarl:6_1|path0:10-14|path0:10-14", "ANT");

        // Snarl 2
        // Path0 goes through once with insertion node 3
        // Path 1 goes through twice with deletion
        // Path 2 goes through twice, with insertion then with deletion 
        truth_fasta.emplace_back(3, ">snarl:2_4|path0:11-12|path0:11-12", "G");

        truth_fasta.emplace_back(4, ">snarl:2_4|path0:11-12|path1:11-14", "N");

        truth_fasta.emplace_back(5, ">snarl:2_4|path0:11-12|path2:11-15", "GN");


        REQUIRE(fasta_equal(output_dir+"/binary_output.fasta", truth_fasta));
    }

    SECTION("Test exact fasta output") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph";
        cmd += " -g " + graph_base + ".hg"
            + " -d " + graph_base + ".dist"
            + " -b " + samples_file
            + " -T exact -r path0 -O fasta";


        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }


        REQUIRE(std::filesystem::exists(output_dir+"/binary_output.fasta"));

        REQUIRE(is_valid_fasta(output_dir+"/binary_output.fasta"));

        std::vector<std::tuple<size_t, std::string, std::string>> truth_fasta;
        truth_fasta.emplace_back(1, ">snarl:6_1|path0:10-14|path1:10-17", "ANTANT");

        truth_fasta.emplace_back(1, ">snarl:6_1|path0:10-14|path2:10-17", "ANTANT");

        truth_fasta.emplace_back(2, ">snarl:6_1|path0:10-14|path0:10-14", "ANT");

        // Snarl 2
        // Path0 goes through once with insertion node 3
        // Path 1 goes through twice with deletion
        // Path 2 goes through twice, with insertion then with deletion 
        truth_fasta.emplace_back(3, ">snarl:2_4|path0:11-12|path0:11-12", "G");

        truth_fasta.emplace_back(4, ">snarl:2_4|path0:11-12|path1:11-14", "N");

        truth_fasta.emplace_back(5, ">snarl:2_4|path0:11-12|path2:11-15", "GN");

        REQUIRE(fasta_equal(output_dir+"/binary_output.fasta", truth_fasta));

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
    const std::string data_path = "../tests/graph_test/";
    const std::string graph_base = "simple_nested_chains";

    SECTION("Test snarls output") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph";

        cmd +=" -g " + data_path + graph_base + ".hg"
            + " -d " + data_path + graph_base + ".dist"
            + " -M 0 -l 0 -s " + output_dir + graph_base + ".snarls.txt" 
            + " -T chi2 -r ref";

        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE(false);
        }

        REQUIRE(std::filesystem::exists(output_dir+graph_base+".snarls.txt"));
        std::ifstream snarlsfile;
        snarlsfile.open(output_dir+graph_base+".snarls.txt");
        REQUIRE(snarlsfile.peek() != std::ifstream::traits_type::eof());

        size_t line_count = 0;
        bool start_counting= false;
        std::string line;
        while (std::getline(snarlsfile, line)) {
            // Only start counting snarls after proper header
            if (start_counting) {
                line_count++;
            }
            if (line == "#SNARLS") {
                start_counting = true;
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

        std::string cmd = "../bin/stoat graph";

        cmd +=" -g " + data_path + graph_base + ".hg"
            + " -d " + data_path + graph_base + ".dist"
            + " -M 0 -l 0 -s " + output_dir + graph_base + ".snarls.txt" 
            + " -T chi2 -r ref -t 4";

        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE(false);
        }


        REQUIRE(std::filesystem::exists(output_dir+graph_base+".snarls.txt"));
        std::ifstream snarlsfile;
        snarlsfile.open(output_dir+graph_base+".snarls.txt");
        REQUIRE(snarlsfile.peek() != std::ifstream::traits_type::eof());

        size_t line_count = 0;
        bool start_counting= false;
        std::string line;
        while (std::getline(snarlsfile, line)) {
            // Only start counting snarls after proper header
            if (start_counting) {
                line_count++;
            }
            if (line == "#SNARLS") {
                start_counting = true;
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

