#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "compare_files_utils.hpp"

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

        cmd +=  " -g " + data_path + "/" + graph_base + ".pg"
            + " -d " + data_path + "/" + graph_base + ".dist"
            + " -b " + data_path + "/phenotype_samples.tsv"
            + " -T chi2 -r ref";

        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
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

    clean_output_dir(output_dir);
}

TEST_CASE("Output simple nested chain", "[graph][bug]") {
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
        REQUIRE(std::filesystem::exists(output_dir+"/top_variant_binary_graph.tsv"));

        std::vector<std::string> truth_lines;
        truth_lines.emplace_back("#CHR\tSTART_POS\tEND_POS\tSNARL\tPATH_LENGTHS\tP_FISHER\tP_CHI2\tP_ADJUSTED\tGROUP_PATHS\tDEPTH");
        truth_lines.emplace_back("path0\t1\t2\t1_4\t1,1\t1.0000\t1.0000\t1.0000\t1:1,1:1\t1");
        truth_lines.emplace_back("path0\t3\t6\t4_8\t0,3\t1.0000\t0.2482\t0.3723\t2:1,0:1\t1");
        truth_lines.emplace_back("path0\t4\t5\t5_7\t0,1\t0.3333\t0.0833\t0.2499\t0:1,2:0\t2");

        REQUIRE(files_equal(output_dir+"/binary_table_graph.tsv", truth_lines));

        truth_lines.clear();
        truth_lines.emplace_back("#CHR\tSTART_POS\tEND_POS\tSNARL\tPATH_LENGTHS\tP_FISHER\tP_CHI2\tP_ADJUSTED\tGROUP_PATHS\tDEPTH");
        REQUIRE(files_equal(output_dir+"/top_variant_binary_graph.tsv", truth_lines));
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
        REQUIRE(std::filesystem::exists(output_dir+"/top_variant_binary_graph.tsv"));

        std::vector<std::string> truth_lines;
        truth_lines.emplace_back("#CHR\tSTART_POS\tEND_POS\tSNARL\tPATH_LENGTHS\tP_FISHER\tP_CHI2\tP_ADJUSTED\tGROUP_PATHS\tDEPTH");
        truth_lines.emplace_back("path0\t4\t5\t5_7\t0,1\tNA\tNA\t1.0000\tNA\t2");

        REQUIRE(files_equal(output_dir+"/binary_table_graph.tsv", truth_lines));

        truth_lines.clear();
        truth_lines.emplace_back("#CHR\tSTART_POS\tEND_POS\tSNARL\tPATH_LENGTHS\tP_FISHER\tP_CHI2\tP_ADJUSTED\tGROUP_PATHS\tDEPTH");
        REQUIRE(files_equal(output_dir+"/top_variant_binary_graph.tsv", truth_lines));
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

        truth_fasta.emplace_back(1, ">snarl:1-4|path0:1-2|path0:1-2", "C");
        truth_fasta.emplace_back(1, ">snarl:1-4|path0:1-2|path1:1-2", "C");
        truth_fasta.emplace_back(2, ">snarl:1-4|path0:1-2|path2:1-2", "C");
        truth_fasta.emplace_back(2, ">snarl:1-4|path0:1-2|path3:1-2", "C");

        truth_fasta.emplace_back(3, ">snarl:4-8|path0:3-6|path0:3-6", "TCA");
        truth_fasta.emplace_back(3, ">snarl:4-8|path0:3-6|path1:3-6", "TA");
        truth_fasta.emplace_back(3, ">snarl:4-8|path0:3-6|path3:3-6", "TA");
        truth_fasta.emplace_back(4, ">snarl:4-8|path0:3-6|path2:3-3", "");

        truth_fasta.emplace_back(5, ">snarl:5-7|path0:4-5|path0:4-5", "C");
        truth_fasta.emplace_back(6, ">snarl:5-7|path0:4-5|path1:4-4", "");
        truth_fasta.emplace_back(6, ">snarl:5-7|path0:4-5|path3:4-4", "");

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

        truth_fasta.emplace_back(2, ">snarl:5-7|path0:4-5|path0:4-5", "C");

        truth_fasta.emplace_back(1, ">snarl:5-7|path0:4-5|path1:4-4", "");
        truth_fasta.emplace_back(1, ">snarl:5-7|path0:4-5|path3:4-4", "");
        REQUIRE(fasta_equal(output_dir+"/binary_output.fasta", truth_fasta));


    }

    clean_output_dir(output_dir);
    fs::remove(samples_file);
}

TEST_CASE("Output loop with snarl", "[graph]") {
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
        REQUIRE(std::filesystem::exists(output_dir+"/top_variant_binary_graph.tsv"));

        std::vector<std::string> truth_lines;
        truth_lines.emplace_back("#CHR\tSTART_POS\tEND_POS\tSNARL\tPATH_LENGTHS\tP_FISHER\tP_CHI2\tP_ADJUSTED\tGROUP_PATHS\tDEPTH");
        truth_lines.emplace_back("path0\t10\t14\t6_1\t3,4\t0.3333\t0.0833\t0.1666\t0:1,2:0\t1");
        truth_lines.emplace_back("path0\t11\t12\t2_4\t0,1\tNA\t0.2231\t0.2231\t0:1,1:0,1:0\t2");

        REQUIRE(files_equal(output_dir+"/binary_table_graph.tsv", truth_lines));

        truth_lines.clear();
        truth_lines.emplace_back("#CHR\tSTART_POS\tEND_POS\tSNARL\tPATH_LENGTHS\tP_FISHER\tP_CHI2\tP_ADJUSTED\tGROUP_PATHS\tDEPTH");
        REQUIRE(files_equal(output_dir+"/top_variant_binary_graph.tsv", truth_lines));
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
        REQUIRE(std::filesystem::exists(output_dir+"/top_variant_binary_graph.tsv"));

        std::vector<std::string> truth_lines;
        truth_lines.emplace_back("#CHR\tSTART_POS\tEND_POS\tSNARL\tPATH_LENGTHS\tP_FISHER\tP_CHI2\tP_ADJUSTED\tGROUP_PATHS\tDEPTH");
        truth_lines.emplace_back("path0\t10\t14\t6_1\t3,4\tNA\tNA\t1.0000\tNA\t1");

        REQUIRE(files_equal(output_dir+"/binary_table_graph.tsv", truth_lines));

        truth_lines.clear();
        truth_lines.emplace_back("#CHR\tSTART_POS\tEND_POS\tSNARL\tPATH_LENGTHS\tP_FISHER\tP_CHI2\tP_ADJUSTED\tGROUP_PATHS\tDEPTH");
        REQUIRE(files_equal(output_dir+"/top_variant_binary_graph.tsv", truth_lines));
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
        truth_fasta.emplace_back(1, ">snarl:6-1|path0:10-14|path0:10-14", "AGCT");

        truth_fasta.emplace_back(2, ">snarl:6-1|path0:10-14|path1:10-16", "ACTACT");
        truth_fasta.emplace_back(2, ">snarl:6-1|path0:10-14|path2:10-17", "ACTAGCT");

        // For snarl 2-4, each path is in a different group, and there needs to be one of each
        truth_fasta.emplace_back(3, ">snarl:2-4|path0:11-12|path0:11-12", "G");

        truth_fasta.emplace_back(4, ">snarl:2-4|path0:11-12|path1:11-11", "");
        truth_fasta.emplace_back(5, ">snarl:2-4|path0:11-12|path1:14-14", "");

        truth_fasta.emplace_back(6, ">snarl:2-4|path0:11-12|path2:11-12", "G");
        truth_fasta.emplace_back(7, ">snarl:2-4|path0:11-12|path2:15-15", "");


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
        truth_fasta.emplace_back(1, ">snarl:6-1|path0:10-14|path1:10-16", "ACTACT");
        truth_fasta.emplace_back(1, ">snarl:6-1|path0:10-14|path2:10-17", "ACTAGCT");

        truth_fasta.emplace_back(2, ">snarl:6-1|path0:10-14|path0:10-14", "AGCT");

        REQUIRE(fasta_equal(output_dir+"/binary_output.fasta", truth_fasta));

    }


    fs::remove(samples_file);
    clean_output_dir(output_dir);
}
