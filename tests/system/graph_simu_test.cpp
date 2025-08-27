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

    std::vector<std::string> samples_of_interest = {"samp_g0_60_h1", "samp_g0_8_h0", "samp_g0_73_h1", "samp_g0_12_h0", "samp_g0_8_h1", "samp_g0_94_h1", "samp_g0_94_h0", "samp_g0_78_h1", "samp_g0_16_h1", "samp_g0_16_h0", "samp_g0_19_h0", "samp_g0_19_h1", "samp_g0_78_h0", "samp_g0_75_h1", "samp_g0_75_h0", "samp_g0_38_h1", "samp_g0_61_h1", "samp_g0_21_h0", "samp_g0_21_h1", "samp_g0_30_h1", "samp_g0_38_h0", "samp_g0_70_h1", "samp_g0_30_h0", "samp_g0_70_h0", "samp_g0_61_h0", "samp_g0_99_h0", "samp_g0_99_h1", "samp_g0_59_h1", "samp_g0_59_h0", "samp_g0_7_h1", "samp_g0_32_h0", "samp_g0_7_h0", "samp_g0_32_h1", "samp_g0_50_h0", "samp_g0_50_h1", "samp_g0_24_h0", "samp_g0_24_h1", "samp_g0_22_h0", "samp_g0_22_h1", "samp_g0_77_h1", "samp_g0_77_h0", "samp_g0_9_h0", "samp_g0_9_h1", "samp_g0_48_h1", "samp_g0_48_h0", "samp_g0_46_h1", "samp_g0_20_h1", "samp_g0_46_h0", "samp_g0_64_h0", "samp_g0_95_h1", "samp_g0_64_h1", "samp_g0_95_h0", "samp_g0_18_h0", "samp_g0_18_h1", "samp_g0_49_h1", "samp_g0_20_h0", "samp_g0_52_h0", "samp_g0_52_h1", "samp_g0_91_h1", "samp_g0_91_h0", "samp_g0_41_h1", "samp_g0_41_h0", "samp_g0_23_h1", "samp_g0_3_h1", "samp_g0_3_h0", "samp_g0_23_h0", "samp_g0_28_h0", "samp_g0_28_h1", "samp_g0_93_h0", "samp_g0_93_h1", "samp_g0_39_h0", "samp_g0_39_h1", "samp_g0_76_h0", "samp_g0_88_h0", "samp_g0_76_h1", "samp_g0_88_h1", "samp_g0_29_h0", "samp_g0_29_h1", "samp_g0_1_h0", "samp_g0_1_h1", "samp_g0_10_h1", "samp_g0_10_h0", "samp_g0_25_h0", "samp_g0_25_h1", "samp_g0_0_h0", "samp_g0_0_h1", "samp_g0_11_h0", "samp_g0_11_h1", "samp_g0_2_h0", "samp_g0_2_h1", "samp_g0_56_h0", "samp_g0_56_h1", "samp_g0_13_h1", "samp_g0_13_h0", "samp_g0_55_h0", "samp_g0_55_h1", "samp_g0_98_h0", "samp_g0_98_h1", "samp_g0_71_h1", "samp_g0_71_h0", "samp_g0_83_h0", "samp_g0_31_h1", "samp_g0_31_h0", "samp_g0_96_h0", "samp_g0_83_h1", "samp_g0_96_h1", "samp_g0_51_h0", "samp_g0_51_h1", "samp_g0_72_h0", "samp_g0_72_h1", "samp_g0_47_h0", "samp_g0_47_h1", "samp_g0_34_h0", "samp_g0_34_h1", "samp_g0_65_h0", "samp_g0_65_h1", "samp_g0_4_h1", "samp_g0_4_h0", "samp_g0_90_h1", "samp_g0_90_h0", "samp_g0_35_h0", "samp_g0_35_h1", "samp_g0_67_h0", "samp_g0_67_h1", "samp_g0_84_h0", "samp_g0_66_h0", "samp_g0_66_h1", "samp_g0_84_h1", "samp_g0_68_h1", "samp_g0_26_h1", "samp_g0_26_h0", "samp_g0_85_h1", "samp_g0_85_h0", "samp_g0_89_h0", "samp_g0_89_h1", "samp_g0_68_h0", "samp_g0_14_h0", "samp_g0_42_h1", "samp_g0_54_h1", "samp_g0_42_h0", "samp_g0_54_h0", "samp_g0_14_h1", "samp_g0_62_h1", "samp_g0_6_h1", "samp_g0_36_h1", "samp_g0_49_h0", "samp_g0_36_h0", "samp_g0_6_h0", "samp_g0_81_h0", "samp_g0_62_h0", "samp_g0_81_h1", "samp_g0_40_h0", "samp_g0_43_h1", "samp_g0_40_h1", "samp_g0_43_h0", "samp_g0_45_h0", "samp_g0_82_h0", "samp_g0_45_h1", "samp_g0_82_h1", "samp_g0_5_h1", "samp_g0_5_h0", "samp_g0_27_h0", "samp_g0_27_h1", "samp_g0_92_h0", "samp_g0_92_h1", "samp_g0_87_h0", "samp_g0_87_h1", "samp_g0_86_h1", "samp_g0_86_h0", "samp_g0_74_h0", "samp_g0_74_h1", "samp_g0_57_h1", "samp_g0_57_h0", "samp_g0_53_h1", "samp_g0_97_h0", "samp_g0_53_h0", "samp_g0_97_h1", "samp_g0_17_h1", "samp_g0_17_h0", "samp_g0_69_h1", "samp_g0_69_h0", "samp_g0_33_h1", "samp_g0_33_h0", "samp_g0_37_h0", "samp_g0_37_h1", "samp_g0_44_h0", "samp_g0_63_h1", "samp_g0_44_h1", "samp_g0_63_h0", "samp_g0_79_h0", "samp_g0_79_h1", "samp_g0_15_h1", "samp_g0_15_h0", "samp_g0_80_h1", "samp_g0_80_h0", "samp_g0_12_h1", "samp_g0_58_h1", "samp_g0_58_h0", "samp_g0_73_h0", "samp_g0_60_h0"};

    SECTION("Test tsv output") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph";

        cmd +=  " -g " + data_path + "/" + graph_base + ".pg"
            + " -d " + data_path + "/" + graph_base + ".dist"
            + " -T chi2 -r ref";

        for (const auto& sample : samples_of_interest) {
            cmd += " -s " + sample;
        }

        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }

        // TODO: Add something that actually checks this
        //bool passed = compare_output_dirs(output_dir, expected_dir);
        //REQUIRE(passed);

    }

    SECTION("Test chi2 fasta output") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph";

        cmd +=  " -g " + data_path + "/" + graph_base + ".pg"
            + " -d " + data_path + "/" + graph_base + ".dist"
            + " -T chi2 -r ref -O fasta";

        for (const auto& sample : samples_of_interest) {
            cmd += " -s " + sample;
        }

        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/associated.fasta"));
        REQUIRE(std::filesystem::exists(output_dir+"/unassociated.fasta"));

        REQUIRE(is_valid_fasta(output_dir+"/associated.fasta"));
        REQUIRE(is_valid_fasta(output_dir+"/unassociated.fasta"));
    }

    SECTION("Test exact fasta output") {


        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph";

        cmd +=  " -g " + data_path + "/" + graph_base + ".pg"
            + " -d " + data_path + "/" + graph_base + ".dist"
            + " -T exact -r ref -O fasta";

        for (const auto& sample : samples_of_interest) {
            cmd += " -s " + sample;
        }

        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/associated.fasta"));
        REQUIRE(std::filesystem::exists(output_dir+"/unassociated.fasta"));

        REQUIRE(is_valid_fasta(output_dir+"/associated.fasta"));
        REQUIRE(is_valid_fasta(output_dir+"/unassociated.fasta"));
    }

    clean_output_dir(output_dir);
}

TEST_CASE("Output simple nested chain", "[graph]") {
    const std::string output_dir = "../output_binary";
    const std::string graph_base = "../tests/graph_test/simple_nested_chain";

    std::vector<std::string> samples_of_interest = {"path1", "path3"};

    SECTION("Test chi2 tsv output") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph";

        cmd +=  " -g " + graph_base + ".hg"
            + " -d " + graph_base + ".dist"
            + " -T chi2 -r path0";

        for (const auto& sample : samples_of_interest) {
            cmd += " -s " + sample;
        }

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
            + " -T exact -r path0";

        for (const auto& sample : samples_of_interest) {
            cmd += " -s " + sample;
        }

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
            + " -T chi2 -r path0 -O fasta";

        for (const auto& sample : samples_of_interest) {
            cmd += " -s " + sample;
        }

        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/associated.fasta"));
        REQUIRE(std::filesystem::exists(output_dir+"/unassociated.fasta"));

        REQUIRE(is_valid_fasta(output_dir+"/associated.fasta"));
        REQUIRE(is_valid_fasta(output_dir+"/unassociated.fasta"));

        std::vector<std::tuple<size_t, std::string, std::string>> truth_fasta;

        // Since we don't know which are associated or not for this test, this should be empty
        REQUIRE(fasta_equal(output_dir+"/unassociated.fasta", truth_fasta));

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
        REQUIRE(fasta_equal(output_dir+"/associated.fasta", truth_fasta));
    }

    SECTION("Test exact fasta output") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph";

        cmd +=  " -g " + graph_base + ".hg"
            + " -d " + graph_base + ".dist"
            + " -T exact -r path0 -O fasta";

        for (const auto& sample : samples_of_interest) {
            cmd += " -s " + sample;
        }

        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/associated.fasta"));
        REQUIRE(std::filesystem::exists(output_dir+"/unassociated.fasta"));

        REQUIRE(is_valid_fasta(output_dir+"/associated.fasta"));
        REQUIRE(is_valid_fasta(output_dir+"/unassociated.fasta"));

        //Test unassociated
        std::vector<std::tuple<size_t, std::string, std::string>> truth_fasta;
        truth_fasta.emplace_back(1, ">snarl:5-7|path0:4-5|path0:4-5", "C");

        REQUIRE(fasta_equal(output_dir+"/unassociated.fasta", truth_fasta));

        // Test associated
        truth_fasta.clear();

        truth_fasta.emplace_back(1, ">snarl:5-7|path0:4-5|path1:4-4", "");
        truth_fasta.emplace_back(1, ">snarl:5-7|path0:4-5|path3:4-4", "");
        REQUIRE(fasta_equal(output_dir+"/associated.fasta", truth_fasta));
    }

    clean_output_dir(output_dir);
}

TEST_CASE("Output loop with snarl", "[graph][bug]") {
    const std::string output_dir = "../output_binary";
    const std::string graph_base = "../tests/graph_test/loop_with_indel";

    std::vector<std::string> samples_of_interest = {"path1", "path2"};

    SECTION("Test chi2 tsv output") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph";
        cmd += " -g " + graph_base + ".hg"
            + " -d " + graph_base + ".dist"
            + " -T chi2 -r path0";

        for (const auto& sample : samples_of_interest) {
            cmd += " -s " + sample;
        }

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
        truth_lines.emplace_back("path0\t11\t12\t2_4\t0,1\t1.0000\t0.2482\t0.2482\t1:1,2:0\t2");

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
            + " -T exact -r path0";

        for (const auto& sample : samples_of_interest) {
            cmd += " -s " + sample;
        }

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
            + " -T chi2 -r path0 -O fasta";

        for (const auto& sample : samples_of_interest) {
            cmd += " -s " + sample;
        }

        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/associated.fasta"));
        REQUIRE(std::filesystem::exists(output_dir+"/unassociated.fasta"));

        REQUIRE(is_valid_fasta(output_dir+"/associated.fasta"));
        REQUIRE(is_valid_fasta(output_dir+"/unassociated.fasta"));

        std::vector<std::tuple<size_t, std::string, std::string>> truth_fasta;
        truth_fasta.emplace_back(1, ">snarl:6-1|path0:10-14|path0:10-14", "AGCT");

        truth_fasta.emplace_back(2, ">snarl:6-1|path0:10-14|path1:10-16", "ACTACT");
        truth_fasta.emplace_back(2, ">snarl:6-1|path0:10-14|path2:10-17", "ACTAGCT");

        truth_fasta.emplace_back(3, ">snarl:2-4|path0:11-12|path0:11-12", "G");
        truth_fasta.emplace_back(3, ">snarl:2-4|path0:11-12|path2:11-11", "");
        truth_fasta.emplace_back(5, ">snarl:2-4|path0:11-12|path2:14-15", "G");

        truth_fasta.emplace_back(4, ">snarl:2-4|path0:11-12|path1:11-11", "");
        truth_fasta.emplace_back(5, ">snarl:2-4|path0:11-12|path1:14-14", "");

        truth_fasta.emplace_back(4, ">snarl:2-4|path0:11-12|path2:11-11", "");
        truth_fasta.emplace_back(5, ">snarl:2-4|path0:11-12|path2:14-15", "G");

        REQUIRE(fasta_equal(output_dir+"/associated.fasta", truth_fasta));
    }

    SECTION("Test exact fasta output") {

        clean_output_dir(output_dir);

        std::string cmd = "../bin/stoat graph";
        cmd += " -g " + graph_base + ".hg"
            + " -d " + graph_base + ".dist"
            + " -T exact -r path0 -O fasta";

        for (const auto& sample : samples_of_interest) {
            cmd += " -s " + sample;
        }

        cmd += " --output " + output_dir;

        std::cout << "Command run : \n" << cmd << std::endl;

        int command_output = std::system(cmd.c_str());
        if (command_output != 0) {
            std::cerr << "Command failed: " << cmd << "\n";
            REQUIRE( false);
        }

        REQUIRE(std::filesystem::exists(output_dir+"/associated.fasta"));
        REQUIRE(std::filesystem::exists(output_dir+"/unassociated.fasta"));

        REQUIRE(is_valid_fasta(output_dir+"/associated.fasta"));
        REQUIRE(is_valid_fasta(output_dir+"/unassociated.fasta"));

        std::vector<std::tuple<size_t, std::string, std::string>> truth_fasta;
        truth_fasta.emplace_back(1, ">snarl:6-1|path0:10-14|path1:10-16", "ACTACT");
        truth_fasta.emplace_back(1, ">snarl:6-1|path0:10-14|path2:10-17", "ACTAGCT");

        REQUIRE(fasta_equal(output_dir+"/associated.fasta", truth_fasta));

        truth_fasta.clear();
        truth_fasta.emplace_back(1, ">snarl:6-1|path0:10-14|path0:10-14", "AGCT");

        REQUIRE(fasta_equal(output_dir+"/unassociated.fasta", truth_fasta));
    }

    clean_output_dir(output_dir);
}
