#include <catch.hpp>
#include "../../src/vcf_parser.hpp"

using namespace stoat_vcf;


class TestVCFParser : VCFParser {
    public: 
    TestVCFParser(bool untangle) :
        VCFParser(untangle) {} 
    using VCFParser::initialize_parser;
    using VCFParser::get_next_chromosome_name;
    using VCFParser::for_each_record_on_chromosome;
    using VCFParser::skip_to_next_chromosome;
    using VCFParser::close_vcf;
    using VCFParser::hap_count;
    using VCFParser::does_sample_have_snarl;
    using VCFParser::get_opposite_snarl_bound;
    using VCFParser::genotypes;
    using VCFParser::snarl_in_to_out;
    using VCFParser::snarl_count;
};

TEST_CASE( "Parse empty vcf", "[vcf_parser]" ) {

    // Write the vcf
    // I didn't actually check that this was a reasonable vcf but the snarls and paths should be right
    std::string vcf_filename = "./test.vcf";
    std::ofstream vcf_out;
    vcf_out.open(vcf_filename);
    vcf_out << "##fileformat=VCFv4.2" << std::endl;
    vcf_out << "##FILTER=<ID=PASS,Description=\"All filters passed\">" << std::endl;
    vcf_out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
    vcf_out << "##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the snarl tree (0=top level)\">" << std::endl;
    vcf_out << "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph\">" << std::endl;
    vcf_out << "##contig=<ID=ref1,length=100>" << std::endl;
    vcf_out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5" << std::endl;
    vcf_out.close();



    SECTION("Make a VCFParser") {
        TestVCFParser parser (true);
        std::vector<std::string> sample_names = parser.initialize_parser (vcf_filename);


        std::string chr = parser.get_next_chromosome_name();
        REQUIRE(chr == "");
    }

    // clean up

    std::string rm_cmd = "rm " + vcf_filename;
    int rm = system(rm_cmd.c_str());

}

TEST_CASE( "Parse vcf simple nested snarl multiple snps", "[vcf_parser]" ) {
    /*      ----------------
          /    3       6    \
         /   /   \   /  \    \
       1 - 2       5  --- 8---9 
             \   /   \   /
               4       7
        repeated 3x
    */

    // Write the vcf
    // I didn't actually check that this was a reasonable vcf but the snarls and paths should be right
    std::string vcf_filename = "./test.vcf";
    std::ofstream vcf_out;
    vcf_out.open(vcf_filename);
    vcf_out << "##fileformat=VCFv4.2" << std::endl;
    vcf_out << "##FILTER=<ID=PASS,Description=\"All filters passed\">" << std::endl;
    vcf_out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
    vcf_out << "##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the snarl tree (0=top level)\">" << std::endl;
    vcf_out << "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph\">" << std::endl;
    vcf_out << "##contig=<ID=ref1,length=100>" << std::endl;
    vcf_out << "##contig=<ID=ref2,length=100>" << std::endl;
    vcf_out << "##contig=<ID=ref3,length=100>" << std::endl;
    vcf_out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5" << std::endl;
    // For the outer snarl, the three alleles are ins, ins, del, so only genotypes 0 or 1 have the ins
    vcf_out << "ref1\t1\t>1>9\tACGTA\tACCTA,A\t60\t.\tLV=0;AT=>1>2>3>5>6>8>9,>1>2>4>5>7>8>9,>1>2>4>5>8>9,>1>6\tGT\t3/0\t2/1\t1/1\t0/3\t1/0" << std::endl;
    vcf_out << "ref1\t2\t>2>5\tCGT\tCCT,A\t60\t.\tLV=1;AT=>2>3>5,>2>4>5\tGT\t0/0\t0/1\t1/1\t0/0\t1/0" << std::endl;
    vcf_out << "ref1\t2\t>5>8\tCGT\tCCT,A\t60\t.\tLV=1;AT=>5>6>8,>5>8\tGT\t0/0\t0/1\t1/1\t0/0\t1/." << std::endl;

    vcf_out << "ref2\t1\t>11>19\tACGTA\tACCTA,A\t60\t.\tLV=0;AT=>11>12>13>15>16>18>19,>11>12>14>15>17>18>19,>11>12>14>15>18>19,>11>16\tGT\t3/0\t2/1\t1/1\t0/3\t1/0" << std::endl;
    vcf_out << "ref2\t2\t>12>15\tCGT\tCCT,A\t60\t.\tLV=1;AT=>12>13>15,>12>14>15\tGT\t0/0\t0/1\t1/1\t0/0\t1/0" << std::endl;
    vcf_out << "ref2\t5\t>15>18\tCGT\tCCT,A\t60\t.\tLV=1;AT=>15>16>18,>15>18\tGT\t0/0\t0/1\t1/1\t0/0\t1/." << std::endl;

    vcf_out << "ref3\t1\t>21>29\tACGTA\tACCTA,A\t60\t.\tLV=0;AT=>21>22>23>25>26>28>29,>21>22>24>25>27>28>29,>21>22>24>25>28>29,>21>26\tGT\t3/0\t2/1\t1/1\t0/3\t1/0" << std::endl;
    vcf_out << "ref3\t2\t>22>25\tCGT\tCCT,A\t60\t.\tLV=1;AT=>22>23>25,>22>24>25\tGT\t0/0\t0/1\t1/1\t0/0\t1/0" << std::endl;
    vcf_out << "ref3\t2\t>25>28\tCGT\tCCT,A\t60\t.\tLV=1;AT=>25>26>28,>25>28\tGT\t0/0\t0/1\t1/1\t0/0\t1/." << std::endl;
    vcf_out.close();



    SECTION("Make a VCFParser without untangling snarls") {
        TestVCFParser parser(false);
        std::vector<std::string> sample_names = parser.initialize_parser (vcf_filename);

        // Check first chr
        std::string chr = parser.get_next_chromosome_name();
        REQUIRE(chr == ("ref1"));
        size_t snarl_num = 0;
        parser.for_each_record_on_chromosome(chr, [&] (const vcf_info_t& vcf_info) {
            if (snarl_num == 0) {
                REQUIRE(vcf_info.lv == 0);
                REQUIRE(vcf_info.paths.size() == 4); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">1>2>3>5>6>8>9");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">1>2>4>5>7>8>9");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[2]) == ">1>2>4>5>8>9");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[3]) == ">1>6");
                REQUIRE(vcf_info.genotype[0] == 3); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 2); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 3); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 1) {
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">2>3>5");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">2>4>5");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 0); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 2) {
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">5>6>8");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">5>8");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 0); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == -1); 
            }
            ++snarl_num;
        });
        REQUIRE(snarl_num == 3);

        // Check second chr
        chr = parser.get_next_chromosome_name();
        REQUIRE(chr == ("ref2"));
        snarl_num = 0;
        parser.for_each_record_on_chromosome(chr, [&] (const vcf_info_t& vcf_info) {
            if (snarl_num == 0) {
                REQUIRE(vcf_info.lv == 0);
                REQUIRE(vcf_info.paths.size() == 4); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">11>12>13>15>16>18>19");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">11>12>14>15>17>18>19");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[2]) == ">11>12>14>15>18>19");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[3]) == ">11>16");
                REQUIRE(vcf_info.genotype[0] == 3); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 2); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 3); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 1) {
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">12>13>15");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">12>14>15");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 0); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 2) {
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">15>16>18");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">15>18");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 0); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == -1); 
            }
            ++snarl_num;
        });
        REQUIRE(snarl_num == 3);

        // Check third chr
        chr = parser.get_next_chromosome_name();
        REQUIRE(chr == ("ref3"));
        snarl_num = 0;
        parser.for_each_record_on_chromosome(chr, [&] (const vcf_info_t& vcf_info) {
            if (snarl_num == 0) {
                REQUIRE(vcf_info.lv == 0);
                REQUIRE(vcf_info.paths.size() == 4); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">21>22>23>25>26>28>29");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">21>22>24>25>27>28>29");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[2]) == ">21>22>24>25>28>29");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[3]) == ">21>26");
                REQUIRE(vcf_info.genotype[0] == 3); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 2); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 3); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 1) {
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">22>23>25");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">22>24>25");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 0); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 2) {
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">25>26>28");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">25>28");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 0); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == -1); 
            }
            ++snarl_num;
        });
        REQUIRE(snarl_num == 3);

        chr = parser.get_next_chromosome_name();
        REQUIRE(chr == "");

        parser.close_vcf();
    }
    SECTION("Skip a chromosome") {
        TestVCFParser parser (false);
        std::vector<std::string> sample_names = parser.initialize_parser (vcf_filename);

        // Check first chr
        std::string chr = parser.get_next_chromosome_name();
        REQUIRE(chr == ("ref1"));
        size_t snarl_num = 0;
        parser.for_each_record_on_chromosome(chr, [&] (const vcf_info_t& vcf_info) {
            if (snarl_num == 0) {
                REQUIRE(vcf_info.lv == 0);
                REQUIRE(vcf_info.paths.size() == 4); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">1>2>3>5>6>8>9");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">1>2>4>5>7>8>9");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[2]) == ">1>2>4>5>8>9");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[3]) == ">1>6");
                REQUIRE(vcf_info.genotype[0] == 3); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 2); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 3); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 1) {
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">2>3>5");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">2>4>5");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 0); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 2) {
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">5>6>8");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">5>8");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 0); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == -1); 
            }
            ++snarl_num;
        });
        REQUIRE(snarl_num == 3);

        // Check second chr
        chr = parser.get_next_chromosome_name();
        REQUIRE(chr == ("ref2"));
        parser.skip_to_next_chromosome(chr);

        // Check third chr
        chr = parser.get_next_chromosome_name();
        REQUIRE(chr == ("ref3"));
        snarl_num = 0;
        parser.for_each_record_on_chromosome(chr, [&] (const vcf_info_t& vcf_info) {
            if (snarl_num == 0) {
                REQUIRE(vcf_info.lv == 0);
                REQUIRE(vcf_info.paths.size() == 4); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">21>22>23>25>26>28>29");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">21>22>24>25>27>28>29");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[2]) == ">21>22>24>25>28>29");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[3]) == ">21>26");
                REQUIRE(vcf_info.genotype[0] == 3); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 2); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 3); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 1) {
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">22>23>25");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">22>24>25");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 0); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 2) {
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">25>26>28");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">25>28");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 0); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == -1); 
            }
            ++snarl_num;
        });
        REQUIRE(snarl_num == 3);

        chr = parser.get_next_chromosome_name();
        REQUIRE(chr == "");

        parser.close_vcf();
    }
    SECTION("Skip last chromosome") {
        TestVCFParser parser (false);
        std::vector<std::string> sample_names = parser.initialize_parser (vcf_filename);

        // Check first chr
        std::string chr = parser.get_next_chromosome_name();
        REQUIRE(chr == ("ref1"));
        size_t snarl_num = 0;
        parser.for_each_record_on_chromosome(chr, [&] (const vcf_info_t& vcf_info) {
            if (snarl_num == 0) {
                REQUIRE(vcf_info.lv == 0);
                REQUIRE(vcf_info.paths.size() == 4); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">1>2>3>5>6>8>9");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">1>2>4>5>7>8>9");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[2]) == ">1>2>4>5>8>9");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[3]) == ">1>6");
                REQUIRE(vcf_info.genotype[0] == 3); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 2); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 3); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 1) {
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">2>3>5");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">2>4>5");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 0); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 2) {
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">5>6>8");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">5>8");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 0); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == -1); 
            }
            ++snarl_num;
        });
        REQUIRE(snarl_num == 3);

        // Check skip chr
        chr = parser.get_next_chromosome_name();
        REQUIRE(chr == ("ref2"));
        parser.skip_to_next_chromosome(chr);

        //skip last chr
        chr = parser.get_next_chromosome_name();
        REQUIRE(chr == ("ref3"));
        parser.skip_to_next_chromosome(chr);

        chr = parser.get_next_chromosome_name();
        REQUIRE(chr == "");

        parser.close_vcf();
    }

    // clean up

    std::string rm_cmd = "rm " + vcf_filename;
    int rm = system(rm_cmd.c_str());

}

TEST_CASE( "Parse pangenie vcf simple nested snarl multiple snps", "[vcf_parser]" ) {
    /*      ----------------
          /    3       6    \
         /   /   \   /  \    \
       1 - 2       5  --- 8---9 
             \   /   \   /
               4       7
        repeated 3x
    */

    // Write the vcf
    // I didn't actually check that this was a reasonable vcf but the snarls and paths should be right
    std::string vcf_filename = "./test.vcf";
    std::ofstream vcf_out;
    vcf_out.open(vcf_filename);
    vcf_out << "##fileformat=VCFv4.2" << std::endl;
    vcf_out << "##FILTER=<ID=PASS,Description=\"All filters passed\">" << std::endl;
    vcf_out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
    vcf_out << "##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the snarl tree (0=top level)\">" << std::endl;
    vcf_out << "##INFO=<ID=ID,Number=A,Type=String,Description=\"Variant IDs per ALT allele.\">" << std::endl;
    vcf_out << "##INFO=<ID=RD,Number=1,Type=String,Description=\"Variant IDs per REF allele.\">" << std::endl;
    vcf_out << "##contig=<ID=ref1,length=100>" << std::endl;
    vcf_out << "##contig=<ID=ref2,length=100>" << std::endl;
    vcf_out << "##contig=<ID=ref3,length=100>" << std::endl;
    vcf_out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5" << std::endl;
    // For the outer snarl, the three alleles are ins, ins, del, so only genotypes 0 or 1 have the ins
    vcf_out << "ref1\t1\t.\tACGTA\tACCTA,A\t60\t.\tLV=0;RD=ref-1-REF->1>2>3>5>6>8>9-3;ID=path0-1-SNP->1>2>4>5>7>8>9-1,path1-1-SNP->1>2>4>5>8>9-3,path2-1-DEL->1>6-1\tGT\t3/0\t2/1\t1/1\t0/3\t1/0" << std::endl;
    vcf_out << "ref1\t2\t.\tCGT\tCCT,A\t60\t.\tLV=1;RD=ref-5-REF->2>3>5-1;ID=path0-3-SNP->2>4>5-2\tGT\t0/0\t0/1\t1/1\t0/0\t1/0" << std::endl;
    vcf_out << "ref1\t2\t.\tCGT\tCCT,A\t60\t.\tLV=1;ID=path0-3-SNP->5>8-3;RD=ref-3-REF->5>6>8-3\tGT\t0/0\t0/1\t1/1\t0/0\t1/." << std::endl;

    vcf_out << "ref2\t1\t.19\tACGTA\tACCTA,A\t60\t.\tLV=0;RD=path0-3-SNP->11>12>13>15>16>18>19;ID=path0-3-SNP->11>12>14>15>17>18>19-3,path2-3-SNP->11>12-3:path2-3-SNP->14>15>18>19-3,path2-3-SNP->11>16-DEL\tGT\t3/0\t2/1\t1/1\t0/3\t1/0" << std::endl;
    vcf_out << "ref2\t2\t.15\tCGT\tCCT,A\t60\t.\tLV=1;RD=path0-3-SNP->12>13>15-1;ID=path0-3-SNP->12>14>15-1\tGT\t0/0\t0/1\t1/1\t0/0\t1/0" << std::endl;
    vcf_out << "ref2\t5\t.18\tCGT\tCCT,A\t60\t.\tLV=1;RD=path0-3-SNP->15>16>18-1;ID=path0-3-SNP->15>18-3\tGT\t0/0\t0/1\t1/1\t0/0\t1/." << std::endl;

    vcf_out << "ref3\t1\t.29\tACGTA\tACCTA,A\t60\t.\tLV=0;RD=path0-3-SNP->21>22-3:path0-3-SNP->23>25>26>28>29-3;ID=path0-3-SNP->21>22>24>25>27>28>29-3,path0-3-SNP->21>22>24>25>28>29-34,path0-3-SNP->21>26-3\tGT\t3/0\t2/1\t1/1\t0/3\t1/0" << std::endl;
    vcf_out << "ref3\t2\t.25\tCGT\tCCT,A\t60\t.\tLV=1;RD=path0-3-SNP->22>23>25-2;ID=path0-3-SNP->22>24>25-1\tGT\t0/0\t0/1\t1/1\t0/0\t1/0" << std::endl;
    vcf_out << "ref3\t2\t.28\tCGT\tCCT,A\t60\t.\tLV=1;RD=path0-3-SNP->25>26>28-1;ID=path0-3-SNP->25>28-5\tGT\t0/0\t0/1\t1/1\t0/0\t1/." << std::endl;
    vcf_out.close();



    SECTION("Make a VCFParser without untangling snarls") {
        TestVCFParser parser(false);
        std::vector<std::string> sample_names = parser.initialize_parser (vcf_filename);

        // Check first chr
        std::string chr = parser.get_next_chromosome_name();
        REQUIRE(chr == ("ref1"));
        size_t snarl_num = 0;
        parser.for_each_record_on_chromosome(chr, [&] (const vcf_info_t& vcf_info) {
            if (snarl_num == 0) {
                REQUIRE(vcf_info.lv == 0);
                REQUIRE(vcf_info.paths.size() == 4); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">1>2>3>5>6>8>9");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">1>2>4>5>7>8>9");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[2]) == ">1>2>4>5>8>9");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[3]) == ">1>6");
                REQUIRE(vcf_info.genotype[0] == 3); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 2); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 3); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 1) {
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">2>3>5");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">2>4>5");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 0); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 2) {
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">5>6>8");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">5>8");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 0); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == -1); 
            }
            ++snarl_num;
        });
        REQUIRE(snarl_num == 3);

        // Check second chr
        chr = parser.get_next_chromosome_name();
        REQUIRE(chr == ("ref2"));
        snarl_num = 0;
        parser.for_each_record_on_chromosome(chr, [&] (const vcf_info_t& vcf_info) {
            if (snarl_num == 0) {
                REQUIRE(vcf_info.lv == 0);
                REQUIRE(vcf_info.paths.size() == 4); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">11>12>13>15>16>18>19");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">11>12>14>15>17>18>19");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[2]) == ">11>12>0>14>15>18>19");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[3]) == ">11>16");
                REQUIRE(vcf_info.genotype[0] == 3); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 2); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 3); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 1) {
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">12>13>15");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">12>14>15");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 0); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 2) {
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">15>16>18");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">15>18");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 0); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == -1); 
            }
            ++snarl_num;
        });
        REQUIRE(snarl_num == 3);

        // Check third chr
        chr = parser.get_next_chromosome_name();
        REQUIRE(chr == ("ref3"));
        snarl_num = 0;
        parser.for_each_record_on_chromosome(chr, [&] (const vcf_info_t& vcf_info) {
            if (snarl_num == 0) {
                REQUIRE(vcf_info.lv == 0);
                REQUIRE(vcf_info.paths.size() == 4); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">21>22>0>23>25>26>28>29");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">21>22>24>25>27>28>29");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[2]) == ">21>22>24>25>28>29");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[3]) == ">21>26");
                REQUIRE(vcf_info.genotype[0] == 3); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 2); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 3); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 1) {
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">22>23>25");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">22>24>25");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 0); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 2) {
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">25>26>28");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">25>28");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 0); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == -1); 
            }
            ++snarl_num;
        });
        REQUIRE(snarl_num == 3);

        chr = parser.get_next_chromosome_name();
        REQUIRE(chr == "");

        parser.close_vcf();
    }
    SECTION("Skip a chromosome") {
        TestVCFParser parser (false);
        std::vector<std::string> sample_names = parser.initialize_parser (vcf_filename);

        // Check first chr
        std::string chr = parser.get_next_chromosome_name();
        REQUIRE(chr == ("ref1"));
        size_t snarl_num = 0;
        parser.for_each_record_on_chromosome(chr, [&] (const vcf_info_t& vcf_info) {
            if (snarl_num == 0) {
                REQUIRE(vcf_info.lv == 0);
                REQUIRE(vcf_info.paths.size() == 4); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">1>2>3>5>6>8>9");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">1>2>4>5>7>8>9");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[2]) == ">1>2>4>5>8>9");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[3]) == ">1>6");
                REQUIRE(vcf_info.genotype[0] == 3); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 2); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 3); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 1) {
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">2>3>5");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">2>4>5");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 0); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 2) {
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">5>6>8");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">5>8");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 0); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == -1); 
            }
            ++snarl_num;
        });
        REQUIRE(snarl_num == 3);

        // Check second chr
        chr = parser.get_next_chromosome_name();
        REQUIRE(chr == ("ref2"));
        parser.skip_to_next_chromosome(chr);

        // Check third chr
        chr = parser.get_next_chromosome_name();
        REQUIRE(chr == ("ref3"));
        snarl_num = 0;
        parser.for_each_record_on_chromosome(chr, [&] (const vcf_info_t& vcf_info) {
            if (snarl_num == 0) {
                REQUIRE(vcf_info.lv == 0);
                REQUIRE(vcf_info.paths.size() == 4); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">21>22>0>23>25>26>28>29");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">21>22>24>25>27>28>29");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[2]) == ">21>22>24>25>28>29");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[3]) == ">21>26");
                REQUIRE(vcf_info.genotype[0] == 3); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 2); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 3); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 1) {
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">22>23>25");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">22>24>25");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 0); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 2) {
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">25>26>28");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">25>28");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 0); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == -1); 
            }
            ++snarl_num;
        });
        REQUIRE(snarl_num == 3);

        chr = parser.get_next_chromosome_name();
        REQUIRE(chr == "");

        parser.close_vcf();
    }
    SECTION("Skip last chromosome") {
        TestVCFParser parser (false);
        std::vector<std::string> sample_names = parser.initialize_parser (vcf_filename);

        // Check first chr
        std::string chr = parser.get_next_chromosome_name();
        REQUIRE(chr == ("ref1"));
        size_t snarl_num = 0;
        parser.for_each_record_on_chromosome(chr, [&] (const vcf_info_t& vcf_info) {
            if (snarl_num == 0) {
                REQUIRE(vcf_info.lv == 0);
                REQUIRE(vcf_info.paths.size() == 4); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">1>2>3>5>6>8>9");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">1>2>4>5>7>8>9");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[2]) == ">1>2>4>5>8>9");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[3]) == ">1>6");
                REQUIRE(vcf_info.genotype[0] == 3); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 2); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 3); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 1) {
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">2>3>5");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">2>4>5");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 0); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 2) {
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">5>6>8");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">5>8");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 0); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == -1); 
            }
            ++snarl_num;
        });
        REQUIRE(snarl_num == 3);

        // Check skip chr
        chr = parser.get_next_chromosome_name();
        REQUIRE(chr == ("ref2"));
        parser.skip_to_next_chromosome(chr);

        //skip last chr
        chr = parser.get_next_chromosome_name();
        REQUIRE(chr == ("ref3"));
        parser.skip_to_next_chromosome(chr);

        chr = parser.get_next_chromosome_name();
        REQUIRE(chr == "");

        parser.close_vcf();
    }

    // clean up

    std::string rm_cmd = "rm " + vcf_filename;
    int rm = system(rm_cmd.c_str());

}


TEST_CASE( "Untangle simple nested snarl", "[vcf_parser]" ) {
    /*      --------
          /    3     \
         /   /   \    \
       1 - 2       5 - 6 
             \   /
               4
    */

    // Write the vcf
    // I didn't actually check that this was a reasonable vcf but the snarls and paths should be right
    std::string vcf_filename = "./test.vcf";
    std::ofstream vcf_out;
    vcf_out.open(vcf_filename);
    vcf_out << "##fileformat=VCFv4.2" << std::endl;
    vcf_out << "##FILTER=<ID=PASS,Description=\"All filters passed\">" << std::endl;
    vcf_out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
    vcf_out << "##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the snarl tree (0=top level)\">" << std::endl;
    vcf_out << "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph\">" << std::endl;
    vcf_out << "##contig=<ID=ref,length=100>" << std::endl;
    vcf_out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5" << std::endl;
    // For the outer snarl, the three alleles are ins, ins, del, so only genotypes 0 or 1 have the ins
    vcf_out << "ref\t1\t>1>6\tACGTA\tACCTA,A\t60\t.\tLV=0;AT=>1>2>3>5>6,>1>2>4>5>6,>1>6\tGT\t2/0\t0/1\t1/1\t0/2\t1/0" << std::endl;
    vcf_out << "ref\t2\t>2>5\tCGT\tCCT,A\t60\t.\tLV=1;AT=>2>3>5,>2>4>5\tGT\t0/0\t0/1\t1/1\t0/0\t1/0" << std::endl;
    vcf_out.close();



    SECTION("Check the untangler") {
        TestVCFParser parser (true);
        std::vector<std::string> sample_names = parser.initialize_parser (vcf_filename);

        parser.for_each_record_on_chromosome("ref", [&](const auto& x ) {});

        REQUIRE(parser.hap_count == 10);
        REQUIRE(parser.snarl_count == 1);

        REQUIRE(parser.snarl_in_to_out.size() == 2);
        REQUIRE(parser.genotypes.size() == 10);

        REQUIRE(parser.get_opposite_snarl_bound(stoat::node_traversal_t(2, false)) == stoat::node_traversal_t(5, false));
        REQUIRE(parser.get_opposite_snarl_bound(stoat::node_traversal_t(5, true)) == stoat::node_traversal_t(2, true));

        // Everything should have the top-level snarl
        REQUIRE(parser.does_sample_have_snarl(0, stoat::node_traversal_t(1, false)));
        REQUIRE(parser.does_sample_have_snarl(1, stoat::node_traversal_t(1, false)));
        REQUIRE(parser.does_sample_have_snarl(2, stoat::node_traversal_t(1, false)));
        REQUIRE(parser.does_sample_have_snarl(3, stoat::node_traversal_t(1, false)));
        REQUIRE(parser.does_sample_have_snarl(4, stoat::node_traversal_t(1, false)));
        REQUIRE(parser.does_sample_have_snarl(5, stoat::node_traversal_t(1, false)));
        REQUIRE(parser.does_sample_have_snarl(6, stoat::node_traversal_t(1, false)));
        REQUIRE(parser.does_sample_have_snarl(7, stoat::node_traversal_t(1, false)));
        REQUIRE(parser.does_sample_have_snarl(8, stoat::node_traversal_t(1, false)));
        REQUIRE(parser.does_sample_have_snarl(9, stoat::node_traversal_t(1, false)));

        REQUIRE(!parser.does_sample_have_snarl(0, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(1, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(2, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(3, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(4, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(5, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(6, stoat::node_traversal_t(2, false)));
        REQUIRE(!parser.does_sample_have_snarl(7, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(8, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(9, stoat::node_traversal_t(2, false)));

        REQUIRE(!parser.does_sample_have_snarl(0, stoat::node_traversal_t(5, true)));
        REQUIRE(parser.does_sample_have_snarl(1, stoat::node_traversal_t(5, true)));
        REQUIRE(parser.does_sample_have_snarl(2, stoat::node_traversal_t(5, true)));
        REQUIRE(parser.does_sample_have_snarl(3, stoat::node_traversal_t(5, true)));
        REQUIRE(parser.does_sample_have_snarl(4, stoat::node_traversal_t(5, true)));
        REQUIRE(parser.does_sample_have_snarl(5, stoat::node_traversal_t(5, true)));
        REQUIRE(parser.does_sample_have_snarl(6, stoat::node_traversal_t(5, true)));
        REQUIRE(!parser.does_sample_have_snarl(7, stoat::node_traversal_t(5, true)));
        REQUIRE(parser.does_sample_have_snarl(8, stoat::node_traversal_t(5, true)));
        REQUIRE(parser.does_sample_have_snarl(9, stoat::node_traversal_t(5, true)));

        parser.close_vcf();

    }

    SECTION("Go through the contents using the untangler") {
        TestVCFParser parser (true);
        std::vector<std::string> sample_names = parser.initialize_parser (vcf_filename);

        // Check first chr
        std::string chr = parser.get_next_chromosome_name();
        REQUIRE(chr == ("ref"));
        size_t snarl_num = 0;
        parser.for_each_record_on_chromosome(chr, [&] (const vcf_info_t& vcf_info) {
            if (snarl_num == 0) {
                // For the outer snarl, should be
                //>1>2>3>5>6   >1>2>4>5>6   >1>6
                // With 2-5 being a nested chain
                //2/0 0/1 1/1 0/2 1/0
                REQUIRE(vcf_info.lv == 0);
                REQUIRE(vcf_info.paths.size() == 3); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">1>2>0>5>6");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">1>2>0>5>6");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[2]) == ">1>6");
                REQUIRE(vcf_info.genotype[0] == 2); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 2); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 1) {
                // For the inner snarl, should be
                //>2>3>5   >2>4>5
                // 0/0  0/1  1/1  0/0  1/0
                // .                .
                // With the . meaning that they weren't present according to the outer snarl 
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">2>3>5");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">2>4>5");
                REQUIRE(vcf_info.genotype[0] == -1); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == -1); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            }
            ++snarl_num;
        });
        REQUIRE(snarl_num == 2);



        parser.close_vcf();

    }

    // clean up

    std::string rm_cmd = "rm " + vcf_filename;
    int rm = system(rm_cmd.c_str());

}

TEST_CASE( "Untangle simple nested snarl multiple snps", "[vcf_parser]" ) {
    /*      ----------------
          /    3       6    \
         /   /   \   /  \    \
       1 - 2       5  --- 8---9 
             \   /   \   /
               4       7
    */

    // Write the vcf
    // I didn't actually check that this was a reasonable vcf but the snarls and paths should be right
    std::string vcf_filename = "./test.vcf";
    std::ofstream vcf_out;
    vcf_out.open(vcf_filename);
    vcf_out << "##fileformat=VCFv4.2" << std::endl;
    vcf_out << "##FILTER=<ID=PASS,Description=\"All filters passed\">" << std::endl;
    vcf_out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
    vcf_out << "##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the snarl tree (0=top level)\">" << std::endl;
    vcf_out << "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph\">" << std::endl;
    vcf_out << "##contig=<ID=ref,length=100>" << std::endl;
    vcf_out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5" << std::endl;
    // For the outer snarl, the three alleles are ins, ins, ins, del, so only genotypes 0, 1, and 2 have the ins
    vcf_out << "ref\t1\t>1>9\tACGTA\tACCTA,A\t60\t.\tLV=0;AT=>1>2>3>5>6>8>9,>1>2>4>5>7>8>9,>1>2>4>5>8>9,>1>6\tGT\t3/0\t2/1\t1/1\t0/3\t1/0" << std::endl;
    vcf_out << "ref\t2\t>2>5\tCGT\tCCT,A\t60\t.\tLV=1;AT=>2>3>5,>2>4>5\tGT\t0/0\t0/1\t1/1\t0/0\t1/0" << std::endl;
    vcf_out << "ref\t2\t>5>8\tCGT\tCCT,A\t60\t.\tLV=1;AT=>5>6>8,>5>8\tGT\t0/0\t0/1\t1/1\t0/0\t1/0" << std::endl;
    vcf_out.close();



    SECTION("Make a VCFParser") {
        TestVCFParser parser (true);
        std::vector<std::string> sample_names = parser.initialize_parser (vcf_filename);

        parser.for_each_record_on_chromosome("ref", [&](const auto& x ) {});

        REQUIRE(parser.hap_count == 10);
        REQUIRE(parser.snarl_count == 2);

        REQUIRE(parser.snarl_in_to_out.size() == 4);
        REQUIRE(parser.genotypes.size() == 20);



        REQUIRE(parser.get_opposite_snarl_bound(stoat::node_traversal_t(2, false)) == stoat::node_traversal_t(5, false));
        REQUIRE(parser.get_opposite_snarl_bound(stoat::node_traversal_t(5, true)) == stoat::node_traversal_t(2, true));
        REQUIRE(parser.get_opposite_snarl_bound(stoat::node_traversal_t(5, false)) == stoat::node_traversal_t(8, false));
        REQUIRE(parser.get_opposite_snarl_bound(stoat::node_traversal_t(8, true)) == stoat::node_traversal_t(5, true));

        REQUIRE(!parser.does_sample_have_snarl(0, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(1, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(2, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(3, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(4, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(5, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(6, stoat::node_traversal_t(2, false)));
        REQUIRE(!parser.does_sample_have_snarl(7, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(8, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(9, stoat::node_traversal_t(2, false)));


        REQUIRE(!parser.does_sample_have_snarl(0, stoat::node_traversal_t(5, true)));
        REQUIRE(parser.does_sample_have_snarl(1, stoat::node_traversal_t(5, true)));
        REQUIRE(parser.does_sample_have_snarl(2, stoat::node_traversal_t(5, true)));
        REQUIRE(parser.does_sample_have_snarl(3, stoat::node_traversal_t(5, true)));
        REQUIRE(parser.does_sample_have_snarl(4, stoat::node_traversal_t(5, true)));
        REQUIRE(parser.does_sample_have_snarl(5, stoat::node_traversal_t(5, true)));
        REQUIRE(parser.does_sample_have_snarl(6, stoat::node_traversal_t(5, true)));
        REQUIRE(!parser.does_sample_have_snarl(7, stoat::node_traversal_t(5, true)));
        REQUIRE(parser.does_sample_have_snarl(8, stoat::node_traversal_t(5, true)));
        REQUIRE(parser.does_sample_have_snarl(9, stoat::node_traversal_t(5, true)));

        REQUIRE(!parser.does_sample_have_snarl(0, stoat::node_traversal_t(5, false)));
        REQUIRE(parser.does_sample_have_snarl(1, stoat::node_traversal_t(5, false)));
        REQUIRE(parser.does_sample_have_snarl(2, stoat::node_traversal_t(5, false)));
        REQUIRE(parser.does_sample_have_snarl(3, stoat::node_traversal_t(5, false)));
        REQUIRE(parser.does_sample_have_snarl(4, stoat::node_traversal_t(5, false)));
        REQUIRE(parser.does_sample_have_snarl(5, stoat::node_traversal_t(5, false)));
        REQUIRE(parser.does_sample_have_snarl(6, stoat::node_traversal_t(5, false)));
        REQUIRE(!parser.does_sample_have_snarl(7, stoat::node_traversal_t(5, false)));
        REQUIRE(parser.does_sample_have_snarl(8, stoat::node_traversal_t(5, false)));
        REQUIRE(parser.does_sample_have_snarl(9, stoat::node_traversal_t(5, false)));

        REQUIRE(!parser.does_sample_have_snarl(0, stoat::node_traversal_t(8, true)));
        REQUIRE(parser.does_sample_have_snarl(1, stoat::node_traversal_t(8, true)));
        REQUIRE(parser.does_sample_have_snarl(2, stoat::node_traversal_t(8, true)));
        REQUIRE(parser.does_sample_have_snarl(3, stoat::node_traversal_t(8, true)));
        REQUIRE(parser.does_sample_have_snarl(4, stoat::node_traversal_t(8, true)));
        REQUIRE(parser.does_sample_have_snarl(5, stoat::node_traversal_t(8, true)));
        REQUIRE(parser.does_sample_have_snarl(6, stoat::node_traversal_t(8, true)));
        REQUIRE(!parser.does_sample_have_snarl(7, stoat::node_traversal_t(8, true)));
        REQUIRE(parser.does_sample_have_snarl(8, stoat::node_traversal_t(8, true)));
        REQUIRE(parser.does_sample_have_snarl(9, stoat::node_traversal_t(8, true)));


        parser.close_vcf();
    }
    SECTION("Go through the contents using the untangler") {
        TestVCFParser parser (true);
        std::vector<std::string> sample_names = parser.initialize_parser (vcf_filename);

        // Check first chr
        std::string chr = parser.get_next_chromosome_name();
        REQUIRE(chr == ("ref"));
        size_t snarl_num = 0;
        parser.for_each_record_on_chromosome(chr, [&] (const vcf_info_t& vcf_info) {
            if (snarl_num == 0) {
                // For the outer snarl, should be
                // >1>2>3>5>6>8>9   >1>2>4>5>7>8>9   >1>2>4>5>8>9   >1>6
                // With 2-5 and 5-8 being nested snarls
                //3/0 2/1 1/1 0/3 1/0
                REQUIRE(vcf_info.lv == 0);
                REQUIRE(vcf_info.paths.size() == 4); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">1>2>0>5>0>8>9");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">1>2>0>5>0>8>9");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[2]) == ">1>2>0>5>0>8>9");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[3]) == ">1>6");
                REQUIRE(vcf_info.genotype[0] == 3); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 2); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 3); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 1) {
                // For the inner snarl, should be
                //>2>3>5   >2>4>5
                // 0/0  0/1  1/1  0/0  1/0
                // .                .
                // With the . meaning that they weren't present according to the outer snarl 
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">2>3>5");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">2>4>5");
                REQUIRE(vcf_info.genotype[0] == -1); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == -1); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 2) {
                // For the inner snarl, should be
                //>5>6>8   >5>8
                // 0/0  0/1  1/1  0/0  1/0
                // .                .
                // With the . meaning that they weren't present according to the outer snarl 
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">5>6>8");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">5>8");
                REQUIRE(vcf_info.genotype[0] == -1); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == -1); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            }
            ++snarl_num;
        });
        REQUIRE(snarl_num == 3);



        parser.close_vcf();

    }

    // clean up

    std::string rm_cmd = "rm " + vcf_filename;
    int rm = system(rm_cmd.c_str());

}
TEST_CASE( "Untangle three nested snarl multiple snps", "[vcf_parser]" ) {
    /*
          -------------
        /   --------   \
       /  /    4     \   \
      /  /   /   \    \   \
     1--2 - 3      6 - 7---8 
             \   /
               5
    */

    // Write the vcf
    // I didn't actually check that this was a reasonable vcf but the snarls and paths should be right
    std::string vcf_filename = "./test.vcf";
    std::ofstream vcf_out;
    vcf_out.open(vcf_filename);
    vcf_out << "##fileformat=VCFv4.2" << std::endl;
    vcf_out << "##FILTER=<ID=PASS,Description=\"All filters passed\">" << std::endl;
    vcf_out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
    vcf_out << "##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the snarl tree (0=top level)\">" << std::endl;
    vcf_out << "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph\">" << std::endl;
    vcf_out << "##contig=<ID=ref,length=100>" << std::endl;
    vcf_out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5" << std::endl;
    // For the outer snarl, the three alleles are ins, ins, del, so only genotypes 0 or 1 have the ins
    vcf_out << "ref\t1\t>1>8\tACGTA\tACCTA,A\t60\t.\tLV=0;AT=>1>2>3>4>6>7>8,>1>2>7>8,>1>8\tGT\t2/0\t1/1\t1/1\t0/2\t1/0" << std::endl;
    // For the middle snarl, the two alleles are ins and del, so the 0's have the insertion
    vcf_out << "ref\t2\t>2>7\tCGT\tCCT,A\t60\t.\tLV=1;AT=>2>3>4>6>7,>2>7\tGT\t0/0\t0/1\t1/1\t0/0\t1/0" << std::endl;
    vcf_out << "ref\t3\t>3>6\tCGT\tCCT,A\t60\t.\tLV=1;AT=>3>4>6,>3>5>6\tGT\t0/1\t1/0\t1/0\t0/1\t1/0" << std::endl;
    vcf_out.close();



    SECTION("Make a VCFParser") {
        TestVCFParser parser (true);
        std::vector<std::string> sample_names = parser.initialize_parser (vcf_filename);

        parser.for_each_record_on_chromosome("ref", [&](const auto& x ) {});

        REQUIRE(parser.hap_count == 10);
        REQUIRE(parser.snarl_count == 2);

        REQUIRE(parser.snarl_in_to_out.size() == 4);
        REQUIRE(parser.genotypes.size() == 20);



        REQUIRE(parser.get_opposite_snarl_bound(stoat::node_traversal_t(2, false)) == stoat::node_traversal_t(7, false));
        REQUIRE(parser.get_opposite_snarl_bound(stoat::node_traversal_t(7, true)) == stoat::node_traversal_t(2, true));
        REQUIRE(parser.get_opposite_snarl_bound(stoat::node_traversal_t(3, false)) == stoat::node_traversal_t(6, false));
        REQUIRE(parser.get_opposite_snarl_bound(stoat::node_traversal_t(6, true)) == stoat::node_traversal_t(3, true));

        REQUIRE(!parser.does_sample_have_snarl(0, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(1, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(2, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(3, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(4, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(5, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(6, stoat::node_traversal_t(2, false)));
        REQUIRE(!parser.does_sample_have_snarl(7, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(8, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(9, stoat::node_traversal_t(2, false)));


        REQUIRE(!parser.does_sample_have_snarl(0, stoat::node_traversal_t(7, true)));
        REQUIRE(parser.does_sample_have_snarl(1, stoat::node_traversal_t(7, true)));
        REQUIRE(parser.does_sample_have_snarl(2, stoat::node_traversal_t(7, true)));
        REQUIRE(parser.does_sample_have_snarl(3, stoat::node_traversal_t(7, true)));
        REQUIRE(parser.does_sample_have_snarl(4, stoat::node_traversal_t(7, true)));
        REQUIRE(parser.does_sample_have_snarl(5, stoat::node_traversal_t(7, true)));
        REQUIRE(parser.does_sample_have_snarl(6, stoat::node_traversal_t(7, true)));
        REQUIRE(!parser.does_sample_have_snarl(7, stoat::node_traversal_t(7, true)));
        REQUIRE(parser.does_sample_have_snarl(8, stoat::node_traversal_t(7, true)));
        REQUIRE(parser.does_sample_have_snarl(9, stoat::node_traversal_t(7, true)));

        REQUIRE(parser.does_sample_have_snarl(0, stoat::node_traversal_t(3, false)));
        REQUIRE(parser.does_sample_have_snarl(1, stoat::node_traversal_t(3, false)));
        REQUIRE(parser.does_sample_have_snarl(2, stoat::node_traversal_t(3, false)));
        REQUIRE(!parser.does_sample_have_snarl(3, stoat::node_traversal_t(3, false)));
        REQUIRE(!parser.does_sample_have_snarl(4, stoat::node_traversal_t(3, false)));
        REQUIRE(!parser.does_sample_have_snarl(5, stoat::node_traversal_t(3, false)));
        REQUIRE(parser.does_sample_have_snarl(6, stoat::node_traversal_t(3, false)));
        REQUIRE(parser.does_sample_have_snarl(7, stoat::node_traversal_t(3, false)));
        REQUIRE(!parser.does_sample_have_snarl(8, stoat::node_traversal_t(3, false)));
        REQUIRE(parser.does_sample_have_snarl(9, stoat::node_traversal_t(3, false)));

        REQUIRE(parser.does_sample_have_snarl(0, stoat::node_traversal_t(6, true)));
        REQUIRE(parser.does_sample_have_snarl(1, stoat::node_traversal_t(6, true)));
        REQUIRE(parser.does_sample_have_snarl(2, stoat::node_traversal_t(6, true)));
        REQUIRE(!parser.does_sample_have_snarl(3, stoat::node_traversal_t(6, true)));
        REQUIRE(!parser.does_sample_have_snarl(4, stoat::node_traversal_t(6, true)));
        REQUIRE(!parser.does_sample_have_snarl(5, stoat::node_traversal_t(6, true)));
        REQUIRE(parser.does_sample_have_snarl(6, stoat::node_traversal_t(6, true)));
        REQUIRE(parser.does_sample_have_snarl(7, stoat::node_traversal_t(6, true)));
        REQUIRE(!parser.does_sample_have_snarl(8, stoat::node_traversal_t(6, true)));
        REQUIRE(parser.does_sample_have_snarl(9, stoat::node_traversal_t(6, true)));

        parser.close_vcf();

    }
    SECTION("Go through the contents using the untangler") {
        TestVCFParser parser (true);
        std::vector<std::string> sample_names = parser.initialize_parser (vcf_filename);

        // Check first chr
        std::string chr = parser.get_next_chromosome_name();
        REQUIRE(chr == ("ref"));
        size_t snarl_num = 0;
        parser.for_each_record_on_chromosome(chr, [&] (const vcf_info_t& vcf_info) {
            if (snarl_num == 0) {
                // For the outer snarl, should be
                // >1>2>3>4>6>7>8  >1>2>7>8  >1>8 
                // With 2-5 being nested snarl
                //2/0  1/1  1/1  0/2  1/0
                REQUIRE(vcf_info.lv == 0);
                REQUIRE(vcf_info.paths.size() == 3); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">1>2>0>7>8");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">1>2>0>7>8");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[2]) == ">1>8");
                REQUIRE(vcf_info.genotype[0] == 2); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 1); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 2); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 1) {
                // For the middle snarl, should be
                //>2>3>4>6>7  >2>7
                // 0/0  0/1  1/1  0/0  1/0
                // .                .
                // With the . meaning that they weren't present according to the outer snarl 
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">2>3>0>6>7");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">2>7");
                REQUIRE(vcf_info.genotype[0] == -1); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == -1); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 2) {
                // For the inner snarl, should be
                //>3>4>6  >3>5>6
                // 0/1  1/0  1/0  0/1  1/0
                //        .  . .       .   
                // With the . meaning that they weren't present according to the outer snarl 
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">3>4>6");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">3>5>6");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 1); 
                REQUIRE(vcf_info.genotype[2] == 1); 
                REQUIRE(vcf_info.genotype[3] == -1); 
                REQUIRE(vcf_info.genotype[4] == -1); 
                REQUIRE(vcf_info.genotype[5] == -1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 1); 
                REQUIRE(vcf_info.genotype[8] == -1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            }
            ++snarl_num;
        });
        REQUIRE(snarl_num == 3);



        parser.close_vcf();

    }

    // clean up

    std::string rm_cmd = "rm " + vcf_filename;
    int rm = system(rm_cmd.c_str());

}
TEST_CASE( "Multiple records from a deletion", "[vcf_parser]" ) {
    /*
          -------------
        /   --------   \
       /  /    4     \   \
      /  /   /   \    \   \
     1--2 - 3      6 - 7---8 
             \   /
               5
    */

    // Write the vcf
    // I didn't actually check that this was a reasonable vcf but the snarls and paths should be right
    std::string vcf_filename = "./test.vcf";
    std::ofstream vcf_out;
    vcf_out.open(vcf_filename);
    vcf_out << "##fileformat=VCFv4.2" << std::endl;
    vcf_out << "##FILTER=<ID=PASS,Description=\"All filters passed\">" << std::endl;
    vcf_out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
    vcf_out << "##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the snarl tree (0=top level)\">" << std::endl;
    vcf_out << "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph\">" << std::endl;
    vcf_out << "##contig=<ID=ref,length=100>" << std::endl;
    vcf_out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5" << std::endl;
    // For the outer snarl, the three alleles are ins, ins, del, so only genotypes 0 or 1 have the ins
    vcf_out << "ref\t1\t>1>8\tACGTA\tACCTA,A\t60\t.\tLV=0;AT=>1>2>3>4>6>7>8,>1>2>7>8\tGT\t./.\t1/1\t1/1\t0/.\t1/." << std::endl;
    vcf_out << "ref\t1\t>1>8\tACGTA\tACCTA,A\t60\t.\tLV=0;AT=>1>2>3>4>6>7>8,>1>8\tGT\t1/0\t./.\t./.\t./1\t./0" << std::endl;
    // For the middle snarl, the two alleles are ins and del, so the 0's have the insertion
    vcf_out << "ref\t2\t>2>7\tCGT\tCCT,A\t60\t.\tLV=1;AT=>2>3>4>6>7,>2>7\tGT\t0/0\t0/1\t1/1\t0/0\t1/0" << std::endl;
    vcf_out << "ref\t3\t>3>6\tCGT\tCCT,A\t60\t.\tLV=1;AT=>3>4>6,>3>5>6\tGT\t0/1\t1/0\t1/0\t0/1\t1/0" << std::endl;
    vcf_out.close();



    SECTION("Check the contents of the utntangler") {
        TestVCFParser parser (true);
        std::vector<std::string> sample_names = parser.initialize_parser (vcf_filename);

        parser.for_each_record_on_chromosome("ref", [&](const auto& x ) {});

        REQUIRE(parser.hap_count == 10);
        REQUIRE(parser.snarl_count == 2);

        REQUIRE(parser.snarl_in_to_out.size() == 4);
        REQUIRE(parser.genotypes.size() == 20);



        REQUIRE(parser.get_opposite_snarl_bound(stoat::node_traversal_t(2, false)) == stoat::node_traversal_t(7, false));
        REQUIRE(parser.get_opposite_snarl_bound(stoat::node_traversal_t(7, true)) == stoat::node_traversal_t(2, true));
        REQUIRE(parser.get_opposite_snarl_bound(stoat::node_traversal_t(3, false)) == stoat::node_traversal_t(6, false));
        REQUIRE(parser.get_opposite_snarl_bound(stoat::node_traversal_t(6, true)) == stoat::node_traversal_t(3, true));

        REQUIRE(!parser.does_sample_have_snarl(0, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(1, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(2, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(3, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(4, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(5, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(6, stoat::node_traversal_t(2, false)));
        REQUIRE(!parser.does_sample_have_snarl(7, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(8, stoat::node_traversal_t(2, false)));
        REQUIRE(parser.does_sample_have_snarl(9, stoat::node_traversal_t(2, false)));


        REQUIRE(!parser.does_sample_have_snarl(0, stoat::node_traversal_t(7, true)));
        REQUIRE(parser.does_sample_have_snarl(1, stoat::node_traversal_t(7, true)));
        REQUIRE(parser.does_sample_have_snarl(2, stoat::node_traversal_t(7, true)));
        REQUIRE(parser.does_sample_have_snarl(3, stoat::node_traversal_t(7, true)));
        REQUIRE(parser.does_sample_have_snarl(4, stoat::node_traversal_t(7, true)));
        REQUIRE(parser.does_sample_have_snarl(5, stoat::node_traversal_t(7, true)));
        REQUIRE(parser.does_sample_have_snarl(6, stoat::node_traversal_t(7, true)));
        REQUIRE(!parser.does_sample_have_snarl(7, stoat::node_traversal_t(7, true)));
        REQUIRE(parser.does_sample_have_snarl(8, stoat::node_traversal_t(7, true)));
        REQUIRE(parser.does_sample_have_snarl(9, stoat::node_traversal_t(7, true)));

        REQUIRE(parser.does_sample_have_snarl(0, stoat::node_traversal_t(3, false)));
        REQUIRE(parser.does_sample_have_snarl(1, stoat::node_traversal_t(3, false)));
        REQUIRE(parser.does_sample_have_snarl(2, stoat::node_traversal_t(3, false)));
        REQUIRE(!parser.does_sample_have_snarl(3, stoat::node_traversal_t(3, false)));
        REQUIRE(!parser.does_sample_have_snarl(4, stoat::node_traversal_t(3, false)));
        REQUIRE(!parser.does_sample_have_snarl(5, stoat::node_traversal_t(3, false)));
        REQUIRE(parser.does_sample_have_snarl(6, stoat::node_traversal_t(3, false)));
        REQUIRE(parser.does_sample_have_snarl(7, stoat::node_traversal_t(3, false)));
        REQUIRE(!parser.does_sample_have_snarl(8, stoat::node_traversal_t(3, false)));
        REQUIRE(parser.does_sample_have_snarl(9, stoat::node_traversal_t(3, false)));

        REQUIRE(parser.does_sample_have_snarl(0, stoat::node_traversal_t(6, true)));
        REQUIRE(parser.does_sample_have_snarl(1, stoat::node_traversal_t(6, true)));
        REQUIRE(parser.does_sample_have_snarl(2, stoat::node_traversal_t(6, true)));
        REQUIRE(!parser.does_sample_have_snarl(3, stoat::node_traversal_t(6, true)));
        REQUIRE(!parser.does_sample_have_snarl(4, stoat::node_traversal_t(6, true)));
        REQUIRE(!parser.does_sample_have_snarl(5, stoat::node_traversal_t(6, true)));
        REQUIRE(parser.does_sample_have_snarl(6, stoat::node_traversal_t(6, true)));
        REQUIRE(parser.does_sample_have_snarl(7, stoat::node_traversal_t(6, true)));
        REQUIRE(!parser.does_sample_have_snarl(8, stoat::node_traversal_t(6, true)));
        REQUIRE(parser.does_sample_have_snarl(9, stoat::node_traversal_t(6, true)));

        parser.close_vcf();

    }
    SECTION("Go through the contents using the untangler") {
        TestVCFParser parser (true);
        std::vector<std::string> sample_names = parser.initialize_parser (vcf_filename);

        // Check first chr
        std::string chr = parser.get_next_chromosome_name();
        REQUIRE(chr == ("ref"));
        size_t snarl_num = 0;
        parser.for_each_record_on_chromosome(chr, [&] (const vcf_info_t& vcf_info) {
            // The first and second records are the same snarl
            if (snarl_num == 0) {
                // For the outer snarl, should be
                // >1>2>3>4>6>7>8  >1>2>7>8 
                // With 2-5 being nested snarl
                //./.  1/1  1/1  0/.  1/.
                REQUIRE(vcf_info.lv == 0);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">1>2>0>7>8");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">1>2>0>7>8");
                REQUIRE(vcf_info.genotype[0] == -1); 
                REQUIRE(vcf_info.genotype[1] == -1); 
                REQUIRE(vcf_info.genotype[2] == 1); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == -1); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == -1); 
            } else if (snarl_num == 1) {
                // For the outer snarl, should be
                // >1>2>3>4>6>7>8  >1>8 
                // With 2-5 being nested snarl
                //1/0  ./.  ./.  ./1  ./0
                REQUIRE(vcf_info.lv == 0);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">1>2>0>7>8");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">1>8");
                REQUIRE(vcf_info.genotype[0] == 1); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == -1); 
                REQUIRE(vcf_info.genotype[3] == -1); 
                REQUIRE(vcf_info.genotype[4] == -1); 
                REQUIRE(vcf_info.genotype[5] == -1); 
                REQUIRE(vcf_info.genotype[6] == -1); 
                REQUIRE(vcf_info.genotype[7] == 1); 
                REQUIRE(vcf_info.genotype[8] == -1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 2) {
                // For the middle snarl, should be
                //>2>3>4>6>7  >2>7
                // 0/0  0/1  1/1  0/0  1/0
                // .                .
                // With the . meaning that they weren't present according to the outer snarl 
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">2>3>0>6>7");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">2>7");
                REQUIRE(vcf_info.genotype[0] == -1); 
                REQUIRE(vcf_info.genotype[1] == 0); 
                REQUIRE(vcf_info.genotype[2] == 0); 
                REQUIRE(vcf_info.genotype[3] == 1); 
                REQUIRE(vcf_info.genotype[4] == 1); 
                REQUIRE(vcf_info.genotype[5] == 1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == -1); 
                REQUIRE(vcf_info.genotype[8] == 1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            } else if (snarl_num == 3) {
                // For the inner snarl, should be
                //>3>4>6  >3>5>6
                // 0/1  1/0  1/0  0/1  1/0
                //        .  . .       .   
                // With the . meaning that they weren't present according to the outer snarl 
                REQUIRE(vcf_info.lv == 1);
                REQUIRE(vcf_info.paths.size() == 2); 
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[0]) == ">3>4>6");
                REQUIRE(path_node_traversal_to_string(vcf_info.paths[1]) == ">3>5>6");
                REQUIRE(vcf_info.genotype[0] == 0); 
                REQUIRE(vcf_info.genotype[1] == 1); 
                REQUIRE(vcf_info.genotype[2] == 1); 
                REQUIRE(vcf_info.genotype[3] == -1); 
                REQUIRE(vcf_info.genotype[4] == -1); 
                REQUIRE(vcf_info.genotype[5] == -1); 
                REQUIRE(vcf_info.genotype[6] == 0); 
                REQUIRE(vcf_info.genotype[7] == 1); 
                REQUIRE(vcf_info.genotype[8] == -1); 
                REQUIRE(vcf_info.genotype[9] == 0); 
            }
            ++snarl_num;
        });
        REQUIRE(snarl_num == 4);



        parser.close_vcf();

    }

    // clean up

    std::string rm_cmd = "rm " + vcf_filename;
    int rm = system(rm_cmd.c_str());

}
