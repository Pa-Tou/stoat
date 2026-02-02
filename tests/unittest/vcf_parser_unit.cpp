#include <catch.hpp>
#include "../../src/vcf_parser.hpp"

using namespace stoat_vcf;


class TestVCFParser : VCFParser {
    public: 
    TestVCFParser() :
        VCFParser() {} 
    using VCFParser::initialize_parser;
    using VCFParser::get_next_chromosome_name;
    using VCFParser::for_each_record_on_chromosome;
    using VCFParser::close_vcf;
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
        TestVCFParser parser;
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



    SECTION("Make a VCFParser") {
        TestVCFParser parser;
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

    // clean up

    std::string rm_cmd = "rm " + vcf_filename;
    int rm = system(rm_cmd.c_str());

}

