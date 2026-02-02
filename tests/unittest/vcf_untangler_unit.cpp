#include <catch.hpp>
#include "../../src/vcf_untangler.hpp"

using namespace stoat_vcf;


class TestVCFUntangler : VCFUntangler {
    public: 
    TestVCFUntangler(const std::string& vcf_filename) :
        VCFUntangler(vcf_filename) {} 
    using VCFUntangler::does_sample_have_snarl;
    using VCFUntangler::get_opposite_snarl_bound;
    using VCFUntangler::fill_in_snarls_for_chromosome;
    using VCFUntangler::genotypes;
    using VCFUntangler::snarl_in_to_out;
    using VCFUntangler::sample_hap_count;
    using VCFUntangler::snarl_count;
};

TEST_CASE( "Simple nested snarl", "[vcf_untangler]" ) {
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



    SECTION("Make a VCFUntangler") {
        TestVCFUntangler untangler (vcf_filename);

        untangler.fill_in_snarls_for_chromosome("ref");

        REQUIRE(untangler.sample_hap_count == 10);
        REQUIRE(untangler.snarl_count == 1);

        REQUIRE(untangler.snarl_in_to_out.size() == 2);
        REQUIRE(untangler.genotypes.size() == 10);

        REQUIRE(untangler.get_opposite_snarl_bound(stoat::node_traversal_t(2, false)) == stoat::node_traversal_t(5, false));
        REQUIRE(untangler.get_opposite_snarl_bound(stoat::node_traversal_t(5, true)) == stoat::node_traversal_t(2, true));

        REQUIRE(!untangler.does_sample_have_snarl(0, stoat::node_traversal_t(2, false)));
        REQUIRE(untangler.does_sample_have_snarl(1, stoat::node_traversal_t(2, false)));
        REQUIRE(untangler.does_sample_have_snarl(2, stoat::node_traversal_t(2, false)));
        REQUIRE(untangler.does_sample_have_snarl(3, stoat::node_traversal_t(2, false)));
        REQUIRE(untangler.does_sample_have_snarl(4, stoat::node_traversal_t(2, false)));
        REQUIRE(untangler.does_sample_have_snarl(5, stoat::node_traversal_t(2, false)));
        REQUIRE(untangler.does_sample_have_snarl(6, stoat::node_traversal_t(2, false)));
        REQUIRE(!untangler.does_sample_have_snarl(7, stoat::node_traversal_t(2, false)));
        REQUIRE(untangler.does_sample_have_snarl(8, stoat::node_traversal_t(2, false)));
        REQUIRE(untangler.does_sample_have_snarl(9, stoat::node_traversal_t(2, false)));

        REQUIRE(!untangler.does_sample_have_snarl(0, stoat::node_traversal_t(5, true)));
        REQUIRE(untangler.does_sample_have_snarl(1, stoat::node_traversal_t(5, true)));
        REQUIRE(untangler.does_sample_have_snarl(2, stoat::node_traversal_t(5, true)));
        REQUIRE(untangler.does_sample_have_snarl(3, stoat::node_traversal_t(5, true)));
        REQUIRE(untangler.does_sample_have_snarl(4, stoat::node_traversal_t(5, true)));
        REQUIRE(untangler.does_sample_have_snarl(5, stoat::node_traversal_t(5, true)));
        REQUIRE(untangler.does_sample_have_snarl(6, stoat::node_traversal_t(5, true)));
        REQUIRE(!untangler.does_sample_have_snarl(7, stoat::node_traversal_t(5, true)));
        REQUIRE(untangler.does_sample_have_snarl(8, stoat::node_traversal_t(5, true)));
        REQUIRE(untangler.does_sample_have_snarl(9, stoat::node_traversal_t(5, true)));

    }

    // clean up

    std::string rm_cmd = "rm " + vcf_filename;
    int rm = system(rm_cmd.c_str());

}

TEST_CASE( "Simple nested snarl multiple snps", "[vcf_untangler]" ) {
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
    // For the outer snarl, the three alleles are ins, ins, del, so only genotypes 0 or 1 have the ins
    vcf_out << "ref\t1\t>1>9\tACGTA\tACCTA,A\t60\t.\tLV=0;AT=>1>2>3>5>6>8>9,>1>2>4>5>7>8>9,>1>2>4>5>8>9,>1>6\tGT\t3/0\t2/1\t1/1\t0/3\t1/0" << std::endl;
    vcf_out << "ref\t2\t>2>5\tCGT\tCCT,A\t60\t.\tLV=1;AT=>2>3>5,>2>4>5\tGT\t0/0\t0/1\t1/1\t0/0\t1/0" << std::endl;
    vcf_out << "ref\t2\t>5>8\tCGT\tCCT,A\t60\t.\tLV=1;AT=>5>6>8,>5>8\tGT\t0/0\t0/1\t1/1\t0/0\t1/0" << std::endl;
    vcf_out.close();



    SECTION("Make a VCFUntangler") {
        TestVCFUntangler untangler (vcf_filename);

        untangler.fill_in_snarls_for_chromosome("ref");

        REQUIRE(untangler.sample_hap_count == 10);
        REQUIRE(untangler.snarl_count == 2);

        REQUIRE(untangler.snarl_in_to_out.size() == 4);
        REQUIRE(untangler.genotypes.size() == 20);



        REQUIRE(untangler.get_opposite_snarl_bound(stoat::node_traversal_t(2, false)) == stoat::node_traversal_t(5, false));
        REQUIRE(untangler.get_opposite_snarl_bound(stoat::node_traversal_t(5, true)) == stoat::node_traversal_t(2, true));
        REQUIRE(untangler.get_opposite_snarl_bound(stoat::node_traversal_t(5, false)) == stoat::node_traversal_t(8, false));
        REQUIRE(untangler.get_opposite_snarl_bound(stoat::node_traversal_t(8, true)) == stoat::node_traversal_t(5, true));

        REQUIRE(!untangler.does_sample_have_snarl(0, stoat::node_traversal_t(2, false)));
        REQUIRE(untangler.does_sample_have_snarl(1, stoat::node_traversal_t(2, false)));
        REQUIRE(untangler.does_sample_have_snarl(2, stoat::node_traversal_t(2, false)));
        REQUIRE(untangler.does_sample_have_snarl(3, stoat::node_traversal_t(2, false)));
        REQUIRE(untangler.does_sample_have_snarl(4, stoat::node_traversal_t(2, false)));
        REQUIRE(untangler.does_sample_have_snarl(5, stoat::node_traversal_t(2, false)));
        REQUIRE(untangler.does_sample_have_snarl(6, stoat::node_traversal_t(2, false)));
        REQUIRE(!untangler.does_sample_have_snarl(7, stoat::node_traversal_t(2, false)));
        REQUIRE(untangler.does_sample_have_snarl(8, stoat::node_traversal_t(2, false)));
        REQUIRE(untangler.does_sample_have_snarl(9, stoat::node_traversal_t(2, false)));


        REQUIRE(!untangler.does_sample_have_snarl(0, stoat::node_traversal_t(5, true)));
        REQUIRE(untangler.does_sample_have_snarl(1, stoat::node_traversal_t(5, true)));
        REQUIRE(untangler.does_sample_have_snarl(2, stoat::node_traversal_t(5, true)));
        REQUIRE(untangler.does_sample_have_snarl(3, stoat::node_traversal_t(5, true)));
        REQUIRE(untangler.does_sample_have_snarl(4, stoat::node_traversal_t(5, true)));
        REQUIRE(untangler.does_sample_have_snarl(5, stoat::node_traversal_t(5, true)));
        REQUIRE(untangler.does_sample_have_snarl(6, stoat::node_traversal_t(5, true)));
        REQUIRE(!untangler.does_sample_have_snarl(7, stoat::node_traversal_t(5, true)));
        REQUIRE(untangler.does_sample_have_snarl(8, stoat::node_traversal_t(5, true)));
        REQUIRE(untangler.does_sample_have_snarl(9, stoat::node_traversal_t(5, true)));

        REQUIRE(!untangler.does_sample_have_snarl(0, stoat::node_traversal_t(5, false)));
        REQUIRE(untangler.does_sample_have_snarl(1, stoat::node_traversal_t(5, false)));
        REQUIRE(untangler.does_sample_have_snarl(2, stoat::node_traversal_t(5, false)));
        REQUIRE(untangler.does_sample_have_snarl(3, stoat::node_traversal_t(5, false)));
        REQUIRE(untangler.does_sample_have_snarl(4, stoat::node_traversal_t(5, false)));
        REQUIRE(untangler.does_sample_have_snarl(5, stoat::node_traversal_t(5, false)));
        REQUIRE(untangler.does_sample_have_snarl(6, stoat::node_traversal_t(5, false)));
        REQUIRE(!untangler.does_sample_have_snarl(7, stoat::node_traversal_t(5, false)));
        REQUIRE(untangler.does_sample_have_snarl(8, stoat::node_traversal_t(5, false)));
        REQUIRE(untangler.does_sample_have_snarl(9, stoat::node_traversal_t(5, false)));

        REQUIRE(!untangler.does_sample_have_snarl(0, stoat::node_traversal_t(8, true)));
        REQUIRE(untangler.does_sample_have_snarl(1, stoat::node_traversal_t(8, true)));
        REQUIRE(untangler.does_sample_have_snarl(2, stoat::node_traversal_t(8, true)));
        REQUIRE(untangler.does_sample_have_snarl(3, stoat::node_traversal_t(8, true)));
        REQUIRE(untangler.does_sample_have_snarl(4, stoat::node_traversal_t(8, true)));
        REQUIRE(untangler.does_sample_have_snarl(5, stoat::node_traversal_t(8, true)));
        REQUIRE(untangler.does_sample_have_snarl(6, stoat::node_traversal_t(8, true)));
        REQUIRE(!untangler.does_sample_have_snarl(7, stoat::node_traversal_t(8, true)));
        REQUIRE(untangler.does_sample_have_snarl(8, stoat::node_traversal_t(8, true)));
        REQUIRE(untangler.does_sample_have_snarl(9, stoat::node_traversal_t(8, true)));

    }

    // clean up

    std::string rm_cmd = "rm " + vcf_filename;
    int rm = system(rm_cmd.c_str());

}
TEST_CASE( "Three nested snarl multiple snps", "[vcf_untangler]" ) {
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
    vcf_out << "ref\t3\t>3>6\tCGT\tCCT,A\t60\t.\tLV=1;AT=>3>4>6,>3>5>6\tGT\t0/0\t0/1\t1/1\t0/0\t1/0" << std::endl;
    vcf_out.close();



    SECTION("Make a VCFUntangler") {
        TestVCFUntangler untangler (vcf_filename);

        untangler.fill_in_snarls_for_chromosome("ref");

        REQUIRE(untangler.sample_hap_count == 10);
        REQUIRE(untangler.snarl_count == 2);

        REQUIRE(untangler.snarl_in_to_out.size() == 4);
        REQUIRE(untangler.genotypes.size() == 20);



        REQUIRE(untangler.get_opposite_snarl_bound(stoat::node_traversal_t(2, false)) == stoat::node_traversal_t(7, false));
        REQUIRE(untangler.get_opposite_snarl_bound(stoat::node_traversal_t(7, true)) == stoat::node_traversal_t(2, true));
        REQUIRE(untangler.get_opposite_snarl_bound(stoat::node_traversal_t(3, false)) == stoat::node_traversal_t(6, false));
        REQUIRE(untangler.get_opposite_snarl_bound(stoat::node_traversal_t(6, true)) == stoat::node_traversal_t(3, true));

        REQUIRE(!untangler.does_sample_have_snarl(0, stoat::node_traversal_t(2, false)));
        REQUIRE(untangler.does_sample_have_snarl(1, stoat::node_traversal_t(2, false)));
        REQUIRE(untangler.does_sample_have_snarl(2, stoat::node_traversal_t(2, false)));
        REQUIRE(untangler.does_sample_have_snarl(3, stoat::node_traversal_t(2, false)));
        REQUIRE(untangler.does_sample_have_snarl(4, stoat::node_traversal_t(2, false)));
        REQUIRE(untangler.does_sample_have_snarl(5, stoat::node_traversal_t(2, false)));
        REQUIRE(untangler.does_sample_have_snarl(6, stoat::node_traversal_t(2, false)));
        REQUIRE(!untangler.does_sample_have_snarl(7, stoat::node_traversal_t(2, false)));
        REQUIRE(untangler.does_sample_have_snarl(8, stoat::node_traversal_t(2, false)));
        REQUIRE(untangler.does_sample_have_snarl(9, stoat::node_traversal_t(2, false)));


        REQUIRE(!untangler.does_sample_have_snarl(0, stoat::node_traversal_t(7, true)));
        REQUIRE(untangler.does_sample_have_snarl(1, stoat::node_traversal_t(7, true)));
        REQUIRE(untangler.does_sample_have_snarl(2, stoat::node_traversal_t(7, true)));
        REQUIRE(untangler.does_sample_have_snarl(3, stoat::node_traversal_t(7, true)));
        REQUIRE(untangler.does_sample_have_snarl(4, stoat::node_traversal_t(7, true)));
        REQUIRE(untangler.does_sample_have_snarl(5, stoat::node_traversal_t(7, true)));
        REQUIRE(untangler.does_sample_have_snarl(6, stoat::node_traversal_t(7, true)));
        REQUIRE(!untangler.does_sample_have_snarl(7, stoat::node_traversal_t(7, true)));
        REQUIRE(untangler.does_sample_have_snarl(8, stoat::node_traversal_t(7, true)));
        REQUIRE(untangler.does_sample_have_snarl(9, stoat::node_traversal_t(7, true)));

        REQUIRE(untangler.does_sample_have_snarl(0, stoat::node_traversal_t(3, false)));
        REQUIRE(untangler.does_sample_have_snarl(1, stoat::node_traversal_t(3, false)));
        REQUIRE(untangler.does_sample_have_snarl(2, stoat::node_traversal_t(3, false)));
        REQUIRE(!untangler.does_sample_have_snarl(3, stoat::node_traversal_t(3, false)));
        REQUIRE(!untangler.does_sample_have_snarl(4, stoat::node_traversal_t(3, false)));
        REQUIRE(!untangler.does_sample_have_snarl(5, stoat::node_traversal_t(3, false)));
        REQUIRE(untangler.does_sample_have_snarl(6, stoat::node_traversal_t(3, false)));
        REQUIRE(untangler.does_sample_have_snarl(7, stoat::node_traversal_t(3, false)));
        REQUIRE(!untangler.does_sample_have_snarl(8, stoat::node_traversal_t(3, false)));
        REQUIRE(untangler.does_sample_have_snarl(9, stoat::node_traversal_t(3, false)));

        REQUIRE(untangler.does_sample_have_snarl(0, stoat::node_traversal_t(6, true)));
        REQUIRE(untangler.does_sample_have_snarl(1, stoat::node_traversal_t(6, true)));
        REQUIRE(untangler.does_sample_have_snarl(2, stoat::node_traversal_t(6, true)));
        REQUIRE(!untangler.does_sample_have_snarl(3, stoat::node_traversal_t(6, true)));
        REQUIRE(!untangler.does_sample_have_snarl(4, stoat::node_traversal_t(6, true)));
        REQUIRE(!untangler.does_sample_have_snarl(5, stoat::node_traversal_t(6, true)));
        REQUIRE(untangler.does_sample_have_snarl(6, stoat::node_traversal_t(6, true)));
        REQUIRE(untangler.does_sample_have_snarl(7, stoat::node_traversal_t(6, true)));
        REQUIRE(!untangler.does_sample_have_snarl(8, stoat::node_traversal_t(6, true)));
        REQUIRE(untangler.does_sample_have_snarl(9, stoat::node_traversal_t(6, true)));

    }

    // clean up

    std::string rm_cmd = "rm " + vcf_filename;
    int rm = system(rm_cmd.c_str());

}
