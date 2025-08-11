#include <catch.hpp>
#include "../../src/snarl_data_t.hpp"
#include "../../src/utils.hpp"

using namespace std;
using namespace bdsg;

using handlegraph::step_handle_t;
using handlegraph::handle_t;
using handlegraph::net_handle_t;

void check_snarl(
    const std::unordered_map<std::string, std::vector<stoat::Snarl_data_t>>& chr_snarl_matrix,
    const std::string& chr,
    const std::string& expected_id_str,
    const std::vector<std::string>& expected_paths,
    size_t expected_start,
    size_t expected_end,
    const std::vector<std::string>& expected_types,
    bdsg::SnarlDistanceIndex& stree) {

    const auto& snarls = chr_snarl_matrix.at(chr);
    const stoat::Snarl_data_t* target = nullptr;

    for (const auto& snarl : snarls) {
        std::string id = stoat::pairToString(stoat::find_snarl_id(stree, snarl.snarl));
        if (id == expected_id_str) {
            target = &snarl;
            break;
        }
    }

    REQUIRE(target != nullptr);  // Ensure the snarl was found
    REQUIRE(stoat::pairToString(stoat::find_snarl_id(stree, target->snarl)) == expected_id_str);
    REQUIRE(stoat::stringToVector<std::string>(stoat::vectorPathToString(target->snarl_paths)) == expected_paths);
    REQUIRE(target->start_positions == expected_start);
    REQUIRE(target->end_positions == expected_end);
    REQUIRE(target->type_variants == expected_types);
}

TEST_CASE("Test simulated case", "[Path]") {

    std::unique_ptr<bdsg::SnarlDistanceIndex> stree;
    std::unique_ptr<bdsg::PackedGraph> pg;
    handlegraph::net_handle_t root;
    std::unique_ptr<handlegraph::PathHandleGraph> path_graph;
    std::unique_ptr<bdsg::PackedPositionOverlay> pp_overlay;

    size_t children_threshold = 50;
    size_t path_length_threshold = 10000;
    std::unordered_set<std::string> ref_chr = {"ref"};
    std::string output_dir = "../tests/graph_test";
    std::string output_snarl_not_analyse = output_dir + "/snarl_not_analyse.tsv";
    std::string output_file = output_dir + "/snarl_analyse.tsv";

    SECTION("simple_snp") {
        std::string pg_path = "../tests/graph_test/simple_snp.pg";
        std::string dist_path = "../tests/graph_test/simple_snp.dist";

        std::tie(stree, pg, root, path_graph, pp_overlay) = stoat::parse_graph_tree(pg_path, dist_path);
        auto snarls = stoat::save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = stoat::loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0);

        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);

        REQUIRE(stoat::pairToString(stoat::find_snarl_id(*stree, snarls_chr["ref"][0].snarl)) == "5_2");
        REQUIRE(stoat::stringToVector<std::string>(stoat::vectorPathToString(snarls_chr["ref"][0].snarl_paths)) == std::vector<std::string>{">2>3>5",">2>4>5"});
        REQUIRE(snarls_chr["ref"][0].start_positions == 8);
        REQUIRE(snarls_chr["ref"][0].end_positions == 10);
        REQUIRE(snarls_chr["ref"][0].type_variants == std::vector<std::string>{"1", "1"});
    }

    SECTION("3th SNP") {
        std::string pg_path = "../tests/graph_test/3th_snp.pg";
        std::string dist_path = "../tests/graph_test/3th_snp.dist";

        std::tie(stree, pg, root, path_graph, pp_overlay) = stoat::parse_graph_tree(pg_path, dist_path);
        auto snarls = stoat::save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = stoat::loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0);

        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);

        REQUIRE(stoat::pairToString(stoat::find_snarl_id(*stree, snarls_chr["ref"][0].snarl)) == "6_2");
        REQUIRE(stoat::stringToVector<std::string>(stoat::vectorPathToString(snarls_chr["ref"][0].snarl_paths)) == std::vector<std::string>{">2>3>6",">2>4>6",">2>5>6"});
        REQUIRE(snarls_chr["ref"][0].start_positions == 8);
        REQUIRE(snarls_chr["ref"][0].end_positions == 10);
        REQUIRE(snarls_chr["ref"][0].type_variants == std::vector<std::string>{"1", "1", "1"});
    }

    SECTION("4th") {
        std::string pg_path = "../tests/graph_test/4th.pg";
        std::string dist_path = "../tests/graph_test/4th.dist";

        std::tie(stree, pg, root, path_graph, pp_overlay) = stoat::parse_graph_tree(pg_path, dist_path);
        auto snarls = stoat::save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = stoat::loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0);

        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);

        REQUIRE(stoat::pairToString(stoat::find_snarl_id(*stree, snarls_chr["ref"][0].snarl)) == "2_7");
        REQUIRE(stoat::stringToVector<std::string>(stoat::vectorPathToString(snarls_chr["ref"][0].snarl_paths)) == std::vector<std::string>{">2>3>5>7", ">2>4>6>7" ,">2>3>6>7"});
        REQUIRE(snarls_chr["ref"][0].start_positions == 8);
        REQUIRE(snarls_chr["ref"][0].end_positions == 13);
        REQUIRE(snarls_chr["ref"][0].type_variants == std::vector<std::string>{"4", "6", "5"});
    }

    SECTION("deletion_snp") {
        std::string pg_path = "../tests/graph_test/deletion_snp.pg";
        std::string dist_path = "../tests/graph_test/deletion_snp.dist";

        std::tie(stree, pg, root, path_graph, pp_overlay) = stoat::parse_graph_tree(pg_path, dist_path);
        auto snarls = stoat::save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = stoat::loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0);

        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);

        REQUIRE(stoat::pairToString(stoat::find_snarl_id(*stree, snarls_chr["ref"][0].snarl)) == "4_2");
        REQUIRE(stoat::stringToVector<std::string>(stoat::vectorPathToString(snarls_chr["ref"][0].snarl_paths)) == std::vector<std::string>{">2>4", ">2>3>4"});
        REQUIRE(snarls_chr["ref"][0].start_positions == 8);
        REQUIRE(snarls_chr["ref"][0].end_positions == 10);
        REQUIRE(snarls_chr["ref"][0].type_variants == std::vector<std::string>{"0", "1"});
    }

    SECTION("insert_deletion") {
        std::string pg_path = "../tests/graph_test/insert_deletion.pg";
        std::string dist_path = "../tests/graph_test/insert_deletion.dist";

        std::tie(stree, pg, root, path_graph, pp_overlay) = stoat::parse_graph_tree(pg_path, dist_path);
        auto snarls = stoat::save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = stoat::loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0);
        
        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);

        REQUIRE(stoat::pairToString(stoat::find_snarl_id(*stree, snarls_chr["ref"][0].snarl)) == "4_2");
        REQUIRE(stoat::stringToVector<std::string>(stoat::vectorPathToString(snarls_chr["ref"][0].snarl_paths)) == std::vector<std::string>{">2>4", ">2>3>4"});
        REQUIRE(snarls_chr["ref"][0].start_positions == 8);
        REQUIRE(snarls_chr["ref"][0].end_positions == 12);
        REQUIRE(snarls_chr["ref"][0].type_variants == std::vector<std::string>{"0", "3"});
    }

    SECTION("insert_snp") {
        std::string pg_path = "../tests/graph_test/insert_snp.pg";
        std::string dist_path = "../tests/graph_test/insert_snp.dist";

        std::tie(stree, pg, root, path_graph, pp_overlay) = stoat::parse_graph_tree(pg_path, dist_path);
        auto snarls = stoat::save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = stoat::loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0);
        
        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);

        REQUIRE(stoat::pairToString(stoat::find_snarl_id(*stree, snarls_chr["ref"][0].snarl)) == "5_2");
        REQUIRE(stoat::stringToVector<std::string>(stoat::vectorPathToString(snarls_chr["ref"][0].snarl_paths)) == std::vector<std::string>{">2>3>5", ">2>4>5"});
        REQUIRE(snarls_chr["ref"][0].start_positions == 8);
        REQUIRE(snarls_chr["ref"][0].end_positions == 10);
        REQUIRE(snarls_chr["ref"][0].type_variants == std::vector<std::string>{"1", "3"});
    }

    SECTION("inversion") {
        std::string pg_path = "../tests/graph_test/inversion.pg";
        std::string dist_path = "../tests/graph_test/inversion.dist";

        std::tie(stree, pg, root, path_graph, pp_overlay) = stoat::parse_graph_tree(pg_path, dist_path);
        auto snarls = stoat::save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = stoat::loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0);
        
        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 2);

        check_snarl(snarls_chr, "ref", "2_6", {">2>6", ">2>3>0>5>6"}, 8, 15, {"0", "6"}, *stree);
        check_snarl(snarls_chr, "ref", "3_5", {">3>4>5", ">3<4>5"}, 9, 12, {"2", "2"}, *stree);
    }

    SECTION("large_del") {
        std::string pg_path = "../tests/graph_test/large_del.pg";
        std::string dist_path = "../tests/graph_test/large_del.dist";

        std::tie(stree, pg, root, path_graph, pp_overlay) = stoat::parse_graph_tree(pg_path, dist_path);
        auto snarls = stoat::save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = stoat::loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0);
        
        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 3);

        check_snarl(snarls_chr, "ref", "2_9", {">2>9", ">2>3>0>8>9"}, 8, 9, {"0", "9/10"}, *stree);
        check_snarl(snarls_chr, "ref", "6_8", {">6>8", ">6>7>8"}, 8, 9, {"0", "1"}, *stree);
        check_snarl(snarls_chr, "ref", "3_6", {">3>5>6", ">3>4>6"}, 8, 9, {"1", "1"}, *stree);
    }

    // SECTION("linear") {
    //     std::string pg_path = "../tests/graph_test/linear.pg";
    //     std::string dist_path = "../tests/graph_test/linear.dist";

    //     std::tie(stree, pg, root, path_graph, pp_overlay) = stoat::parse_graph_tree(pg_path, dist_path);
    //     auto snarls = stoat::save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
    //     auto snarls_chr = stoat::loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0);
    //     REQUIRE(snarls_chr.size() == 0);
    // }

    SECTION("loop_simple") {
        std::string pg_path = "../tests/graph_test/loop_simple.pg";
        std::string dist_path = "../tests/graph_test/loop_simple.dist";

        std::tie(stree, pg, root, path_graph, pp_overlay) = stoat::parse_graph_tree(pg_path, dist_path);
        auto snarls = stoat::save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = stoat::loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 1);

        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);

        REQUIRE(stoat::pairToString(stoat::find_snarl_id(*stree, snarls_chr["ref"][0].snarl)) == "5_2");
        REQUIRE(stoat::stringToVector<std::string>(stoat::vectorPathToString(snarls_chr["ref"][0].snarl_paths)) == std::vector<std::string>{">2>3>5", ">2>3>3>5", ">2>4>5"});
        REQUIRE(snarls_chr["ref"][0].start_positions == 8);
        REQUIRE(snarls_chr["ref"][0].end_positions == 10);
        REQUIRE(snarls_chr["ref"][0].type_variants == std::vector<std::string>{"1", "2", "2"});
    }

    SECTION("loop") {
        std::string pg_path = "../tests/graph_test/loop.pg";
        std::string dist_path = "../tests/graph_test/loop.dist";

        std::tie(stree, pg, root, path_graph, pp_overlay) = stoat::parse_graph_tree(pg_path, dist_path);
        auto snarls = stoat::save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = stoat::loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 2);

        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);

        REQUIRE(stoat::pairToString(stoat::find_snarl_id(*stree, snarls_chr["ref"][0].snarl)) == "5_1");
        REQUIRE(stoat::stringToVector<std::string>(stoat::vectorPathToString(snarls_chr["ref"][0].snarl_paths)) == std::vector<std::string>{">1>2>3>5", ">1>2>3>2>3>5", ">1>2>3>2>3>2>3>5", ">1>2>3>2>3>2>4>5", ">1>2>3>2>4>5", ">1>2>4>5"});
        REQUIRE(snarls_chr["ref"][0].start_positions == 4);
        REQUIRE(snarls_chr["ref"][0].end_positions == 10);
        REQUIRE(snarls_chr["ref"][0].type_variants == std::vector<std::string>{"5", "10", "15", "16", "11", "6"});
    }

    SECTION("loop_double") {
        std::string pg_path = "../tests/graph_test/loop_double.pg";
        std::string dist_path = "../tests/graph_test/loop_double.dist";

        std::tie(stree, pg, root, path_graph, pp_overlay) = stoat::parse_graph_tree(pg_path, dist_path);
        auto snarls = stoat::save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = stoat::loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 2);

        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);

        REQUIRE(stoat::pairToString(stoat::find_snarl_id(*stree, snarls_chr["ref"][0].snarl)) == "7_2");
        REQUIRE(stoat::stringToVector<std::string>(stoat::vectorPathToString(snarls_chr["ref"][0].snarl_paths)) == std::vector<std::string>{">2>3>4>5>7", ">2>3>4>5>3>4>5>7", ">2>3>4>5>3>4>5>3>4>5>7", ">2>3>4>5>3>4>3>4>5>7", ">2>3>4>3>4>5>7", ">2>3>4>3>4>5>3>4>5>7", ">2>3>4>3>4>3>4>5>7", ">2>6>7"});
        REQUIRE(snarls_chr["ref"][0].start_positions == 8);
        REQUIRE(snarls_chr["ref"][0].end_positions == 11);
        REQUIRE(snarls_chr["ref"][0].type_variants == std::vector<std::string>{"10", "20", "30", "28", "18", "28", "26", "2"});
    }

    SECTION("loop_plus") {
        std::string pg_path = "../tests/graph_test/loop_plus.pg";
        std::string dist_path = "../tests/graph_test/loop_plus.dist";

        std::tie(stree, pg, root, path_graph, pp_overlay) = stoat::parse_graph_tree(pg_path, dist_path);
        auto snarls = stoat::save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = stoat::loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 1);

        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 2);

        check_snarl(snarls_chr, "ref", "8_2", {">2>3>0>6>8", ">2>3>0>6>3>0>6>8", ">2>7>8"}, 8, 10, {"3/4", "6/8", "1"}, *stree);
        check_snarl(snarls_chr, "ref", "3_6", {">3>5>6", ">3>4>6"}, 8, 10, {"2", "1"}, *stree);
    }

    SECTION("repetition") {
        std::string pg_path = "../tests/graph_test/repetition.pg";
        std::string dist_path = "../tests/graph_test/repetition.dist";

        std::tie(stree, pg, root, path_graph, pp_overlay) = stoat::parse_graph_tree(pg_path, dist_path);
        auto snarls = stoat::save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = stoat::loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0);

        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);

        REQUIRE(stoat::pairToString(stoat::find_snarl_id(*stree, snarls_chr["ref"][0].snarl)) == "6_2");
        REQUIRE(stoat::stringToVector<std::string>(stoat::vectorPathToString(snarls_chr["ref"][0].snarl_paths)) == std::vector<std::string>{">2>6", ">2>3>6", ">2>3>4>6", ">2>3>4>5>6"});
        REQUIRE(snarls_chr["ref"][0].start_positions == 8);
        REQUIRE(snarls_chr["ref"][0].end_positions == 9);
        REQUIRE(snarls_chr["ref"][0].type_variants == std::vector<std::string>{"0", "3", "6", "9"});
    }

    SECTION("complex_ins") {
        std::string pg_path = "../tests/graph_test/complex_ins.pg";
        std::string dist_path = "../tests/graph_test/complex_ins.dist";

        std::tie(stree, pg, root, path_graph, pp_overlay) = stoat::parse_graph_tree(pg_path, dist_path);
        auto snarls = stoat::save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = stoat::loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0);

        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);

        REQUIRE(stoat::pairToString(stoat::find_snarl_id(*stree, snarls_chr["ref"][0].snarl)) == "8_2");
        REQUIRE(stoat::stringToVector<std::string>(stoat::vectorPathToString(snarls_chr["ref"][0].snarl_paths)) == std::vector<std::string>{">2>8", ">2>3>4>6>8",">2>3>5>6>8",">2>3>5>7>8",">2>7>8"});
        REQUIRE(snarls_chr["ref"][0].start_positions == 8);
        REQUIRE(snarls_chr["ref"][0].end_positions == 10);
        REQUIRE(snarls_chr["ref"][0].type_variants == std::vector<std::string>{"0","3","3","3","1"});
    }

    SECTION("snp_and_nested_snp") {
        std::string pg_path = "../tests/graph_test/snp_and_nested_snp.pg";
        std::string dist_path = "../tests/graph_test/snp_and_nested_snp.dist";

        std::tie(stree, pg, root, path_graph, pp_overlay) = stoat::parse_graph_tree(pg_path, dist_path);
        auto snarls = stoat::save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = stoat::loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0);

        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 2);

        check_snarl(snarls_chr, "ref", "8_2", {">2>3>0>6>8", ">2>7>8"}, 8, 10, {"3/4", "1"}, *stree);
        check_snarl(snarls_chr, "ref", "6_3", {">3>4>6", ">3>5>6"}, 8, 10, {"1", "2"}, *stree);        
    }

    SECTION("nested_plus") {
        std::string pg_path = "../tests/graph_test/nested_plus.pg";
        std::string dist_path = "../tests/graph_test/nested_plus.dist";

        std::tie(stree, pg, root, path_graph, pp_overlay) = stoat::parse_graph_tree(pg_path, dist_path);
        auto snarls = stoat::save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = stoat::loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0);

        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 2);

        check_snarl(snarls_chr, "ref", "2_8", {">2>8", ">2>3>0>6>7>8", ">2>3>0>6>8"}, 8, 13, {"0", "5", "4"}, *stree);
        check_snarl(snarls_chr, "ref", "3_6", {">3>5>6", ">3>4>6"}, 9, 12, {"2", "2"}, *stree);
    }
}
