#include <catch.hpp>
#include "../../src/types_and_structs.hpp"
#include "../../src/utils.hpp"

using namespace std;
using namespace bdsg;

using handlegraph::step_handle_t;
using handlegraph::handle_t;
using handlegraph::net_handle_t;


TEST_CASE("Test input format", "[Path]") {

    std::unique_ptr<bdsg::SnarlDistanceIndex> stree;
    std::unique_ptr<handlegraph::PathHandleGraph> graph;
    // handlegraph::net_handle_t root;

    SECTION("Packedgraph (pg)") {
        std::string graph_path = "../tests/graph_test/simple_snp.pg";
        std::string dist_path = "../tests/graph_test/simple_snp.dist";
        std::tie(stree, graph) = stoat::load_graph_tree(graph_path, dist_path);
        REQUIRE(stree != nullptr);
        REQUIRE(graph != nullptr);
    }

    SECTION("Hashgraph (hg)") {
        std::string graph_path = "../tests/graph_test/simple_snp.hg";
        std::string dist_path = "../tests/graph_test/simple_snp.dist";
        std::tie(stree, graph) = stoat::load_graph_tree(graph_path, dist_path);
        REQUIRE(stree != nullptr);
        REQUIRE(graph != nullptr);
    }

    // SECTION("XG format (xg)") {
    //     std::string graph_path = "../tests/graph_test/simple_snp.hg";
    //     std::string dist_path = "../tests/graph_test/simple_snp.dist";
    //     std::tie(stree, graph) = stoat::load_graph_tree(graph_path, dist_path);
    //     REQUIRE(stree != nullptr);
    //     REQUIRE(graph != nullptr);
    // }

    SECTION("GBZ graph (gbz)") {
        std::string graph_path = "../tests/graph_test/simple_snp.gbz";
        std::string dist_path = "../tests/graph_test/simple_snp.dist";
        std::tie(stree, graph) = stoat::load_graph_tree(graph_path, dist_path);
        REQUIRE(stree != nullptr);
        REQUIRE(graph != nullptr);
    }
}

