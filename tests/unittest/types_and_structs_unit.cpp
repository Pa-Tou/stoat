#include <catch.hpp>
#include <vg/io/vpkg.hpp>
#include "../../src/io/register_io.hpp"
#include "../../src/types_and_structs.hpp"
#include "../../src/utils.hpp"

using namespace std;
using namespace bdsg;

using handlegraph::step_handle_t;
using handlegraph::handle_t;
using handlegraph::net_handle_t;

std::tuple<std::unique_ptr<bdsg::SnarlDistanceIndex>,std::unique_ptr<handlegraph::PathHandleGraph>> load_graph_tree(const std::string& graph_file, const std::string& dist_file) {

    // Load the snarl tree and graph
    // Tell the IO library about libvg types.
    if (!stoat::io::register_libvg_io()) {
        throw std::runtime_error("error[stoat vgio]: Could not register libvg types with libvgio");
    }

    std::unique_ptr<handlegraph::PathHandleGraph> graph = vg::io::VPKG::load_one<handlegraph::PathHandleGraph>(graph_file);
    std::unique_ptr<bdsg::SnarlDistanceIndex> distance_index = std::make_unique<bdsg::SnarlDistanceIndex>();
    distance_index->deserialize(dist_file);
    return std::make_tuple(std::move(distance_index), std::move(graph));
}

TEST_CASE("Test input format", "[Path]") {

    std::unique_ptr<bdsg::SnarlDistanceIndex> stree;
    std::unique_ptr<handlegraph::PathHandleGraph> graph;
    // handlegraph::net_handle_t root;

    SECTION("Packedgraph (pg)") {
        std::string graph_path = "../tests/test_data/test_graphs/simple_snp.pg";
        std::string dist_path = "../tests/test_data/test_graphs/simple_snp.dist";
        std::tie(stree, graph) = load_graph_tree(graph_path, dist_path);
        REQUIRE(stree != nullptr);
        REQUIRE(graph != nullptr);
    }

    SECTION("Hashgraph (hg)") {
        std::string graph_path = "../tests/test_data/test_graphs/simple_snp.hg";
        std::string dist_path = "../tests/test_data/test_graphs/simple_snp.dist";
        std::tie(stree, graph) = load_graph_tree(graph_path, dist_path);
        REQUIRE(stree != nullptr);
        REQUIRE(graph != nullptr);
    }

    // SECTION("XG format (xg)") {
    //     std::string graph_path = "../tests/test_data/test_graphs/simple_snp.hg";
    //     std::string dist_path = "../tests/test_data/test_graphs/simple_snp.dist";
    //     std::tie(stree, graph) = stoat::load_graph_tree(graph_path, dist_path);
    //     REQUIRE(stree != nullptr);
    //     REQUIRE(graph != nullptr);
    // }

    SECTION("GBZ graph (gbz)") {
        std::string graph_path = "../tests/test_data/test_graphs/simple_snp.gbz";
        std::string dist_path = "../tests/test_data/test_graphs/simple_snp.dist";
        std::tie(stree, graph) = load_graph_tree(graph_path, dist_path);
        REQUIRE(stree != nullptr);
        REQUIRE(graph != nullptr);
    }
}

