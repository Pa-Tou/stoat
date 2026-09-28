#include <catch.hpp>
#include <bdsg/hash_graph.hpp>
#include <bdsg/overlays/overlay_helper.hpp>
#include "../../src/snarl_traversals.hpp"
#include "../../src/log.hpp"
#include "../../src/gbzgraph.hpp"
#include <vg/io/vpkg.hpp>
#include "../../src/io/register_io.hpp"

using namespace stoat;

TEST_CASE( "snarl traversals nested bubbles",
          "[snarl_traversal]" ) {

    /*
                       5
                     /   \
            1       4 ----6    8
          /   \   /         \ / \
        0       3  ----------7---9
          \   /
            2

   */

    //bdsg::HashGraph graph;

    //std::vector<std::string> sequences = { "C", "C", "C", "A", "T", "C", "A", "C", "A", "A"};

    //std::vector<handlegraph::handle_t> nodes;
    //for (auto& seq : sequences) {
    //    nodes.emplace_back(graph.create_handle(seq));
    //}

    //graph.create_edge(nodes[0], nodes[1]);
    //graph.create_edge(nodes[0], nodes[2]);
    //graph.create_edge(nodes[1], nodes[3]);
    //graph.create_edge(nodes[2], nodes[3]);
    //graph.create_edge(nodes[3], nodes[4]);
    //graph.create_edge(nodes[3], nodes[7]);
    //graph.create_edge(nodes[4], nodes[5]);
    //graph.create_edge(nodes[4], nodes[6]);
    //graph.create_edge(nodes[5], nodes[6]);
    //graph.create_edge(nodes[6], nodes[7]);
    //graph.create_edge(nodes[7], nodes[8]);
    //graph.create_edge(nodes[7], nodes[9]);
    //graph.create_edge(nodes[8], nodes[9]);

    //// TODO one of these should really be the reference but idk how to add reference paths to a graph
    //std::vector<std::vector<std::size_t>> paths_seqs = { {0, 1, 3, 4, 5, 6, 7}, {0, 1, 3, 4, 6, 7}, {0, 2, 3, 7}, {0, 2, 3, 4, 6, 7}};
    //std::vector<handlegraph::path_handle_t> paths;

    //for (int path_i = 0 ; path_i < paths_seqs.size() ; path_i++) {
    //    paths.emplace_back(graph.create_path_handle("path"+std::to_string(path_i)));
    //    for (size_t node_i : paths_seqs[path_i]) {
    //        graph.append_step(paths.back(), nodes[node_i]);
    //    }
    //}

    //// vg isn't included so the distance index can only be built from the command line
    //graph.serialize("../tests/test_data/test_graphs/simple_nested_chain.hg");
    //int built = system("vg index -j ../tests/test_data/test_graphs/simple_nested_chain.dist ../tests/test_data/test_graphs/simple_nested_chain.hg"); 
    //
    ////Change sense of paths
    //built = system("vg convert --hap-locus path0 --new-sample path0 ../tests/test_data/test_graphs/simple_nested_chain.hg >../tests/test_data/test_graphs/simple_nested_chain1.hg"); 
    //built = system("vg convert --hap-locus path1 --new-sample path1 ../tests/test_data/test_graphs/simple_nested_chain1.hg >../tests/test_data/test_graphs/simple_nested_chain2.hg"); 
    //built = system("vg convert --ref-sample path0 ../tests/test_data/test_graphs/simple_nested_chain2.hg | vg convert -a - >../tests/test_data/test_graphs/simple_nested_chain3.hg"); 
    //built = system("mv ../tests/test_data/test_graphs/simple_nested_chain3.hg ../tests/test_data/test_graphs/simple_nested_chain.hg"); 
    //built = system("rm ../tests/test_data/test_graphs/simple_nested_chain1.hg"); 
    //built = system("rm ../tests/test_data/test_graphs/simple_nested_chain2.hg"); 
    //built = system("vg gbwt -x ../tests/test_data/test_graphs/simple_nested_chain.hg -E --gbz-format -g ../tests/test_data/test_graphs/simple_nested_chain.gbz "); 



    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/test_data/test_graphs/simple_nested_chain.dist");

    bdsg::HashGraph graph;
    graph.deserialize("../tests/test_data/test_graphs/simple_nested_chain.hg");

    std::vector<handlegraph::path_handle_t> paths;

    paths.emplace_back(graph.get_path_handle("path0#0#path0"));
    paths.emplace_back(graph.get_path_handle("path1#0#path1#0"));
    paths.emplace_back(graph.get_path_handle("path2"));
    paths.emplace_back(graph.get_path_handle("path3"));

    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&graph);

    GBZGraph gbz;
    std::ifstream instream;
    instream.open("../tests/test_data/test_graphs/simple_nested_chain.gbz");
    gbz.gbz.simple_sds_load(instream);
    instream.close();

    gbwt::GBWT* gbwt = &gbz.gbz.index;

    gbwt::FastLocate r_index;
    std::ifstream r_instream;
    r_instream.open("../tests/test_data/test_graphs/simple_nested_chain.ri");
    r_index.load(r_instream);
    r_instream.close();
    r_index.setGBWT(*gbwt);

    handlegraph::net_handle_t snarl1 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(2)));
    handlegraph::net_handle_t snarl2 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(5)));
    handlegraph::net_handle_t snarl3 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(6)));
    handlegraph::net_handle_t snarl4 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(9)));
    handlegraph::net_handle_t root_chain = distance_index.get_parent(snarl1);
    handlegraph::net_handle_t nested_chain = distance_index.get_parent(snarl3);

    std::vector<stoat::sample_hap_t> all_samples({stoat::sample_hap_t(*path_graph, paths[0]),
                                                   stoat::sample_hap_t(*path_graph, paths[1]),
                                                   stoat::sample_hap_t(*path_graph, paths[2]),
                                                   stoat::sample_hap_t(*path_graph, paths[3])});



    SECTION("get_all_walks_through_snarl") {
        SECTION("snarl1") {

            std::vector<stoat::PathTraversal> walks;
            get_all_walks_through_snarl(*path_graph, distance_index, snarl1, walks);
            REQUIRE(walks.size() == 2); 
            REQUIRE(walks[0].get_path().size() == 3);
            REQUIRE(walks[1].get_path().size() == 3);
            REQUIRE((walks[0].get_path()[1].get_node_id() == 2 || walks[1].get_path()[1].get_node_id() == 2));
            REQUIRE((walks[0].get_path()[1].get_node_id() == 3 || walks[1].get_path()[1].get_node_id() == 3));
            REQUIRE((walks[0].get_path()[1].get_node_id() != walks[1].get_path()[1].get_node_id()));
        }
        SECTION("snarl2") {

            std::vector<stoat::PathTraversal> walks;
            get_all_walks_through_snarl(*path_graph, distance_index, snarl2, walks);
            REQUIRE(walks.size() == 2); 
            REQUIRE((walks[0].get_path().size() == 2 || walks[0].get_path().size() == 5));
            REQUIRE((walks[1].get_path().size() == 2 || walks[1].get_path().size() == 5));
            REQUIRE(walks[0].get_path().size() != walks[1].get_path().size());
            const auto& del_path = walks[0].get_path().size() == 2 ? walks[0] : walks[1];
            const auto& ins_path = walks[0].get_path().size() == 2 ? walks[1] : walks[0];

            REQUIRE((del_path.to_string() == ">4>8" || del_path.to_string() == "<8<4"));
            REQUIRE((ins_path.to_string() == ">4>5>0>7>8" || ins_path.to_string() == "<8<7>0<5<4"));
        }
        SECTION("snarl3") {

            std::vector<stoat::PathTraversal> walks;
            get_all_walks_through_snarl(*path_graph, distance_index, snarl3, walks);
            REQUIRE(walks.size() == 2); 
            REQUIRE((walks[0].get_path().size() == 2 || walks[0].get_path().size() == 3));
            REQUIRE((walks[1].get_path().size() == 2 || walks[1].get_path().size() == 3));
            REQUIRE(walks[0].get_path().size() != walks[1].get_path().size());
            const auto& del_path = walks[0].get_path().size() == 2 ? walks[0] : walks[1];
            const auto& ins_path = walks[0].get_path().size() == 2 ? walks[1] : walks[0];

            REQUIRE((del_path.to_string() == ">5>7" || del_path.to_string() == "<7<5"));
            REQUIRE((ins_path.to_string() == ">5>6>7" || ins_path.to_string() == "<7<6<5"));
        }
        SECTION("snarl4") {

            std::vector<stoat::PathTraversal> walks;
            get_all_walks_through_snarl(*path_graph, distance_index, snarl4, walks);
            REQUIRE(walks.size() == 2); 
            REQUIRE((walks[0].get_path().size() == 2 || walks[0].get_path().size() == 3));
            REQUIRE((walks[1].get_path().size() == 2 || walks[1].get_path().size() == 3));
            REQUIRE(walks[0].get_path().size() != walks[1].get_path().size());
            const auto& del_path = walks[0].get_path().size() == 2 ? walks[0] : walks[1];
            const auto& ins_path = walks[0].get_path().size() == 2 ? walks[1] : walks[0];

            REQUIRE((del_path.to_string() == ">8>10" || del_path.to_string() == "<10<8"));
            REQUIRE((ins_path.to_string() == ">8>9>10" || ins_path.to_string() == "<10<9<8"));
        }

    }

    SECTION("get_haplotype_walks_through_snarl gbz") {
        SECTION("snarl1") {

            std::vector<stoat::PathTraversal> walks;
            get_haplotype_walks_through_snarl(*path_graph, *gbwt, r_index, distance_index, snarl1, walks);
            REQUIRE(walks.size() == 2); 
            REQUIRE(walks[0].get_path().size() == 3);
            REQUIRE(walks[1].get_path().size() == 3);
            REQUIRE((walks[0].get_path()[1].get_node_id() == 2 || walks[1].get_path()[1].get_node_id() == 2));
            REQUIRE((walks[0].get_path()[1].get_node_id() == 3 || walks[1].get_path()[1].get_node_id() == 3));
            REQUIRE((walks[0].get_path()[1].get_node_id() != walks[1].get_path()[1].get_node_id()));
        }
        SECTION("snarl2") {

            std::vector<stoat::PathTraversal> walks;
            get_haplotype_walks_through_snarl(*path_graph, *gbwt, r_index, distance_index, snarl2, walks);
            REQUIRE(walks.size() == 2); 
            REQUIRE((walks[0].get_path().size() == 2 || walks[0].get_path().size() == 5));
            REQUIRE((walks[1].get_path().size() == 2 || walks[1].get_path().size() == 5));
            REQUIRE(walks[0].get_path().size() != walks[1].get_path().size());
            const auto& del_path = walks[0].get_path().size() == 2 ? walks[0] : walks[1];
            const auto& ins_path = walks[0].get_path().size() == 2 ? walks[1] : walks[0];

            REQUIRE((del_path.to_string() == ">4>8" || del_path.to_string() == "<8<4"));
            REQUIRE((ins_path.to_string() == ">4>5>0>7>8" || ins_path.to_string() == "<8<7>0<5<4"));
        }
        SECTION("snarl3") {

            std::vector<stoat::PathTraversal> walks;
            get_haplotype_walks_through_snarl(*path_graph, *gbwt, r_index, distance_index, snarl3, walks);
            REQUIRE(walks.size() == 2); 
            REQUIRE((walks[0].get_path().size() == 2 || walks[0].get_path().size() == 3));
            REQUIRE((walks[1].get_path().size() == 2 || walks[1].get_path().size() == 3));
            REQUIRE(walks[0].get_path().size() != walks[1].get_path().size());
            const auto& del_path = walks[0].get_path().size() == 2 ? walks[0] : walks[1];
            const auto& ins_path = walks[0].get_path().size() == 2 ? walks[1] : walks[0];

            REQUIRE((del_path.to_string() == ">5>7" || del_path.to_string() == "<7<5"));
            REQUIRE((ins_path.to_string() == ">5>6>7" || ins_path.to_string() == "<7<6<5"));
        }
        SECTION("snarl4") {

            std::vector<stoat::PathTraversal> walks;
            get_haplotype_walks_through_snarl(*path_graph, *gbwt, r_index, distance_index, snarl4, walks);
            REQUIRE(walks.size() == 0); 
        }

    }
}



TEST_CASE( "traversal finder looping snarl", "[snarl_traversal]" ) {

    /*

             --------
            |   2    |
            \ / \    /
        0 ---1---3--4----5

    */

    bdsg::HashGraph graph;
    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/test_data/test_graphs/loop_with_indel.dist");

    graph.deserialize("../tests/test_data/test_graphs/loop_with_indel.hg");
    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&graph);

    GBZGraph gbz;
    std::ifstream instream;
    instream.open("../tests/test_data/test_graphs/loop_with_indel.gbz");
    gbz.gbz.simple_sds_load(instream);
    instream.close();

    gbwt::GBWT* gbwt = &gbz.gbz.index;

    gbwt::FastLocate r_index;
    std::ifstream r_instream;
    r_instream.open("../tests/test_data/test_graphs/loop_with_indel.ri");
    r_index.load(r_instream);
    r_instream.close();
    r_index.setGBWT(*gbwt);


    // Nested snarl
    handlegraph::net_handle_t snarl2 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(3)));
    // Duplication snarl
    handlegraph::net_handle_t snarl1 = distance_index.get_parent(distance_index.get_parent(snarl2));
    handlegraph::net_handle_t root_chain = distance_index.get_parent(snarl1);


    SECTION("get_all_walks_through_snarl") {
        SECTION("snarl1") {

            std::vector<stoat::PathTraversal> walks;
            get_all_walks_through_snarl(*path_graph, distance_index, snarl1, walks, 1);
            REQUIRE(walks.size() == 3); 
            REQUIRE((walks[0].get_path().size() == 5 || walks[0].get_path().size() == 8 || walks[0].get_path().size() == 11));
            REQUIRE((walks[1].get_path().size() == 5 || walks[1].get_path().size() == 8 || walks[1].get_path().size() == 11));
            REQUIRE(walks[0].get_path().size() != walks[1].get_path().size());
            REQUIRE(walks[0].get_path().size() != walks[2].get_path().size());
            REQUIRE(walks[1].get_path().size() != walks[2].get_path().size());
            const auto& del_path = walks[0].get_path().size() == 5 ? walks[0] : (walks[1].get_path().size() == 5 ? walks[1] : walks[2]);
            const auto& dup1_path = walks[0].get_path().size() == 8 ? walks[0] : (walks[1].get_path().size() == 8 ? walks[1] : walks[2]);
            const auto& dup2_path = walks[0].get_path().size() == 11 ? walks[0] : (walks[1].get_path().size() == 11 ? walks[1] : walks[2]);

            REQUIRE((del_path.to_string() == ">1>2>0>5>6" || del_path.to_string() == "<6<5>0<2<1"));
            REQUIRE((dup1_path.to_string() == ">1>2>0>5>2>0>5>6" || dup1_path.to_string() == "<6<5>0<2<5>0<2<1"));
            REQUIRE((dup2_path.to_string() == ">1>2>0>5>2>0>5>2>0>5>6" || dup2_path.to_string() == "<6<5>0<2<5>0<2<5>0<2<1"));
        }
        SECTION("snarl2") {

            std::vector<stoat::PathTraversal> walks;
            get_all_walks_through_snarl(*path_graph, distance_index, snarl2, walks);
            REQUIRE(walks.size() == 2); 
            REQUIRE((walks[0].get_path().size() == 2 || walks[0].get_path().size() == 3));
            REQUIRE((walks[1].get_path().size() == 2 || walks[1].get_path().size() == 3));
            REQUIRE(walks[0].get_path().size() != walks[1].get_path().size());
            const auto& del_path = walks[0].get_path().size() == 2 ? walks[0] : walks[1];
            const auto& ins_path = walks[0].get_path().size() == 2 ? walks[1] : walks[0];

            REQUIRE((del_path.to_string() == ">2>4" || del_path.to_string() == "<4<2"));
            REQUIRE((ins_path.to_string() == ">2>3>4" || ins_path.to_string() == "<4<3<2"));
        }
    }
    SECTION("get_haplotype_walks_through_snarl") {
        SECTION("snarl1") {

            std::vector<stoat::PathTraversal> walks;
            get_haplotype_walks_through_snarl(*path_graph, *gbwt, r_index, distance_index, snarl1, walks);
            REQUIRE(walks.size() == 2); 
            REQUIRE((walks[0].get_path().size() == 5 || walks[0].get_path().size() == 8));
            REQUIRE((walks[1].get_path().size() == 5 || walks[1].get_path().size() == 8));
            REQUIRE(walks[0].get_path().size() != walks[1].get_path().size());
            const auto& del_path = walks[0].get_path().size() == 5 ? walks[0] : walks[1];
            const auto& dup_path = walks[0].get_path().size() == 5 ? walks[1] : walks[0];

            REQUIRE((del_path.to_string() == ">1>2>0>5>6" || del_path.to_string() == "<6<5>0<2<1"));
            REQUIRE((dup_path.to_string() == ">1>2>0>5>2>0>5>6" || dup_path.to_string() == "<6<5>0<2<5>0<2<1"));
        }
        SECTION("snarl2") {

            std::vector<stoat::PathTraversal> walks;
            get_haplotype_walks_through_snarl(*path_graph, *gbwt, r_index, distance_index, snarl2, walks);
            REQUIRE(walks.size() == 2); 
            REQUIRE((walks[0].get_path().size() == 2 || walks[0].get_path().size() == 3));
            REQUIRE((walks[1].get_path().size() == 2 || walks[1].get_path().size() == 3));
            REQUIRE(walks[0].get_path().size() != walks[1].get_path().size());
            const auto& del_path = walks[0].get_path().size() == 2 ? walks[0] : walks[1];
            const auto& ins_path = walks[0].get_path().size() == 2 ? walks[1] : walks[0];

            REQUIRE((del_path.to_string() == ">2>4" || del_path.to_string() == "<4<2"));
            REQUIRE((ins_path.to_string() == ">2>3>4" || ins_path.to_string() == "<4<3<2"));
        }
    }


}


