#include <catch.hpp>
#include <bdsg/hash_graph.hpp>
#include <bdsg/overlays/overlay_helper.hpp>
#include "../../src/path_partitioner.hpp"
#include "../../src/log.hpp"
#include "../../src/gbzgraph.hpp"
#include <vg/io/vpkg.hpp>
#include "../../src/io/register_io.hpp"

using namespace stoat_graph;

TEST_CASE( "Path partitioner nested bubbles",
          "[path_partitioner]" ) {

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



    SECTION("partition_embedded_paths_in_snarl") {
        // This isn't really a good test because all the snarls are regular

        // Should be {0,1} and {2,3}
        std::vector<size_t> alleles_per_sample1 = partition_embedded_paths_in_snarl(*path_graph, distance_index, snarl1, all_samples);
        REQUIRE(alleles_per_sample1.size() == all_samples.size());
        REQUIRE(alleles_per_sample1[0] == alleles_per_sample1[1]);
        REQUIRE(alleles_per_sample1[0] != alleles_per_sample1[2]);
        REQUIRE(alleles_per_sample1[2] == alleles_per_sample1[3]);
        REQUIRE(alleles_per_sample1[0] != std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample1[2] != std::numeric_limits<size_t>::max());

        // Should be {0,1,3} and {2}
        std::vector<size_t> alleles_per_sample2 = partition_embedded_paths_in_snarl(*path_graph, distance_index, snarl2, all_samples);
        REQUIRE(alleles_per_sample2.size() == all_samples.size());
        REQUIRE(alleles_per_sample2[0] == alleles_per_sample2[1]);
        REQUIRE(alleles_per_sample2[0] == alleles_per_sample2[3]);
        REQUIRE(alleles_per_sample2[2] != alleles_per_sample2[3]);
        REQUIRE(alleles_per_sample2[0] != std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample2[3] != std::numeric_limits<size_t>::max());

        // Should be {0}, {1,3}. 2 didn't go through this snarl
        std::vector<size_t> alleles_per_sample3 = partition_embedded_paths_in_snarl(*path_graph, distance_index, snarl3, all_samples);
        REQUIRE(alleles_per_sample3.size() == all_samples.size());
        REQUIRE(alleles_per_sample3[0] != alleles_per_sample3[1]);
        REQUIRE(alleles_per_sample3[1] == alleles_per_sample3[3]);
        REQUIRE(alleles_per_sample3[2] == std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample3[0] != std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample3[1] != std::numeric_limits<size_t>::max());
    }
}

TEST_CASE( "Path partitioner nested bubbles gbz",
          "[path_partitioner]" ) {

    /*
                       5
                     /   \
            1       4 ----6    8
          /   \   /         \ / \
        0       3  ----------7---9
          \   /
            2

   */
    //built = system("vg gbwt -x ../tests/test_data/test_graphs/simple_nested_chain.hg -E --gbz-format -g ../tests/test_data/test_graphs/simple_nested_chain.gbz "); 


    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/test_data/test_graphs/simple_nested_chain.dist");


    //std::unique_ptr<GBZGraph> gbz = vg::io::VPKG::load_one<GBZGraph>("../tests/test_data/test_graphs/simple_nested_chain.gbz");
    GBZGraph gbz;
    std::ifstream instream;
    instream.open("../tests/test_data/test_graphs/simple_nested_chain.gbz");
    gbz.gbz.simple_sds_load(instream);
    instream.close();

    gbwt::GBWT* gbwt = &gbz.gbz.index;

    //handlegraph::PathHandleGraph* graph = gbz.get();




    std::vector<handlegraph::path_handle_t> paths;

    paths.emplace_back(gbz.get_path_handle("path0#0#path0"));
    paths.emplace_back(gbz.get_path_handle("path1#0#path1#0"));
    paths.emplace_back(gbz.get_path_handle("path2"));
    paths.emplace_back(gbz.get_path_handle("path3"));

    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&gbz);


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



    SECTION("partition_embedded_paths_in_snarl") {
        // This isn't really a good test because all the snarls are regular

        // Should be {0,1} and {2,3}
        
        std::vector<PathTraversal> paths_per_allele1;
        std::vector<size_t> alleles_per_sample1 = partition_embedded_paths_in_snarl_with_gbwt(*path_graph, *gbwt, distance_index, snarl1,
                                                         all_samples, paths_per_allele1);

        REQUIRE(alleles_per_sample1.size() == all_samples.size());
        REQUIRE(alleles_per_sample1[0] == alleles_per_sample1[1]);
        REQUIRE(alleles_per_sample1[0] != alleles_per_sample1[2]);
        REQUIRE(alleles_per_sample1[2] == alleles_per_sample1[3]);
        REQUIRE(alleles_per_sample1[0] != std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample1[2] != std::numeric_limits<size_t>::max());

        // Should be {0,1,3} and {2}
        
        std::vector<PathTraversal> paths_per_allele2;
        std::vector<size_t> alleles_per_sample2 = partition_embedded_paths_in_snarl_with_gbwt(*path_graph, *gbwt, distance_index, snarl2,
                                                         all_samples, paths_per_allele2);
        REQUIRE(alleles_per_sample2.size() == all_samples.size());
        REQUIRE(alleles_per_sample2[0] == alleles_per_sample2[1]);
        REQUIRE(alleles_per_sample2[0] == alleles_per_sample2[3]);
        REQUIRE(alleles_per_sample2[2] != alleles_per_sample2[3]);
        REQUIRE(alleles_per_sample2[0] != std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample2[3] != std::numeric_limits<size_t>::max());

        // Should be {0}, {1,3}. 2 didn't go through this snarl
        std::vector<PathTraversal> paths_per_allele3;
        std::vector<size_t> alleles_per_sample3 = partition_embedded_paths_in_snarl_with_gbwt(*path_graph, *gbwt, distance_index, snarl3,
                                                         all_samples, paths_per_allele3);
        REQUIRE(alleles_per_sample3.size() == all_samples.size());
        REQUIRE(alleles_per_sample3[0] != alleles_per_sample3[1]);
        REQUIRE(alleles_per_sample3[1] == alleles_per_sample3[3]);
        REQUIRE(alleles_per_sample3[2] == std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample3[0] != std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample3[1] != std::numeric_limits<size_t>::max());
    }
}

TEST_CASE( "Path partitioner multiple nested bubbles",
          "[path_partitioner]" ) {

    /*
                       5
                     /   \
            1       4 ----6    8
          /   \   /         \ / \
        0       3  ----------7---9
          \   /
            2
    this graph is duplicated 4x

   */

    //bdsg::HashGraph graph;

    //std::vector<std::string> sequences = { "C", "C", "C", "A", "T", "C", "A", "C", "A", "A", 
    //                                       "C", "C", "C", "A", "T", "C", "A", "C", "A", "A",  
    //                                       "C", "C", "C", "A", "T", "C", "A", "C", "A", "A",
    //                                       "C", "C", "C", "A", "T", "C", "A", "C", "A", "A", 
    //                                       "C", "C", "C", "A", "T", "C", "A", "C", "A", "A", 
    //                                       "C", "C", "C", "A", "T", "C", "A", "C", "A", "A"};

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

    //graph.create_edge(nodes[9], nodes[10]);

    //graph.create_edge(nodes[10], nodes[11]);
    //graph.create_edge(nodes[10], nodes[12]);
    //graph.create_edge(nodes[11], nodes[13]);
    //graph.create_edge(nodes[12], nodes[13]);
    //graph.create_edge(nodes[13], nodes[14]);
    //graph.create_edge(nodes[13], nodes[17]);
    //graph.create_edge(nodes[14], nodes[15]);
    //graph.create_edge(nodes[14], nodes[16]);
    //graph.create_edge(nodes[15], nodes[16]);
    //graph.create_edge(nodes[16], nodes[17]);
    //graph.create_edge(nodes[17], nodes[18]);
    //graph.create_edge(nodes[17], nodes[19]);
    //graph.create_edge(nodes[18], nodes[19]);


    //graph.create_edge(nodes[20], nodes[21]);
    //graph.create_edge(nodes[20], nodes[22]);
    //graph.create_edge(nodes[21], nodes[23]);
    //graph.create_edge(nodes[22], nodes[23]);
    //graph.create_edge(nodes[23], nodes[24]);
    //graph.create_edge(nodes[23], nodes[27]);
    //graph.create_edge(nodes[24], nodes[25]);
    //graph.create_edge(nodes[24], nodes[26]);
    //graph.create_edge(nodes[25], nodes[26]);
    //graph.create_edge(nodes[26], nodes[27]);
    //graph.create_edge(nodes[27], nodes[28]);
    //graph.create_edge(nodes[27], nodes[29]);
    //graph.create_edge(nodes[28], nodes[29]);

    //graph.create_edge(nodes[29], nodes[30]);

    //graph.create_edge(nodes[30], nodes[31]);
    //graph.create_edge(nodes[30], nodes[32]);
    //graph.create_edge(nodes[31], nodes[33]);
    //graph.create_edge(nodes[32], nodes[33]);
    //graph.create_edge(nodes[33], nodes[34]);
    //graph.create_edge(nodes[33], nodes[37]);
    //graph.create_edge(nodes[34], nodes[35]);
    //graph.create_edge(nodes[34], nodes[36]);
    //graph.create_edge(nodes[35], nodes[36]);
    //graph.create_edge(nodes[36], nodes[37]);
    //graph.create_edge(nodes[37], nodes[38]);
    //graph.create_edge(nodes[37], nodes[39]);
    //graph.create_edge(nodes[38], nodes[39]);

    //graph.create_edge(nodes[39], nodes[40]);

    //graph.create_edge(nodes[40], nodes[41]);
    //graph.create_edge(nodes[40], nodes[42]);
    //graph.create_edge(nodes[41], nodes[43]);
    //graph.create_edge(nodes[42], nodes[43]);
    //graph.create_edge(nodes[43], nodes[44]);
    //graph.create_edge(nodes[43], nodes[47]);
    //graph.create_edge(nodes[44], nodes[45]);
    //graph.create_edge(nodes[44], nodes[46]);
    //graph.create_edge(nodes[45], nodes[46]);
    //graph.create_edge(nodes[46], nodes[47]);
    //graph.create_edge(nodes[47], nodes[48]);
    //graph.create_edge(nodes[47], nodes[49]);
    //graph.create_edge(nodes[48], nodes[49]);

    //graph.create_edge(nodes[49], nodes[50]);

    //graph.create_edge(nodes[50], nodes[51]);
    //graph.create_edge(nodes[50], nodes[52]);
    //graph.create_edge(nodes[51], nodes[53]);
    //graph.create_edge(nodes[52], nodes[53]);
    //graph.create_edge(nodes[53], nodes[54]);
    //graph.create_edge(nodes[53], nodes[57]);
    //graph.create_edge(nodes[54], nodes[55]);
    //graph.create_edge(nodes[54], nodes[56]);
    //graph.create_edge(nodes[55], nodes[56]);
    //graph.create_edge(nodes[56], nodes[57]);
    //graph.create_edge(nodes[57], nodes[58]);
    //graph.create_edge(nodes[57], nodes[59]);
    //graph.create_edge(nodes[58], nodes[59]);

    //// TODO one of these should really be the reference but idk how to add reference paths to a graph
    //std::vector<std::vector<std::size_t>> paths_seqs = { {0, 1, 3, 4, 5, 6, 7}, {0, 1, 3, 4, 6, 7}, {0, 2, 3, 7}, {0, 2, 3, 4, 6, 7},
    //                                                      {10, 11, 13, 14, 15, 16, 17}, {10, 11, 13, 14, 16, 17}, {10, 12, 13, 17}, {10, 12, 13, 14, 16, 17},
    //                                                      {20, 21, 23, 24, 25, 26, 27}, {20, 21, 23, 24, 26, 27}, {20, 22, 23, 27}, {20, 22, 23, 24, 26, 27},
    //                                                      {30, 31, 33, 34, 35, 36, 37}, {30, 31, 33, 34, 36, 37}, {30, 32, 33, 37}, {30, 32, 33, 34, 36, 37},
    //                                                      {39, 40, 41, 43, 44, 45, 46, 47}, {40, 41, 43, 44, 46, 47}, {40, 42, 43, 47}, {40, 42, 43, 44, 46, 47},
    //                                                      {49, 50, 51, 53, 54, 55, 56, 57}, {50, 51, 53, 54, 56, 57}, {50, 52, 53, 57}, {50, 52, 53, 54, 56, 57}};
    //std::vector<handlegraph::path_handle_t> paths;

    //for (int path_i = 0 ; path_i < paths_seqs.size() ; path_i++) {
    //    paths.emplace_back(graph.create_path_handle("path"+std::to_string(path_i)+"#0#0#0"));
    //    for (size_t node_i : paths_seqs[path_i]) {
    //        graph.append_step(paths.back(), nodes[node_i]);
    //    }
    //}

    //// vg isn't included so the distance index can only be built from the command line
    //graph.serialize("../tests/test_data/test_graphs/simple_nested_chains.hg");
    //int built = system("vg index -j ../tests/test_data/test_graphs/simple_nested_chains.dist ../tests/test_data/test_graphs/simple_nested_chains.hg"); 
    //   bdsg::SnarlDistanceIndex distance_index;
    //distance_index.deserialize("../tests/test_data/test_graphs/simple_nested_chains.dist");

   // bdsg::HashGraph graph;
   // graph.deserialize("../tests/test_data/test_graphs/simple_nested_chains.hg");

    //bdsg::PathPositionOverlayHelper overlay_helper;
    //auto path_graph = overlay_helper.apply(&graph);

}

TEST_CASE( "Path partitioner nested bubbles distanceless index",
          "[path_partitioner]" ) {

    /*
                       5
                     /   \
            1       4 ----6    8
          /   \   /         \ / \
        0       3  ----------7---9
          \   /
            2

   */

    //int built = system("vg index --snarl-limit 0 -j ../tests/test_data/test_graphs/simple_nested_chain.nodist.dist ../tests/test_data/test_graphs/simple_nested_chain.hg"); 
    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/test_data/test_graphs/simple_nested_chain.nodist.dist");

    bdsg::HashGraph graph;
    graph.deserialize("../tests/test_data/test_graphs/simple_nested_chain.hg");

    std::vector<handlegraph::path_handle_t> paths;

    paths.emplace_back(graph.get_path_handle("path0#0#path0"));
    paths.emplace_back(graph.get_path_handle("path1#0#path1#0"));
    paths.emplace_back(graph.get_path_handle("path2"));
    paths.emplace_back(graph.get_path_handle("path3"));

    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&graph);


    handlegraph::net_handle_t snarl1 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(2)));
    handlegraph::net_handle_t snarl2 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(5)));
    handlegraph::net_handle_t snarl3 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(6)));
    handlegraph::net_handle_t snarl4 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(9)));
    handlegraph::net_handle_t root_chain = distance_index.get_parent(snarl1);
    handlegraph::net_handle_t nested_chain = distance_index.get_parent(snarl3);

    std::vector<stoat::sample_hap_t> all_samples ({stoat::sample_hap_t(*path_graph, paths[0]),
                                         stoat::sample_hap_t(*path_graph, paths[1]),
                                         stoat::sample_hap_t(*path_graph, paths[2]),
                                         stoat::sample_hap_t(*path_graph, paths[3])});

    SECTION("partition_embedded_paths_in_snarl") {
        // This isn't really a good test because all the snarls are regular

        // Should be {0,1} and {2,3}
        std::vector<size_t> alleles_per_sample1 = partition_embedded_paths_in_snarl(*path_graph, distance_index, snarl1, all_samples);
        REQUIRE(alleles_per_sample1.size() == all_samples.size());
        REQUIRE(alleles_per_sample1[0] == alleles_per_sample1[1]);
        REQUIRE(alleles_per_sample1[1] != alleles_per_sample1[2]);
        REQUIRE(alleles_per_sample1[2] == alleles_per_sample1[3]);
        REQUIRE(alleles_per_sample1[0] != std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample1[2] != std::numeric_limits<size_t>::max());


        // Should be {0,1,3} and {2}
        std::vector<size_t> alleles_per_sample2 = partition_embedded_paths_in_snarl(*path_graph, distance_index, snarl2, all_samples);
        REQUIRE(alleles_per_sample2.size() == all_samples.size());
        REQUIRE(alleles_per_sample2[0] == alleles_per_sample2[1]);
        REQUIRE(alleles_per_sample2[0] == alleles_per_sample2[3]);
        REQUIRE(alleles_per_sample2[2] != alleles_per_sample2[3]);
        REQUIRE(alleles_per_sample2[0] != std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample2[2] != std::numeric_limits<size_t>::max());

        // Should be {0}, {1,3}. 2 didn't go through this snarl
        std::vector<size_t> alleles_per_sample3 = partition_embedded_paths_in_snarl(*path_graph, distance_index, snarl3, all_samples);
        REQUIRE(alleles_per_sample3.size() == all_samples.size());
        REQUIRE(alleles_per_sample3[0] != alleles_per_sample3[1]);
        REQUIRE(alleles_per_sample3[1] == alleles_per_sample3[3]);
        REQUIRE(alleles_per_sample3[2] == std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample3[0] != std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample3[1] != std::numeric_limits<size_t>::max());
    }

}

TEST_CASE( "Path partitioner finder looping snarl", "[path_partitioner]" ) {

    /*

             --------
            |   2    |
            \ / \    /
        0 ---1---3--4----5

    */

    bdsg::HashGraph graph;

    //std::vector<std::string> sequences = {"AAAAAAAAAA", "A", "G", "C", "T",  "AAAAAAAAA"};

    //std::vector<handlegraph::handle_t> nodes;
    //for (auto& seq : sequences) {
    //    nodes.emplace_back(graph.create_handle(seq));
    //}

    //graph.create_edge(nodes[0], nodes[1]);
    //graph.create_edge(nodes[1], nodes[2]);
    //graph.create_edge(nodes[1], nodes[3]);
    //graph.create_edge(nodes[2], nodes[3]);
    //graph.create_edge(nodes[3], nodes[4]);
    //graph.create_edge(nodes[4], nodes[1]);
    //graph.create_edge(nodes[4], nodes[5]);


    //// Paths 0 and 2 take the insertion, but paths 1 and 2 take the duplication, and the deletion
    //std::vector<std::vector<std::size_t>> path_seqs = { {0, 1, 2, 3, 4, 5}, {0, 1, 3, 4, 1, 3, 4, 5}, {0, 1, 2, 3, 4, 1, 3, 4, 5}};
    //std::vector<handlegraph::path_handle_t> paths;

    //for (int path_i = 0 ; path_i < path_seqs.size() ; path_i++) {
    //    paths.emplace_back(graph.create_path_handle("path"+std::to_string(path_i)));
    //    for (size_t node_i : path_seqs[path_i]) {
    //        graph.append_step(paths.back(), nodes[node_i]);
    //    }
    //}

    //// vg isn't included so the distance index can only be built from the command line
    //graph.serialize("../tests/test_data/test_graphs/loop_with_indel.hg");
    //int built = system("vg index -j ../tests/test_data/test_graphs/loop_with_indel.dist ../tests/test_data/test_graphs/loop_with_indel.hg"); 

    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/test_data/test_graphs/loop_with_indel.dist");

    graph.deserialize("../tests/test_data/test_graphs/loop_with_indel.hg");
    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&graph);

    std::vector<handlegraph::path_handle_t> paths;

    for (int path_i = 0 ; path_i < 3 ; path_i++) {
        paths.emplace_back(graph.get_path_handle("path"+std::to_string(path_i)));
    }

    // Nested snarl
    handlegraph::net_handle_t snarl2 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(3)));
    // Duplication snarl
    handlegraph::net_handle_t snarl1 = distance_index.get_parent(distance_index.get_parent(snarl2));
    handlegraph::net_handle_t root_chain = distance_index.get_parent(snarl1);


    std::vector<stoat::sample_hap_t> all_samples ({stoat::sample_hap_t(*path_graph, paths[0]),
                                         stoat::sample_hap_t(*path_graph, paths[1]),
                                         stoat::sample_hap_t(*path_graph, paths[2])});

    SECTION("partition_embedded_paths_in_snarl") {
        // This isn't really a good test because all the snarls are regular

        // Should be {0} and {1,2}
        std::vector<size_t> alleles_per_sample1 = partition_embedded_paths_in_snarl(*path_graph, distance_index, snarl1, all_samples);
        REQUIRE(alleles_per_sample1.size() == all_samples.size());
        REQUIRE(alleles_per_sample1[0] != alleles_per_sample1[1]);
        REQUIRE(alleles_per_sample1[1] == alleles_per_sample1[2]);
        REQUIRE(alleles_per_sample1[0] != std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample1[1] != std::numeric_limits<size_t>::max());

        // Should be {0}, {1} and {2}

        std::vector<size_t> alleles_per_sample2 = partition_embedded_paths_in_snarl(*path_graph, distance_index, snarl2, all_samples);
        REQUIRE(alleles_per_sample2.size() == all_samples.size());
        REQUIRE(alleles_per_sample2[0] != alleles_per_sample2[1]);
        REQUIRE(alleles_per_sample2[0] != alleles_per_sample2[2]);
        REQUIRE(alleles_per_sample2[1] != alleles_per_sample2[2]);
        REQUIRE(alleles_per_sample2[0] != std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample2[1] != std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample2[2] != std::numeric_limits<size_t>::max());
    }


}

TEST_CASE( "Path partitioner finder looping snarl gbz", "[path_partitioner]" ) {

    /*

             --------
            |   2    |
            \ / \    /
        0 ---1---3--4----5

    */


    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/test_data/test_graphs/loop_with_indel.dist");

    GBZGraph gbz;
    std::ifstream instream;
    instream.open("../tests/test_data/test_graphs/loop_with_indel.gbz");
    gbz.gbz.simple_sds_load(instream);
    instream.close();

    gbwt::GBWT* gbwt = &gbz.gbz.index;

    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&gbz);

    std::vector<handlegraph::path_handle_t> paths;

    for (int path_i = 0 ; path_i < 3 ; path_i++) {
        paths.emplace_back(gbz.get_path_handle("path"+std::to_string(path_i)));
    }

    // Nested snarl
    handlegraph::net_handle_t snarl2 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(3)));
    // Duplication snarl
    handlegraph::net_handle_t snarl1 = distance_index.get_parent(distance_index.get_parent(snarl2));
    handlegraph::net_handle_t root_chain = distance_index.get_parent(snarl1);


    std::vector<stoat::sample_hap_t> all_samples ({stoat::sample_hap_t(*path_graph, paths[0]),
                                         stoat::sample_hap_t(*path_graph, paths[1]),
                                         stoat::sample_hap_t(*path_graph, paths[2])});

    SECTION("partition_embedded_paths_in_snarl") {
        // This isn't really a good test because all the snarls are regular

        // Should be {0} and {1,2}

        std::vector<PathTraversal> paths_per_allele1;
        std::vector<size_t> alleles_per_sample1 = partition_embedded_paths_in_snarl_with_gbwt(*path_graph, *gbwt, distance_index, snarl1,
                                                         all_samples, paths_per_allele1);
        REQUIRE(alleles_per_sample1.size() == all_samples.size());
        REQUIRE(alleles_per_sample1[0] != alleles_per_sample1[1]);
        REQUIRE(alleles_per_sample1[1] == alleles_per_sample1[2]);
        REQUIRE(alleles_per_sample1[0] != std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample1[1] != std::numeric_limits<size_t>::max());

        // Should be {0}, {1} and {2}

        std::vector<PathTraversal> paths_per_allele2;
        std::vector<size_t> alleles_per_sample2 = partition_embedded_paths_in_snarl_with_gbwt(*path_graph, *gbwt, distance_index, snarl2,
                                                         all_samples, paths_per_allele2);
        REQUIRE(alleles_per_sample2.size() == all_samples.size());
        REQUIRE(alleles_per_sample2[0] != alleles_per_sample2[1]);
        REQUIRE(alleles_per_sample2[0] != alleles_per_sample2[2]);
        REQUIRE(alleles_per_sample2[1] != alleles_per_sample2[2]);
        REQUIRE(alleles_per_sample2[0] != std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample2[1] != std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample2[2] != std::numeric_limits<size_t>::max());
    }


}

TEST_CASE( "Path partitioner finder looping snarl with fragments", "[path_partitioner]" ) {
    // Different fragments with the same sample and haplotype count as separate paths

    /*

             --------
            |   2    |
            \ / \    /
        0 ---1---3--4----5

    */

    bdsg::HashGraph graph;

    std::vector<std::string> sequences = {"AAAAAAAAAA", "A", "G", "C", "T",  "AAAAAAAAA"};

    std::vector<handlegraph::handle_t> nodes;
    for (auto& seq : sequences) {
        nodes.emplace_back(graph.create_handle(seq));
    }

    graph.create_edge(nodes[0], nodes[1]);
    graph.create_edge(nodes[1], nodes[2]);
    graph.create_edge(nodes[1], nodes[3]);
    graph.create_edge(nodes[2], nodes[3]);
    graph.create_edge(nodes[3], nodes[4]);
    graph.create_edge(nodes[4], nodes[1]);
    graph.create_edge(nodes[4], nodes[5]);

    // Add path 0 that goes through the loop twice
    // The paths are given the same name and haplotype but different loci (I think)
    std::vector<std::vector<std::size_t>> path_seqs = { {0, 1, 2, 3}, {4, 1, 2, 3, 4}};
    std::vector<handlegraph::path_handle_t> paths;

    for (int path_i = 0 ; path_i < path_seqs.size() ; path_i++) {
        paths.emplace_back(graph.create_path_handle("path0#0#"+std::to_string(path_i)+"#0"));
        for (size_t node_i : path_seqs[path_i]) {
            graph.append_step(paths.back(), nodes[node_i]);
        }
    }

    // Add path 1 that goes through three times
    path_seqs = {{0, 1, 2, 3}, {4, 1, 2, 3, 4}, {4, 1, 2, 3, 4, 5}};

    for (int path_i = 0 ; path_i < path_seqs.size() ; path_i++) {
        paths.emplace_back(graph.create_path_handle("path1#0#"+std::to_string(path_i)+"#0"));
        for (size_t node_i : path_seqs[path_i]) {
            graph.append_step(paths.back(), nodes[node_i]);
        }
    }
    // Add path 2 that goes through four times
    path_seqs = { {0, 1, 2, 3}, {4, 1, 2, 3, 4}, {3,4, 1, 2, 3, 4}, {4, 1, 2, 3, 4, 5}};

    for (int path_i = 0 ; path_i < path_seqs.size() ; path_i++) {
        paths.emplace_back(graph.create_path_handle("path2#0#"+std::to_string(path_i)+"#0"));
        for (size_t node_i : path_seqs[path_i]) {
            graph.append_step(paths.back(), nodes[node_i]);
        }
    }

    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/test_data/test_graphs/loop_with_indel.dist");

    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&graph);

    // Nested snarl
    handlegraph::net_handle_t snarl2 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(3)));
    // Duplication snarl
    handlegraph::net_handle_t snarl1 = distance_index.get_parent(distance_index.get_parent(snarl2));
    handlegraph::net_handle_t root_chain = distance_index.get_parent(snarl1);

    std::vector<stoat::sample_hap_t> all_samples;
    for (const auto& path : paths) {
        all_samples.emplace_back(*path_graph, path);
    }


    SECTION("partition_embedded_paths_in_snarl") {

        std::vector<size_t> alleles_per_sample2 = partition_embedded_paths_in_snarl(*path_graph, distance_index, snarl2, all_samples);
        REQUIRE(alleles_per_sample2.size() == all_samples.size());
        REQUIRE(alleles_per_sample2[0] == alleles_per_sample2[1]);
        REQUIRE(alleles_per_sample2[0] == alleles_per_sample2[2]);
        REQUIRE(alleles_per_sample2[0] == alleles_per_sample2[3]);
        REQUIRE(alleles_per_sample2[0] == alleles_per_sample2[4]);
        REQUIRE(alleles_per_sample2[0] == alleles_per_sample2[5]);
        REQUIRE(alleles_per_sample2[0] == alleles_per_sample2[6]);
        REQUIRE(alleles_per_sample2[0] == alleles_per_sample2[7]);
        REQUIRE(alleles_per_sample2[0] == alleles_per_sample2[8]);
        REQUIRE(alleles_per_sample2[0] != std::numeric_limits<size_t>::max());

    }

}
TEST_CASE( "Path partitioner finder looping snarl with fragments gbz", "[path_partitioner][bug]" ) {
    // Different fragments with the same sample and haplotype count as separate paths

    /*

             --------
            |   2    |
            \ / \    /
        0 ---1---3--4----5

    */

    //bdsg::HashGraph graph;

    //std::vector<std::string> sequences = {"AAAAAAAAAA", "A", "G", "C", "T",  "AAAAAAAAA"};

    //std::vector<handlegraph::handle_t> nodes;
    //for (auto& seq : sequences) {
    //    nodes.emplace_back(graph.create_handle(seq));
    //}

    //graph.create_edge(nodes[0], nodes[1]);
    //graph.create_edge(nodes[1], nodes[2]);
    //graph.create_edge(nodes[1], nodes[3]);
    //graph.create_edge(nodes[2], nodes[3]);
    //graph.create_edge(nodes[3], nodes[4]);
    //graph.create_edge(nodes[4], nodes[1]);
    //graph.create_edge(nodes[4], nodes[5]);

    //// Add path 0 that goes through the loop twice
    //// The paths are given the same name and haplotype but different loci (I think)
    //std::vector<std::vector<std::size_t>> path_seqs = { {0, 1, 2, 3}, {4, 1, 2, 3, 4}};
    //std::vector<handlegraph::path_handle_t> paths;

    //for (int path_i = 0 ; path_i < path_seqs.size() ; path_i++) {
    //    paths.emplace_back(graph.create_path_handle("path0#0#"+std::to_string(path_i)+"#0"));
    //    for (size_t node_i : path_seqs[path_i]) {
    //        graph.append_step(paths.back(), nodes[node_i]);
    //    }
    //}

    //// Add path 1 that goes through three times
    //path_seqs = {{0, 1, 2, 3}, {4, 1, 2, 3, 4}, {4, 1, 2, 3, 4, 5}};

    //for (int path_i = 0 ; path_i < path_seqs.size() ; path_i++) {
    //    paths.emplace_back(graph.create_path_handle("path1#0#"+std::to_string(path_i)+"#0"));
    //    for (size_t node_i : path_seqs[path_i]) {
    //        graph.append_step(paths.back(), nodes[node_i]);
    //    }
    //}
    //// Add path 2 that goes through four times
    //path_seqs = { {0, 1, 2, 3}, {4, 1, 2, 3, 4}, {3,4, 1, 2, 3, 4}, {4, 1, 2, 3, 4, 5}};

    //for (int path_i = 0 ; path_i < path_seqs.size() ; path_i++) {
    //    paths.emplace_back(graph.create_path_handle("path2#0#"+std::to_string(path_i)+"#0"));
    //    for (size_t node_i : path_seqs[path_i]) {
    //        graph.append_step(paths.back(), nodes[node_i]);
    //    }
    //}
    //graph.serialize("../tests/test_data/test_graphs/loop_with_indel_fragmented.hg");
    //int built = system("vg gbwt -x ../tests/test_data/test_graphs/split_paths.hg -E --gbz-format -g ../tests/test_data/test_graphs/loop_with_indel_fragmented.gbz "); 

    GBZGraph gbz;
    std::ifstream instream;
    instream.open("../tests/test_data/test_graphs/loop_with_indel_fragmented.gbz");
    gbz.gbz.simple_sds_load(instream);
    instream.close();

    gbwt::GBWT* gbwt = &gbz.gbz.index;

    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&gbz);

    std::vector<handlegraph::path_handle_t> paths;


    for (int path_i = 0 ; path_i < 2 ; path_i++) {
        paths.emplace_back(gbz.get_path_handle("path0#0#"+std::to_string(path_i)+"#0"));
    }
    for (int path_i = 0 ; path_i < 3 ; path_i++) {
        paths.emplace_back(gbz.get_path_handle("path1#0#"+std::to_string(path_i)+"#0"));
    }

    for (int path_i = 0 ; path_i < 4 ; path_i++) {
        paths.emplace_back(gbz.get_path_handle("path2#0#"+std::to_string(path_i)+"#0"));
    }

    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/test_data/test_graphs/loop_with_indel.dist");


    // Nested snarl
    handlegraph::net_handle_t snarl2 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(3)));
    // Duplication snarl
    handlegraph::net_handle_t snarl1 = distance_index.get_parent(distance_index.get_parent(snarl2));
    handlegraph::net_handle_t root_chain = distance_index.get_parent(snarl1);

    std::vector<stoat::sample_hap_t> all_samples;
    for (const auto& path : paths) {
        all_samples.emplace_back(*path_graph, path);
    }


    SECTION("partition_embedded_paths_in_snarl") {

        std::vector<PathTraversal> paths_per_allele2;
        std::vector<size_t> alleles_per_sample2 = partition_embedded_paths_in_snarl_with_gbwt(*path_graph, *gbwt, distance_index, snarl2,
                                                         all_samples, paths_per_allele2);
        REQUIRE(alleles_per_sample2.size() == all_samples.size());
        REQUIRE(alleles_per_sample2[0] == alleles_per_sample2[1]);
        REQUIRE(alleles_per_sample2[0] == alleles_per_sample2[2]);
        REQUIRE(alleles_per_sample2[0] == alleles_per_sample2[3]);
        REQUIRE(alleles_per_sample2[0] == alleles_per_sample2[4]);
        REQUIRE(alleles_per_sample2[0] == alleles_per_sample2[5]);
        REQUIRE(alleles_per_sample2[0] == alleles_per_sample2[6]);
        REQUIRE(alleles_per_sample2[0] == alleles_per_sample2[7]);
        REQUIRE(alleles_per_sample2[0] == alleles_per_sample2[8]);
        REQUIRE(alleles_per_sample2[0] != std::numeric_limits<size_t>::max());

    }

}
TEST_CASE( "Path partitioner finder bubble with three nodes",
          "[path_partitioner]" ) {

    /*
           1    
         /   \
        0--2--4
         \   /
           3

    */

    bdsg::HashGraph graph;

    //std::vector<std::string> sequences = {"AAAAAAAAAA", "A", "G", "C",  "AAAAAAAAA"};

    //std::vector<handlegraph::handle_t> nodes;
    //for (auto& seq : sequences) {
    //    nodes.emplace_back(graph.create_handle(seq));
    //}

    //graph.create_edge(nodes[0], nodes[1]);
    //graph.create_edge(nodes[0], nodes[2]);
    //graph.create_edge(nodes[0], nodes[3]);
    //graph.create_edge(nodes[1], nodes[4]);
    //graph.create_edge(nodes[2], nodes[4]);
    //graph.create_edge(nodes[3], nodes[4]);


    //// Two paths go through node 2, path 2 is associated
    //std::vector<std::vector<std::size_t>> path_seqs = { {0, 1, 4}, {0, 1, 4}, {0, 2, 4}, {0, 3, 4}};
    //std::vector<handlegraph::path_handle_t> paths;

    //for (int path_i = 0 ; path_i < path_seqs.size() ; path_i++) {
    //    paths.emplace_back(graph.create_path_handle("path"+std::to_string(path_i)));
    //    for (size_t node_i : path_seqs[path_i]) {
    //        graph.append_step(paths.back(), nodes[node_i]);
    //    }
    //}

    //// vg isn't included so the distance index can only be built from the command line
    //graph.serialize("../tests/test_data/test_graphs/simple_bubble.hg");
    //int built = system("vg index -j ../tests/test_data/test_graphs/simple_bubble.dist ../tests/test_data/test_graphs/simple_bubble.hg"); 

    graph.deserialize("../tests/test_data/test_graphs/simple_bubble.hg");
    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/test_data/test_graphs/simple_bubble.dist");

    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&graph);


    handlegraph::net_handle_t snarl = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(3)));
    std::vector<handlegraph::path_handle_t> paths;

    for (int path_i = 0 ; path_i < 4 ; path_i++) {
        paths.emplace_back(graph.get_path_handle("path"+std::to_string(path_i)));
    }

    std::vector<stoat::sample_hap_t> all_samples ({stoat::sample_hap_t(*path_graph, paths[0]),
                                         stoat::sample_hap_t(*path_graph, paths[1]),
                                         stoat::sample_hap_t(*path_graph, paths[2]),
                                         stoat::sample_hap_t(*path_graph, paths[3])});


    SECTION("partition_embedded_paths_in_snarl") {
        // This isn't really a good test because all the snarls are regular

        // Should be {0,1} {2} {3}

        std::vector<size_t> alleles_per_sample1 = partition_embedded_paths_in_snarl(*path_graph, distance_index, snarl, all_samples);
        REQUIRE(alleles_per_sample1.size() == all_samples.size());
        REQUIRE(alleles_per_sample1[0] == alleles_per_sample1[1]);
        REQUIRE(alleles_per_sample1[0] != alleles_per_sample1[2]);
        REQUIRE(alleles_per_sample1[0] != alleles_per_sample1[3]);
        REQUIRE(alleles_per_sample1[2] != alleles_per_sample1[3]);
        REQUIRE(alleles_per_sample1[0] != std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample1[2] != std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample1[3] != std::numeric_limits<size_t>::max());
    }

}
TEST_CASE( "Path partitioner finder looping snarl same edges different order ", "[path_partitioner]" ) {

    /*

             --------
            |   2    |
            \ / \    /
        0 ---1---3--4----5

    */

//stoat::Logger::instance().setLevel(stoat::LogLevel::Trace);
    bdsg::HashGraph graph;

    //std::vector<std::string> sequences = {"AAAAAAAAAA", "A", "G", "C", "T",  "AAAAAAAAA"};

    //std::vector<handlegraph::handle_t> nodes;
    //for (auto& seq : sequences) {
    //    nodes.emplace_back(graph.create_handle(seq));
    //}

    //graph.create_edge(nodes[0], nodes[1]);
    //graph.create_edge(nodes[1], nodes[2]);
    //graph.create_edge(nodes[1], nodes[3]);
    //graph.create_edge(nodes[2], nodes[3]);
    //graph.create_edge(nodes[3], nodes[4]);
    //graph.create_edge(nodes[4], nodes[1]);
    //graph.create_edge(nodes[4], nodes[5]);


    //// path 0 takes the deletion then the insertion, path 1 takes the insertion then the deletion
    //std::vector<std::vector<std::size_t>> path_seqs = { {0, 1, 3, 4, 1, 2, 3, 4, 5}, {0, 1, 2, 3, 4, 1, 3, 4, 5}};
    //std::vector<handlegraph::path_handle_t> paths;

    //for (int path_i = 0 ; path_i < path_seqs.size() ; path_i++) {
    //    paths.emplace_back(graph.create_path_handle("path"+std::to_string(path_i)));
    //    for (size_t node_i : path_seqs[path_i]) {
    //        graph.append_step(paths.back(), nodes[node_i]);
    //    }
    //}

    //// vg isn't included so the distance index can only be built from the command line
    //graph.serialize("../tests/test_data/test_graphs/loop_with_indel_two_paths.hg");
    //int built = system("vg index -j ../tests/test_data/test_graphs/loop_with_indel_two_paths.dist ../tests/test_data/test_graphs/loop_with_indel_two_paths.hg"); 

    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/test_data/test_graphs/loop_with_indel_two_paths.dist");

    graph.deserialize("../tests/test_data/test_graphs/loop_with_indel_two_paths.hg");
    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&graph);

    std::vector<handlegraph::path_handle_t> paths;

    for (int path_i = 0 ; path_i < 2 ; path_i++) {
        paths.emplace_back(graph.get_path_handle("path"+std::to_string(path_i)));
    }

    // Nested snarl
    handlegraph::net_handle_t snarl2 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(3)));
    // Duplication snarl
    handlegraph::net_handle_t snarl1 = distance_index.get_parent(distance_index.get_parent(snarl2));
    handlegraph::net_handle_t root_chain = distance_index.get_parent(snarl1);


    std::vector<stoat::sample_hap_t> all_samples ({stoat::sample_hap_t(*path_graph, paths[0]),
                                                stoat::sample_hap_t(*path_graph, paths[1])});


    SECTION("partition_embedded_paths_in_snarl") {
        // This isn't really a good test because all the snarls are regular

        // Outer snarl, hould be {0, 1}

        std::vector<size_t> alleles_per_sample1 = partition_embedded_paths_in_snarl(*path_graph, distance_index, snarl1, all_samples);
        REQUIRE(alleles_per_sample1.size() == all_samples.size());
        REQUIRE(alleles_per_sample1[0] == alleles_per_sample1[1]);
        REQUIRE(alleles_per_sample1[0] != std::numeric_limits<size_t>::max());

        // Inner snarl, should be {0} and {1}
        std::vector<size_t> alleles_per_sample2 = partition_embedded_paths_in_snarl(*path_graph, distance_index, snarl2, all_samples);
        REQUIRE(alleles_per_sample1.size() == all_samples.size());
        REQUIRE(alleles_per_sample2[0] != alleles_per_sample2[1]);
        REQUIRE(alleles_per_sample2[0] != std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample2[1] != std::numeric_limits<size_t>::max());
    }

}
TEST_CASE( "Path partitioner bubble with three nodes",
          "[path_partitioner]" ) {

    /*
           1    
         /   \
        0--2--4
         \   /
           3

    */

    bdsg::HashGraph graph;

    //std::vector<std::string> sequences = {"AAAAAAAAAA", "A", "G", "C",  "AAAAAAAAA"};

    //std::vector<handlegraph::handle_t> nodes;
    //for (auto& seq : sequences) {
    //    nodes.emplace_back(graph.create_handle(seq));
    //}

    //graph.create_edge(nodes[0], nodes[1]);
    //graph.create_edge(nodes[0], nodes[2]);
    //graph.create_edge(nodes[0], nodes[3]);
    //graph.create_edge(nodes[1], nodes[4]);
    //graph.create_edge(nodes[2], nodes[4]);
    //graph.create_edge(nodes[3], nodes[4]);


    //// Two paths go through node 2, path 2 is associated
    //std::vector<std::vector<std::size_t>> path_seqs = { {0, 1, 4}, {0, 1, 4}, {0, 2, 4}, {0, 3, 4}};
    //std::vector<handlegraph::path_handle_t> paths;

    //for (int path_i = 0 ; path_i < path_seqs.size() ; path_i++) {
    //    paths.emplace_back(graph.create_path_handle("path"+std::to_string(path_i)));
    //    for (size_t node_i : path_seqs[path_i]) {
    //        graph.append_step(paths.back(), nodes[node_i]);
    //    }
    //}

    //// vg isn't included so the distance index can only be built from the command line
    //graph.serialize("../tests/test_data/test_graphs/simple_bubble.hg");
    //int built = system("vg index -j ../tests/test_data/test_graphs/simple_bubble.dist ../tests/test_data/test_graphs/simple_bubble.hg"); 

    graph.deserialize("../tests/test_data/test_graphs/simple_bubble.hg");
    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/test_data/test_graphs/simple_bubble.dist");

    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&graph);


    handlegraph::net_handle_t snarl = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(3)));
    std::vector<handlegraph::path_handle_t> paths;

    for (int path_i = 0 ; path_i < 4 ; path_i++) {
        paths.emplace_back(graph.get_path_handle("path"+std::to_string(path_i)));
    }

    std::vector<stoat::sample_hap_t> all_samples ({stoat::sample_hap_t(*path_graph, paths[0]),
                                         stoat::sample_hap_t(*path_graph, paths[1]),
                                         stoat::sample_hap_t(*path_graph, paths[2]),
                                         stoat::sample_hap_t(*path_graph, paths[3])});


    SECTION("partition_embedded_paths_in_snarl") {

        // Should be {0,1} {2} {3}
        std::vector<size_t> alleles_per_sample1 = partition_embedded_paths_in_snarl(*path_graph, distance_index, snarl, all_samples);
        REQUIRE(alleles_per_sample1.size() == all_samples.size());
        REQUIRE(alleles_per_sample1[0] == alleles_per_sample1[1]);
        REQUIRE(alleles_per_sample1[0] != alleles_per_sample1[2]);
        REQUIRE(alleles_per_sample1[0] != alleles_per_sample1[3]);
        REQUIRE(alleles_per_sample1[2] != alleles_per_sample1[3]);
        REQUIRE(alleles_per_sample1[0] != std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample1[2] != std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample1[3] != std::numeric_limits<size_t>::max());
    }

}
TEST_CASE( "Path partitioner nested bubbles with path fragments",
          "[path_partitioner]" ) {

    /*
                       5
                     /   \
            1       4 ----6    8
          /   \   /         \ / \
        0       3  ----------7---9
          \   /
            2

   */

   //This uses the simple_nested_chain distance index but reubilds the graph with different paths 

    bdsg::HashGraph graph;

    std::vector<std::string> sequences = { "C", "C", "C", "A", "T", "C", "A", "C", "A", "A"};

    std::vector<handlegraph::handle_t> nodes;
    for (auto& seq : sequences) {
        nodes.emplace_back(graph.create_handle(seq));
    }

    graph.create_edge(nodes[0], nodes[1]);
    graph.create_edge(nodes[0], nodes[2]);
    graph.create_edge(nodes[1], nodes[3]);
    graph.create_edge(nodes[2], nodes[3]);
    graph.create_edge(nodes[3], nodes[4]);
    graph.create_edge(nodes[3], nodes[7]);
    graph.create_edge(nodes[4], nodes[5]);
    graph.create_edge(nodes[4], nodes[6]);
    graph.create_edge(nodes[5], nodes[6]);
    graph.create_edge(nodes[6], nodes[7]);
    graph.create_edge(nodes[7], nodes[8]);
    graph.create_edge(nodes[7], nodes[9]);
    graph.create_edge(nodes[8], nodes[9]);

    std::vector<std::vector<std::size_t>> paths_seqs = { {0, 1, 3, 7, 8, 9},  {0, 2, 3, 4, 6, 7}};
    std::vector<handlegraph::path_handle_t> paths;

    // Reference taking insertion
    paths.emplace_back(graph.create_path_handle("path0#0#path0"));
    for (size_t node_i : paths_seqs[1]) {
        graph.append_step(paths.back(), nodes[node_i]);
    }

    // Path 1, hap0, two fragments (loci?) going through the deletion 
    paths.emplace_back(graph.create_path_handle("path1#0#path1_0#0"));
    for (size_t node_i : paths_seqs[0]) {
        graph.append_step(paths.back(), nodes[node_i]);
    }
    paths.emplace_back(graph.create_path_handle("path1#0#path1_1#0"));
    for (size_t node_i : paths_seqs[0]) {
        graph.append_step(paths.back(), nodes[node_i]);
    }

    // Path 2, hap0, two fragments (loci?) going through the deletion 
    paths.emplace_back(graph.create_path_handle("path2#0#path2_0#0"));
    for (size_t node_i : paths_seqs[1]) {
        graph.append_step(paths.back(), nodes[node_i]);
    }
    paths.emplace_back(graph.create_path_handle("path2#0#path2_1#0"));
    for (size_t node_i : paths_seqs[0]) {
        graph.append_step(paths.back(), nodes[node_i]);
    }


    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/test_data/test_graphs/simple_nested_chain.dist");


    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&graph);


    handlegraph::net_handle_t snarl1 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(2)));
    handlegraph::net_handle_t snarl2 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(5)));
    handlegraph::net_handle_t snarl3 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(6)));
    handlegraph::net_handle_t snarl4 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(9)));
    handlegraph::net_handle_t root_chain = distance_index.get_parent(snarl1);
    handlegraph::net_handle_t nested_chain = distance_index.get_parent(snarl3);

    std::vector<stoat::sample_hap_t> all_samples ({stoat::sample_hap_t(*path_graph, paths[0]),
                                         stoat::sample_hap_t(*path_graph, paths[1]),
                                         stoat::sample_hap_t(*path_graph, paths[2]),
                                         stoat::sample_hap_t(*path_graph, paths[3]),
                                         stoat::sample_hap_t(*path_graph, paths[4])});



    SECTION("Snarl with multiple fragments") {

        // Should be {0, 3}, {1,2,4}

        std::vector<size_t> alleles_per_sample1 = partition_embedded_paths_in_snarl(*path_graph, distance_index, snarl2, all_samples);
        REQUIRE(alleles_per_sample1.size() == all_samples.size());
        REQUIRE(alleles_per_sample1[0] == alleles_per_sample1[3]);
        REQUIRE(alleles_per_sample1[0] != alleles_per_sample1[1]);
        REQUIRE(alleles_per_sample1[1] == alleles_per_sample1[2]);
        REQUIRE(alleles_per_sample1[1] == alleles_per_sample1[4]);
        REQUIRE(alleles_per_sample1[0] != std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample1[1] != std::numeric_limits<size_t>::max());
    }
}
TEST_CASE( "Path partitioner doesn't go through snarl bounds",
          "[path_partitioner]" ) {

    /*
                       5
                     /   \
            1       4 ----6    8
          /   \   /         \ / \
        0       3  ----------7---9
          \   /
            2

   */

   //This uses the simple_nested_chain distance index but reubilds the graph with different paths 

    bdsg::HashGraph graph;

    std::vector<std::string> sequences = { "C", "C", "C", "A", "T", "C", "A", "C", "A", "A"};

    std::vector<handlegraph::handle_t> nodes;
    for (auto& seq : sequences) {
        nodes.emplace_back(graph.create_handle(seq));
    }

    graph.create_edge(nodes[0], nodes[1]);
    graph.create_edge(nodes[0], nodes[2]);
    graph.create_edge(nodes[1], nodes[3]);
    graph.create_edge(nodes[2], nodes[3]);
    graph.create_edge(nodes[3], nodes[4]);
    graph.create_edge(nodes[3], nodes[7]);
    graph.create_edge(nodes[4], nodes[5]);
    graph.create_edge(nodes[4], nodes[6]);
    graph.create_edge(nodes[5], nodes[6]);
    graph.create_edge(nodes[6], nodes[7]);
    graph.create_edge(nodes[7], nodes[8]);
    graph.create_edge(nodes[7], nodes[9]);
    graph.create_edge(nodes[8], nodes[9]);

    std::vector<std::vector<std::size_t>> paths_seqs = { {0, 1, 3, 7, 8, 9},  {4, 6},  {3, 4, 6},  {4, 6, 7}};
    std::vector<handlegraph::path_handle_t> paths;

    // Reference taking insertion
    paths.emplace_back(graph.create_path_handle("path0#0#path0"));
    for (size_t node_i : paths_seqs[0]) {
        graph.append_step(paths.back(), nodes[node_i]);
    }

    // Path 1, hap0, two fragments (loci?) going through the deletion 
    paths.emplace_back(graph.create_path_handle("path1#0#0#0"));
    for (size_t node_i : paths_seqs[1]) {
        graph.append_step(paths.back(), nodes[node_i]);
    }
    paths.emplace_back(graph.create_path_handle("path2#0#0#0"));
    for (size_t node_i : paths_seqs[2]) {
        graph.append_step(paths.back(), nodes[node_i]);
    }

    // Path 2, hap0, two fragments (loci?) going through the deletion 
    paths.emplace_back(graph.create_path_handle("path3#0#0#0"));
    for (size_t node_i : paths_seqs[3]) {
        graph.append_step(paths.back(), nodes[node_i]);
    }


    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/test_data/test_graphs/simple_nested_chain.dist");


    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&graph);


    handlegraph::net_handle_t snarl1 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(2)));
    handlegraph::net_handle_t snarl2 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(5)));
    handlegraph::net_handle_t snarl3 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(6)));
    handlegraph::net_handle_t snarl4 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(9)));
    handlegraph::net_handle_t root_chain = distance_index.get_parent(snarl1);
    handlegraph::net_handle_t nested_chain = distance_index.get_parent(snarl3);

    std::vector<stoat::sample_hap_t> all_samples ({stoat::sample_hap_t(*path_graph, paths[0]),
                                         stoat::sample_hap_t(*path_graph, paths[1]),
                                         stoat::sample_hap_t(*path_graph, paths[2]),
                                         stoat::sample_hap_t(*path_graph, paths[3])});



    SECTION("Snarl with multiple fragments") {

        // Should be {0}

        std::vector<size_t> alleles_per_sample1 = partition_embedded_paths_in_snarl(*path_graph, distance_index, snarl2, all_samples);
        REQUIRE(alleles_per_sample1.size() == all_samples.size());
        REQUIRE(alleles_per_sample1[0] != std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample1[1] == std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample1[2] == std::numeric_limits<size_t>::max());
        REQUIRE(alleles_per_sample1[3] == std::numeric_limits<size_t>::max());
    }
}
TEST_CASE( "Path partitioner rejoining paths",
          "[path_partitioner]" ) {

    /*
                 
                 
                      4    
                    /   \
        -----------3      6-----------
       /         /  \   /   \          \
      1         /     5      \          8
       \       /              \        /
        ------2 ----------------7------

   */

    //bdsg::HashGraph graph;

    //std::vector<std::string> sequences = { "CCCCCCCCCCCC", "C", "T", "C", "A", "G", "T", "AAAAAAAAAAAAAAAAAAAAA"};

    //std::vector<handlegraph::handle_t> nodes;
    //for (auto& seq : sequences) {
    //    nodes.emplace_back(graph.create_handle(seq));
    //}

    //// Add the -1 to make them 1-offset because it's easier
    //graph.create_edge(nodes[1-1], nodes[2-1]);
    //graph.create_edge(nodes[1-1], nodes[3-1]);
    //graph.create_edge(nodes[2-1], nodes[3-1]);
    //graph.create_edge(nodes[2-1], nodes[7-1]);
    //graph.create_edge(nodes[3-1], nodes[4-1]);
    //graph.create_edge(nodes[3-1], nodes[5-1]);
    //graph.create_edge(nodes[4-1], nodes[6-1]);
    //graph.create_edge(nodes[5-1], nodes[6-1]);
    //graph.create_edge(nodes[6-1], nodes[7-1]);
    //graph.create_edge(nodes[6-1], nodes[8-1]);
    //graph.create_edge(nodes[7-1], nodes[8-1]);

    //// These paths are the same going into the snarl, then the split in the snarl, and must be rejoined after the snarl but split right after
    //std::vector<std::vector<std::size_t>> paths_seqs = {{0, 1, 2, 4, 5, 6, 7}, {0, 1, 2, 3, 5, 7}};
    //std::vector<handlegraph::path_handle_t> paths;

    //for (int path_i = 0 ; path_i < paths_seqs.size() ; path_i++) {
    //    paths.emplace_back(graph.create_path_handle("path"+std::to_string(path_i)));
    //    for (size_t node_i : paths_seqs[path_i]) {
    //        graph.append_step(paths.back(), nodes[node_i]);
    //    }
    //}

    //// vg isn't included so the distance index can only be built from the command line
    //graph.serialize("../tests/test_data/test_graphs/split_paths.hg");
    //int built = system("vg index -j ../tests/test_data/test_graphs/split_paths.dist ../tests/test_data/test_graphs/split_paths.hg"); 
    //built = system("vg gbwt -x ../tests/test_data/test_graphs/split_paths.hg -E --gbz-format -g ../tests/test_data/test_graphs/split_paths.gbz "); 



    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/test_data/test_graphs/split_paths.dist");

    GBZGraph gbz;
    std::ifstream instream;
    instream.open("../tests/test_data/test_graphs/split_paths.gbz");
    gbz.gbz.simple_sds_load(instream);
    instream.close();

    gbwt::GBWT* gbwt = &gbz.gbz.index;

    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&gbz);

    std::vector<handlegraph::path_handle_t> paths;
    //paths.clear();

    paths.emplace_back(gbz.get_path_handle("path0"));
    paths.emplace_back(gbz.get_path_handle("path1"));
    paths.emplace_back(gbz.get_path_handle("path2"));
    paths.emplace_back(gbz.get_path_handle("path3"));



    handlegraph::net_handle_t snarl1 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(2)));

    std::vector<stoat::sample_hap_t> all_samples({stoat::sample_hap_t(*path_graph, paths[0]),
                                                  stoat::sample_hap_t(*path_graph, paths[1])});

    SECTION("partition_embedded_paths_in_snarl") {

        // Should be {0,1} and {2,3}
        std::vector<size_t> alleles_per_sample1 = partition_embedded_paths_in_snarl(*path_graph, distance_index, snarl1, all_samples);
        REQUIRE(alleles_per_sample1.size() == all_samples.size());
        REQUIRE(alleles_per_sample1[0] != alleles_per_sample1[1]);
    }

    SECTION("partition_embedded_paths_in_snarl_with_gbwt") {

        // Should be {0,1} and {2,3}
        std::vector<PathTraversal> paths_per_allele1;
        std::vector<size_t> alleles_per_sample1 = partition_embedded_paths_in_snarl_with_gbwt(*path_graph, *gbwt, distance_index, snarl1,
                                                         all_samples, paths_per_allele1);
        REQUIRE(alleles_per_sample1.size() == all_samples.size());
        REQUIRE(alleles_per_sample1[0] != alleles_per_sample1[1]);
    }
}

