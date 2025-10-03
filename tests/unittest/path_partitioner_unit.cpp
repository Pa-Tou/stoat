#include <catch.hpp>
#include <bdsg/hash_graph.hpp>
#include <bdsg/overlays/overlay_helper.hpp>
#include "../../src/partitioner.hpp"
#include "../../src/log.hpp"

namespace stoat_graph{


class TestSnarlTraverserAndPathPartitioner : SnarlTraverserAndPathPartitioner {
    public: 
    TestSnarlTraverserAndPathPartitioner(const std::set<stoat::sample_hap_t>& all_sample_haplotypes, const bdsg::SnarlDistanceIndex* distance_index, const std::string& reference_sample, bool     save_partitions) :
        SnarlTraverserAndPathPartitioner(all_sample_haplotypes, distance_index, reference_sample, save_partitions) {} 
    using SnarlTraverserAndPathPartitioner::partition_samples_in_snarl;
    using SnarlTraverserAndPathPartitioner::get_walk_sets;
    using SnarlTraverserAndPathPartitioner::for_each_snarl_partition;
    using SnarlTraverserAndPathPartitioner::serialize;
    using SnarlTraverserAndPathPartitioner::deserialize;
};

TEST_CASE( "Path partitioner finder one node", "[path_partitioner]" ) {


    bdsg::HashGraph graph;
        
    //handlegraph::handle_t n1 = graph.create_handle("GCAAACAGATT");

    //handlegraph::path_handle_t path = graph.create_path_handle("path");
    //graph.append_step(path, n1);

    // vg isn't included so the distance index can only be built from the command line
    //graph.serialize("../tests/graph_test/one_node.hg");
    //int built = system("vg index -j ../tests/graph_test/one_node.dist ../tests/graph_test/one_node.hg"); 
    bdsg::SnarlDistanceIndex distance_index;
    graph.deserialize("../tests/graph_test/one_node.hg");
    distance_index.deserialize("../tests/graph_test/one_node.dist");

    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&graph);

    std::set<stoat::sample_hap_t> all_samples ({stoat::get_sample_and_haplotype(*path_graph,graph.get_path_handle("path"))});

    SECTION("Make partitioner") {
        // There isn't much to do with one node so just make sure we can run the constructor without crashing
        TestSnarlTraverserAndPathPartitioner af(all_samples, &distance_index, "path0", false);
    }
    SECTION("Serialize partitioner") {
        // There isn't much to do with one node so just make sure we can run the constructor without crashing
        TestSnarlTraverserAndPathPartitioner af(all_samples, &distance_index, "path0", true);
        af.serialize("./test.snarl_partitions.txt");
        
        TestSnarlTraverserAndPathPartitioner af_loaded(all_samples, nullptr, "path0", false);
        af_loaded.deserialize("./test.snarl_partitions.txt", *path_graph);

        int rm = system("rm ./test.snarl_partitions.txt"); 
    }

}

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
    //graph.serialize("../tests/graph_test/simple_nested_chain.hg");
    //int built = system("vg index -j ../tests/graph_test/simple_nested_chain.dist ../tests/graph_test/simple_nested_chain.hg"); 
    //
    ////Change sense of paths
    //built = system("vg convert --hap-locus path0 --new-sample path0 ../tests/graph_test/simple_nested_chain.hg >../tests/graph_test/simple_nested_chain1.hg"); 
    //built = system("vg convert --hap-locus path1 --new-sample path1 ../tests/graph_test/simple_nested_chain1.hg >../tests/graph_test/simple_nested_chain2.hg"); 
    //built = system("vg convert --ref-sample path0 ../tests/graph_test/simple_nested_chain2.hg | vg convert -a - >../tests/graph_test/simple_nested_chain3.hg"); 
    //built = system("mv ../tests/graph_test/simple_nested_chain3.hg ../tests/graph_test/simple_nested_chain.hg"); 
    //built = system("rm ../tests/graph_test/simple_nested_chain1.hg"); 
    //built = system("rm ../tests/graph_test/simple_nested_chain2.hg"); 
    //built = system("vg gbwt -x ../tests/graph_test/simple_nested_chain.hg -E --gbz-format -g ../tests/graph_test/simple_nested_chain.gbz "); 



    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/graph_test/simple_nested_chain.dist");

    bdsg::HashGraph graph;
    graph.deserialize("../tests/graph_test/simple_nested_chain.hg");

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

    // snarl3 should be associated
    std::set<std::string> samples ({"path1", "path3"});
    std::set<stoat::sample_hap_t> all_samples ({stoat::get_sample_and_haplotype(*path_graph, paths[0]),
                                         stoat::get_sample_and_haplotype(*path_graph, paths[1]),
                                         stoat::get_sample_and_haplotype(*path_graph, paths[2]),
                                         stoat::get_sample_and_haplotype(*path_graph, paths[3])});

    TestSnarlTraverserAndPathPartitioner af(all_samples, &distance_index, "path0", false);


    SECTION("get_walk_set") {
        // This isn't really a good test because all the snarls are regular

        // Should be {0,1} and {2,3}
        std::vector<std::set<stoat::sample_hap_t>> walks1 = af.get_walk_sets(*path_graph, snarl1, false);
        REQUIRE(walks1.size() == 2);
        for ( const auto& walk_set : walks1) {
            REQUIRE(walk_set.size() == 2);
            REQUIRE( ((walk_set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[0]), stoat::get_sample_and_haplotype(*path_graph, paths[1])})) || 
                     (walk_set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[2]), stoat::get_sample_and_haplotype(*path_graph, paths[3])}))));
        }

        // Should be {0,1,3} and {2}
        std::vector<std::set<stoat::sample_hap_t>> walks2 = af.get_walk_sets(*path_graph, snarl2, false);
        REQUIRE(walks2.size() == 2);
        for ( const auto& set : walks2) {
            REQUIRE(((set.size() == 3) || (set.size() == 1)));
            REQUIRE( ((set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[0]), stoat::get_sample_and_haplotype(*path_graph, paths[1]), stoat::get_sample_and_haplotype(*path_graph, paths[3])})) || 
                      (set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[2])}))));
        }

        // Should be {0}, {1,3}. 2 didn't go through this snarl
        std::vector<std::set<stoat::sample_hap_t>> walks3 = af.get_walk_sets(*path_graph, snarl3, false);
        REQUIRE(walks3.size() == 2);
        for ( const auto& set : walks3) {
            REQUIRE(((set.size() == 2) || (set.size() == 1)));
            REQUIRE( ((set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[0])}) ) ||
                      (set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[1]), stoat::get_sample_and_haplotype(*path_graph, paths[3])}))));
        }
    }
    SECTION("get start edge sets") {
        // Should be {0,1} and {2,3}
        std::vector<std::set<stoat::sample_hap_t>> edges1 = af.get_walk_sets(*path_graph, snarl1, true);
        REQUIRE(edges1.size() == 2);
        for ( const auto& set : edges1) {
            REQUIRE(set.size() == 2);
            REQUIRE( ((set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[0]), stoat::get_sample_and_haplotype(*path_graph, paths[1])})) || 
                     (set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[2]), stoat::get_sample_and_haplotype(*path_graph, paths[3])}))));
        }

        // Should be {0,1,3} and {2}
        std::vector<std::set<stoat::sample_hap_t>> edges2 = af.get_walk_sets(*path_graph, snarl2, true);
        REQUIRE(edges2.size() == 2);
        for ( const auto& set : edges2) {
            REQUIRE(((set.size() == 3) || (set.size() == 1)));
            REQUIRE( ((set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[0]), stoat::get_sample_and_haplotype(*path_graph, paths[1]), stoat::get_sample_and_haplotype(*path_graph, paths[3])})) || 
                      (set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[2])}))));
        }

        // Should be {0} and {1,3}
        std::vector<std::set<stoat::sample_hap_t>> edges3 = af.get_walk_sets(*path_graph, snarl3, true);
        REQUIRE(edges3.size() == 2);
        for ( const auto& set : edges3) {
            REQUIRE(((set.size() == 2) || (set.size() == 1)));
            REQUIRE( ((set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[0])}) ) ||
                      (set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[1]), stoat::get_sample_and_haplotype(*path_graph, paths[3])}))));
        }
    }
    SECTION ("Traverse snarls" ) {
        // Get all the snarl partitions
        std::vector<stoat::snarl_partition_t> partitions;
        af.for_each_snarl_partition(*path_graph, [&](const bdsg::SnarlDistanceIndex& dist_index, const handlegraph::net_handle_t& net) {return true;}, [&] (const stoat::snarl_partition_t& snarl_info) {
            partitions.emplace_back(snarl_info);
        });
        REQUIRE(partitions.size() == 4);

        // Make the set of true partitions of snarls
        std::vector<stoat::snarl_partition_t>truth_partitions;
        std::vector<std::set<stoat::sample_hap_t>> temp_partitions;
        temp_partitions.emplace_back();
        temp_partitions.back().emplace("path0", 0);
        temp_partitions.back().emplace("path1", 0);
        temp_partitions.emplace_back();
        temp_partitions.back().emplace("path2", std::numeric_limits<size_t>::max());
        temp_partitions.back().emplace("path3", std::numeric_limits<size_t>::max());
        truth_partitions.emplace_back(snarl1, graph.get_handle(1, false), graph.get_handle(4, true), std::make_pair<size_t, size_t> (1,4), 1, 2, 1, 1, 1, "path0#0#path0", temp_partitions); 

        temp_partitions.clear();
        temp_partitions.emplace_back();
        temp_partitions.back().emplace("path2", std::numeric_limits<size_t>::max());
        temp_partitions.emplace_back();
        temp_partitions.back().emplace("path0", 0);
        temp_partitions.back().emplace("path1", 0);
        temp_partitions.back().emplace("path3", std::numeric_limits<size_t>::max());
        truth_partitions.emplace_back(snarl2, graph.get_handle(4, false), graph.get_handle(8, true), std::make_pair<size_t, size_t> (4,8), 3, 6, 1, 0, 3, "path0#0#path0", temp_partitions); 

        temp_partitions.clear();
        temp_partitions.emplace_back();
        temp_partitions.back().emplace("path0", 0);
        temp_partitions.emplace_back();
        temp_partitions.back().emplace("path1", 0);
        temp_partitions.back().emplace("path3", std::numeric_limits<size_t>::max());
        truth_partitions.emplace_back(snarl3, graph.get_handle(5, false), graph.get_handle(7, true), std::make_pair<size_t, size_t> (5,7), 4, 5, 2, 0, 1, "path0#0#path0", temp_partitions); 

        temp_partitions.clear();
        truth_partitions.emplace_back(snarl4, graph.get_handle(8, false), graph.get_handle(10, true), std::make_pair<size_t, size_t> (8, 10), 0, 0, 1, 0, 1, "NA", temp_partitions); 

        for (const stoat::snarl_partition_t& test_partition : partitions) {
            stoat::snarl_partition_t& truth_partition = truth_partitions[0];
            if (distance_index.start_end_traversal_of(test_partition.snarl) == snarl1) {
                truth_partition = truth_partitions[0];
            } else if (distance_index.start_end_traversal_of(test_partition.snarl) == snarl2) {
                truth_partition = truth_partitions[1];
            } else if (distance_index.start_end_traversal_of(test_partition.snarl) == snarl3) {
                truth_partition = truth_partitions[2];
            } else {
                REQUIRE(distance_index.start_end_traversal_of(test_partition.snarl) == snarl4);
                truth_partition = truth_partitions[3];
           }

           REQUIRE(handlegraph::as_integer(test_partition.snarl) == handlegraph::as_integer(truth_partition.snarl));
           REQUIRE(test_partition.start_handle == truth_partition.start_handle);
           REQUIRE(test_partition.end_handle == truth_partition.end_handle);
           REQUIRE(test_partition.ref_path == truth_partition.ref_path);
           REQUIRE(test_partition.start_positions == truth_partition.start_positions);
           REQUIRE(test_partition.end_positions == truth_partition.end_positions);
           REQUIRE(test_partition.depth == truth_partition.depth);
           REQUIRE(test_partition.min_length == truth_partition.min_length);
           REQUIRE(test_partition.max_length == truth_partition.max_length);
           
           REQUIRE(test_partition.partitions.size() == truth_partition.partitions.size());
           if (truth_partition.partitions.size() != 0) {

                REQUIRE((test_partition.partitions[0] == truth_partition.partitions[0] || 
                         test_partition.partitions[0] == truth_partition.partitions[1]));
                REQUIRE((test_partition.partitions[1] == truth_partition.partitions[0] || 
                         test_partition.partitions[1] == truth_partition.partitions[1]));
           }
        }
    }
    SECTION("Serialize partitioner") {
        TestSnarlTraverserAndPathPartitioner af_serialized(all_samples, &distance_index, "path0", true);


        // Get all the snarl partitions
        std::vector<stoat::snarl_partition_t> serialized_partitions;
        af_serialized.for_each_snarl_partition(*path_graph, [&](const bdsg::SnarlDistanceIndex& dist_index, const handlegraph::net_handle_t& net) {return true;}, [&] (const stoat::snarl_partition_t& snarl_info) {
            serialized_partitions.emplace_back(snarl_info);
        });

        //Serialize it
        af_serialized.serialize("./test.snarl_partitions.txt");


        TestSnarlTraverserAndPathPartitioner af_deserialized(all_samples, nullptr, "path0", false);
        af_deserialized.deserialize("./test.snarl_partitions.txt", *path_graph);

        std::vector<stoat::snarl_partition_t> deserialized_partitions;
        af_deserialized.for_each_snarl_partition(*path_graph, [&](const bdsg::SnarlDistanceIndex& dist_index, const handlegraph::net_handle_t& net) {return true;}, [&] (const stoat::snarl_partition_t& snarl_info) {
            deserialized_partitions.emplace_back(snarl_info);
        });

        REQUIRE(serialized_partitions.size() == deserialized_partitions.size());

        // The serialized and deserialized snarls should be the same
        // The order isn't required to be the same but it will be
        for (size_t i = 0 ; i < serialized_partitions.size() ; i++ ) {
            REQUIRE(serialized_partitions[i].partitions == deserialized_partitions[i].partitions);
            REQUIRE(serialized_partitions[i].start_positions == deserialized_partitions[i].start_positions);
            REQUIRE(serialized_partitions[i].end_positions == deserialized_partitions[i].end_positions);
            REQUIRE(serialized_partitions[i].start_handle == deserialized_partitions[i].start_handle);
            REQUIRE(serialized_partitions[i].end_handle == deserialized_partitions[i].end_handle);
            REQUIRE(serialized_partitions[i].min_length == deserialized_partitions[i].min_length);
            REQUIRE(serialized_partitions[i].max_length == deserialized_partitions[i].max_length);
            REQUIRE(serialized_partitions[i].depth == deserialized_partitions[i].depth);
            REQUIRE(serialized_partitions[i].ref_path == deserialized_partitions[i].ref_path);
        }

        int rm = system("rm ./test.snarl_partitions.txt"); 
    }

}
TEST_CASE( "Path partitioner nested bubbles distanceless index",
          "[path_partitioner][bug]" ) {

    /*
                       5
                     /   \
            1       4 ----6    8
          /   \   /         \ / \
        0       3  ----------7---9
          \   /
            2

   */

    //int built = system("vg index --snarl-limit 0 -j ../tests/graph_test/simple_nested_chain.nodist.dist ../tests/graph_test/simple_nested_chain.hg"); 
    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/graph_test/simple_nested_chain.nodist.dist");

    bdsg::HashGraph graph;
    graph.deserialize("../tests/graph_test/simple_nested_chain.hg");

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

    // snarl3 should be associated
    std::set<std::string> samples ({"path1", "path3"});
    std::set<stoat::sample_hap_t> all_samples ({stoat::get_sample_and_haplotype(*path_graph, paths[0]),
                                         stoat::get_sample_and_haplotype(*path_graph, paths[1]),
                                         stoat::get_sample_and_haplotype(*path_graph, paths[2]),
                                         stoat::get_sample_and_haplotype(*path_graph, paths[3])});

    TestSnarlTraverserAndPathPartitioner af(all_samples, &distance_index, "path0", false);


    SECTION("get_walk_set") {
        // This isn't really a good test because all the snarls are regular

        // Should be {0,1} and {2,3}
        std::vector<std::set<stoat::sample_hap_t>> walks1 = af.get_walk_sets(*path_graph, snarl1, false);
        REQUIRE(walks1.size() == 2);
        for ( const auto& walk_set : walks1) {
            REQUIRE(walk_set.size() == 2);
            REQUIRE( ((walk_set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[0]), stoat::get_sample_and_haplotype(*path_graph, paths[1])})) || 
                     (walk_set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[2]), stoat::get_sample_and_haplotype(*path_graph, paths[3])}))));
        }

        // Should be {0,1,3} and {2}
        std::vector<std::set<stoat::sample_hap_t>> walks2 = af.get_walk_sets(*path_graph, snarl2, false);
        REQUIRE(walks2.size() == 2);
        for ( const auto& set : walks2) {
            REQUIRE(((set.size() == 3) || (set.size() == 1)));
            REQUIRE( ((set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[0]), stoat::get_sample_and_haplotype(*path_graph, paths[1]), stoat::get_sample_and_haplotype(*path_graph, paths[3])})) || 
                      (set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[2])}))));
        }
        std::vector<std::set<stoat::sample_hap_t>> partitions2 = af.partition_samples_in_snarl(*path_graph, snarl2);
        REQUIRE(partitions2.size() == 2);
        for ( const auto& set : partitions2) {
            REQUIRE(((set.size() == 3) || (set.size() == 1)));
        }

        // Should be {0}, {1,3}. 2 didn't go through this snarl
        std::vector<std::set<stoat::sample_hap_t>> walks3 = af.get_walk_sets(*path_graph, snarl3, false);
        REQUIRE(walks3.size() == 2);
        for ( const auto& set : walks3) {
            REQUIRE(((set.size() == 2) || (set.size() == 1)));
            REQUIRE( ((set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[0])}) ) ||
                      (set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[1]), stoat::get_sample_and_haplotype(*path_graph, paths[3])}))));
        }
    }
    SECTION("get start edge sets") {
        // Should be {0,1} and {2,3}
        std::vector<std::set<stoat::sample_hap_t>> edges1 = af.get_walk_sets(*path_graph, snarl1, true);
        REQUIRE(edges1.size() == 2);
        for ( const auto& set : edges1) {
            REQUIRE(set.size() == 2);
            REQUIRE( ((set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[0]), stoat::get_sample_and_haplotype(*path_graph, paths[1])})) || 
                     (set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[2]), stoat::get_sample_and_haplotype(*path_graph, paths[3])}))));
        }

        // Should be {0,1,3} and {2}
        std::vector<std::set<stoat::sample_hap_t>> edges2 = af.get_walk_sets(*path_graph, snarl2, true);
        REQUIRE(edges2.size() == 2);
        for ( const auto& set : edges2) {
            REQUIRE(((set.size() == 3) || (set.size() == 1)));
            REQUIRE( ((set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[0]), stoat::get_sample_and_haplotype(*path_graph, paths[1]), stoat::get_sample_and_haplotype(*path_graph, paths[3])})) || 
                      (set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[2])}))));
        }

        // Should be {0} and {1,3}
        std::vector<std::set<stoat::sample_hap_t>> edges3 = af.get_walk_sets(*path_graph, snarl3, true);
        REQUIRE(edges3.size() == 2);
        for ( const auto& set : edges3) {
            REQUIRE(((set.size() == 2) || (set.size() == 1)));
            REQUIRE( ((set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[0])}) ) ||
                      (set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[1]), stoat::get_sample_and_haplotype(*path_graph, paths[3])}))));
        }
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
    //graph.serialize("../tests/graph_test/loop_with_indel.hg");
    //int built = system("vg index -j ../tests/graph_test/loop_with_indel.dist ../tests/graph_test/loop_with_indel.hg"); 

    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/graph_test/loop_with_indel.dist");

    graph.deserialize("../tests/graph_test/loop_with_indel.hg");
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


    std::set<std::string> samples ({"path1", "path2"});
    std::set<stoat::sample_hap_t> all_samples ({stoat::get_sample_and_haplotype(*path_graph, paths[0]),
                                         stoat::get_sample_and_haplotype(*path_graph, paths[1]),
                                         stoat::get_sample_and_haplotype(*path_graph, paths[2])});
    TestSnarlTraverserAndPathPartitioner af(all_samples, &distance_index, "path0", false);

    SECTION("get_walk_set") {
        // This isn't really a good test because all the snarls are regular

        // Should be {0} and {1,2}
        std::vector<std::set<stoat::sample_hap_t>> walks1 = af.get_walk_sets(*path_graph, snarl1, false);
        REQUIRE(walks1.size() == 2);
        for ( const auto& set : walks1) {
            REQUIRE( ((set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[1]), stoat::get_sample_and_haplotype(*path_graph, paths[2])})) || 
                     (set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[0])}))));
        }

        // Should be {0}, {1} and {2}
        std::vector<std::set<stoat::sample_hap_t>> walks2 = af.get_walk_sets(*path_graph, snarl2, false);
        REQUIRE(walks2.size() == 3);
    }
    SECTION("only check start edge") {

        // Should be {0,2} and {1,2}
        std::vector<std::set<stoat::sample_hap_t>> edges2 = af.get_walk_sets(*path_graph, snarl2, true);
        REQUIRE(edges2.size() == 3);
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
    //graph.serialize("../tests/graph_test/simple_bubble.hg");
    //int built = system("vg index -j ../tests/graph_test/simple_bubble.dist ../tests/graph_test/simple_bubble.hg"); 

    graph.deserialize("../tests/graph_test/simple_bubble.hg");
    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/graph_test/simple_bubble.dist");

    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&graph);


    handlegraph::net_handle_t snarl = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(3)));
    std::vector<handlegraph::path_handle_t> paths;

    for (int path_i = 0 ; path_i < 4 ; path_i++) {
        paths.emplace_back(graph.get_path_handle("path"+std::to_string(path_i)));
    }


    // This file is meant to test the base association finder but since it is technically an interface with some implementations,
    // build the path version and only test the base functions
    std::set<std::string> samples ({"path2"});
    std::set<stoat::sample_hap_t> all_samples ({stoat::get_sample_and_haplotype(*path_graph, paths[0]),
                                         stoat::get_sample_and_haplotype(*path_graph, paths[1]),
                                         stoat::get_sample_and_haplotype(*path_graph, paths[2]),
                                         stoat::get_sample_and_haplotype(*path_graph, paths[3])});

    TestSnarlTraverserAndPathPartitioner af(all_samples, &distance_index, "path0", false);

    SECTION("get_walk_set") {
        // This isn't really a good test because all the snarls are regular

        // Should be {0,1} {2} {3}
        std::vector<std::set<stoat::sample_hap_t>> walks1 = af.get_walk_sets(*path_graph, snarl, false);
        REQUIRE(walks1.size() == 3);
        for ( const auto& set : walks1) {
            REQUIRE( ((set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[0]), stoat::get_sample_and_haplotype(*path_graph, paths[1])})) || 
                     (set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[2])})) || 
                     (set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[3])}))));
        }
    }
    SECTION("only check start edge") {

        // Should be {0,1} {2} {3}
        std::vector<std::set<stoat::sample_hap_t>> walks1 = af.get_walk_sets(*path_graph, snarl, false);
        REQUIRE(walks1.size() == 3);
        for ( const auto& set : walks1) {
            REQUIRE( ((set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[0]), stoat::get_sample_and_haplotype(*path_graph, paths[1])})) || 
                     (set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[2])})) || 
                     (set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[3])}))));
        }
    }
}
TEST_CASE( "Path partitioner finder looping snarl same edges different order ", "[path_partitioner][bug]" ) {

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
    //graph.serialize("../tests/graph_test/loop_with_indel_two_paths.hg");
    //int built = system("vg index -j ../tests/graph_test/loop_with_indel_two_paths.dist ../tests/graph_test/loop_with_indel_two_paths.hg"); 

    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/graph_test/loop_with_indel_two_paths.dist");

    graph.deserialize("../tests/graph_test/loop_with_indel_two_paths.hg");
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


    // This file is meant to test the base association finder but since it is technically an interface with some implementations,
    // build the path version and only test the base functions
    std::set<std::string> samples ({"path0", "path1"});
    std::set<stoat::sample_hap_t> all_samples ({stoat::get_sample_and_haplotype(*path_graph, paths[0]),
                                                stoat::get_sample_and_haplotype(*path_graph, paths[1])});
    TestSnarlTraverserAndPathPartitioner af(all_samples, &distance_index, "path0", false);


    SECTION("get_walk_set") {
        // This isn't really a good test because all the snarls are regular

        // Outer snarl, hould be {0, 1}
        std::vector<std::set<stoat::sample_hap_t>> walks1 = af.get_walk_sets(*path_graph, snarl1, false);
        REQUIRE(walks1.size() == 1);
        REQUIRE( (walks1[0] == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[0]), stoat::get_sample_and_haplotype(*path_graph, paths[1])})));

        // Inner snarl, should be {0} and {1}
        std::vector<std::set<stoat::sample_hap_t>> walks2 = af.get_walk_sets(*path_graph, snarl2, false);
        REQUIRE(walks2.size() == 2);
        for ( const auto& set : walks2) {
            REQUIRE( ((set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[0])})) || 
                      (set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[1])}))));
        }
    }

}
TEST_CASE( "Path association finder bubble with three nodes",
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
    //graph.serialize("../tests/graph_test/simple_bubble.hg");
    //int built = system("vg index -j ../tests/graph_test/simple_bubble.dist ../tests/graph_test/simple_bubble.hg"); 

    graph.deserialize("../tests/graph_test/simple_bubble.hg");
    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/graph_test/simple_bubble.dist");

    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&graph);


    handlegraph::net_handle_t snarl = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(3)));
    std::vector<handlegraph::path_handle_t> paths;

    for (int path_i = 0 ; path_i < 4 ; path_i++) {
        paths.emplace_back(graph.get_path_handle("path"+std::to_string(path_i)));
    }


    // This file is meant to test the base association finder but since it is technically an interface with some implementations,
    // build the path version and only test the base functions
    std::set<std::string> samples ({"path2"});
    std::set<stoat::sample_hap_t> all_samples ({stoat::get_sample_and_haplotype(*path_graph, paths[0]),
                                         stoat::get_sample_and_haplotype(*path_graph, paths[1]),
                                         stoat::get_sample_and_haplotype(*path_graph, paths[2]),
                                         stoat::get_sample_and_haplotype(*path_graph, paths[3])});

    TestSnarlTraverserAndPathPartitioner af(all_samples, &distance_index, "path0", false);

    SECTION("get_walk_set") {
        // This isn't really a good test because all the snarls are regular

        // Should be {0,1} {2} {3}
        std::vector<std::set<stoat::sample_hap_t>> walks1 = af.get_walk_sets(*path_graph, snarl, false);
        REQUIRE(walks1.size() == 3);
        for ( const auto& set : walks1) {
            REQUIRE( ((set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[0]), stoat::get_sample_and_haplotype(*path_graph, paths[1])})) || 
                     (set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[2])})) || 
                     (set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[3])}))));
        }
    }
    SECTION("only check start edge") {

        // Should be {0,1} {2} {3}
        std::vector<std::set<stoat::sample_hap_t>> walks1 = af.get_walk_sets(*path_graph, snarl, false);
        REQUIRE(walks1.size() == 3);
        for ( const auto& set : walks1) {
            REQUIRE( ((set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[0]), stoat::get_sample_and_haplotype(*path_graph, paths[1])})) || 
                     (set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[2])})) || 
                     (set == std::set<stoat::sample_hap_t> ({stoat::get_sample_and_haplotype(*path_graph, paths[3])}))));
        }
    }
}
}
