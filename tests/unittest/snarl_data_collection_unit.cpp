#include <catch.hpp>
#include <bdsg/hash_graph.hpp>
#include <bdsg/overlays/overlay_helper.hpp>
#include "../../src/snarl_data_collection.hpp"
#include "../../src/log.hpp"

namespace stoat{


class TestSnarlDataCollection : SnarlDataCollection {
    public: 
    TestSnarlDataCollection(size_t allele_size_limit, size_t snarl_child_limit, size_t walk_cycle_limit, size_t walk_steps_limit) :
        SnarlDataCollection(allele_size_limit, snarl_child_limit, walk_cycle_limit, walk_steps_limit) {} 
    using SnarlDataCollection::fill_in_snarl_info;
    using SnarlDataCollection::add_snarl_partitions;
    using SnarlDataCollection::for_each_snarl;
    using SnarlDataCollection::write_snarl_data_collection;
    using SnarlDataCollection::load_snarl_data_collection;
};

TEST_CASE( "Snarl collection one node", "[snarl_collection]" ) {


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

    SECTION("Make and fill in snarl collection") {
        // There isn't much to do with one node so just make sure we can run the constructor without crashing
        TestSnarlDataCollection snarl_collection(1,1,1,1);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, false, [&](const handlegraph::net_handle_t& net) { 
            std::vector<std::set<sample_hap_t>> partitions;
            return partitions;
        }, false, "", false, cout);

    }
    SECTION("Serialize partitioner") {
        // There isn't much to do with one node so just make sure we can run the constructor without crashing

        TestSnarlDataCollection snarl_collection(1,1,1,1);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, false, [&](const handlegraph::net_handle_t& net) { 
            std::vector<std::set<sample_hap_t>> partitions;
            return partitions;
        }, false, "", false, cout);

        std::string test_file = "./test.snarl_collection.txt";

        ofstream outstream;
        outstream.open(test_file);
        snarl_collection.write_snarl_data_collection(outstream);
        outstream.close();
        
        TestSnarlDataCollection snarl_collection_loaded(1,1,1,1);
        ifstream instream;
        instream.open(test_file);
        snarl_collection_loaded.load_snarl_data_collection(instream);
        instream.close();

        std::string rm_cmd = "rm " + test_file;

        int rm = system(rm_cmd.c_str()); 
    }

}

TEST_CASE( "Snarl collection nested bubbles",
          "[snarl_collection]" ) {

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


    SECTION("Make and fill in snarl collection with no data") {
        // Don't get the partitions or anything else
        TestSnarlDataCollection snarl_collection(1,10,1,10);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, false, [&](const handlegraph::net_handle_t& net) { 
            std::vector<std::set<sample_hap_t>> partitions;
            return partitions;
        }, false, "", false, cout);

        // Check that we got all snarls and that we got the correct snarls
        size_t snarl_count = 0;

        std::vector<Path_traversal_t> snarl_walks;
        std::vector<std::set<sample_hap_t>> partitions;
        std::vector<std::string> sequences;

        stoat::snarl_info_t snarl1 (stoat::Node_traversal_t(1, false), // start (or end) node 
                                    stoat::Node_traversal_t(4, true), // end (or start) node 
                                    "path0#0#path0",                   // reference path 
                                    1,                                 // start offset
                                    2,                                 // end offset
                                    1,                                 // depth
                                    "",                                // variant type (allele length counts)
                                    snarl_walks,                       // walks through the snarl (not checking this)
                                    partitions,                        // set of samples per walk (not checking this)
                                    sequences                         // sequences per walk (not checking this)
                                    );
        stoat::snarl_info_t snarl2 (stoat::Node_traversal_t(4, false), // start (or end) node 
                                    stoat::Node_traversal_t(8, true), // end (or start) node 
                                    "path0#0#path0",                   // reference path 
                                    3,                                 // start offset
                                    6,                                 // end offset
                                    1,                                 // depth
                                    "",                                // variant type (allele length counts)
                                    snarl_walks,                       // walks through the snarl (not checking this)
                                    partitions,                        // set of samples per walk (not checking this)
                                    sequences                          // sequences per walk (not checking this)
                                    );
        stoat::snarl_info_t snarl3 (stoat::Node_traversal_t(5, false), // start (or end) node 
                                    stoat::Node_traversal_t(7, true), // end (or start) node 
                                    "path0#0#path0",                   // reference path 
                                    4,                                 // start offset
                                    5,                                 // end offset
                                    2,                                 // depth
                                    "",                                // variant type (allele length counts)
                                    snarl_walks,                       // walks through the snarl (not checking this)
                                    partitions,                        // set of samples per walk (not checking this)
                                    sequences                          // sequences per walk (not checking this)
                                    );
        stoat::snarl_info_t snarl4 (stoat::Node_traversal_t(8, false), // start (or end) node 
                                    stoat::Node_traversal_t(10, true), // end (or start) node 
                                    "NA",                              // reference path 
                                    0,                                 // start offset
                                    0,                                 // end offset
                                    1,                                 // depth
                                    "",                                // variant type (allele length counts)
                                    snarl_walks,                       // walks through the snarl (not checking this)
                                    partitions,                        // set of samples per walk (not checking this)
                                    sequences                          // sequences per walk (not checking this)
                                    );
        snarl_collection.for_each_snarl([&](const snarl_info_t& snarl_info) {
            snarl_count++;

            if ((snarl_info.start_node == snarl1.start_node && snarl_info.end_node == snarl1.end_node) ||
                (snarl_info.start_node == snarl1.end_node && snarl_info.end_node == snarl1.start_node)) {
                REQUIRE(snarl_info.ref_path == snarl1.ref_path);
                REQUIRE(snarl_info.start_position == snarl1.start_position);
                REQUIRE(snarl_info.end_position == snarl1.end_position);
                REQUIRE(snarl_info.depth == snarl1.depth);
            } else if ((snarl_info.start_node == snarl2.start_node && snarl_info.end_node == snarl2.end_node) ||
                (snarl_info.start_node == snarl2.end_node && snarl_info.end_node == snarl2.start_node)) {
                REQUIRE(snarl_info.ref_path == snarl2.ref_path);
                REQUIRE(snarl_info.start_position == snarl2.start_position);
                REQUIRE(snarl_info.end_position == snarl2.end_position);
                REQUIRE(snarl_info.depth == snarl2.depth);
            }  else if ((snarl_info.start_node == snarl3.start_node && snarl_info.end_node == snarl3.end_node) ||
                (snarl_info.start_node == snarl3.end_node && snarl_info.end_node == snarl3.start_node)) {
                REQUIRE(snarl_info.ref_path == snarl3.ref_path);
                REQUIRE(snarl_info.start_position == snarl3.start_position);
                REQUIRE(snarl_info.end_position == snarl3.end_position);
                REQUIRE(snarl_info.depth == snarl3.depth);
            } else if ((snarl_info.start_node == snarl4.start_node && snarl_info.end_node == snarl4.end_node) ||
                (snarl_info.start_node == snarl4.end_node && snarl_info.end_node == snarl4.start_node)) {
                REQUIRE(snarl_info.ref_path == snarl4.ref_path);
                REQUIRE(snarl_info.start_position == snarl4.start_position);
                REQUIRE(snarl_info.end_position == snarl4.end_position);
                REQUIRE(snarl_info.depth == snarl4.depth);
            } else {
                REQUIRE(false);
            }
        });

        REQUIRE(snarl_count == 4);
    }
    SECTION("Make and fill in snarl collection with just the reference in partitions") {
        // Make the partitions but only the reference sample

        sample_hap_t ref_sample;

        path_graph->for_each_path_handle([&] (handlegraph::path_handle_t path) {
            if (graph.get_sense(path) == handlegraph::PathSense::REFERENCE) {
                ref_sample = stoat::get_sample_and_haplotype(*path_graph, path);
            }
            return true;
        });


        TestSnarlDataCollection snarl_collection(1,10,1,10);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, true, [&](const handlegraph::net_handle_t& net) { 
            std::vector<std::set<sample_hap_t>> partitions;
            partitions.emplace_back();
            partitions.back().emplace(ref_sample);
            return partitions;
        }, true, "", false, cout);

        // Check that we got all snarls and that we got the correct snarls
        size_t snarl_count = 0;

        std::vector<std::set<sample_hap_t>> partitions;
        partitions.emplace_back();
        partitions.back().emplace(ref_sample);
        std::vector<std::string> sequences;


        // Get just the walk for the reference since this only tests the partitions containing only the reference
        std::vector<Path_traversal_t> snarl_walks1;
        snarl_walks1.emplace_back();
        snarl_walks1.back().add_node_traversal_t(Node_traversal_t(2, false));
        stoat::snarl_info_t snarl1 (stoat::Node_traversal_t(1, false), // start (or end) node 
                                    stoat::Node_traversal_t(4, true), // end (or start) node 
                                    "path0#0#path0",                   // reference path 
                                    1,                                 // start offset
                                    2,                                 // end offset
                                    1,                                 // depth
                                    "1",                                // variant type (allele length counts)
                                    snarl_walks1,                       // walks through the snarl 
                                    partitions,                        // set of samples per walk 
                                    sequences                         // sequences per walk 
                                    );
        std::vector<Path_traversal_t> snarl_walks2;
        snarl_walks2.emplace_back();
        snarl_walks2.back().add_node_traversal_t(Node_traversal_t(5, false));
        snarl_walks2.back().add_node_traversal_t(Node_traversal_t(0, false));
        snarl_walks2.back().add_node_traversal_t(Node_traversal_t(7, false));
        stoat::snarl_info_t snarl2 (stoat::Node_traversal_t(4, false), // start (or end) node 
                                    stoat::Node_traversal_t(8, true), // end (or start) node 
                                    "path0#0#path0",                   // reference path 
                                    3,                                 // start offset
                                    6,                                 // end offset
                                    1,                                 // depth
                                    "2,3",                                // variant type (allele length counts)
                                    snarl_walks2,                       // walks through the snarl 
                                    partitions,                        // set of samples per walk 
                                    sequences                          // sequences per walk 
                                    );
        std::vector<Path_traversal_t> snarl_walks3;
        snarl_walks3.emplace_back();
        snarl_walks3.back().add_node_traversal_t(Node_traversal_t(6, false));
        stoat::snarl_info_t snarl3 (stoat::Node_traversal_t(5, false), // start (or end) node 
                                    stoat::Node_traversal_t(7, true), // end (or start) node 
                                    "path0#0#path0",                   // reference path 
                                    4,                                 // start offset
                                    5,                                 // end offset
                                    2,                                 // depth
                                    "1",                                // variant type (allele length counts)
                                    snarl_walks3,                       // walks through the snarl 
                                    partitions,                        // set of samples per walk 
                                    sequences                          // sequences per walk 
                                    );

        std::vector<Path_traversal_t> snarl_walks4;
        snarl_walks4.emplace_back();
        stoat::snarl_info_t snarl4 (stoat::Node_traversal_t(8, false), // start (or end) node 
                                    stoat::Node_traversal_t(10, true), // end (or start) node 
                                    "NA",                              // reference path 
                                    0,                                 // start offset
                                    0,                                 // end offset
                                    1,                                 // depth
                                    "",                                // variant type (allele length counts)
                                    snarl_walks4,                       // walks through the snarl 
                                    partitions,                        // set of samples per walk 
                                    sequences                          // sequences per walk 
                                    );
        snarl_collection.for_each_snarl([&](const snarl_info_t& snarl_info) {
            snarl_count++;

            if ((snarl_info.start_node == snarl1.start_node && snarl_info.end_node == snarl1.end_node) ||
                (snarl_info.start_node == snarl1.end_node && snarl_info.end_node == snarl1.start_node)) {
                REQUIRE(snarl_info.ref_path == snarl1.ref_path);
                REQUIRE(snarl_info.start_position == snarl1.start_position);
                REQUIRE(snarl_info.end_position == snarl1.end_position);
                REQUIRE(snarl_info.depth == snarl1.depth);
                REQUIRE(snarl_info.variant_type == snarl1.variant_type);
                REQUIRE(snarl_info.partitions == snarl1.partitions);
                REQUIRE(snarl_info.snarl_walks.size() == 1);
                REQUIRE(snarl_info.snarl_walks[0].get_paths().size() == 1);
                REQUIRE(snarl_info.snarl_walks[0].get_paths()[0] == snarl_walks1[0].get_paths()[0]);
                REQUIRE(snarl_info.sequences.size() == 1);
                REQUIRE(snarl_info.sequences[0] == "C");
            } else if ((snarl_info.start_node == snarl2.start_node && snarl_info.end_node == snarl2.end_node) ||
                (snarl_info.start_node == snarl2.end_node && snarl_info.end_node == snarl2.start_node)) {
                REQUIRE(snarl_info.ref_path == snarl2.ref_path);
                REQUIRE(snarl_info.start_position == snarl2.start_position);
                REQUIRE(snarl_info.end_position == snarl2.end_position);
                REQUIRE(snarl_info.depth == snarl2.depth);
                REQUIRE(snarl_info.variant_type == snarl2.variant_type);
                REQUIRE(snarl_info.partitions == snarl2.partitions);
                REQUIRE(snarl_info.snarl_walks.size() == 1);
                REQUIRE(snarl_info.snarl_walks[0].get_paths().size() == 3);
                REQUIRE(snarl_info.snarl_walks[0].get_paths()[0] == snarl_walks2[0].get_paths()[0]);
                REQUIRE(snarl_info.snarl_walks[0].get_paths()[1] == snarl_walks2[0].get_paths()[1]);
                REQUIRE(snarl_info.snarl_walks[0].get_paths()[2] == snarl_walks2[0].get_paths()[2]);
                REQUIRE(snarl_info.sequences.size() == 1);
                REQUIRE(snarl_info.sequences[0] == "TA");
            }  else if ((snarl_info.start_node == snarl3.start_node && snarl_info.end_node == snarl3.end_node) ||
                (snarl_info.start_node == snarl3.end_node && snarl_info.end_node == snarl3.start_node)) {
                REQUIRE(snarl_info.ref_path == snarl3.ref_path);
                REQUIRE(snarl_info.start_position == snarl3.start_position);
                REQUIRE(snarl_info.end_position == snarl3.end_position);
                REQUIRE(snarl_info.depth == snarl3.depth);
                REQUIRE(snarl_info.variant_type == snarl3.variant_type);
                REQUIRE(snarl_info.partitions == snarl3.partitions);
                REQUIRE(snarl_info.snarl_walks.size() == 1);
                REQUIRE(snarl_info.snarl_walks[0].get_paths().size() == 1);
                REQUIRE(snarl_info.snarl_walks[0].get_paths()[0] == snarl_walks3[0].get_paths()[0]);
                REQUIRE(snarl_info.sequences.size() == 1);
                REQUIRE(snarl_info.sequences[0] == "C");
            } else if ((snarl_info.start_node == snarl4.start_node && snarl_info.end_node == snarl4.end_node) ||
                (snarl_info.start_node == snarl4.end_node && snarl_info.end_node == snarl4.start_node)) {
                REQUIRE(snarl_info.ref_path == snarl4.ref_path);
                REQUIRE(snarl_info.start_position == snarl4.start_position);
                REQUIRE(snarl_info.end_position == snarl4.end_position);
                REQUIRE(snarl_info.depth == snarl4.depth);
                REQUIRE(snarl_info.variant_type == snarl4.variant_type);
                REQUIRE(snarl_info.partitions == snarl4.partitions);
                REQUIRE(snarl_info.snarl_walks.size() == 0);
                REQUIRE(snarl_info.sequences.size() == 0);
            } else {
                REQUIRE(false);
            }
        });

        REQUIRE(snarl_count == 4);
    }



}

}
