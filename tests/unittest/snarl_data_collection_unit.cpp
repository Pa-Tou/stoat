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



    SECTION("Make and fill in snarl collection with no partitions") {
        // Don't get the partitions or anything else
        TestSnarlDataCollection snarl_collection(1,10,1,10);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, false, [&](const handlegraph::net_handle_t& net) { 
            std::vector<std::set<sample_hap_t>> partitions;
            return partitions;
        }, true, "", false, cout);

        // Check that we got all snarls and that we got the correct snarls
        size_t snarl_count = 0;

        snarl_collection.for_each_snarl([&](const snarl_info_t& snarl_info) {
            snarl_count++;

            if ((snarl_info.start_node == Node_traversal_t(1, false) && snarl_info.end_node == stoat::Node_traversal_t(4, true)) ||
                (snarl_info.start_node == stoat::Node_traversal_t(4, true) && snarl_info.end_node == Node_traversal_t(1, false))) {
                // First snarl

                REQUIRE(snarl_info.ref_path == "path0#0#path0");
                REQUIRE(snarl_info.start_position == 1);
                REQUIRE(snarl_info.end_position == 2);
                REQUIRE(snarl_info.depth == 1);
                REQUIRE(snarl_info.snarl_walks.size() == 2);
                REQUIRE(snarl_info.snarl_walks[0].get_paths().size() == 3);
                REQUIRE(snarl_info.snarl_walks[1].get_paths().size() == 3);

                // The paths should be 1-2-4 and 1-3-4
                REQUIRE(snarl_info.snarl_walks[0].get_paths()[0] == Node_traversal_t(1, false));
                REQUIRE(snarl_info.snarl_walks[1].get_paths()[0] == Node_traversal_t(1, false));
                REQUIRE(((snarl_info.snarl_walks[0].get_paths()[1] == Node_traversal_t(2, false) && snarl_info.snarl_walks[1].get_paths()[1] == Node_traversal_t(3, false))
                     || (snarl_info.snarl_walks[0].get_paths()[1] == Node_traversal_t(3, false) && snarl_info.snarl_walks[1].get_paths()[1] == Node_traversal_t(2, false))));

                REQUIRE(snarl_info.snarl_walks[0].get_paths()[2] == Node_traversal_t(4, false));
                REQUIRE(snarl_info.snarl_walks[1].get_paths()[2] == Node_traversal_t(4, false));

                REQUIRE(snarl_info.variant_type == "1,1");

                // The sequences should both be CCA, oops
                REQUIRE(snarl_info.sequences.size() == 2);
                REQUIRE(snarl_info.sequences[0] == "CCA");
                REQUIRE(snarl_info.sequences[1] == "CCA");

            } else if ((snarl_info.start_node == stoat::Node_traversal_t(4, false) && snarl_info.end_node == stoat::Node_traversal_t(8, true)) ||
                (snarl_info.start_node == stoat::Node_traversal_t(8, true) && snarl_info.end_node == stoat::Node_traversal_t(4, false))) {
                // Snarl 4-8
                REQUIRE(snarl_info.ref_path == "path0#0#path0");
                REQUIRE(snarl_info.start_position == 3);
                REQUIRE(snarl_info.end_position == 6);
                REQUIRE(snarl_info.depth == 1);

                REQUIRE(snarl_info.snarl_walks.size() == 2);
                size_t insertion_index;
                size_t deletion_index;
                if (snarl_info.snarl_walks[0].get_paths().size() == 2) {
                    // First is deletion, second is insertion
                    REQUIRE(snarl_info.snarl_walks[1].get_paths().size() == 5);
                    deletion_index = 0;
                    insertion_index = 1;
                    REQUIRE(snarl_info.variant_type == "0,2/3");
                } else {
                    // First is insertion, second is deletion
                    REQUIRE(snarl_info.snarl_walks[0].get_paths().size() == 5);
                    REQUIRE(snarl_info.snarl_walks[1].get_paths().size() == 2);
                    deletion_index = 1;
                    insertion_index = 0;
                    REQUIRE(snarl_info.variant_type == "2/3,0");
                }

                // The paths should be 4-8 and 4-5-0-7-8
                REQUIRE(snarl_info.snarl_walks[deletion_index].get_paths()[0] == Node_traversal_t(4, false));
                REQUIRE(snarl_info.snarl_walks[deletion_index].get_paths()[1] == Node_traversal_t(8, false));

                REQUIRE(snarl_info.snarl_walks[insertion_index].get_paths()[0] == Node_traversal_t(4, false));
                REQUIRE(snarl_info.snarl_walks[insertion_index].get_paths()[1] == Node_traversal_t(5, false));
                REQUIRE(snarl_info.snarl_walks[insertion_index].get_paths()[2] == Node_traversal_t(0, false));
                REQUIRE(snarl_info.snarl_walks[insertion_index].get_paths()[3] == Node_traversal_t(7, false));
                REQUIRE(snarl_info.snarl_walks[insertion_index].get_paths()[4] == Node_traversal_t(8, false));


                // The sequences 
                REQUIRE(snarl_info.sequences.size() == 2);
                REQUIRE(snarl_info.sequences[insertion_index] == "ATAC");
                REQUIRE(snarl_info.sequences[deletion_index] == "AC");

            }  else if ((snarl_info.start_node == stoat::Node_traversal_t(5, false) && snarl_info.end_node == stoat::Node_traversal_t(7, true)) ||
                (snarl_info.start_node == stoat::Node_traversal_t(7, true) && snarl_info.end_node == stoat::Node_traversal_t(5, false))) {
                REQUIRE(snarl_info.ref_path == "path0#0#path0");
                REQUIRE(snarl_info.start_position == 4);
                REQUIRE(snarl_info.end_position == 5);
                REQUIRE(snarl_info.depth == 2);

                REQUIRE(snarl_info.snarl_walks.size() == 2);

                size_t insertion_index;
                size_t deletion_index;
                if (snarl_info.snarl_walks[0].get_paths().size() == 2) {
                    // First is deletion, second is insertion
                    REQUIRE(snarl_info.snarl_walks[1].get_paths().size() == 3);
                    deletion_index = 0;
                    insertion_index = 1;
                    REQUIRE(snarl_info.variant_type == "0,1");
                } else {
                    // First is insertion, second is deletion
                    REQUIRE(snarl_info.snarl_walks[0].get_paths().size() == 3);
                    REQUIRE(snarl_info.snarl_walks[1].get_paths().size() == 2);
                    deletion_index = 1;
                    insertion_index = 0;
                    REQUIRE(snarl_info.variant_type == "1,0");
                }

                // The paths should be 5-7 and 5-6-7
                REQUIRE(snarl_info.snarl_walks[deletion_index].get_paths()[0] == Node_traversal_t(5, false));
                REQUIRE(snarl_info.snarl_walks[deletion_index].get_paths()[1] == Node_traversal_t(7, false));

                REQUIRE(snarl_info.snarl_walks[insertion_index].get_paths()[0] == Node_traversal_t(5, false));
                REQUIRE(snarl_info.snarl_walks[insertion_index].get_paths()[1] == Node_traversal_t(6, false));
                REQUIRE(snarl_info.snarl_walks[insertion_index].get_paths()[2] == Node_traversal_t(7, false));

                // The sequences 
                REQUIRE(snarl_info.sequences.size() == 2);
                REQUIRE(snarl_info.sequences[insertion_index] == "TCA");
                REQUIRE(snarl_info.sequences[deletion_index] == "TA");

            } else if ((snarl_info.start_node == stoat::Node_traversal_t(8, false) && snarl_info.end_node == stoat::Node_traversal_t(10, true)) ||
                (snarl_info.start_node == stoat::Node_traversal_t(10, true) && snarl_info.end_node == stoat::Node_traversal_t(8, false))) {
                REQUIRE(snarl_info.ref_path == "NA");
                REQUIRE(snarl_info.start_position == 0);
                REQUIRE(snarl_info.end_position == 0);
                REQUIRE(snarl_info.depth == 1);
                REQUIRE(snarl_info.snarl_walks.size() == 2);

                size_t insertion_index;
                size_t deletion_index;
                if (snarl_info.snarl_walks[0].get_paths().size() == 2) {
                    // First is deletion, second is insertion
                    REQUIRE(snarl_info.snarl_walks[1].get_paths().size() == 3);
                    deletion_index = 0;
                    insertion_index = 1;
                    REQUIRE(snarl_info.variant_type == "0,1");
                } else {
                    // First is insertion, second is deletion
                    REQUIRE(snarl_info.snarl_walks[0].get_paths().size() == 3);
                    REQUIRE(snarl_info.snarl_walks[1].get_paths().size() == 2);
                    deletion_index = 1;
                    insertion_index = 0;
                    REQUIRE(snarl_info.variant_type == "1,0");
                }

                // The paths should be 8-10 and 8-9-10
                REQUIRE(snarl_info.snarl_walks[deletion_index].get_paths()[0] == Node_traversal_t(8, false));
                REQUIRE(snarl_info.snarl_walks[deletion_index].get_paths()[1] == Node_traversal_t(10, false));

                REQUIRE(snarl_info.snarl_walks[insertion_index].get_paths()[0] == Node_traversal_t(8, false));
                REQUIRE(snarl_info.snarl_walks[insertion_index].get_paths()[1] == Node_traversal_t(9, false));
                REQUIRE(snarl_info.snarl_walks[insertion_index].get_paths()[2] == Node_traversal_t(10, false));

                // The sequences 
                REQUIRE(snarl_info.sequences.size() == 2);
                REQUIRE(snarl_info.sequences[insertion_index] == "CAA");
                REQUIRE(snarl_info.sequences[deletion_index] == "CA");

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

        std::string reference_path = "path0#0#path0";

        // Get just the walk for the reference since this only tests the partitions containing only the reference
        std::vector<std::set<sample_hap_t>> partitions;
        partitions.emplace_back();
        partitions.back().emplace(ref_sample);
        std::vector<std::string> sequences;



        snarl_collection.for_each_snarl([&](const snarl_info_t& snarl_info) {
            snarl_count++;

            if ((snarl_info.start_node == Node_traversal_t(1, false) && snarl_info.end_node == stoat::Node_traversal_t(4, true)) ||
                (snarl_info.start_node == stoat::Node_traversal_t(4, true) && snarl_info.end_node == Node_traversal_t(1, false))) {
                // First snarl

                REQUIRE(snarl_info.ref_path == "path0#0#path0");
                REQUIRE(snarl_info.start_position == 1);
                REQUIRE(snarl_info.end_position == 2);
                REQUIRE(snarl_info.depth == 1);
                REQUIRE(snarl_info.snarl_walks.size() == 1);
                REQUIRE(snarl_info.snarl_walks[0].get_paths().size() == 3);

                // The paths should be 1-2-4 and 1-3-4
                REQUIRE(snarl_info.snarl_walks[0].get_paths()[0] == Node_traversal_t(1, false));
                REQUIRE(snarl_info.snarl_walks[0].get_paths()[1] == Node_traversal_t(2, false));
                REQUIRE(snarl_info.snarl_walks[0].get_paths()[2] == Node_traversal_t(4, false));

                REQUIRE(snarl_info.variant_type == "1");

                // The sequence should be CCA
                REQUIRE(snarl_info.sequences.size() == 1);
                REQUIRE(snarl_info.sequences[0] == "CCA");

                REQUIRE(snarl_info.partitions == partitions);
            } else if ((snarl_info.start_node == stoat::Node_traversal_t(4, false) && snarl_info.end_node == stoat::Node_traversal_t(8, true)) ||
                (snarl_info.start_node == stoat::Node_traversal_t(8, true) && snarl_info.end_node == stoat::Node_traversal_t(4, false))) {
                // Snarl 4-8
                REQUIRE(snarl_info.ref_path == "path0#0#path0");
                REQUIRE(snarl_info.start_position == 3);
                REQUIRE(snarl_info.end_position == 6);
                REQUIRE(snarl_info.depth == 1);

                REQUIRE(snarl_info.snarl_walks.size() == 1);
                REQUIRE(snarl_info.snarl_walks[0].get_paths().size() == 5);
                REQUIRE(snarl_info.variant_type == "2/3");

                // The paths should be and 4-5-0-7-8

                REQUIRE(snarl_info.snarl_walks[0].get_paths()[0] == Node_traversal_t(4, false));
                REQUIRE(snarl_info.snarl_walks[0].get_paths()[1] == Node_traversal_t(5, false));
                REQUIRE(snarl_info.snarl_walks[0].get_paths()[2] == Node_traversal_t(0, false));
                REQUIRE(snarl_info.snarl_walks[0].get_paths()[3] == Node_traversal_t(7, false));
                REQUIRE(snarl_info.snarl_walks[0].get_paths()[4] == Node_traversal_t(8, false));


                // The sequences 
                REQUIRE(snarl_info.sequences.size() == 1);
                REQUIRE(snarl_info.sequences[0] == "ATAC");

                REQUIRE(snarl_info.partitions == partitions);

            }  else if ((snarl_info.start_node == stoat::Node_traversal_t(5, false) && snarl_info.end_node == stoat::Node_traversal_t(7, true)) ||
                (snarl_info.start_node == stoat::Node_traversal_t(7, true) && snarl_info.end_node == stoat::Node_traversal_t(5, false))) {
                REQUIRE(snarl_info.ref_path == "path0#0#path0");
                REQUIRE(snarl_info.start_position == 4);
                REQUIRE(snarl_info.end_position == 5);
                REQUIRE(snarl_info.depth == 2);

                REQUIRE(snarl_info.snarl_walks.size() == 1);

                REQUIRE(snarl_info.snarl_walks[0].get_paths().size() == 3);
                REQUIRE(snarl_info.variant_type == "1");

                // The paths should be 5-6-7

                REQUIRE(snarl_info.snarl_walks[0].get_paths()[0] == Node_traversal_t(5, false));
                REQUIRE(snarl_info.snarl_walks[0].get_paths()[1] == Node_traversal_t(6, false));
                REQUIRE(snarl_info.snarl_walks[0].get_paths()[2] == Node_traversal_t(7, false));

                // The sequences 
                REQUIRE(snarl_info.sequences.size() == 1);
                REQUIRE(snarl_info.sequences[0] == "TCA");

                REQUIRE(snarl_info.partitions == partitions);

            } else if ((snarl_info.start_node == stoat::Node_traversal_t(8, false) && snarl_info.end_node == stoat::Node_traversal_t(10, true)) ||
                (snarl_info.start_node == stoat::Node_traversal_t(10, true) && snarl_info.end_node == stoat::Node_traversal_t(8, false))) {
                REQUIRE(snarl_info.ref_path == "NA");
                REQUIRE(snarl_info.start_position == 0);
                REQUIRE(snarl_info.end_position == 0);
                REQUIRE(snarl_info.depth == 1);
                REQUIRE(snarl_info.snarl_walks.size() == 0);

                // The sequences 
                REQUIRE(snarl_info.sequences.size() == 0);


            } else {
                REQUIRE(false);
            }
        });

        REQUIRE(snarl_count == 4);
    }



}

}
