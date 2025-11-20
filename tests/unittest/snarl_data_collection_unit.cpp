#include <catch.hpp>
#include <bdsg/hash_graph.hpp>
#include <bdsg/overlays/overlay_helper.hpp>
#include "../../src/snarl_data_collection.hpp"
#include "../../src/log.hpp"

using namespace stoat;


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
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, 
            true, // partitions before walks but it doesn't matter here
            false, // don't get walks 
            [&](const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                                                         const net_handle_t& snarl, const snarl_info_t& snarl_data,
                                                         std::vector<std::vector<handlegraph::net_handle_t>>& walks) { 
                return;
            }, false, // don't get partitions 
            [&](const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                                                         const net_handle_t& snarl, const snarl_info_t& snarl_data,
                                                         std::vector<std::set<sample_hap_t>>& partitions) { 
                return;
            },
            false, // don't get sequences
            "", false);

    }
    SECTION("Serialize partitioner") {
        // There isn't much to do with one node so just make sure we can run the constructor without crashing

        TestSnarlDataCollection snarl_collection(1,1,1,1);

        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, 
            true, // partitions before walks but it doesn't matter here
            false, // don't get walks 
            [&](const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                                                         const net_handle_t& snarl, const snarl_info_t& snarl_data,
                                                         std::vector<std::vector<handlegraph::net_handle_t>>& walks) { 
                return;
            }, false, // don't get partitions 
            [&](const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                                                         const net_handle_t& snarl, const snarl_info_t& snarl_data,
                                                         std::vector<std::set<sample_hap_t>>& partitions) { 
                return;
            },
            false, // don't get sequences
            "", false);

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

    // If the walks were found without partitions then there is a walk for the last snarl
    auto check_collection = [&] (const TestSnarlDataCollection& snarl_collection, bool check_walks, bool check_partitions, bool check_sequences, bool get_all_walks) {

        // Get the partitions again so we can check them. the order might be different though

        // Make partitions for each snarl
        std::vector<stoat::sample_hap_t> sample_haps;
        for (const auto& path : paths) {
            sample_haps.emplace_back(stoat::get_sample_and_haplotype(*path_graph, path));
        }

        std::vector<std::vector<std::set<sample_hap_t>>> partitions_by_snarl;
        // snarl 1-4
        //[0,1], [2,3]
        partitions_by_snarl.emplace_back();
        partitions_by_snarl.back().emplace_back();
        partitions_by_snarl.back().back().emplace(sample_haps[0]);
        partitions_by_snarl.back().back().emplace(sample_haps[1]);
        partitions_by_snarl.back().emplace_back();
        partitions_by_snarl.back().back().emplace(sample_haps[2]);
        partitions_by_snarl.back().back().emplace(sample_haps[3]);

        // snarl 4-8
        //[0,1,3], [2]
        partitions_by_snarl.emplace_back();
        partitions_by_snarl.back().emplace_back();
        partitions_by_snarl.back().back().emplace(sample_haps[0]);
        partitions_by_snarl.back().back().emplace(sample_haps[1]);
        partitions_by_snarl.back().back().emplace(sample_haps[3]);
        partitions_by_snarl.back().emplace_back();
        partitions_by_snarl.back().back().emplace(sample_haps[2]);

        // snarl 5-7
        //[0], [1,3]
        partitions_by_snarl.emplace_back();
        partitions_by_snarl.back().emplace_back();
        partitions_by_snarl.back().back().emplace(sample_haps[0]);
        partitions_by_snarl.back().emplace_back();
        partitions_by_snarl.back().back().emplace(sample_haps[1]);
        partitions_by_snarl.back().back().emplace(sample_haps[3]);

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
                if (check_walks) {
                    REQUIRE(snarl_info.snarl_walks.size() == 2);
                    REQUIRE(snarl_info.snarl_walks[0].get_paths().size() == 3);
                    REQUIRE(snarl_info.snarl_walks[1].get_paths().size() == 3);
                }

                // The paths should be 1-2-4 and 1-3-4
                if (check_partitions) {
                    REQUIRE(snarl_info.partitions.size() == 2);
                    REQUIRE(snarl_info.partitions[0].size() == 2);
                    REQUIRE(snarl_info.partitions[1].size() == 2);
                    REQUIRE(((snarl_info.partitions[0] == partitions_by_snarl[0][0] && snarl_info.partitions[1] == partitions_by_snarl[0][1])
                         || (snarl_info.partitions[1] == partitions_by_snarl[0][0] && snarl_info.partitions[0] == partitions_by_snarl[0][1])));
                }
                if (check_walks) {
                    REQUIRE(snarl_info.snarl_walks[0].get_paths()[0] == Node_traversal_t(1, false));
                    REQUIRE(snarl_info.snarl_walks[1].get_paths()[0] == Node_traversal_t(1, false));

                    // The second node is either 2 or 3
                    REQUIRE((snarl_info.snarl_walks[0].get_paths()[1] == Node_traversal_t(2, false) ||
                             snarl_info.snarl_walks[0].get_paths()[1] == Node_traversal_t(3, false)));

                    REQUIRE((snarl_info.snarl_walks[1].get_paths()[1] == Node_traversal_t(2, false) ||
                             snarl_info.snarl_walks[1].get_paths()[1] == Node_traversal_t(3, false)));

                    REQUIRE(!(snarl_info.snarl_walks[0].get_paths()[1] == snarl_info.snarl_walks[1].get_paths()[1]));

                    REQUIRE(snarl_info.snarl_walks[0].get_paths()[2] == Node_traversal_t(4, false));
                    REQUIRE(snarl_info.snarl_walks[1].get_paths()[2] == Node_traversal_t(4, false));
                    REQUIRE(snarl_info.variant_type == "1,1");
                }

                // Check that the walks and partitions match
                if (check_walks && check_partitions) {

                    if (snarl_info.snarl_walks[0].get_paths()[1] == Node_traversal_t(2, false)) {
                        REQUIRE(snarl_info.snarl_walks[1].get_paths()[1] == Node_traversal_t(3, false));
                        if (check_partitions) {
                            REQUIRE(snarl_info.partitions[0] == partitions_by_snarl[0][0]);
                            REQUIRE(snarl_info.partitions[1] == partitions_by_snarl[0][1]);
                        }
                    } else {
                        REQUIRE(snarl_info.snarl_walks[0].get_paths()[1] == Node_traversal_t(3, false));
                        REQUIRE(snarl_info.snarl_walks[1].get_paths()[1] == Node_traversal_t(2, false));
                        if (check_partitions) {
                            REQUIRE(snarl_info.partitions[0] == partitions_by_snarl[0][1]);
                            REQUIRE(snarl_info.partitions[1] == partitions_by_snarl[0][0]);
                        }
                    }
                }


                if (check_sequences) {
                    // The sequences should both be CCA, oops
                    REQUIRE(snarl_info.sequences.size() == 2);
                    REQUIRE(snarl_info.sequences[0] == "CCA");
                    REQUIRE(snarl_info.sequences[1] == "CCA");
                }

            } else if ((snarl_info.start_node == stoat::Node_traversal_t(4, false) && snarl_info.end_node == stoat::Node_traversal_t(8, true)) ||
                (snarl_info.start_node == stoat::Node_traversal_t(8, true) && snarl_info.end_node == stoat::Node_traversal_t(4, false))) {
                // Snarl 4-8
                REQUIRE(snarl_info.ref_path == "path0#0#path0");
                REQUIRE(snarl_info.start_position == 3);
                REQUIRE(snarl_info.end_position == 6);
                REQUIRE(snarl_info.depth == 1);

                size_t insertion_index;
                size_t deletion_index;
                if (check_partitions) {
                    REQUIRE(snarl_info.partitions.size() == 2);
                }
                if (check_walks) {
                    REQUIRE(snarl_info.snarl_walks.size() == 2);
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

                    if (check_sequences) {
                        // The sequences 
                        REQUIRE(snarl_info.sequences.size() == 2);
                        REQUIRE(snarl_info.sequences[insertion_index] == "ATAC");
                        REQUIRE(snarl_info.sequences[deletion_index] == "AC");
                    }
                    if (check_partitions) {
                        REQUIRE(snarl_info.partitions.size() == 2);
                        REQUIRE(snarl_info.partitions[insertion_index] == partitions_by_snarl[1][0]);
                        REQUIRE(snarl_info.partitions[deletion_index] == partitions_by_snarl[1][1]); 
                    }
                } else if (check_partitions) {
                    // If we want to check the partitions but don't have the walks
                    REQUIRE(snarl_info.partitions.size() == 2);
                    REQUIRE(((snarl_info.partitions[0] == partitions_by_snarl[1][0] && snarl_info.partitions[1] == partitions_by_snarl[1][1])
                         || (snarl_info.partitions[1] == partitions_by_snarl[1][0] && snarl_info.partitions[0] == partitions_by_snarl[1][1]))); 
                }


            }  else if ((snarl_info.start_node == stoat::Node_traversal_t(5, false) && snarl_info.end_node == stoat::Node_traversal_t(7, true)) ||
                (snarl_info.start_node == stoat::Node_traversal_t(7, true) && snarl_info.end_node == stoat::Node_traversal_t(5, false))) {
                REQUIRE(snarl_info.ref_path == "path0#0#path0");
                REQUIRE(snarl_info.start_position == 4);
                REQUIRE(snarl_info.end_position == 5);
                REQUIRE(snarl_info.depth == 2);


                if (check_walks) {
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
                    if (check_sequences) {
                        // The sequences 
                        REQUIRE(snarl_info.sequences.size() == 2);
                        REQUIRE(snarl_info.sequences[insertion_index] == "TCA");
                        REQUIRE(snarl_info.sequences[deletion_index] == "TA");
                    }
                    if (check_partitions) {
                        REQUIRE(snarl_info.partitions.size() == 2);
                        REQUIRE(snarl_info.partitions[insertion_index] == partitions_by_snarl[2][0]);
                        REQUIRE(snarl_info.partitions[deletion_index] == partitions_by_snarl[2][1]); 
                    }
                } else if (check_partitions) {
                    // If we want to check the partitions but don't have the walks
                    REQUIRE(snarl_info.partitions.size() == 2);
                    REQUIRE(((snarl_info.partitions[0] == partitions_by_snarl[2][0] && snarl_info.partitions[1] == partitions_by_snarl[2][1])
                         || (snarl_info.partitions[1] == partitions_by_snarl[2][0] && snarl_info.partitions[0] == partitions_by_snarl[2][1]))); 
                }



            } else if ((snarl_info.start_node == stoat::Node_traversal_t(8, false) && snarl_info.end_node == stoat::Node_traversal_t(10, true)) ||
                (snarl_info.start_node == stoat::Node_traversal_t(10, true) && snarl_info.end_node == stoat::Node_traversal_t(8, false))) {
                REQUIRE(snarl_info.ref_path == "NA");
                REQUIRE(snarl_info.start_position == 0);
                REQUIRE(snarl_info.end_position == 0);
                REQUIRE(snarl_info.depth == 1);
                if (get_all_walks && check_walks) {
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

                    if (check_sequences) {
                        // The sequences 
                        REQUIRE(snarl_info.sequences.size() == 2);
                        REQUIRE(snarl_info.sequences[insertion_index] == "CAA");
                        REQUIRE(snarl_info.sequences[deletion_index] == "CA");
                    }

                    if (check_partitions) {
                        REQUIRE(snarl_info.partitions.size() == 0);
                    }
                }

            } else {
                REQUIRE(false);
            }
        });

        REQUIRE(snarl_count == 4);
    };


    SECTION("Make and fill in snarl collection with no walks, partitions, or sequences") {

        TestSnarlDataCollection snarl_collection(1,10,1,10);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, 
            true, // partitions before walks but it doesn't matter here
            false, // don't get walks 
            [&](const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                                                         const net_handle_t& snarl, const snarl_info_t& snarl_data,
                                                         std::vector<std::vector<handlegraph::net_handle_t>>& walks) { 
                return;
            }, false, // don't get partitions 
            [&](const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                                                         const net_handle_t& snarl, const snarl_info_t& snarl_data,
                                                         std::vector<std::set<sample_hap_t>>& partitions) { 
                return;
            },
            false, // don't get sequences
            "", false);


        check_collection(snarl_collection, false, false, false, true);
    }
    SECTION("Make and fill in snarl collection with walks and sequences but no partitions") {

        TestSnarlDataCollection snarl_collection(1,10,1,10);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, 
            true, // partitions before walks but it doesn't matter here
            true, // get walks 
            SnarlDataCollection::get_all_walks_through_snarl, 
            false, // don't get partitions 
            [&](const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                                                         const net_handle_t& snarl, const snarl_info_t& snarl_data,
                                                         std::vector<std::set<sample_hap_t>>& partitions) { 
                return;
            },
            true, // get sequences
            "", false);


        check_collection(snarl_collection, true, false, true, true);
    }

    SECTION("Make and fill in snarl collection with walks, sequences, and no partitions, then fill in partitions later") {

        // Don't get the partitions or anything else
        TestSnarlDataCollection snarl_collection(1,10,1,10);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, 
            false, //walks before partitions but it doesn't matter here
            true, // get walks 
            SnarlDataCollection::get_all_walks_through_snarl, 
            false, // don't get partitions 
            [&](const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                                                         const net_handle_t& snarl, const snarl_info_t& snarl_data,
                                                         std::vector<std::set<sample_hap_t>>& partitions) { 
                return;
            },
            true, // get sequences
            "", false);


        check_collection(snarl_collection, true, false, true, true);



        // Make partitions for each snarl
        std::vector<stoat::sample_hap_t> sample_haps;
        for (const auto& path : paths) {
            sample_haps.emplace_back(stoat::get_sample_and_haplotype(*path_graph, path));
        }


        std::unordered_map<stoat::sample_hap_t, size_t> sample_haplotype_to_index;
        snarl_collection.add_snarl_partitions(sample_haplotype_to_index, [&](const snarl_info_t& snarl_info) {

            //Add partitions to follow the walks
            std::vector<std::set<sample_hap_t>> partitions;
            if (snarl_info.start_node == Node_traversal_t(1, false) || snarl_info.start_node == stoat::Node_traversal_t(4, true)) {
                // snarl 1-4
                //[0,1], [2,3]
                if (snarl_info.snarl_walks[0].get_paths()[1] == Node_traversal_t(2, false)) {
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[0]);
                    partitions.back().emplace(sample_haps[1]);
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[2]);
                    partitions.back().emplace(sample_haps[3]);
                } else {
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[2]);
                    partitions.back().emplace(sample_haps[3]);
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[0]);
                    partitions.back().emplace(sample_haps[1]);
                }
            } else if (snarl_info.start_node == stoat::Node_traversal_t(4, false) || snarl_info.start_node == stoat::Node_traversal_t(8, true)) {
                // snarl 4-8
                //[0,1,3], [2]
                if (snarl_info.snarl_walks[0].get_paths().size() == 5) {
                    //insertion first
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[0]);
                    partitions.back().emplace(sample_haps[1]);
                    partitions.back().emplace(sample_haps[3]);
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[2]);
                } else {
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[2]);
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[0]);
                    partitions.back().emplace(sample_haps[1]);
                    partitions.back().emplace(sample_haps[3]);
                }
            } else if (snarl_info.start_node == stoat::Node_traversal_t(5, false) || snarl_info.start_node == stoat::Node_traversal_t(7, true)) {

                // snarl 5-7
                //[0], [1,3]
                if (snarl_info.snarl_walks[0].get_paths().size() == 3) {
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[0]);
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[1]);
                    partitions.back().emplace(sample_haps[3]);
                } else {
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[1]);
                    partitions.back().emplace(sample_haps[3]);
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[0]);
                }
            } else {
                // For the last one there are no partitions because there are no paths
            }
            return partitions;
        });

        check_collection(snarl_collection, true, true, true, true);

    }
    SECTION("Make and fill in snarl collection with no partitions, then fill in partitions later by chromosome") {

        // Don't get the partitions or anything else
        TestSnarlDataCollection snarl_collection(1,10,1,10);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, 
            false, // walks before partitions but it doesn't matter here
            true, // get walks 
            SnarlDataCollection::get_all_walks_through_snarl,
            false, // don't get partitions 
            [&](const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                                                         const net_handle_t& snarl, const snarl_info_t& snarl_data,
                                                         std::vector<std::set<sample_hap_t>>& partitions) { 
                return;
            },
            true, // get sequences
            "", false);

        // Make partitions for each snarl
        std::vector<stoat::sample_hap_t> sample_haps;
        for (const auto& path : paths) {
            sample_haps.emplace_back(stoat::get_sample_and_haplotype(*path_graph, path));
        }


        // First fill it in with a non-existent path, which should do nothing
        std::unordered_map<stoat::sample_hap_t, size_t> sample_haplotype_to_index;
        snarl_collection.add_snarl_partitions(sample_haplotype_to_index, [&](const snarl_info_t& snarl_info) {

            //Add partitions to follow the walks
            std::vector<std::set<sample_hap_t>> partitions;
            if (snarl_info.start_node == Node_traversal_t(1, false) || snarl_info.start_node == stoat::Node_traversal_t(4, true)) {
                // snarl 1-4
                //[0,1], [2,3]
                if (snarl_info.snarl_walks[0].get_paths()[1] == Node_traversal_t(2, false)) {
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[0]);
                    partitions.back().emplace(sample_haps[1]);
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[2]);
                    partitions.back().emplace(sample_haps[3]);
                } else {
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[2]);
                    partitions.back().emplace(sample_haps[3]);
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[0]);
                    partitions.back().emplace(sample_haps[1]);
                }
            } else if (snarl_info.start_node == stoat::Node_traversal_t(4, false) || snarl_info.start_node == stoat::Node_traversal_t(8, true)) {
                // snarl 4-8
                //[0,1,3], [2]
                if (snarl_info.snarl_walks[0].get_paths().size() == 5) {
                    //insertion first
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[0]);
                    partitions.back().emplace(sample_haps[1]);
                    partitions.back().emplace(sample_haps[3]);
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[2]);
                } else {
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[2]);
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[0]);
                    partitions.back().emplace(sample_haps[1]);
                    partitions.back().emplace(sample_haps[3]);
                }
            } else if (snarl_info.start_node == stoat::Node_traversal_t(5, false) || snarl_info.start_node == stoat::Node_traversal_t(7, true)) {

                // snarl 5-7
                //[0], [1,3]
                if (snarl_info.snarl_walks[0].get_paths().size() == 3) {
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[0]);
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[1]);
                    partitions.back().emplace(sample_haps[3]);
                } else {
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[1]);
                    partitions.back().emplace(sample_haps[3]);
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[0]);
                }
            } else {
                // For the last one there are no partitions because there are no paths
            }
            return partitions;
        }, "empty_path");

        REQUIRE(sample_haplotype_to_index.size() == 0);

        snarl_collection.for_each_snarl([&](const snarl_info_t& snarl_info) {
            // Make sure the partitions haven't been filled in
            REQUIRE(snarl_info.partitions.size() == 0);
        });

        // Now fill it in with the paths
        snarl_collection.add_snarl_partitions(sample_haplotype_to_index, [&](const snarl_info_t& snarl_info) {

            //Add partitions to follow the walks
            std::vector<std::set<sample_hap_t>> partitions;
            if (snarl_info.start_node == Node_traversal_t(1, false) || snarl_info.start_node == stoat::Node_traversal_t(4, true)) {
                // snarl 1-4
                //[0,1], [2,3]
                if (snarl_info.snarl_walks[0].get_paths()[1] == Node_traversal_t(2, false)) {
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[0]);
                    partitions.back().emplace(sample_haps[1]);
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[2]);
                    partitions.back().emplace(sample_haps[3]);
                } else {
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[2]);
                    partitions.back().emplace(sample_haps[3]);
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[0]);
                    partitions.back().emplace(sample_haps[1]);
                }
            } else if (snarl_info.start_node == stoat::Node_traversal_t(4, false) || snarl_info.start_node == stoat::Node_traversal_t(8, true)) {
                // snarl 4-8
                //[0,1,3], [2]
                if (snarl_info.snarl_walks[0].get_paths().size() == 5) {
                    //insertion first
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[0]);
                    partitions.back().emplace(sample_haps[1]);
                    partitions.back().emplace(sample_haps[3]);
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[2]);
                } else {
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[2]);
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[0]);
                    partitions.back().emplace(sample_haps[1]);
                    partitions.back().emplace(sample_haps[3]);
                }
            } else if (snarl_info.start_node == stoat::Node_traversal_t(5, false) || snarl_info.start_node == stoat::Node_traversal_t(7, true)) {

                // snarl 5-7
                //[0], [1,3]
                if (snarl_info.snarl_walks[0].get_paths().size() == 3) {
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[0]);
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[1]);
                    partitions.back().emplace(sample_haps[3]);
                } else {
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[1]);
                    partitions.back().emplace(sample_haps[3]);
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[0]);
                }
            } else {
                // For the last one there are no partitions because there are no paths
            }
            return partitions;
        }, "path0#0#path0");

        check_collection(snarl_collection, true, true, true, true);

    }
    SECTION("Make and fill in snarl collection with walks, partitions, and sequences") {

        // Make partitions for each snarl
        std::vector<stoat::sample_hap_t> sample_haps;
        for (const auto& path : paths) {
            sample_haps.emplace_back(stoat::get_sample_and_haplotype(*path_graph, path));
        }


        TestSnarlDataCollection snarl_collection(1,10,1,10);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, 
            true, // partitions before walks but it doesn't matter here
            true, // get walks 
            SnarlDataCollection::get_walks_from_partitions,
            true, // get partitions 
            [&](const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                                                         const net_handle_t& snarl, const snarl_info_t& snarl_data,
                                                         std::vector<std::set<sample_hap_t>>& partitions) { 
                if (snarl == distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(2)))) {
                    // snarl 1-4
                    //[0,1], [2,3]
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[0]);
                    partitions.back().emplace(sample_haps[1]);
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[2]);
                    partitions.back().emplace(sample_haps[3]);
                } else if (snarl == distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(5)))) {
                    // snarl 4-8
                    //[0,1,3], [2]
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[0]);
                    partitions.back().emplace(sample_haps[1]);
                    partitions.back().emplace(sample_haps[3]);
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[2]);
                } else if (snarl == distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(6)))) {
                    // snarl 5-7
                    //[0], [1,3]
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[0]);
                    partitions.emplace_back();
                    partitions.back().emplace(sample_haps[1]);
                    partitions.back().emplace(sample_haps[3]);
                } else {
                    // For the last one there are no partitions because there are no paths
                }
                return;
            },
            true, // get sequences
            "", false);

        check_collection(snarl_collection, true, true, true, false);

        SECTION("Serialize it") {

            std::string test_file = "./test_snarls.txt";
            ofstream outstream;
            outstream.open(test_file);
            snarl_collection.write_snarl_data_collection(outstream);
            outstream.close();

            TestSnarlDataCollection loaded_snarl_collection(1,10,1,10);
            ifstream instream;
            instream.open(test_file);
            loaded_snarl_collection.load_snarl_data_collection(instream);
            instream.close();

            check_collection(loaded_snarl_collection, true, true, true, false);

            std::string rm_cmd = "rm " + test_file;

            int rm = system(rm_cmd.c_str()); 
        }
    }

}
TEST_CASE( "Snarl collection multiple connected components",
          "[snarl_collection]" ) {

    /*
                       5
                     /   \
            1       4 ----6    8
          /   \   /         \ / \
        0       3  ----------7---9
          \   /
            2
    this graph is duplicated 6x

   */

    bdsg::HashGraph graph;

    std::vector<std::string> sequences = { "C", "C", "C", "A", "T", "C", "A", "C", "A", "A", 
                                           "C", "C", "C", "A", "T", "C", "A", "C", "A", "A",  
                                           "C", "C", "C", "A", "T", "C", "A", "C", "A", "A",  
                                           "C", "C", "C", "A", "T", "C", "A", "C", "A", "A",  
                                           "C", "C", "C", "A", "T", "C", "A", "C", "A", "A", 
                                           "C", "C", "C", "A", "T", "C", "A", "C", "A", "A"};

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

    graph.create_edge(nodes[9], nodes[10]);

    graph.create_edge(nodes[10], nodes[11]);
    graph.create_edge(nodes[10], nodes[12]);
    graph.create_edge(nodes[11], nodes[13]);
    graph.create_edge(nodes[12], nodes[13]);
    graph.create_edge(nodes[13], nodes[14]);
    graph.create_edge(nodes[13], nodes[17]);
    graph.create_edge(nodes[14], nodes[15]);
    graph.create_edge(nodes[14], nodes[16]);
    graph.create_edge(nodes[15], nodes[16]);
    graph.create_edge(nodes[16], nodes[17]);
    graph.create_edge(nodes[17], nodes[18]);
    graph.create_edge(nodes[17], nodes[19]);
    graph.create_edge(nodes[18], nodes[19]);


    graph.create_edge(nodes[20], nodes[21]);
    graph.create_edge(nodes[20], nodes[22]);
    graph.create_edge(nodes[21], nodes[23]);
    graph.create_edge(nodes[22], nodes[23]);
    graph.create_edge(nodes[23], nodes[24]);
    graph.create_edge(nodes[23], nodes[27]);
    graph.create_edge(nodes[24], nodes[25]);
    graph.create_edge(nodes[24], nodes[26]);
    graph.create_edge(nodes[25], nodes[26]);
    graph.create_edge(nodes[26], nodes[27]);
    graph.create_edge(nodes[27], nodes[28]);
    graph.create_edge(nodes[27], nodes[29]);
    graph.create_edge(nodes[28], nodes[29]);

    graph.create_edge(nodes[29], nodes[30]);

    graph.create_edge(nodes[30], nodes[31]);
    graph.create_edge(nodes[30], nodes[32]);
    graph.create_edge(nodes[31], nodes[33]);
    graph.create_edge(nodes[32], nodes[33]);
    graph.create_edge(nodes[33], nodes[34]);
    graph.create_edge(nodes[33], nodes[37]);
    graph.create_edge(nodes[34], nodes[35]);
    graph.create_edge(nodes[34], nodes[36]);
    graph.create_edge(nodes[35], nodes[36]);
    graph.create_edge(nodes[36], nodes[37]);
    graph.create_edge(nodes[37], nodes[38]);
    graph.create_edge(nodes[37], nodes[39]);
    graph.create_edge(nodes[38], nodes[39]);

    graph.create_edge(nodes[39], nodes[40]);
    
    graph.create_edge(nodes[40], nodes[41]);
    graph.create_edge(nodes[40], nodes[42]);
    graph.create_edge(nodes[41], nodes[43]);
    graph.create_edge(nodes[42], nodes[43]);
    graph.create_edge(nodes[43], nodes[44]);
    graph.create_edge(nodes[43], nodes[47]);
    graph.create_edge(nodes[44], nodes[45]);
    graph.create_edge(nodes[44], nodes[46]);
    graph.create_edge(nodes[45], nodes[46]);
    graph.create_edge(nodes[46], nodes[47]);
    graph.create_edge(nodes[47], nodes[48]);
    graph.create_edge(nodes[47], nodes[49]);
    graph.create_edge(nodes[48], nodes[49]);
    
    graph.create_edge(nodes[49], nodes[50]);
    
    graph.create_edge(nodes[50], nodes[51]);
    graph.create_edge(nodes[50], nodes[52]);
    graph.create_edge(nodes[51], nodes[53]);
    graph.create_edge(nodes[52], nodes[53]);
    graph.create_edge(nodes[53], nodes[54]);
    graph.create_edge(nodes[53], nodes[57]);
    graph.create_edge(nodes[54], nodes[55]);
    graph.create_edge(nodes[54], nodes[56]);
    graph.create_edge(nodes[55], nodes[56]);
    graph.create_edge(nodes[56], nodes[57]);
    graph.create_edge(nodes[57], nodes[58]);
    graph.create_edge(nodes[57], nodes[59]);
    graph.create_edge(nodes[58], nodes[59]);
    
    // TODO one of these should really be the reference but idk how to add reference paths to a graph
    std::vector<std::vector<std::size_t>> paths_seqs = { {0, 1, 3, 4, 5, 6, 7,9}, {0, 1, 3, 4, 6, 7}, {0, 2, 3, 7}, {0, 2, 3, 4, 6, 7},
                                                          {10, 11, 13, 14, 15, 16, 17, 19}, {10, 11, 13, 14, 16, 17}, {10, 12, 13, 17}, {10, 12, 13, 14, 16, 17},
                                                          {20, 21, 23, 24, 25, 26, 27, 29}, {20, 21, 23, 24, 26, 27}, {20, 22, 23, 27}, {20, 22, 23, 24, 26, 27},
                                                          {30, 31, 33, 34, 35, 36, 37, 39}, {30, 31, 33, 34, 36, 37}, {30, 32, 33, 37}, {30, 32, 33, 34, 36, 37},
                                                          {39, 40, 41, 43, 44, 45, 46, 47, 49}, {40, 41, 43, 44, 46, 47}, {40, 42, 43, 47}, {40, 42, 43, 44, 46, 47},
                                                          {49, 50, 51, 53, 54, 55, 56, 57, 59}, {50, 51, 53, 54, 56, 57}, {50, 52, 53, 57}, {50, 52, 53, 54, 56, 57}};

    std::vector<handlegraph::path_handle_t> paths;

    for (int path_i = 0 ; path_i < paths_seqs.size() ; path_i++) {
        //idk where the chromosme is supposed to go either
        size_t chr = path_i / 4;
        size_t sample_num = path_i % 4; 
        if (sample_num == 0) {
            paths.emplace_back(graph.create_path_handle("path"+std::to_string(sample_num)+"#" + std::to_string(chr) + "#path0"));
        } else {
            paths.emplace_back(graph.create_path_handle("path"+std::to_string(sample_num)+"#" + std::to_string(chr) + "#0#0"));
        }
        for (size_t node_i : paths_seqs[path_i]) {
            graph.append_step(paths.back(), nodes[node_i]);
        }
    }

    //// vg isn't included so the distance index can only be built from the command line
    //graph.serialize("../tests/graph_test/simple_nested_chains.hg");
    //int built = system("vg index -j ../tests/graph_test/simple_nested_chains.dist ../tests/graph_test/simple_nested_chains.hg"); 

    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/graph_test/simple_nested_chains.dist");

   // bdsg::HashGraph graph;
   // graph.deserialize("../tests/graph_test/simple_nested_chains.hg");

    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&graph);

    SECTION("Make and fill in snarl collection with no partitions, then fill in partitions later by chromosome") {

        // Don't get the partitions or anything else
        TestSnarlDataCollection snarl_collection(1,10,1,10);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, 
            true, // partitions before walks but it doesn't matter here
            true, // get walks 
            SnarlDataCollection::get_all_walks_through_snarl,
            false, // don't get partitions 
            [&](const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                                                         const net_handle_t& snarl, const snarl_info_t& snarl_data,
                                                         std::vector<std::set<sample_hap_t>>& partitions) { 
                return;
            },
            true, // get sequences
            "", false);




        // Make partitions for each snarl
        std::vector<stoat::sample_hap_t> sample_haps;
        for (const auto& path : paths) {
            sample_haps.emplace_back(stoat::get_sample_and_haplotype(*path_graph, path));
        }


        // First fill it in with a non-existent path, which should do nothing
        std::unordered_map<stoat::sample_hap_t, size_t> sample_haplotype_to_index;
        snarl_collection.add_snarl_partitions(sample_haplotype_to_index, [&](const snarl_info_t& snarl_info) {

            //Add partitions to follow the walks
            std::vector<std::set<sample_hap_t>> partitions;
            if (snarl_info.start_node == Node_traversal_t(1, false) || snarl_info.start_node == stoat::Node_traversal_t(4, true) || 
                snarl_info.start_node == stoat::Node_traversal_t(4, false) || snarl_info.start_node == stoat::Node_traversal_t(8, true) ||
                snarl_info.start_node == stoat::Node_traversal_t(5, false) || snarl_info.start_node == stoat::Node_traversal_t(7, true) || 
                snarl_info.start_node == stoat::Node_traversal_t(8, false) || snarl_info.start_node == stoat::Node_traversal_t(10, true)) {
                // First connected component
                partitions.back().emplace(sample_haps[0]);
            } else if (snarl_info.start_node == Node_traversal_t(11, false) || snarl_info.start_node == stoat::Node_traversal_t(14, true) || 
                snarl_info.start_node == stoat::Node_traversal_t(14, false) || snarl_info.start_node == stoat::Node_traversal_t(18, true) ||
                snarl_info.start_node == stoat::Node_traversal_t(15, false) || snarl_info.start_node == stoat::Node_traversal_t(17, true) || 
                snarl_info.start_node == stoat::Node_traversal_t(18, false) || snarl_info.start_node == stoat::Node_traversal_t(20, true)) {
                // Second connected component
                partitions.back().emplace(sample_haps[4]);
            } else if (snarl_info.start_node == Node_traversal_t(21, false) || snarl_info.start_node == stoat::Node_traversal_t(24, true) || 
                snarl_info.start_node == stoat::Node_traversal_t(24, false) || snarl_info.start_node == stoat::Node_traversal_t(28, true) ||
                snarl_info.start_node == stoat::Node_traversal_t(25, false) || snarl_info.start_node == stoat::Node_traversal_t(27, true) || 
                snarl_info.start_node == stoat::Node_traversal_t(28, false) || snarl_info.start_node == stoat::Node_traversal_t(30, true)) {
                // Third connected component
                partitions.back().emplace(sample_haps[8]);
            } else {
                // I'm not going to bother testing anything else
            }
            return partitions;
        }, "empty_path");
        REQUIRE(sample_haplotype_to_index.size() == 0);

        // Make sure the partitions haven't been filled in
        snarl_collection.for_each_snarl([&](const snarl_info_t& snarl_info) {
            REQUIRE(snarl_info.partitions.size() == 0);
        });

        // Now fill it in with just the reference path for just the first chromosome
        snarl_collection.add_snarl_partitions(sample_haplotype_to_index, [&](const snarl_info_t& snarl_info) {

            //Add partitions to follow the walks
            std::vector<std::set<sample_hap_t>> partitions;
            if (snarl_info.start_node == Node_traversal_t(1, false) || snarl_info.start_node == stoat::Node_traversal_t(4, true) || 
                snarl_info.start_node == stoat::Node_traversal_t(4, false) || snarl_info.start_node == stoat::Node_traversal_t(8, true) ||
                snarl_info.start_node == stoat::Node_traversal_t(5, false) || snarl_info.start_node == stoat::Node_traversal_t(7, true) || 
                snarl_info.start_node == stoat::Node_traversal_t(8, false) || snarl_info.start_node == stoat::Node_traversal_t(10, true)) {
                // First connected component
                partitions.emplace_back();
                partitions.back().emplace(sample_haps[0]);
            } else if (snarl_info.start_node == Node_traversal_t(11, false) || snarl_info.start_node == stoat::Node_traversal_t(14, true) || 
                snarl_info.start_node == stoat::Node_traversal_t(14, false) || snarl_info.start_node == stoat::Node_traversal_t(18, true) ||
                snarl_info.start_node == stoat::Node_traversal_t(15, false) || snarl_info.start_node == stoat::Node_traversal_t(17, true) || 
                snarl_info.start_node == stoat::Node_traversal_t(18, false) || snarl_info.start_node == stoat::Node_traversal_t(20, true)) {
                // Second connected component
                partitions.emplace_back();
                partitions.back().emplace(sample_haps[4]);
            } else if (snarl_info.start_node == Node_traversal_t(21, false) || snarl_info.start_node == stoat::Node_traversal_t(24, true) || 
                snarl_info.start_node == stoat::Node_traversal_t(24, false) || snarl_info.start_node == stoat::Node_traversal_t(28, true) ||
                snarl_info.start_node == stoat::Node_traversal_t(25, false) || snarl_info.start_node == stoat::Node_traversal_t(27, true) || 
                snarl_info.start_node == stoat::Node_traversal_t(28, false) || snarl_info.start_node == stoat::Node_traversal_t(30, true)) {
                // Third connected component
                partitions.emplace_back();
                partitions.back().emplace(sample_haps[8]);
            } else {
                // I'm not going to bother testing anything else
            }
            return partitions;
        }, "path0#0#path0");

        // Make sure that just the partitions on the first chromosome been filled in
        snarl_collection.for_each_snarl([&](const snarl_info_t& snarl_info) {
            if (snarl_info.start_node == Node_traversal_t(1, false) || snarl_info.start_node == stoat::Node_traversal_t(4, true) || 
                snarl_info.start_node == stoat::Node_traversal_t(4, false) || snarl_info.start_node == stoat::Node_traversal_t(8, true) ||
                snarl_info.start_node == stoat::Node_traversal_t(5, false) || snarl_info.start_node == stoat::Node_traversal_t(7, true) || 
                snarl_info.start_node == stoat::Node_traversal_t(8, false) || snarl_info.start_node == stoat::Node_traversal_t(10, true)) {
                REQUIRE(snarl_info.partitions.size() == 1);
            } else {
                REQUIRE(snarl_info.partitions.size() == 0);
            }
        });

        // Now fill it in with just the reference path for just the second chromosome
        snarl_collection.add_snarl_partitions(sample_haplotype_to_index, [&](const snarl_info_t& snarl_info) {

            //Add partitions to follow the walks
            std::vector<std::set<sample_hap_t>> partitions;
            if (snarl_info.start_node == Node_traversal_t(1, false) || snarl_info.start_node == stoat::Node_traversal_t(4, true) || 
                snarl_info.start_node == stoat::Node_traversal_t(4, false) || snarl_info.start_node == stoat::Node_traversal_t(8, true) ||
                snarl_info.start_node == stoat::Node_traversal_t(5, false) || snarl_info.start_node == stoat::Node_traversal_t(7, true) || 
                snarl_info.start_node == stoat::Node_traversal_t(8, false) || snarl_info.start_node == stoat::Node_traversal_t(10, true)) {
                // First connected component
                partitions.emplace_back();
                partitions.back().emplace(sample_haps[0]);
            } else if (snarl_info.start_node == Node_traversal_t(11, false) || snarl_info.start_node == stoat::Node_traversal_t(14, true) || 
                snarl_info.start_node == stoat::Node_traversal_t(14, false) || snarl_info.start_node == stoat::Node_traversal_t(18, true) ||
                snarl_info.start_node == stoat::Node_traversal_t(15, false) || snarl_info.start_node == stoat::Node_traversal_t(17, true) || 
                snarl_info.start_node == stoat::Node_traversal_t(18, false) || snarl_info.start_node == stoat::Node_traversal_t(20, true)) {
                // Second connected component
                partitions.emplace_back();
                partitions.back().emplace(sample_haps[4]);
            } else if (snarl_info.start_node == Node_traversal_t(21, false) || snarl_info.start_node == stoat::Node_traversal_t(24, true) || 
                snarl_info.start_node == stoat::Node_traversal_t(24, false) || snarl_info.start_node == stoat::Node_traversal_t(28, true) ||
                snarl_info.start_node == stoat::Node_traversal_t(25, false) || snarl_info.start_node == stoat::Node_traversal_t(27, true) || 
                snarl_info.start_node == stoat::Node_traversal_t(28, false) || snarl_info.start_node == stoat::Node_traversal_t(30, true)) {
                // Third connected component
                partitions.emplace_back();
                partitions.back().emplace(sample_haps[8]);
            } else {
                // I'm not going to bother testing anything else
            }
            return partitions;
        }, "path0#1#path0");

        // Make sure that just the partitions on the first chromosome been filled in
        snarl_collection.for_each_snarl([&](const snarl_info_t& snarl_info) {
            if (snarl_info.start_node == Node_traversal_t(1, false) || snarl_info.start_node == stoat::Node_traversal_t(4, true) || 
                snarl_info.start_node == stoat::Node_traversal_t(4, false) || snarl_info.start_node == stoat::Node_traversal_t(8, true) ||
                snarl_info.start_node == stoat::Node_traversal_t(5, false) || snarl_info.start_node == stoat::Node_traversal_t(7, true) || 
                snarl_info.start_node == stoat::Node_traversal_t(8, false) || snarl_info.start_node == stoat::Node_traversal_t(10, true) ||
                snarl_info.start_node == Node_traversal_t(11, false) || snarl_info.start_node == stoat::Node_traversal_t(14, true) || 
                snarl_info.start_node == stoat::Node_traversal_t(14, false) || snarl_info.start_node == stoat::Node_traversal_t(18, true) ||
                snarl_info.start_node == stoat::Node_traversal_t(15, false) || snarl_info.start_node == stoat::Node_traversal_t(17, true) || 
                snarl_info.start_node == stoat::Node_traversal_t(18, false) || snarl_info.start_node == stoat::Node_traversal_t(20, true)) {
                REQUIRE(snarl_info.partitions.size() == 1);
            } else {
                REQUIRE(snarl_info.partitions.size() == 0);
            }
        });

    }
}

