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

    // If the walks were found without partitions then there is a walk for the last snarl
    auto check_collection = [&] (const TestSnarlDataCollection& snarl_collection, bool check_partitions, bool check_sequences, bool get_all_walks) {

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
                REQUIRE(snarl_info.snarl_walks.size() == 2);
                REQUIRE(snarl_info.snarl_walks[0].get_paths().size() == 3);
                REQUIRE(snarl_info.snarl_walks[1].get_paths().size() == 3);

                // The paths should be 1-2-4 and 1-3-4
                if (check_partitions) {
                    REQUIRE(snarl_info.partitions.size() == 2);
                }
                REQUIRE(snarl_info.snarl_walks[0].get_paths()[0] == Node_traversal_t(1, false));
                REQUIRE(snarl_info.snarl_walks[1].get_paths()[0] == Node_traversal_t(1, false));
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

                REQUIRE(snarl_info.snarl_walks[0].get_paths()[2] == Node_traversal_t(4, false));
                REQUIRE(snarl_info.snarl_walks[1].get_paths()[2] == Node_traversal_t(4, false));

                REQUIRE(snarl_info.variant_type == "1,1");

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

                REQUIRE(snarl_info.snarl_walks.size() == 2);
                size_t insertion_index;
                size_t deletion_index;
                if (check_partitions) {
                    REQUIRE(snarl_info.partitions.size() == 2);
                }
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

            } else if ((snarl_info.start_node == stoat::Node_traversal_t(8, false) && snarl_info.end_node == stoat::Node_traversal_t(10, true)) ||
                (snarl_info.start_node == stoat::Node_traversal_t(10, true) && snarl_info.end_node == stoat::Node_traversal_t(8, false))) {
                REQUIRE(snarl_info.ref_path == "NA");
                REQUIRE(snarl_info.start_position == 0);
                REQUIRE(snarl_info.end_position == 0);
                REQUIRE(snarl_info.depth == 1);
                if (get_all_walks) {
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


    SECTION("Make and fill in snarl collection with no partitions") {
        // Don't get the partitions or anything else
        TestSnarlDataCollection snarl_collection(1,10,1,10);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, false, [&](const handlegraph::net_handle_t& net) { 
            std::vector<std::set<sample_hap_t>> partitions;
            return partitions;
        }, true, "", false, cout);

        check_collection(snarl_collection, false, true, true);
    }

    SECTION("Make and fill in snarl collection with no partitions, then fill in partitions later") {

        // Don't get the partitions or anything else
        TestSnarlDataCollection snarl_collection(1,10,1,10);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, false, [&](const handlegraph::net_handle_t& net) { 
            std::vector<std::set<sample_hap_t>> partitions;
            return partitions;
        }, true, "", false, cout);


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

        check_collection(snarl_collection, true, true, true);

    }
    SECTION("Make and fill in snarl collection with partitions") {

        // Make partitions for each snarl
        std::vector<stoat::sample_hap_t> sample_haps;
        for (const auto& path : paths) {
            sample_haps.emplace_back(stoat::get_sample_and_haplotype(*path_graph, path));
        }


        TestSnarlDataCollection snarl_collection(1,10,1,10);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, true, [&](const handlegraph::net_handle_t& net) { 
            std::vector<std::set<sample_hap_t>> partitions;
            if (net == distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(2)))) {
                // snarl 1-4
                //[0,1], [2,3]
                partitions.emplace_back();
                partitions.back().emplace(sample_haps[0]);
                partitions.back().emplace(sample_haps[1]);
                partitions.emplace_back();
                partitions.back().emplace(sample_haps[2]);
                partitions.back().emplace(sample_haps[3]);
            } else if (net == distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(5)))) {
                // snarl 4-8
                //[0,1,3], [2]
                partitions.emplace_back();
                partitions.back().emplace(sample_haps[0]);
                partitions.back().emplace(sample_haps[1]);
                partitions.back().emplace(sample_haps[3]);
                partitions.emplace_back();
                partitions.back().emplace(sample_haps[2]);
            } else if (net == distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(6)))) {
                // snarl 5-7
                //[0], [1,3]
                partitions.emplace_back();
                partitions.back().emplace(sample_haps[0]);
                partitions.emplace_back();
                partitions.back().emplace(sample_haps[1]);
                partitions.back().emplace(sample_haps[3]);
            } else {
                // For the last one there are no partitions because there are no paths
                std::vector<std::set<sample_hap_t>> partitions;
                return partitions;
            }
            return partitions;
        }, true, "", false, cout);

        check_collection(snarl_collection, true, true, false);

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

            check_collection(loaded_snarl_collection, true, true, false);

            std::string rm_cmd = "rm " + test_file;

            int rm = system(rm_cmd.c_str()); 
        }
    }

}

}
