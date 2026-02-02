#include <catch.hpp>
#include <bdsg/hash_graph.hpp>
#include <bdsg/overlays/overlay_helper.hpp>
#include "../../src/snarl_data_collection.hpp"
#include "../../src/log.hpp"

using namespace stoat;


class TestSnarlDataCollection : public SnarlDataCollection {
    public: 
    TestSnarlDataCollection(size_t allele_size_limit, size_t snarl_child_limit, size_t walk_steps_limit) :
        SnarlDataCollection(allele_size_limit, snarl_child_limit, walk_steps_limit) {} 
    using SnarlDataCollection::fill_in_snarl_info;
    using SnarlDataCollection::add_alleles_by_sample;
    using SnarlDataCollection::for_each_snarl;
    using SnarlDataCollection::write_snarl_data_collection;
    using SnarlDataCollection::load_snarl_data_collection;
    using SnarlDataCollection::is_equivalent;
};

TEST_CASE( "Snarl collection one node", "[snarl_collection]" ) {


    bdsg::HashGraph graph;
        
    //handlegraph::handle_t n1 = graph.create_handle("GCAAACAGATT");

    //handlegraph::path_handle_t path = graph.create_path_handle("path");
    //graph.append_step(path, n1);

    // vg isn't included so the distance index can only be built from the command line
    //graph.serialize("../tests/test_data/test_graphs/one_node.hg");
    //int built = system("vg index -j ../tests/test_data/test_graphs/one_node.dist ../tests/test_data/test_graphs/one_node.hg"); 
    bdsg::SnarlDistanceIndex distance_index;
    graph.deserialize("../tests/test_data/test_graphs/one_node.hg");
    distance_index.deserialize("../tests/test_data/test_graphs/one_node.dist");

    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&graph);

    std::vector<stoat::sample_hap_t> all_samples ({stoat::sample_hap_t(*path_graph,graph.get_path_handle("path"))});
    std::unordered_map<std::string, size_t> sample_to_index;

    SECTION("Make and fill in snarl collection") {
        // There isn't much to do with one node so just make sure we can run the constructor without crashing
        TestSnarlDataCollection snarl_collection(1,1,1);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, all_samples, 
            true, // alleles before walks but it doesn't matter here
            false, // don't get walks 
            [&](const net_handle_t& snarl, const snarl_info_t& snarl_data,std::vector<stoat::PathTraversal>& walks) { 
                return;
            }, false, // don't get alleles 
            [&](const net_handle_t& snarl, const snarl_info_t& snarl_data, const std::vector<sample_hap_t>& samples) { 
                return std::vector<size_t>();
            },
            false, // don't get sequences
            std::unordered_set<std::string>(), false);

    }
    SECTION("Serialize snarl collection") {
        // There isn't much to do with one node so just make sure we can run the constructor without crashing

        TestSnarlDataCollection snarl_collection(1,1,1);

        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, all_samples, 
            true, // alleles before walks but it doesn't matter here
            false, // don't get walks 
            [&](const net_handle_t& snarl, const snarl_info_t& snarl_data,std::vector<stoat::PathTraversal>& walks) { 
                return;
            }, false, // don't get alleles 
            [&](const net_handle_t& snarl, const snarl_info_t& snarl_data, const std::vector<sample_hap_t>& samples) { 
                return std::vector<size_t>();
            },
            false, // don't get sequences
            std::unordered_set<std::string>(), false);

        std::string test_file = "./test.snarl_collection.txt";

        std::ofstream outstream;
        outstream.open(test_file);
        snarl_collection.write_snarl_data_collection(outstream);
        outstream.close();
        
        TestSnarlDataCollection snarl_collection_loaded(1,1,1);
        snarl_collection_loaded.load_snarl_data_collection(test_file);

        REQUIRE(SnarlDataCollection::is_equivalent(snarl_collection, snarl_collection_loaded));

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

    // If the walks were found without allele then there is a walk for the last snarl
    auto check_collection = [&] (const TestSnarlDataCollection& snarl_collection, bool check_walks, bool check_alleles, bool check_sequences, bool get_all_walks) {

        // Get the alleles again so we can check them. the order might be different though

        // Make alleles for each snarl
        std::vector<stoat::sample_hap_t> sample_haps;
        for (const auto& path : paths) {
            sample_haps.emplace_back(stoat::sample_hap_t(*path_graph, path));
        }


        // Check that we got all snarls and that we got the correct snarls
        size_t snarl_count = 0;

        snarl_collection.for_each_snarl([&](const snarl_info_t& snarl_info) {
            snarl_count++;

            if ((snarl_info.start_node == node_traversal_t(1, false) && snarl_info.end_node == stoat::node_traversal_t(4, true)) ||
                (snarl_info.start_node == stoat::node_traversal_t(4, true) && snarl_info.end_node == node_traversal_t(1, false))) {
                // First snarl

                REQUIRE(snarl_info.ref_path == "path0#0#path0");
                REQUIRE(snarl_info.start_position == 1);
                REQUIRE(snarl_info.end_position == 2);
                REQUIRE(snarl_info.depth == 1);

                // The paths should be 1-2-4 and 1-3-4
                if (check_alleles) {

                    REQUIRE(snarl_info.alleles_by_sample.allele_count == 2);
                    //[0,1], [2,3]
                    std::string genotype0 =  snarl_info.genotypes.get_genotype_as_string("path0");
                    std::string genotype1 =  snarl_info.genotypes.get_genotype_as_string("path1");
                    std::string genotype2 =  snarl_info.genotypes.get_genotype_as_string("path2");
                    std::string genotype3 =  snarl_info.genotypes.get_genotype_as_string("path3");
                    REQUIRE(genotype0 == genotype1);
                    REQUIRE(genotype1 != genotype2);
                    REQUIRE(genotype2 == genotype3);

                    // Each sample has just one count of one
                    REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path0", 0) != snarl_info.genotypes.get_count_for_sample_and_allele("path0", 1));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path0", 0) == 0 || snarl_info.genotypes.get_count_for_sample_and_allele("path0", 1) == 0 ));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path0", 0) == 1 || snarl_info.genotypes.get_count_for_sample_and_allele("path0", 1) == 1 ));

                    REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path1", 0) != snarl_info.genotypes.get_count_for_sample_and_allele("path1", 1));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path1", 0) == 0 || snarl_info.genotypes.get_count_for_sample_and_allele("path1", 1) == 0 ));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path1", 0) == 1 || snarl_info.genotypes.get_count_for_sample_and_allele("path1", 1) == 1 ));

                    REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path2", 0) != snarl_info.genotypes.get_count_for_sample_and_allele("path2", 1));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path2", 0) == 0 || snarl_info.genotypes.get_count_for_sample_and_allele("path2", 1) == 0 ));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path2", 0) == 1 || snarl_info.genotypes.get_count_for_sample_and_allele("path2", 1) == 1 ));

                    REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path3", 0) != snarl_info.genotypes.get_count_for_sample_and_allele("path3", 1));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path3", 0) == 0 || snarl_info.genotypes.get_count_for_sample_and_allele("path3", 1) == 0 ));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path3", 0) == 1 || snarl_info.genotypes.get_count_for_sample_and_allele("path3", 1) == 1 ));
                }
                if (check_walks) {
                    REQUIRE(snarl_info.walks_by_allele.size() == 2);
                    REQUIRE(snarl_info.walks_by_allele[0].get_path().size() == 3);
                    REQUIRE(snarl_info.walks_by_allele[1].get_path().size() == 3);

                    REQUIRE(snarl_info.walks_by_allele[0].get_path()[0] == node_traversal_t(1, false));
                    REQUIRE(snarl_info.walks_by_allele[1].get_path()[0] == node_traversal_t(1, false));

                    // The second node is either 2 or 3
                    REQUIRE((snarl_info.walks_by_allele[0].get_path()[1] == node_traversal_t(2, false) ||
                             snarl_info.walks_by_allele[0].get_path()[1] == node_traversal_t(3, false)));

                    REQUIRE((snarl_info.walks_by_allele[1].get_path()[1] == node_traversal_t(2, false) ||
                             snarl_info.walks_by_allele[1].get_path()[1] == node_traversal_t(3, false)));

                    REQUIRE(!(snarl_info.walks_by_allele[0].get_path()[1] == snarl_info.walks_by_allele[1].get_path()[1]));

                    REQUIRE(snarl_info.walks_by_allele[0].get_path()[2] == node_traversal_t(4, false));
                    REQUIRE(snarl_info.walks_by_allele[1].get_path()[2] == node_traversal_t(4, false));
                    REQUIRE(stoat::vectorPathToString(snarl_info.walks_by_allele, true) == "1,1");
                }

                // Check that the walks and alleles match
                if (check_walks && check_alleles) {

                    if (snarl_info.walks_by_allele[0].get_path()[1] == node_traversal_t(2, false)) {
                        REQUIRE(snarl_info.walks_by_allele[1].get_path()[1] == node_traversal_t(3, false));
                        // If allele 0 is 2 and allele 1 is 3
                        if (check_alleles) {
                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path0", 0) == 1);
                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path0", 1) == 0);

                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path1", 0) == 1);
                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path1", 1) == 0);

                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path2", 0) == 0);
                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path2", 1) == 1);

                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path3", 0) == 0);
                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path3", 1) == 1);

                        }
                    } else {
                        // If allele 0 is 3 and allele 1 is 1
                        REQUIRE(snarl_info.walks_by_allele[0].get_path()[1] == node_traversal_t(3, false));
                        REQUIRE(snarl_info.walks_by_allele[1].get_path()[1] == node_traversal_t(2, false));
                        if (check_alleles) {
                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path0", 0) == 0);
                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path0", 1) == 1);

                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path1", 0) == 0);
                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path1", 1) == 1);

                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path2", 0) == 1);
                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path2", 1) == 0);

                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path3", 0) == 1);
                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path3", 1) == 0);
                        }
                    }
                }


                if (check_sequences) {
                    // The sequences should both be CCA, oops
                    REQUIRE(snarl_info.sequences_by_allele.size() == 2);
                    REQUIRE(snarl_info.sequences_by_allele[0] == "C");
                    REQUIRE(snarl_info.sequences_by_allele[1] == "C");
                }

            } else if ((snarl_info.start_node == stoat::node_traversal_t(4, false) && snarl_info.end_node == stoat::node_traversal_t(8, true)) ||
                (snarl_info.start_node == stoat::node_traversal_t(8, true) && snarl_info.end_node == stoat::node_traversal_t(4, false))) {
                // Snarl 4-8
                REQUIRE(snarl_info.ref_path == "path0#0#path0");
                REQUIRE(snarl_info.start_position == 3);
                REQUIRE(snarl_info.end_position == 6);
                REQUIRE(snarl_info.depth == 1);

                size_t insertion_index;
                size_t deletion_index;
                if (check_alleles) {
                    //[0,1,3], [2]
                    REQUIRE(snarl_info.alleles_by_sample.allele_count == 2);
                    std::string genotype0 =  snarl_info.genotypes.get_genotype_as_string("path0");
                    std::string genotype1 =  snarl_info.genotypes.get_genotype_as_string("path1");
                    std::string genotype2 =  snarl_info.genotypes.get_genotype_as_string("path2");
                    std::string genotype3 =  snarl_info.genotypes.get_genotype_as_string("path3");
                    REQUIRE(genotype0 == genotype1);
                    REQUIRE(genotype1 == genotype3);
                    REQUIRE(genotype2 != genotype3);

                    // Each sample has just one count of one
                    REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path0", 0) != snarl_info.genotypes.get_count_for_sample_and_allele("path0", 1));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path0", 0) == 0 || snarl_info.genotypes.get_count_for_sample_and_allele("path0", 1) == 0 ));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path0", 0) == 1 || snarl_info.genotypes.get_count_for_sample_and_allele("path0", 1) == 1 ));

                    REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path1", 0) != snarl_info.genotypes.get_count_for_sample_and_allele("path1", 1));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path1", 0) == 0 || snarl_info.genotypes.get_count_for_sample_and_allele("path1", 1) == 0 ));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path1", 0) == 1 || snarl_info.genotypes.get_count_for_sample_and_allele("path1", 1) == 1 ));

                    REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path2", 0) != snarl_info.genotypes.get_count_for_sample_and_allele("path2", 1));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path2", 0) == 0 || snarl_info.genotypes.get_count_for_sample_and_allele("path2", 1) == 0 ));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path2", 0) == 1 || snarl_info.genotypes.get_count_for_sample_and_allele("path2", 1) == 1 ));

                    REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path3", 0) != snarl_info.genotypes.get_count_for_sample_and_allele("path3", 1));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path3", 0) == 0 || snarl_info.genotypes.get_count_for_sample_and_allele("path3", 1) == 0 ));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path3", 0) == 1 || snarl_info.genotypes.get_count_for_sample_and_allele("path3", 1) == 1 ));
                }
                if (check_walks) {
                    REQUIRE(snarl_info.walks_by_allele.size() == 2);
                    if (snarl_info.walks_by_allele[0].get_path().size() == 2) {
                        // First is deletion, second is insertion
                        REQUIRE(snarl_info.walks_by_allele[1].get_path().size() == 5);
                        deletion_index = 0;
                        insertion_index = 1;
                        REQUIRE(stoat::vectorPathToString(snarl_info.walks_by_allele, true) == "0,2/3");
                    } else {
                        // First is insertion, second is deletion
                        REQUIRE(snarl_info.walks_by_allele[0].get_path().size() == 5);
                        REQUIRE(snarl_info.walks_by_allele[1].get_path().size() == 2);
                        deletion_index = 1;
                        insertion_index = 0;
                        REQUIRE(stoat::vectorPathToString(snarl_info.walks_by_allele, true) == "2/3,0");
                    }
               

                    // The paths should be 4-8 and 4-5-0-7-8
                    REQUIRE(snarl_info.walks_by_allele[deletion_index].get_path()[0] == node_traversal_t(4, false));
                    REQUIRE(snarl_info.walks_by_allele[deletion_index].get_path()[1] == node_traversal_t(8, false));

                    REQUIRE(snarl_info.walks_by_allele[insertion_index].get_path()[0] == node_traversal_t(4, false));
                    REQUIRE(snarl_info.walks_by_allele[insertion_index].get_path()[1] == node_traversal_t(5, false));
                    REQUIRE(snarl_info.walks_by_allele[insertion_index].get_path()[2].get_node_id() == 0);
                    REQUIRE(snarl_info.walks_by_allele[insertion_index].get_path()[3] == node_traversal_t(7, false));
                    REQUIRE(snarl_info.walks_by_allele[insertion_index].get_path()[4] == node_traversal_t(8, false));

                    if (check_sequences) {
                        // The sequences 
                        REQUIRE(snarl_info.sequences_by_allele.size() == 2);
                        REQUIRE(snarl_info.sequences_by_allele[insertion_index] == "TNA");
                        REQUIRE(snarl_info.sequences_by_allele[deletion_index] == "");
                    }
                    if (check_alleles) {
                        // 0,1,3 take the insertion, 2 takes the deletion
                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path0", insertion_index) == 1);
                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path0", deletion_index) == 0);

                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path1", insertion_index) == 1);
                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path1", deletion_index) == 0);

                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path2", insertion_index) == 0);
                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path2", deletion_index) == 1);

                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path3", insertion_index) == 1);
                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path3", deletion_index) == 0);

                    }
                }


            }  else if ((snarl_info.start_node == stoat::node_traversal_t(5, false) && snarl_info.end_node == stoat::node_traversal_t(7, true)) ||
                (snarl_info.start_node == stoat::node_traversal_t(7, true) && snarl_info.end_node == stoat::node_traversal_t(5, false))) {
                REQUIRE(snarl_info.ref_path == "path0#0#path0");
                REQUIRE(snarl_info.start_position == 4);
                REQUIRE(snarl_info.end_position == 5);
                REQUIRE(snarl_info.depth == 2);


                if (check_walks) {
                    REQUIRE(snarl_info.walks_by_allele.size() == 2);
                    size_t insertion_index;
                    size_t deletion_index;
                    if (snarl_info.walks_by_allele[0].get_path().size() == 2) {
                        // First is deletion, second is insertion
                        REQUIRE(snarl_info.walks_by_allele[1].get_path().size() == 3);
                        deletion_index = 0;
                        insertion_index = 1;
                        REQUIRE(stoat::vectorPathToString(snarl_info.walks_by_allele, true) == "0,1");
                    } else {
                        // First is insertion, second is deletion
                        REQUIRE(snarl_info.walks_by_allele[0].get_path().size() == 3);
                        REQUIRE(snarl_info.walks_by_allele[1].get_path().size() == 2);
                        deletion_index = 1;
                        insertion_index = 0;
                        REQUIRE(stoat::vectorPathToString(snarl_info.walks_by_allele, true) == "1,0");
                    }

                    // The paths should be 5-7 and 5-6-7
                    REQUIRE(snarl_info.walks_by_allele[deletion_index].get_path()[0] == node_traversal_t(5, false));
                    REQUIRE(snarl_info.walks_by_allele[deletion_index].get_path()[1] == node_traversal_t(7, false));

                    REQUIRE(snarl_info.walks_by_allele[insertion_index].get_path()[0] == node_traversal_t(5, false));
                    REQUIRE(snarl_info.walks_by_allele[insertion_index].get_path()[1] == node_traversal_t(6, false));
                    REQUIRE(snarl_info.walks_by_allele[insertion_index].get_path()[2] == node_traversal_t(7, false));
                    if (check_sequences) {
                        // The sequences 
                        REQUIRE(snarl_info.sequences_by_allele.size() == 2);
                        REQUIRE(snarl_info.sequences_by_allele[insertion_index] == "C");
                        REQUIRE(snarl_info.sequences_by_allele[deletion_index] == "");
                    }
                    if (check_alleles) {
                        REQUIRE(snarl_info.alleles_by_sample.allele_count == 2);
                        // 0 takes the insertion, 1,3 take the deletion 
                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path0", insertion_index) == 1);
                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path0", deletion_index) == 0);

                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path1", insertion_index) == 0);
                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path1", deletion_index) == 1);

                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path2", insertion_index) == 0);
                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path2", deletion_index) == 0);

                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path3", insertion_index) == 0);
                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path3", deletion_index) == 1);
                    }
                } else if (check_alleles) {
                    REQUIRE(snarl_info.alleles_by_sample.allele_count == 2);
                    // If we want to check the alleles but don't have the walks
                    //[0], [1,3]
                    std::string genotype0 =  snarl_info.genotypes.get_genotype_as_string("path0");
                    std::string genotype1 =  snarl_info.genotypes.get_genotype_as_string("path1");
                    std::string genotype2 =  snarl_info.genotypes.get_genotype_as_string("path2");
                    std::string genotype3 =  snarl_info.genotypes.get_genotype_as_string("path3");
                    REQUIRE(genotype0 != genotype1);
                    REQUIRE(genotype0 != genotype2);
                    REQUIRE(genotype1 == genotype3);
                    REQUIRE(genotype2 == "0,0");

                    REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path0", 0) != snarl_info.genotypes.get_count_for_sample_and_allele("path0", 1));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path0", 0) == 0 || snarl_info.genotypes.get_count_for_sample_and_allele("path0", 1) == 0 ));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path0", 0) == 1 || snarl_info.genotypes.get_count_for_sample_and_allele("path0", 1) == 1 ));

                    REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path1", 0) != snarl_info.genotypes.get_count_for_sample_and_allele("path1", 1));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path1", 0) == 0 || snarl_info.genotypes.get_count_for_sample_and_allele("path1", 1) == 0 ));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path1", 0) == 1 || snarl_info.genotypes.get_count_for_sample_and_allele("path1", 1) == 1 ));

                    REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path2", 0) == 0);
                    REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path2", 1) == 0 );

                    REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path3", 0) != snarl_info.genotypes.get_count_for_sample_and_allele("path3", 1));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path3", 0) == 0 || snarl_info.genotypes.get_count_for_sample_and_allele("path3", 1) == 0 ));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path3", 0) == 1 || snarl_info.genotypes.get_count_for_sample_and_allele("path3", 1) == 1 ));
                }



            } else if ((snarl_info.start_node == stoat::node_traversal_t(8, false) && snarl_info.end_node == stoat::node_traversal_t(10, true)) ||
                (snarl_info.start_node == stoat::node_traversal_t(10, true) && snarl_info.end_node == stoat::node_traversal_t(8, false))) {
                REQUIRE(snarl_info.ref_path == "NA");
                REQUIRE(snarl_info.start_position == 0);
                REQUIRE(snarl_info.end_position == 0);
                REQUIRE(snarl_info.depth == 1);
                if (get_all_walks && check_walks) {
                    REQUIRE(snarl_info.walks_by_allele.size() == 2);

                    size_t insertion_index;
                    size_t deletion_index;
                    if (snarl_info.walks_by_allele[0].get_path().size() == 2) {
                        // First is deletion, second is insertion
                        REQUIRE(snarl_info.walks_by_allele[1].get_path().size() == 3);
                        deletion_index = 0;
                        insertion_index = 1;
                        REQUIRE(stoat::vectorPathToString(snarl_info.walks_by_allele, true) == "0,1");
                    } else {
                        // First is insertion, second is deletion
                        REQUIRE(snarl_info.walks_by_allele[0].get_path().size() == 3);
                        REQUIRE(snarl_info.walks_by_allele[1].get_path().size() == 2);
                        deletion_index = 1;
                        insertion_index = 0;
                        REQUIRE(stoat::vectorPathToString(snarl_info.walks_by_allele, true) == "1,0");
                    }

                    // The paths should be 8-10 and 8-9-10
                    REQUIRE(snarl_info.walks_by_allele[deletion_index].get_path()[0] == node_traversal_t(8, false));
                    REQUIRE(snarl_info.walks_by_allele[deletion_index].get_path()[1] == node_traversal_t(10, false));

                    REQUIRE(snarl_info.walks_by_allele[insertion_index].get_path()[0] == node_traversal_t(8, false));
                    REQUIRE(snarl_info.walks_by_allele[insertion_index].get_path()[1] == node_traversal_t(9, false));
                    REQUIRE(snarl_info.walks_by_allele[insertion_index].get_path()[2] == node_traversal_t(10, false));

                    if (check_sequences) {
                        // The sequences 
                        REQUIRE(snarl_info.sequences_by_allele.size() == 2);
                        REQUIRE(snarl_info.sequences_by_allele[insertion_index] == "A");
                        REQUIRE(snarl_info.sequences_by_allele[deletion_index] == "");
                    }

                    if (check_alleles) {
                    }
                }

            } else {
                REQUIRE(false);
            }
        });

        REQUIRE(snarl_count == 4);
    };

    auto get_alleles_per_snarl = [&](const snarl_info_t& snarl_info, const std::vector<stoat::sample_hap_t>& haplotypes) {

        std::vector<size_t> alleles_by_snarl;

        if (snarl_info.start_node == node_traversal_t(1, false) || snarl_info.start_node == stoat::node_traversal_t(4, true)) {
            // snarl 1-4
            //[0,1], [2,3]
            if (snarl_info.walks_by_allele.size() == 0 || snarl_info.walks_by_allele[0].get_path()[1] == node_traversal_t(2, false)) {
                alleles_by_snarl = std::vector<size_t>({0,0,1,1});
            } else {
                alleles_by_snarl = std::vector<size_t>({1,1,0,0});
            }
        } else if (snarl_info.start_node == stoat::node_traversal_t(4, false) || snarl_info.start_node == stoat::node_traversal_t(8, true)) {
            // snarl 4-8
            //[0,1,3], [2]
            if (snarl_info.walks_by_allele.size() == 0 || snarl_info.walks_by_allele[0].get_path().size() == 5) {
                //insertion first
                alleles_by_snarl = std::vector<size_t>({0,0,1,0});
            } else {
                alleles_by_snarl = std::vector<size_t>({1,1,0,1});
            }
        } else if (snarl_info.start_node == stoat::node_traversal_t(5, false) || snarl_info.start_node == stoat::node_traversal_t(7, true)) {
    
            // snarl 5-7
            //[0], [1,3]
            if (snarl_info.walks_by_allele.size() == 0 || snarl_info.walks_by_allele[0].get_path().size() == 3) {

                alleles_by_snarl = std::vector<size_t>({0,1,std::numeric_limits<size_t>::max(),1});
            } else {
                alleles_by_snarl = std::vector<size_t>({1,0,std::numeric_limits<size_t>::max(),0});
            }
        } else {
            // For the last one there are no alleles because there are no paths
        }
        return alleles_by_snarl;
    };
    std::unordered_map<std::string, size_t> sample_to_index;
    sample_to_index.emplace("path0", 0);
    sample_to_index.emplace("path1", 1);
    sample_to_index.emplace("path2", 2);
    sample_to_index.emplace("path3", 3);

    std::vector<stoat::sample_hap_t> all_samples ({stoat::sample_hap_t(*path_graph, paths[0]),
                                         stoat::sample_hap_t(*path_graph, paths[1]),
                                         stoat::sample_hap_t(*path_graph, paths[2]),
                                         stoat::sample_hap_t(*path_graph, paths[3])});


    SECTION("Make and fill in snarl collection with no walks, alleles, or sequences") {

        TestSnarlDataCollection snarl_collection(1,10,10);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, all_samples, 
            true, // alleles before walks but it doesn't matter here
            false, // don't get walks 
            [&](const net_handle_t& snarl, const snarl_info_t& snarl_data, std::vector<stoat::PathTraversal>& walks) { 
                return;
            }, false, // don't get alleles 
            [&](const net_handle_t& snarl, const snarl_info_t& snarl_data, const std::vector<sample_hap_t>& sample_haps) { 
                return std::vector<size_t>();
            },
            false, // don't get sequences
            std::unordered_set<std::string>(), false);


        check_collection(snarl_collection, false, false, false, true);
    }
    SECTION("Make and fill in snarl collection with walks and sequences but no allele") {

        TestSnarlDataCollection snarl_collection(1,10,10);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, all_samples, 
            true, // allele before walks but it doesn't matter here
            true, // get walks 
            [&](const net_handle_t& snarl, const snarl_info_t& snarl_data, std::vector<stoat::PathTraversal>& walks) {
                return SnarlDataCollection::get_all_walks_through_snarl(*path_graph, distance_index, snarl, snarl_data, walks, 1);
            },
            false, // don't get allele 
            [&](const net_handle_t& snarl, const snarl_info_t& snarl_data, const std::vector<sample_hap_t>& sample_haps) { 
                return std::vector<size_t>();
            },
            true, // get sequences
            std::unordered_set<std::string>(), false);


        check_collection(snarl_collection, true, false, true, true);
    }

    SECTION("Make and fill in snarl collection with walks, sequences, and no alleles, then fill in alleles later") {

        // Don't get the alleles or anything else
        TestSnarlDataCollection snarl_collection(1,10,10);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, all_samples, 
            false, //walks before alleles but it doesn't matter here
            true, // get walks 
            [&](const net_handle_t& snarl, const snarl_info_t& snarl_data, std::vector<stoat::PathTraversal>& walks) {
                return SnarlDataCollection::get_all_walks_through_snarl(*path_graph, distance_index, snarl, snarl_data, walks, 1);
            },
            false, // don't get alleles 
            [&](const net_handle_t& snarl, const snarl_info_t& snarl_data,const std::vector<sample_hap_t>& samples) { 
                return std::vector<size_t>();
            },
            true, // get sequences
            std::unordered_set<std::string>(), false);


        check_collection(snarl_collection, true, false, true, true);




        std::unordered_map<stoat::sample_hap_t, size_t> sample_haplotype_to_index;
        snarl_collection.add_alleles_by_sample(get_alleles_per_snarl, "");

        check_collection(snarl_collection, true, true, true, true);

    }
    SECTION("Make and fill in snarl collection with no alleles, then fill in alleles later by chromosome") {

        // Don't get the alleles or anything else
        TestSnarlDataCollection snarl_collection(1,10,10);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, all_samples, 
            false, // walks before alleles but it doesn't matter here
            true, // get walks 
            [&](const net_handle_t& snarl, const snarl_info_t& snarl_data, std::vector<stoat::PathTraversal>& walks) {
                return SnarlDataCollection::get_all_walks_through_snarl(*path_graph, distance_index, snarl, snarl_data, walks, 0);
            },
            false, // don't get alleles 
            [&](const net_handle_t& snarl, const snarl_info_t& snarl_data,const std::vector<sample_hap_t>& samples) { 
                return std::vector<size_t>();
            },
            true, // get sequences
            std::unordered_set<std::string>(), false);

        // First fill it in with a non-existent path, which should do nothing
        std::unordered_map<stoat::sample_hap_t, size_t> sample_haplotype_to_index;
        snarl_collection.add_alleles_by_sample(get_alleles_per_snarl, "empty_path");

        REQUIRE(sample_haplotype_to_index.size() == 0);

        snarl_collection.for_each_snarl([&](const snarl_info_t& snarl_info) {
            // Make sure the alleles haven't been filled in
            REQUIRE(snarl_info.genotypes.get_allele_count() == 0);
        });

        // Now fill it in with the paths
        snarl_collection.add_alleles_by_sample(get_alleles_per_snarl, "path0#0#path0");

        check_collection(snarl_collection, true, true, true, true);

    }
    SECTION("Make and fill in snarl collection with walks, alleles, and sequences") {



        TestSnarlDataCollection snarl_collection(1,10,10);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, all_samples, 
            true, // alleles before walks
            true, // get walks 
            [&](const net_handle_t& snarl, const snarl_info_t& snarl_data, std::vector<PathTraversal>& walks) {
                SnarlDataCollection::get_walks_from_alleles(*path_graph, distance_index, snarl, snarl_data, walks);
                return;
            },
            true, // get alleles 
            [&](const net_handle_t& snarl, const snarl_info_t& snarl_data, const std::vector<sample_hap_t>& samples) { 
                std::vector<size_t> alleles = get_alleles_per_snarl(snarl_data, samples);
                return alleles;
            },
            true, // get sequences
            std::unordered_set<std::string>(), false);

        check_collection(snarl_collection, true, true, true, false);

        SECTION("Serialize it") {

            std::string test_file = "./test_snarls.txt";
            std::ofstream outstream;
            outstream.open(test_file);
            snarl_collection.write_snarl_data_collection(outstream);
            outstream.close();

            TestSnarlDataCollection loaded_snarl_collection(1,10,10);
            loaded_snarl_collection.load_snarl_data_collection(test_file);

            check_collection(loaded_snarl_collection, true, true, true, false);

            REQUIRE(SnarlDataCollection::is_equivalent(snarl_collection, loaded_snarl_collection));

            std::string rm_cmd = "rm " + test_file;

            int rm = system(rm_cmd.c_str()); 
        }
        SECTION("Serialize it and load one line at a time") {
            // This uses the same code as the normal loader so just test that it doesn't crash and gets the right number of snarls

            std::string test_file = "./test_snarls.txt";
            std::ofstream outstream;
            outstream.open(test_file);
            snarl_collection.write_snarl_data_collection(outstream);
            outstream.close();

            TestSnarlDataCollection loaded_snarl_collection(1,10,10);
            std::ifstream instream;
            instream.open(test_file);
            size_t loaded_snarl_count = 0;
            loaded_snarl_collection.for_each_snarl_in_file(instream, [&](const snarl_info_t& snarl_info) {
                ++loaded_snarl_count;
            });
            instream.close();

            REQUIRE(loaded_snarl_count == 4);

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
    std::vector<stoat::sample_hap_t> all_samples;

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
        all_samples.emplace_back(stoat::sample_hap_t(graph, paths.back()));
    }

    //// vg isn't included so the distance index can only be built from the command line
    //graph.serialize("../tests/test_data/test_graphs/simple_nested_chains.hg");
    //int built = system("vg index -j ../tests/test_data/test_graphs/simple_nested_chains.dist ../tests/test_data/test_graphs/simple_nested_chains.hg"); 

    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/test_data/test_graphs/simple_nested_chains.dist");

   // bdsg::HashGraph graph;
   // graph.deserialize("../tests/test_data/test_graphs/simple_nested_chains.hg");

    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&graph);

    std::unordered_map<std::string, size_t> sample_to_index;
    sample_to_index.emplace("path0", 0);
    sample_to_index.emplace("path1", 1);
    sample_to_index.emplace("path2", 2);
    sample_to_index.emplace("path3", 3);

    

    SECTION("Make and fill in snarl collection with no alleles, then fill in alleles later by chromosome") {

        // Don't get the alleles or anything else
        TestSnarlDataCollection snarl_collection(1,10,10);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, all_samples, 
            true, // alleles before walks but it doesn't matter here
            true, // get walks  
            [&](const net_handle_t& snarl, const snarl_info_t& snarl_data, std::vector<stoat::PathTraversal>& walks) {
                return SnarlDataCollection::get_all_walks_through_snarl(*path_graph, distance_index, snarl, snarl_data, walks, 1);
            },
            false, // don't get alleles 
            [&](const net_handle_t& snarl, const snarl_info_t& snarl_data, const std::vector<sample_hap_t>& samples) { 
                return std::vector<size_t>();
            },
            true, // get sequences
            std::unordered_set<std::string>(), false);




        std::vector<stoat::sample_hap_t> sample_haps;
        for (const auto& path : paths) {
            sample_haps.emplace_back(stoat::sample_hap_t(*path_graph, path));
        }

        // First fill it in with a non-existent path, which should do nothing
        std::unordered_map<stoat::sample_hap_t, size_t> sample_haplotype_to_index;
        snarl_collection.add_alleles_by_sample([&](const snarl_info_t& snarl_info, const std::vector<stoat::sample_hap_t>& all_sample_haplotypes) {

            //Add alleles to follow the walks
            std::vector<size_t> alleles;
            if (snarl_info.start_node == node_traversal_t(1, false) || snarl_info.start_node == stoat::node_traversal_t(4, true) || 
                snarl_info.start_node == stoat::node_traversal_t(4, false) || snarl_info.start_node == stoat::node_traversal_t(8, true) ||
                snarl_info.start_node == stoat::node_traversal_t(5, false) || snarl_info.start_node == stoat::node_traversal_t(7, true) || 
                snarl_info.start_node == stoat::node_traversal_t(8, false) || snarl_info.start_node == stoat::node_traversal_t(10, true)) {
                // First connected component
                alleles = std::vector<size_t>( {1,1,1,1, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0} );
            } else if (snarl_info.start_node == node_traversal_t(11, false) || snarl_info.start_node == stoat::node_traversal_t(14, true) || 
                snarl_info.start_node == stoat::node_traversal_t(14, false) || snarl_info.start_node == stoat::node_traversal_t(18, true) ||
                snarl_info.start_node == stoat::node_traversal_t(15, false) || snarl_info.start_node == stoat::node_traversal_t(17, true) || 
                snarl_info.start_node == stoat::node_traversal_t(18, false) || snarl_info.start_node == stoat::node_traversal_t(20, true)) {
                // Second connected component
                alleles = std::vector<size_t>({0,0,0,0, 1,1,1,1, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0} );
            } else if (snarl_info.start_node == node_traversal_t(21, false) || snarl_info.start_node == stoat::node_traversal_t(24, true) || 
                snarl_info.start_node == stoat::node_traversal_t(24, false) || snarl_info.start_node == stoat::node_traversal_t(28, true) ||
                snarl_info.start_node == stoat::node_traversal_t(25, false) || snarl_info.start_node == stoat::node_traversal_t(27, true) || 
                snarl_info.start_node == stoat::node_traversal_t(28, false) || snarl_info.start_node == stoat::node_traversal_t(30, true)) {
                // Third connected component
                alleles = std::vector<size_t>({0,0,0,0, 0,0,0,0, 1,1,1,1, 0,0,0,0, 0,0,0,0, 0,0,0,0} );
            } else {
                // I'm not going to bother testing anything else
                alleles = std::vector<size_t>({0,0,0,0, 0,0,0,0, 0,0,0,0, 1,1,1,1, 0,0,0,0, 0,0,0,0} );
            }
            return alleles;
        }, "empty_path");
        REQUIRE(sample_haplotype_to_index.size() == 0);

        // Make sure the alleles haven't been filled in
        snarl_collection.for_each_snarl([&](const snarl_info_t& snarl_info) {
            REQUIRE(snarl_info.genotypes.get_allele_count() == 0);
        });

        // Now fill it in with just the reference path for just the first chromosome
        snarl_collection.add_alleles_by_sample([&](const snarl_info_t& snarl_info, const std::vector<stoat::sample_hap_t>& all_sample_haplotypes) { 
            //Add alleles to follow the walks
            std::vector<size_t> alleles;
            if (snarl_info.start_node == node_traversal_t(1, false) || snarl_info.start_node == stoat::node_traversal_t(4, true) || 
                snarl_info.start_node == stoat::node_traversal_t(4, false) || snarl_info.start_node == stoat::node_traversal_t(8, true) ||
                snarl_info.start_node == stoat::node_traversal_t(5, false) || snarl_info.start_node == stoat::node_traversal_t(7, true) || 
                snarl_info.start_node == stoat::node_traversal_t(8, false) || snarl_info.start_node == stoat::node_traversal_t(10, true)) {
                // First connected component
                alleles = std::vector<size_t>({1,1,1,1, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0} );
            } else if (snarl_info.start_node == node_traversal_t(11, false) || snarl_info.start_node == stoat::node_traversal_t(14, true) || 
                snarl_info.start_node == stoat::node_traversal_t(14, false) || snarl_info.start_node == stoat::node_traversal_t(18, true) ||
                snarl_info.start_node == stoat::node_traversal_t(15, false) || snarl_info.start_node == stoat::node_traversal_t(17, true) || 
                snarl_info.start_node == stoat::node_traversal_t(18, false) || snarl_info.start_node == stoat::node_traversal_t(20, true)) {
                // Second connected component
                alleles = std::vector<size_t>({0,0,0,0, 1,1,1,1, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0} );
            } else if (snarl_info.start_node == node_traversal_t(21, false) || snarl_info.start_node == stoat::node_traversal_t(24, true) || 
                snarl_info.start_node == stoat::node_traversal_t(24, false) || snarl_info.start_node == stoat::node_traversal_t(28, true) ||
                snarl_info.start_node == stoat::node_traversal_t(25, false) || snarl_info.start_node == stoat::node_traversal_t(27, true) || 
                snarl_info.start_node == stoat::node_traversal_t(28, false) || snarl_info.start_node == stoat::node_traversal_t(30, true)) {
                // Third connected component
                alleles = std::vector<size_t>({0,0,0,0, 0,0,0,0, 1,1,1,1, 0,0,0,0, 0,0,0,0, 0,0,0,0} );
            } else {
                // I'm not going to bother testing anything else
                alleles= std::vector<size_t>({0,0,0,0, 0,0,0,0, 0,0,0,0, 1,1,1,1, 0,0,0,0, 0,0,0,0} );
            }
            return alleles;
        }, "path0#0#path0");

        // Make sure that just the alleles on the first chromosome been filled in
        snarl_collection.for_each_snarl([&](const snarl_info_t& snarl_info) {
            if (snarl_info.start_node == node_traversal_t(1, false) || snarl_info.start_node == stoat::node_traversal_t(4, true) || 
                snarl_info.start_node == stoat::node_traversal_t(4, false) || snarl_info.start_node == stoat::node_traversal_t(8, true) ||
                snarl_info.start_node == stoat::node_traversal_t(5, false) || snarl_info.start_node == stoat::node_traversal_t(7, true) || 
                snarl_info.start_node == stoat::node_traversal_t(8, false) || snarl_info.start_node == stoat::node_traversal_t(10, true)) {
                REQUIRE(snarl_info.genotypes.get_allele_count() >= 0);
            } else {
                REQUIRE(snarl_info.genotypes.get_allele_count() == 0);
            }
        });

        // Now fill it in with just the reference path for just the second chromosome
        snarl_collection.add_alleles_by_sample([&](const snarl_info_t& snarl_info, const std::vector<stoat::sample_hap_t>& all_sample_haplotypes) { 
            //Add alleles to follow the walks
            std::vector<size_t> alleles;
            if (snarl_info.start_node == node_traversal_t(1, false) || snarl_info.start_node == stoat::node_traversal_t(4, true) || 
                snarl_info.start_node == stoat::node_traversal_t(4, false) || snarl_info.start_node == stoat::node_traversal_t(8, true) ||
                snarl_info.start_node == stoat::node_traversal_t(5, false) || snarl_info.start_node == stoat::node_traversal_t(7, true) || 
                snarl_info.start_node == stoat::node_traversal_t(8, false) || snarl_info.start_node == stoat::node_traversal_t(10, true)) {
                // First connected component
                alleles = std::vector<size_t>({1,1,1,1, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0} );
            } else if (snarl_info.start_node == node_traversal_t(11, false) || snarl_info.start_node == stoat::node_traversal_t(14, true) || 
                snarl_info.start_node == stoat::node_traversal_t(14, false) || snarl_info.start_node == stoat::node_traversal_t(18, true) ||
                snarl_info.start_node == stoat::node_traversal_t(15, false) || snarl_info.start_node == stoat::node_traversal_t(17, true) || 
                snarl_info.start_node == stoat::node_traversal_t(18, false) || snarl_info.start_node == stoat::node_traversal_t(20, true)) {
                // Second connected component
                alleles = std::vector<size_t>({0,0,0,0, 1,1,1,1, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0} );
            } else if (snarl_info.start_node == node_traversal_t(21, false) || snarl_info.start_node == stoat::node_traversal_t(24, true) || 
                snarl_info.start_node == stoat::node_traversal_t(24, false) || snarl_info.start_node == stoat::node_traversal_t(28, true) ||
                snarl_info.start_node == stoat::node_traversal_t(25, false) || snarl_info.start_node == stoat::node_traversal_t(27, true) || 
                snarl_info.start_node == stoat::node_traversal_t(28, false) || snarl_info.start_node == stoat::node_traversal_t(30, true)) {
                // Third connected component
                alleles = std::vector<size_t>({0,0,0,0, 0,0,0,0, 1,1,1,1, 0,0,0,0, 0,0,0,0, 0,0,0,0} );
            } else {
                // I'm not going to bother testing anything else
                alleles = std::vector<size_t>({0,0,0,0, 0,0,0,0, 0,0,0,0, 1,1,1,1, 0,0,0,0, 0,0,0,0} );
            }
            return alleles;
        }, "path0#1#path0");

        // Make sure that just the alleles on the first chromosome been filled in
        snarl_collection.for_each_snarl([&](const snarl_info_t& snarl_info) {
            if (snarl_info.start_node == node_traversal_t(1, false) || snarl_info.start_node == stoat::node_traversal_t(4, true) || 
                snarl_info.start_node == stoat::node_traversal_t(4, false) || snarl_info.start_node == stoat::node_traversal_t(8, true) ||
                snarl_info.start_node == stoat::node_traversal_t(5, false) || snarl_info.start_node == stoat::node_traversal_t(7, true) || 
                snarl_info.start_node == stoat::node_traversal_t(8, false) || snarl_info.start_node == stoat::node_traversal_t(10, true) ||
                snarl_info.start_node == node_traversal_t(11, false) || snarl_info.start_node == stoat::node_traversal_t(14, true) || 
                snarl_info.start_node == stoat::node_traversal_t(14, false) || snarl_info.start_node == stoat::node_traversal_t(18, true) ||
                snarl_info.start_node == stoat::node_traversal_t(15, false) || snarl_info.start_node == stoat::node_traversal_t(17, true) || 
                snarl_info.start_node == stoat::node_traversal_t(18, false) || snarl_info.start_node == stoat::node_traversal_t(20, true)) {
                REQUIRE(snarl_info.genotypes.get_allele_count() >= 0);
            } else {
                REQUIRE(snarl_info.genotypes.get_allele_count() == 0);
            }
        });

    }
}

TEST_CASE( "snarl collection looping snarl", "[snarl_collection]" ) {

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


    // Paths 0 and 2 take the insertion, but paths 1 and 2 take the duplication, and the deletion
    std::vector<std::vector<std::size_t>> path_seqs = { {0, 1, 2, 3, 4, 5}, {0, 1, 3, 4, 1, 3, 4, 5}, {0, 1, 2, 3, 4, 1, 3, 4, 5}};
    std::vector<handlegraph::path_handle_t> paths;

    std::vector<stoat::sample_hap_t> all_samples;
    for (int path_i = 0 ; path_i < path_seqs.size() ; path_i++) {
        if (path_i == 0) {
            // First path is the reference
            paths.emplace_back(graph.create_path_handle("path"+std::to_string(path_i)+"#0#path0"));
        } else {
            paths.emplace_back(graph.create_path_handle("path"+std::to_string(path_i)+"#0#0#0"));
        }
        for (size_t node_i : path_seqs[path_i]) {
            graph.append_step(paths.back(), nodes[node_i]);
        }
        all_samples.emplace_back(stoat::sample_hap_t(graph,paths.back()));
    }

    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&graph);

    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("../tests/test_data/test_graphs/loop_with_indel.dist");



    // Nested snarl
    handlegraph::net_handle_t snarl2 = distance_index.get_parent(distance_index.get_parent(distance_index.get_node_net_handle(3)));
    // Duplication snarl
    handlegraph::net_handle_t snarl1 = distance_index.get_parent(distance_index.get_parent(snarl2));
    handlegraph::net_handle_t root_chain = distance_index.get_parent(snarl1);


    // Function to check if the collection is valid
    // When we find the walks through the snarl first and then the sample sets, the walks will not include duplications and we will not know how to deal with
    // samples that take multiple walks. So don't check this if samples_first is false
    auto check_collection = [&] (const TestSnarlDataCollection& snarl_collection, bool check_walks, bool check_alleles, bool check_sequences, bool get_all_walks, 
                                 bool samples_first) {

        // There is an inner snarl (2-4, indel) and an outer snarl (1-6, duplication)
        // Path0 takes the insertion but not the duplication
        // Path1 takes the deletion twice
        // Path2 takes the insertion then the deletion 

        // Get the alleles again so we can check them. the order might be different though

        // snarl 4-8, indel
        //[0], [1], [2]

        // Check that we got all snarls and that we got the correct snarls
        size_t snarl_count = 0;

        snarl_collection.for_each_snarl([&](const snarl_info_t& snarl_info) {
            snarl_count++;

            if ((snarl_info.start_node == node_traversal_t(1, false) && snarl_info.end_node == stoat::node_traversal_t(6, true)) ||
                (snarl_info.start_node == stoat::node_traversal_t(6, true) && snarl_info.end_node == node_traversal_t(1, false))) {
                // Outer duplication snarl

                REQUIRE(snarl_info.ref_path == "path0#0#path0");
                REQUIRE(snarl_info.start_position == 10);
                REQUIRE(snarl_info.end_position == 14);
                REQUIRE(snarl_info.depth == 1);

                if (check_alleles) {
                    // snarl 1-6, duplication
                    //[0], [1,2]

                    std::string genotype0 =  snarl_info.genotypes.get_genotype_as_string("path0");
                    std::string genotype1 =  snarl_info.genotypes.get_genotype_as_string("path1");
                    std::string genotype2 =  snarl_info.genotypes.get_genotype_as_string("path2");
                    REQUIRE(genotype0 != genotype1);
                    REQUIRE(genotype1 == genotype2);

                    // Each sample has just one count of one
                    REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path0", 0) != snarl_info.genotypes.get_count_for_sample_and_allele("path0", 1));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path0", 0) == 0 || snarl_info.genotypes.get_count_for_sample_and_allele("path0", 1) == 0 ));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path0", 0) == 1 || snarl_info.genotypes.get_count_for_sample_and_allele("path0", 1) == 1 ));

                    REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path1", 0) != snarl_info.genotypes.get_count_for_sample_and_allele("path1", 1));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path1", 0) == 0 || snarl_info.genotypes.get_count_for_sample_and_allele("path1", 1) == 0 ));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path1", 0) == 1 || snarl_info.genotypes.get_count_for_sample_and_allele("path1", 1) == 1 ));

                    REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path2", 0) != snarl_info.genotypes.get_count_for_sample_and_allele("path2", 1));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path2", 0) == 0 || snarl_info.genotypes.get_count_for_sample_and_allele("path2", 1) == 0 ));
                    REQUIRE((snarl_info.genotypes.get_count_for_sample_and_allele("path2", 0) == 1 || snarl_info.genotypes.get_count_for_sample_and_allele("path2", 1) == 1 ));

                }
                if (check_walks) {
                    // There should be one allele with one duplication and one with none, but we don't know which is which yet
                    // Since the child is a chain, it is added to the walk as the boundary nodes and a fake node 0 for the chain
                    REQUIRE(snarl_info.walks_by_allele.size() == 2); 
                    bool walk_0_is_nodup = snarl_info.walks_by_allele[0].get_path().size() == 5 &&
                                        snarl_info.walks_by_allele[0].get_path()[0] == node_traversal_t(1, false) &&
                                        snarl_info.walks_by_allele[0].get_path()[1] == node_traversal_t(2, false) &&
                                        snarl_info.walks_by_allele[0].get_path()[2].get_node_id() == 0 && //TODO: Change this to be  == node_traversal_t(0, false)
                                        snarl_info.walks_by_allele[0].get_path()[3] == node_traversal_t(5, false) &&
                                        snarl_info.walks_by_allele[0].get_path()[4] == node_traversal_t(6, false);
                    bool walk_0_is_dup = snarl_info.walks_by_allele[0].get_path().size() == 8 &&
                                        snarl_info.walks_by_allele[0].get_path()[0] == node_traversal_t(1, false) &&
                                        snarl_info.walks_by_allele[0].get_path()[1] == node_traversal_t(2, false) &&
                                        snarl_info.walks_by_allele[0].get_path()[2].get_node_id() == 0 &&
                                        snarl_info.walks_by_allele[0].get_path()[3] == node_traversal_t(5, false) &&
                                        snarl_info.walks_by_allele[0].get_path()[4] == node_traversal_t(2, false) &&
                                        snarl_info.walks_by_allele[0].get_path()[5].get_node_id() == 0 &&
                                        snarl_info.walks_by_allele[0].get_path()[6] == node_traversal_t(5, false) &&
                                        snarl_info.walks_by_allele[0].get_path()[7] == node_traversal_t(6, false);

                    bool walk_1_is_nodup = snarl_info.walks_by_allele[1].get_path().size() == 5 &&
                                        snarl_info.walks_by_allele[1].get_path()[0] == node_traversal_t(1, false) &&
                                        snarl_info.walks_by_allele[1].get_path()[1] == node_traversal_t(2, false) &&
                                        snarl_info.walks_by_allele[1].get_path()[2].get_node_id() == 0 &&
                                        snarl_info.walks_by_allele[1].get_path()[3] == node_traversal_t(5, false) &&
                                        snarl_info.walks_by_allele[1].get_path()[4] == node_traversal_t(6, false);
                    bool walk_1_is_dup =  snarl_info.walks_by_allele[1].get_path().size() == 8 &&
                                        snarl_info.walks_by_allele[1].get_path()[0] == node_traversal_t(1, false) &&
                                        snarl_info.walks_by_allele[1].get_path()[1] == node_traversal_t(2, false) &&
                                        snarl_info.walks_by_allele[1].get_path()[2].get_node_id() == 0 &&
                                        snarl_info.walks_by_allele[1].get_path()[3] == node_traversal_t(5, false) &&
                                        snarl_info.walks_by_allele[1].get_path()[4] == node_traversal_t(2, false) &&
                                        snarl_info.walks_by_allele[1].get_path()[5].get_node_id() == 0 &&
                                        snarl_info.walks_by_allele[1].get_path()[6] == node_traversal_t(5, false) &&
                                        snarl_info.walks_by_allele[1].get_path()[7] == node_traversal_t(6, false);
                    REQUIRE(((walk_0_is_dup && walk_1_is_nodup) || (walk_0_is_nodup && walk_1_is_dup)));

                    REQUIRE((stoat::vectorPathToString(snarl_info.walks_by_allele, true) == "3/4,6/8" || stoat::vectorPathToString(snarl_info.walks_by_allele, true) == "6/8,3/4"));
                }

                // Check that the walks and sample_sets match
                if (check_walks && check_alleles) {
                    if (snarl_info.walks_by_allele[0].get_path().size() == 5) {
                        // First allele is no duplication ([0])

                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path0", 0) == 1);
                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path0", 1) == 0);

                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path1", 0) == 0);
                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path1", 1) == 1);

                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path2", 0) == 0);
                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path2", 1) == 1);
                        
                        if (check_sequences) {
                            REQUIRE(snarl_info.sequences_by_allele[0] == "ANT");
                            REQUIRE(snarl_info.sequences_by_allele[1] == "ANTANT");
                        }
                    } else {
                        // First allele is duplication ([1,2])

                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path0", 0) == 0);
                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path0", 1) == 1);

                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path1", 0) == 1);
                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path1", 1) == 0);

                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path2", 0) == 1);
                        REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path2", 1) == 0);
                        if (check_sequences) {
                            REQUIRE(snarl_info.sequences_by_allele[0] == "ANTANT");
                            REQUIRE(snarl_info.sequences_by_allele[1] == "ANT");
                        }
                    }
                }


                if (check_sequences) {
                    REQUIRE(snarl_info.sequences_by_allele.size() == 2);
                    REQUIRE(((snarl_info.sequences_by_allele[0] == "ANT" && snarl_info.sequences_by_allele[1] == "ANTANT") ||
                             (snarl_info.sequences_by_allele[1] == "ANT" && snarl_info.sequences_by_allele[0] == "ANTANT")));
                }

            } else if ((snarl_info.start_node == stoat::node_traversal_t(2, false) && snarl_info.end_node == stoat::node_traversal_t(4, true)) ||
                (snarl_info.start_node == stoat::node_traversal_t(4, true) && snarl_info.end_node == stoat::node_traversal_t(2, false))) {
                // Inner indel snarl

                REQUIRE(snarl_info.ref_path == "path0#0#path0");
                REQUIRE(snarl_info.start_position == 11);
                REQUIRE(snarl_info.end_position == 12);
                REQUIRE(snarl_info.depth == 2);

                if (check_alleles) {
                    if (samples_first) {
                        // If we got the sample sets before the walks, then there are three sample sets. Otherwise we haven't decided what the right answer is so don't check
                        std::string genotype0 =  snarl_info.genotypes.get_genotype_as_string("path0");
                        std::string genotype1 =  snarl_info.genotypes.get_genotype_as_string("path1");
                        std::string genotype2 =  snarl_info.genotypes.get_genotype_as_string("path2");
                        REQUIRE(genotype0 != genotype1);
                        REQUIRE(genotype1 != genotype2);
                        REQUIRE(genotype0 != genotype2);
                    }
                }
                if (check_walks) {
                    if (samples_first) {
                        // Each path takes a different walk
                        // figure out which allele takes which walk
                        REQUIRE(snarl_info.walks_by_allele.size() == 3);

                        for (size_t allele_i = 0 ; allele_i < 3 ; allele_i++) {
                            if (snarl_info.walks_by_allele[allele_i].get_path().size() == 3) {
                                // This is walk 0 taking the insertion once
                                REQUIRE(snarl_info.walks_by_allele[allele_i].get_path()[0] == node_traversal_t(2, false));
                                REQUIRE(snarl_info.walks_by_allele[allele_i].get_path()[1] == node_traversal_t(3, false));
                                REQUIRE(snarl_info.walks_by_allele[allele_i].get_path()[2] == node_traversal_t(4, false));
                            } else if (snarl_info.walks_by_allele[allele_i].get_path().size() == 5) {
                                // This is walk1 that takes the deletion twice
                                REQUIRE(snarl_info.walks_by_allele[allele_i].get_path()[0] == node_traversal_t(2, false));
                                REQUIRE(snarl_info.walks_by_allele[allele_i].get_path()[1] == node_traversal_t(4, false));
                                REQUIRE(snarl_info.walks_by_allele[allele_i].get_path()[2] == node_traversal_t(0, true));
                                REQUIRE(snarl_info.walks_by_allele[allele_i].get_path()[3] == node_traversal_t(2, false));
                                REQUIRE(snarl_info.walks_by_allele[allele_i].get_path()[4] == node_traversal_t(4, false));
                            } else if (snarl_info.walks_by_allele[allele_i].get_path().size() == 6) {
                                // This is walk2 that takes insertion then deletion
                                REQUIRE(snarl_info.walks_by_allele[allele_i].get_path()[0] == node_traversal_t(2, false));
                                REQUIRE(snarl_info.walks_by_allele[allele_i].get_path()[1] == node_traversal_t(3, false));
                                REQUIRE(snarl_info.walks_by_allele[allele_i].get_path()[2] == node_traversal_t(4, false));
                                REQUIRE(snarl_info.walks_by_allele[allele_i].get_path()[3] == node_traversal_t(0, true));
                                REQUIRE(snarl_info.walks_by_allele[allele_i].get_path()[4] == node_traversal_t(2, false));
                                REQUIRE(snarl_info.walks_by_allele[allele_i].get_path()[5] == node_traversal_t(4, false));
                            }
                        }

                        REQUIRE((stoat::vectorPathToString(snarl_info.walks_by_allele, true) == "0,1,1" || stoat::vectorPathToString(snarl_info.walks_by_allele, true) == "1,0,1" || stoat::vectorPathToString(snarl_info.walks_by_allele, true) == "1,1,0"));
                    } else {
                        // There are only two walks, the insertion and the deletion
                        // figure out which allele takes which walk
                        REQUIRE(snarl_info.walks_by_allele.size() == 2);

                        for (size_t allele_i = 0 ; allele_i < 2 ; allele_i++) {
                            if (snarl_info.walks_by_allele[allele_i].get_path().size() == 2) {
                                // deletion
                                REQUIRE(snarl_info.walks_by_allele[allele_i].get_path()[0] == node_traversal_t(2, false));
                                REQUIRE(snarl_info.walks_by_allele[allele_i].get_path()[1] == node_traversal_t(4, false));
                            } else if (snarl_info.walks_by_allele[allele_i].get_path().size() == 3) {
                                // insertion
                                REQUIRE(snarl_info.walks_by_allele[allele_i].get_path()[0] == node_traversal_t(2, false));
                                REQUIRE(snarl_info.walks_by_allele[allele_i].get_path()[1] == node_traversal_t(3, false));
                                REQUIRE(snarl_info.walks_by_allele[allele_i].get_path()[2] == node_traversal_t(4, false));
                            }
                        }

                        REQUIRE((stoat::vectorPathToString(snarl_info.walks_by_allele, true) == "0,1" || stoat::vectorPathToString(snarl_info.walks_by_allele, true) == "1,0"));
                    }
                }

                // Check that the walks and alleles match
                // Don't do this if we got the walks first
                if (samples_first && check_walks && check_alleles) {
                    for (size_t allele_i = 0 ; allele_i < 3 ; allele_i++) {
                        if (snarl_info.walks_by_allele[allele_i].get_path().size() == 3) {
                            // This is walk 0 taking the insertion once

                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path0", allele_i) == 1);
                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path1", allele_i) == 0);
                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path2", allele_i) == 0);

                            if (check_sequences) {
                                REQUIRE(snarl_info.sequences_by_allele[allele_i] == "G");
                            }
                        } else if (snarl_info.walks_by_allele[allele_i].get_path().size() == 5) {
                            // This is walk1 that takes the deletion twice

                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path0", allele_i) == 0);
                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path1", allele_i) == 1);
                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path2", allele_i) == 0);

                            if (check_sequences) {
                                REQUIRE(snarl_info.sequences_by_allele[allele_i] == "N");
                            }
                        } else if (snarl_info.walks_by_allele[allele_i].get_path().size() == 6) {
                            // This is walk2 that takes insertion then deletion

                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path0", allele_i) == 0);
                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path1", allele_i) == 0);
                            REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path2", allele_i) == 1);
                            if (check_sequences) {
                                REQUIRE(snarl_info.sequences_by_allele[allele_i] == "GN");
                            }
                        }
                    }
                }


                if (check_sequences) {

                    if (samples_first) {
                        REQUIRE(snarl_info.sequences_by_allele.size() == 3);
                        for (size_t allele_i = 0 ; allele_i < 3 ; allele_i++) {
                            REQUIRE((snarl_info.sequences_by_allele[allele_i] == "G" || 
                                     snarl_info.sequences_by_allele[allele_i] == "N" || 
                                     snarl_info.sequences_by_allele[allele_i] == "GN"));
                        }
                    } else {
                        REQUIRE(snarl_info.sequences_by_allele.size() == 2);
                        REQUIRE(((snarl_info.sequences_by_allele[0] == "" && snarl_info.sequences_by_allele[1] == "G") ||
                                 (snarl_info.sequences_by_allele[1] == "" && snarl_info.sequences_by_allele[0] == "G")));
                    }
                }

            } else {
                REQUIRE(false);
            }
        });

        REQUIRE(snarl_count == 2);
    };

    auto get_alleles_per_snarl = [&](const snarl_info_t& snarl_info, const std::vector<stoat::sample_hap_t>& haplotypes) {

        std::vector<size_t> alleles_by_snarl;

        //Add alleles to follow the walks
        if ((snarl_info.start_node == node_traversal_t(1, false) && snarl_info.end_node == stoat::node_traversal_t(6, true)) ||
            (snarl_info.start_node == stoat::node_traversal_t(6, true) && snarl_info.end_node == node_traversal_t(1, false))) {
            // Outer duplication snarl

            if (snarl_info.walks_by_allele.size() == 0 || snarl_info.walks_by_allele[0].get_path().size() == 5) {
                // If we aren't doing it by the walks or if the first walk is the non-dulpication
                alleles_by_snarl = std::vector<size_t>({0,1,1});
            } else {
                alleles_by_snarl = std::vector<size_t>({1,0,0});
            }
        } else if ((snarl_info.start_node == stoat::node_traversal_t(2, false) && snarl_info.end_node == stoat::node_traversal_t(4, true)) ||
                (snarl_info.start_node == stoat::node_traversal_t(4, true) && snarl_info.end_node == stoat::node_traversal_t(2, false))) {
            // snarl 4-8
            // inner indel snarl
            if (snarl_info.walks_by_allele.size() == 0) {
                //if we don't care about matching the order of the sample sets with the order of walks
                alleles_by_snarl = std::vector<size_t>({0,1,2});
            } else {
                alleles_by_snarl = std::vector<size_t>({0,0,0});
                for (size_t allele_i = 0 ; allele_i < 3 ; allele_i ++ ) {
                    if (snarl_info.walks_by_allele[allele_i].get_path().size() == 3) { 
                        alleles_by_snarl[0] = allele_i;
                    } else if (snarl_info.walks_by_allele[allele_i].get_path().size() == 4) { 
                        alleles_by_snarl[1] = allele_i;
                    } else if (snarl_info.walks_by_allele[allele_i].get_path().size() == 5) {
                        alleles_by_snarl[2] = allele_i;
                    }
                }
            }
        }
        return alleles_by_snarl;
    };

    std::unordered_map<std::string, size_t> sample_to_index;
    sample_to_index.emplace("path0", 0);
    sample_to_index.emplace("path1", 1);
    sample_to_index.emplace("path2", 2);

    SECTION("Make and fill in snarl collection with no walks, alleles, or sequences") {

        TestSnarlDataCollection snarl_collection(1,10,10);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, all_samples, 
            true, // alleles before walks but it doesn't matter here
            false, // don't get walks 
            [&](const net_handle_t& snarl, const snarl_info_t& snarl_data,std::vector<stoat::PathTraversal>& walks) {
                return;
            }, false, // don't get alleles 
            [&](const net_handle_t& snarl, const snarl_info_t& snarl_data,const std::vector<sample_hap_t>& samples) { 
                return std::vector<size_t>();
            },
            false, // don't get sequences
            std::unordered_set<std::string>(), false);


        check_collection(snarl_collection, false, false, false, true, true);
    }
    SECTION("Make and fill in snarl collection with walks, alleles, and sequences") {



        TestSnarlDataCollection snarl_collection(1,10,10);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, all_samples, 
            true, // alleles before walks but it doesn't matter here
            true, // get walks 
            [&](const net_handle_t& snarl, const snarl_info_t& snarl_data, std::vector<PathTraversal>& walks) {
                SnarlDataCollection::get_walks_from_alleles(*path_graph, distance_index, snarl, snarl_data, walks);
                return;
            },
            true, // get alleles 
            [&](const net_handle_t& snarl, const snarl_info_t& snarl_data, const std::vector<sample_hap_t>& samples) { 
                return  get_alleles_per_snarl(snarl_data, samples);
            },
            true, // get sequences
            std::unordered_set<std::string>(), false);

        check_collection(snarl_collection, true, true, true, false, true);

        SECTION("Serialize it") {

            std::string test_file = "./test_snarls.txt";
            std::ofstream outstream;
            outstream.open(test_file);
            snarl_collection.write_snarl_data_collection(outstream);
            outstream.close();

            TestSnarlDataCollection loaded_snarl_collection(1,10,10);
            loaded_snarl_collection.load_snarl_data_collection(test_file);

            check_collection(loaded_snarl_collection, true, true, true, false, true);

            REQUIRE(SnarlDataCollection::is_equivalent(snarl_collection, loaded_snarl_collection));

            std::string rm_cmd = "rm " + test_file;

            int rm = system(rm_cmd.c_str()); 
        }
    }
    SECTION("Make and fill in snarl collection with no alleles, then fill in alleles later by chromosome") {

        // Don't get the alleles or anything else
        TestSnarlDataCollection snarl_collection(1,10,10);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, all_samples, 
            false, // walks before alleles but it doesn't matter here
            true, // get walks  
            [&](const net_handle_t& snarl, const snarl_info_t& snarl_data, std::vector<stoat::PathTraversal>& walks) {
                return SnarlDataCollection::get_all_walks_through_snarl(*path_graph, distance_index, snarl, snarl_data, walks, 0);
            },
            false, // don't get alleles 
            [&](const net_handle_t& snarl, const snarl_info_t& snarl_data,const std::vector<sample_hap_t>& samples) { 
                return std::vector<size_t>();
            },
            true, // get sequences
            std::unordered_set<std::string>(), false);

        // First fill it in with a non-existent path, which should do nothing
        std::unordered_map<stoat::sample_hap_t, size_t> sample_haplotype_to_index;
        snarl_collection.add_alleles_by_sample(get_alleles_per_snarl, "empty_path");

        check_collection(snarl_collection, true, false, true, true, false);

        REQUIRE(sample_haplotype_to_index.size() == 0);

        snarl_collection.for_each_snarl([&](const snarl_info_t& snarl_info) {
            // Make sure the alleles haven't been filled in
            REQUIRE(snarl_info.genotypes.get_allele_count() == 0);
        });

        // Now fill it in with the paths
        snarl_collection.add_alleles_by_sample(get_alleles_per_snarl, "path0#0#path0");

        check_collection(snarl_collection, true, true, true, true, false);

    }
}
TEST_CASE( "Snarl collection nested bubbles with path fragments",
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

   //This uses the simple_nested_chain distance index but rebuilds the graph with different paths 

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

    // Function to check if the collection is valid
    auto check_collection = [&] (const TestSnarlDataCollection& snarl_collection) {

        // The point of this test is to check what happens when there are multiple fragments
        // So just check snarl 3-7

        // Get the alleles again so we can check them. the order might be different though

        // Make alleles for each snarl
        std::vector<stoat::sample_hap_t> sample_haps;
        for (const auto& path : paths) {
            sample_haps.emplace_back(stoat::sample_hap_t(*path_graph, path));
        }

        // Check that we got all snarls and that we got the correct snarls
        size_t snarl_count = 0;

        snarl_collection.for_each_snarl([&](const snarl_info_t& snarl_info) {
            snarl_count++;

            if ((snarl_info.start_node == node_traversal_t(4, false) && snarl_info.end_node == stoat::node_traversal_t(8, true)) ||
                (snarl_info.start_node == stoat::node_traversal_t(8, true) && snarl_info.end_node == node_traversal_t(4, false))) {
                // Outer duplication snarl

                REQUIRE(snarl_info.ref_path == "path0#0#path0");
                REQUIRE(snarl_info.start_position == 3);
                REQUIRE(snarl_info.end_position == 5);
                REQUIRE(snarl_info.depth == 1);


                // Since I'm giving it the order of the sets, it should be ordered

                // The genotypes are correct
                // Allele 0 is insertion, 1 is deletion
                // ins, del/del, del/ins
                REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path0", 0) == 1);
                REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path0", 1) == 0);

                REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path1", 0) == 0);
                REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path1", 1) == 2);

                REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path2", 0) == 1);
                REQUIRE(snarl_info.genotypes.get_count_for_sample_and_allele("path2", 1) == 1);

                // first path takes the insertion 
                REQUIRE(snarl_info.walks_by_allele[0].get_path().size() == 5);
                REQUIRE((snarl_info.walks_by_allele[0].get_path()[0] == node_traversal_t(4, false) &&
                         snarl_info.walks_by_allele[0].get_path()[1] == node_traversal_t(5, false) &&
                         snarl_info.walks_by_allele[0].get_path()[2] == node_traversal_t(0, false) &&
                         snarl_info.walks_by_allele[0].get_path()[3] == node_traversal_t(7, false) &&
                         snarl_info.walks_by_allele[0].get_path()[4] == node_traversal_t(8, false)));

                REQUIRE(snarl_info.sequences_by_allele[0] == "TNA");

                REQUIRE(snarl_info.walks_by_allele[1].get_path().size() == 2);
                REQUIRE((snarl_info.walks_by_allele[1].get_path()[0] == node_traversal_t(4, false) &&
                         snarl_info.walks_by_allele[1].get_path()[1] == node_traversal_t(8, false)));

                REQUIRE(snarl_info.sequences_by_allele[1] == "");

            }
        });

    };

    auto get_alleles_per_snarl = [&](const snarl_info_t& snarl_info, const std::vector<stoat::sample_hap_t>& haplotypes) {

        std::vector<size_t> alleles_by_snarl;

        if (snarl_info.start_node == node_traversal_t(1, false) || snarl_info.start_node == stoat::node_traversal_t(4, true)) {
            // snarl 1-4
            alleles_by_snarl = std::vector<size_t>({0,1,1,0,1});
        } else if (snarl_info.start_node == stoat::node_traversal_t(4, false) || snarl_info.start_node == stoat::node_traversal_t(8, true)) {
            // snarl 4-8
            alleles_by_snarl = std::vector<size_t>({0,1,1,0,1});
        } else if (snarl_info.start_node == stoat::node_traversal_t(5, false) || snarl_info.start_node == stoat::node_traversal_t(7, true)) {
    
            // snarl 5-7
            alleles_by_snarl = std::vector<size_t>({0,std::numeric_limits<size_t>::max(),std::numeric_limits<size_t>::max(),0,std::numeric_limits<size_t>::max()});
        } else {
            // For the last one there are no alleles because there are no paths
            alleles_by_snarl = std::vector<size_t>({std::numeric_limits<size_t>::max(),std::numeric_limits<size_t>::max(),std::numeric_limits<size_t>::max(),std::numeric_limits<size_t>::max(),std::numeric_limits<size_t>::max()});
        }


        return alleles_by_snarl;
    };

    std::unordered_map<std::string, size_t> sample_to_index;
    sample_to_index.emplace("path0", 0);
    sample_to_index.emplace("path1", 1);
    sample_to_index.emplace("path2", 2);

    SECTION("Make and fill in snarl collection with walks, alleles, and sequences") {

        TestSnarlDataCollection snarl_collection(1,10,10);
        snarl_collection.fill_in_snarl_info(*path_graph, distance_index, all_samples, 
            true, // alleles before walks but it doesn't matter here
            true, // get walks 
            [&](const net_handle_t& snarl, const snarl_info_t& snarl_data, std::vector<PathTraversal>& walks) {
                return SnarlDataCollection::get_walks_from_alleles(*path_graph, distance_index, snarl, snarl_data, walks);
            },
            true, // get alleles 
            [&](const net_handle_t& snarl, const snarl_info_t& snarl_data, const std::vector<sample_hap_t>& samples) { 
                return get_alleles_per_snarl(snarl_data, samples);
            },
            true, // get sequences
            std::unordered_set<std::string>(), false);

        check_collection(snarl_collection);

    }
}
