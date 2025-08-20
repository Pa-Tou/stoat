#include "graph_path_association_finder.hpp"
#include "utils.hpp"
#include "writer.hpp"
#include "binary_table.hpp"

//#define DEBUG_ASSOCIATION_FINDER

namespace stoat_graph {

AssociationFinder::AssociationFinder(const handlegraph::PathPositionHandleGraph& graph, 
                                     const bdsg::SnarlDistanceIndex& distance_index,
                                     std::shared_ptr<Partitioner> partitioner,
                                     const std::set<std::string>& samples_of_interest, 
                                     const std::string& reference_sample,
                                     const std::string& test_method,
                                     const std::string& output_format,
                                     size_t allele_size_limit,
                                     std::ostream& out_associated,
                                     std::ostream& out_unassociated) :
    graph(graph), 
    distance_index(distance_index), 
    partitioner(std::move(partitioner)),
    samples_of_interest(samples_of_interest), 
    reference_sample(reference_sample),
    test_method(test_method),
    output_format(output_format),
    allele_size_limit(allele_size_limit),
    out_associated(out_associated),
    out_unassociated(out_unassociated),
    check_distances(distance_index.has_distances())
    {}

void AssociationFinder::test_snarls() const {

    //TODO: Make this general
    // If the file output has a header, write it
    if (output_format == "tsv") {
        stoat::write_binary_header(out_associated);
    }

    std::vector<handlegraph::net_handle_t> chains;
    chains.reserve(graph.get_node_count()/100);
    handlegraph::net_handle_t root = distance_index.get_root();
    distance_index.for_each_child(root, [&] (handlegraph::net_handle_t chain) {
        chains.emplace_back(chain);
        return true;
    });

    FisherKhi2 fisher_chi2_tester;
    while (!chains.empty()) {
        handlegraph::net_handle_t chain = chains.back();
        chains.pop_back();

        distance_index.for_each_child(chain, [&] (handlegraph::net_handle_t snarl) {

            //TODO: For now it's fine to check is_eligible here because it's only checking size and we don't want to look at small chains anyway
            if (distance_index.is_snarl(snarl) && snarl_is_eligible(snarl) ) {
                #ifdef DEBUG_ASSOCIATION_FINDER
                    cerr << "Test snarl " << distance_index.net_handle_as_string(snarl) << endl;
                #endif

                // Should we write this?
                bool write_output = false;

                // the strings we are going to output
                string group_paths = "NA";
                string fastfisher_p_value = "NA";
                string chi2_p_value = "NA";
                // Get the path lengths, except since we don't know the lengths of the alleles, it's just the min and max length of the snarl
                std::stringstream ss;
                if (check_distances) {
                    ss << distance_index.minimum_length(snarl) << "/" << distance_index.maximum_length(snarl);
                } else {
                    ss << "NA";
                }
                string path_lengths = ss.str();

                // Each set represents a partition of samples that takes the same path through the snarl's netgraph
                std::vector<std::set<std::string>> sample_partitions = partitioner->partition_samples_in_snarl(graph, distance_index, snarl);

                // Do we test nested snarls? Don't test snarls that are already flagged as significant
                bool test_nested_snarls = true;

                #ifdef DEBUG_ASSOCIATION_FINDER
                    cerr << "\tTRUTH" << endl;
                    for (const std::string& sample : samples_of_interest) {
                        cerr << "\t\t" << sample << endl;
                    }

                    for (const std::set<std::string>& partition : sample_partitions) {
                        cerr << "\tPARTITION" << endl;
                        for (const std::string& sample : partition) {
                            cerr << "\t\t" << sample << endl;
                        }
                    }
                #endif

                if (sample_partitions.size() > 1) {

                    // If we are writing a fasta, then pick one sample from each partition to write
                    std::unordered_map<std::string, bool> samples_to_write;

                    if (test_method == "exact") {



                        for (const std::set<std::string>& partition : sample_partitions) {
                            if (partition == samples_of_interest) {

                                // For the exact test, since we already know the result of the test, write only those snarls that pass the test
                                write_output = true;
                                // Don't look for nested snarls
                                test_nested_snarls = false;
                                if (output_format == "fasta") {
                                    samples_to_write[*partition.begin()] = true;
                                } else {
                                    break;
                                }
                            } else if (output_format == "fasta") {
                                samples_to_write[*partition.begin()] = false;
                            }
                        }

                    } else {

                        // If we are using a real statistical test, then always write the output because the BH correction will need all the p-values
                        // TODO: This could do what pangwas was doing to keep track of only good p-values instead of writing everything
                        write_output = true;

                        // Fill in the genotypes. Each item in these vectors is an allele (path/sample partition)
                        std::vector<size_t> genotype_associated(sample_partitions.size(), 0);
                        std::vector<size_t> genotype_unassociated(sample_partitions.size(), 0);
                        for (size_t i = 0 ; i < sample_partitions.size() ; i++) {
                            const std::set<std::string> sample_set = sample_partitions[i];
                            for (const std::string sample : sample_set) {
                                if (samples_of_interest.count(sample) != 0) {
                                    genotype_associated[i]++;
                                } else {
                                    genotype_unassociated[i]++;
                                }
                            }
                        }

                        //Get a bunch of strings that get used for the output
                        // TODO: This function should probably be part of the output function

                        //Get a bunch of strings that get used for the output
                        group_paths = stoat_vcf::format_group_paths(genotype_associated, genotype_unassociated);
 
                        // Run the statistical test
                        std::tie(chi2_p_value, fastfisher_p_value) = fisher_chi2_tester.fisher_khi2(genotype_associated, genotype_unassociated);

                        if (output_format == "fasta") {
                            // Figure out which samples we want to write
                            // Since we don't know which partition is actually associated, just write everything to one file
                            for (const std::set<std::string>& partition : sample_partitions) {
                                samples_to_write[*partition.begin()] = true;
                            }
                        }

                    }
                
                    if (write_output) {
                        if (output_format == "tsv") {

                            string chr = "NA"; 
                            // TODO: Maybe I sould keep the snarls as snarl_data_t's? 
                            // TODO: get the type properly
                            stoat::Snarl_data_t snarl_data_s(snarl, graph, distance_index);

                            // Get the offsets of the start and end nodes along the reference
                            std::vector<stoat::path_range_t> ranges = stoat::get_coordinates_of_snarl(graph, distance_index, snarl, true, reference_sample, false);
                            if (ranges.size() != 0) {
                                std::tie(chr, snarl_data_s.start_positions, snarl_data_s.end_positions) = get_name_and_offsets_of_snarl_path_range(graph, distance_index, ranges.front());
                            }

                            # pragma omp critical (out_associated) 
                            {
                                // Leave adjusted p-value blank, to be filled in later
                                stoat::write_binary(out_associated, chr, snarl_data_s, path_lengths, fastfisher_p_value, chi2_p_value, "",  group_paths);
                            }
                        } else if (output_format == "fasta") {

                            # pragma omp critical (out_associated) 
                            {
                                stoat::write_fasta(out_associated, out_unassociated, graph, distance_index, snarl, samples_to_write, reference_sample);
                            }
                        }
                    }
                }

                if (test_nested_snarls) { 
                    // Add the child chains to the stack
                    distance_index.for_each_child(snarl, [&] (handlegraph::net_handle_t child) {

                        chains.emplace_back(child);
                        return true;
                    });
                }
            }
            return true;
        });
    }
}

bool AssociationFinder::snarl_is_eligible(const handlegraph::net_handle_t& snarl) const {
    if (!check_distances) {
        // If the distance index doesn't let us check distances, just return true
        return true;
    } else {
        return distance_index.maximum_length(snarl) >= allele_size_limit;
    }
}

} //end pangwas namespace
