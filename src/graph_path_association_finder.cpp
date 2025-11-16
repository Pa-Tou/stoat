#include "graph_path_association_finder.hpp"
#include "utils.hpp"
#include "writer.hpp"
#include "binary_table.hpp"

//#define DEBUG_GRAPH_SNARLS

using namespace stoat;
namespace stoat_graph {

AssociationFinder::AssociationFinder(const handlegraph::PathPositionHandleGraph& graph, 
                                     std::shared_ptr<SnarlTraverserAndPartitioner> partitioner,
                                     const std::pair<std::set<std::string>, std::set<std::string>>& sample_sets, 
                                     const std::string& reference_sample,
                                     const std::string& test_method,
                                     double maf_threshold,
                                     size_t min_individuals,
                                     const std::string& output_format,
                                     std::ostream& out_associated) :
    graph(graph), 
    partitioner(std::move(partitioner)),
    sample_sets(sample_sets),
    reference_sample(reference_sample),
    test_method(test_method),
    maf_threshold(maf_threshold),
    min_individuals(min_individuals),
    output_format(output_format),
    out_associated(out_associated)
    {}

void AssociationFinder::test_snarls() const {

    //TODO: Make this general
    // If the file output has a header, write it
    if (output_format == "tsv") {
        stoat::write_binary_header(out_associated);
    }

    stoat::FisherKhi2 fisher_chi2_tester;

    // Go through all snarls, which have already been partitioned by the partitioner
    partitioner->for_each_snarl_partition(graph, 
    [&] (const stoat::snarl_partition_t& snarl_info) {


        // Each set represents a partition of samples that takes the same path through the snarl's netgraph
        const std::vector<std::set<sample_hap_t>>& sample_partitions = snarl_info.partitions;

        // Do we test nested snarls? Don't test snarls that are already flagged as significant
        // This can only be used for the exact test, so anything nested will not pass anyways
        bool test_nested_snarls = true;

        #ifdef DEBUG_GRAPH_SNARLS
        cerr << "\tTRUTH 1" << endl;
        for (const auto& sample : sample_sets.first) {
            cerr << "\t\t" << sample << endl;
        }
        cerr << "\tTRUTH 2" << endl;
        for (const auto& sample : sample_sets.second) {
            cerr << "\t\t" << sample << endl;
        }

        for (const std::set<sample_hap_t>& partition : sample_partitions) {
            cerr << "\tPARTITION" << endl;
            for (const sample_hap_t& sample : partition) {
                cerr << "\t\t" << sample.sample << endl;
            }
        }
        #endif

        if (sample_partitions.size() > 1) {


            // Should we write this snarl immediately?
            bool write_output = false;

            // the strings we are going to output
            string group_paths = "NA";
            string fastfisher_p_value = "NA";
            string chi2_p_value = "NA";
            // Get the path lengths, except since we don't know the lengths of the alleles, it's just the min and max length of the snarl
            std::stringstream ss;
            if (snarl_info.min_length == std::numeric_limits<size_t>::max() && snarl_info.max_length == std::numeric_limits<size_t>::max()) {
                ss << "NA,NA";
            } else {
                ss << snarl_info.min_length << "," << snarl_info.max_length;
            }

            string path_lengths = ss.str();

            // If we are writing a fasta, then pick one sample from each partition to write
            std::unordered_map<std::string, bool> samples_to_write;

            if (test_method == "exact") {

                for (const std::set<sample_hap_t>& sample_hap_partition : sample_partitions) {

                    // Make a set of just the sample names
                    // TODO: This isn't super efficient
                    std::set<std::string> partition;
                    for (const sample_hap_t& sample : sample_hap_partition) {
                        partition.emplace(sample.sample);
                    }

                    // If one partition exactly matches one group we want, then all other partitions combined (including things not in the snarl) will match
                    // the other.
                    // TODO: This could be better but I don't think it's worth working on it yet
                    if (partition == sample_sets.first || partition == sample_sets.second) {

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
                // Else we are doing chi2

                // Fill in the genotypes. Each item in these vectors is an allele (path/sample partition)
                std::vector<size_t> genotype_associated(sample_partitions.size(), 0);
                std::vector<size_t> genotype_unassociated(sample_partitions.size(), 0);

                // How many individuals/samples are included?
                std::unordered_set<std::string> seen_samples;
                // Put the sample partitions into genotype vectors (one entry in the vector is one allele/path/genotype)
                for (size_t i = 0 ; i < sample_partitions.size() ; i++) {
                    const std::set<sample_hap_t>& sample_set = sample_partitions[i];
                    for (const sample_hap_t sample : sample_set) {
                        if (sample_sets.first.count(sample.sample) == 1) {
                            genotype_associated[i]++;
                        } else if (sample_sets.second.count(sample.sample) == 1) {
                            genotype_unassociated[i]++;
                        }
                        seen_samples.insert(sample.sample);
                    }
                }
                if (!stoat::filter_binary_table(genotype_associated, genotype_unassociated, seen_samples.size(), min_individuals, maf_threshold)) {
                    // If we are using a real statistical test, then always write the output because the BH correction will need all the p-values
                    // TODO: This could do what pangwas was doing to keep track of only good p-values instead of writing everything
                    write_output = true;

                    //Get a bunch of strings that get used for the output
                    // TODO: This function should probably be part of the output function

                    //Get a bunch of strings that get used for the output
                    group_paths = stoat_vcf::format_group_paths(genotype_associated, genotype_unassociated);
 
                    // Run the statistical test
                    std::tie(chi2_p_value, fastfisher_p_value) = fisher_chi2_tester.fisher_khi2(genotype_associated, genotype_unassociated);

                    if (output_format == "fasta") {
                        // Figure out which samples we want to write
                        // Since we don't know which partition is actually associated, just write everything to one file
                        for (const std::set<sample_hap_t>& partition : sample_partitions) {
                            samples_to_write[partition.begin()->sample] = true;
                        }
                    }
                }

            }
        
            if (write_output) {
                if (output_format == "tsv") {


                    # pragma omp critical (out_associated) 
                    {
                        // Leave adjusted p-value blank, to be filled in later
                        stoat::write_binary(out_associated, snarl_info.ref_path, snarl_info, path_lengths, fastfisher_p_value, chi2_p_value, group_paths);
                    }
                } else if (output_format == "fasta") {

                    # pragma omp critical (out_associated) 
                    {
                        stoat::write_fasta(out_associated, graph, snarl_info, samples_to_write, reference_sample);
                    }
                }
            }
        }
    });
}


} //end pangwas namespace
