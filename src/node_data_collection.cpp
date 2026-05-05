#include <fstream>
#include <filesystem>

#include "node_data_collection.hpp"
#include "node_data_collection.hpp"
#include "matrix.hpp"
#include "utils.hpp"

//#define DEBUG_NODE_DATA_COLLECTION

namespace stoat {

// Constructor
NodeDataCollection::NodeDataCollection() {}

NodeDataCollection::NodeDataCollection(std::unordered_map<std::string, size_t> sample_to_index) :
                    sample_to_index(sample_to_index) {}

// This goes through all the nodes and fills in the data
void NodeDataCollection::fill_in_node_info(const handlegraph::PathPositionHandleGraph& graph, 
                                             const bdsg::SnarlDistanceIndex& distance_index,
                                             const std::vector<stoat::sample_hap_t>& sample_haplotypes,
                                             bool find_alleles_first,
                                             bool walks_requested,
                                             const std::function<void(const net_handle_t& node, const node_info_t& node_data,
                                                                      std::vector<stoat::PathTraversal>& walks)>& find_walks,
                                             bool alleles_requested,
                                             const std::function<std::vector<size_t>(const net_handle_t& node, const node_info_t& node_data,
                                                                                     const std::vector<stoat::sample_hap_t>& all_sample_haplotypes)>& find_alleles_by_sample,
                                             bool sequence_requested,
                                             const std::unordered_set<std::string>& reference_samples, 
                                             bool check_distances,
                                             stoat::Writer& out_writer, 
                                             bool keep_nodes) {

    // If we are going to write the nodes, then we are going to write all the nodes to a temporary file, then write the header (which isn't done until we find
    // all the nodes since there could be new references added), then copy the temporary file after the header
    std::string out_filename = out_writer.get_file_path();
    std::string out_temp_filename = out_filename + ".temp";
    std::shared_ptr<BgzWriter> temp_writer;

    // log info initialization
    number_node_analyzed = 0;

    if (out_filename != "") {

        // Make sure that the temporary file we write doesn't already exist 
        // Since the actual file name is given by the user (probably) it can be overwritten (probably)
        while (std::filesystem::exists(out_temp_filename)) {
            out_temp_filename += "1";
        }
        out_temp_filename += ".gz";
        temp_writer.reset(new BgzWriter(out_temp_filename));
    }

    // Make a copy of the samples. Not a reference since it has to be able to be loaded from a file and stored
    all_sample_haplotypes = sample_haplotypes;
    size_t sample_index = 0;

    // Fill in sample_to_index
    for (const sample_hap_t& sample_hap : all_sample_haplotypes) {
        if (!sample_to_index.count(sample_hap.sample)) {
            sample_to_index.emplace(sample_hap.sample, sample_index++);
        }
    }
    
    // Keep track of which references we've seen and their index in reference_names
    std::unordered_map<std::string, size_t> reference_name_to_index;

    // Get all the references in reference_samples first
    for (const std::string& ref_name : reference_samples) {
        reference_names.emplace_back(ref_name);
        reference_name_to_index.emplace(ref_name, reference_names.size() - 1);
    }

    // Go through the contents of chains in parallel
    // Everything touching chains needs to be in an omp critical block so they don't collide. 
    #pragma omp parallel shared(reference_names, all_sample_haplotypes, number_node_analyzed, all_node_data, out_writer)
    {
        // The actual while loop is run on a single thread
        #pragma omp single
        {
            // write header to file
            write_node_data_collection_header(out_writer);
            #pragma omp task
            {
                // Iterate over handle_t
                graph.for_each_handle([&](const handle_t& handle) {

                    size_t node_id = graph.get_id(handle);
                    bool is_reverse = graph.get_is_reverse(handle);
                    node_traversal_t node(node_id, is_reverse);
                    number_node_analyzed++;

                    // node sequence
                    std::string node_sequence = graph.get_sequence(handle);

                    // Compute depth
                    net_handle_t net = distance_index.get_net(handle, &graph);
                    size_t depth = distance_index.get_depth(net) - 1;
                    path_handle_t path;
                    std::string ref_name;
                    size_t position;

                    // Iterate over each step handle
                    graph.for_each_step_on_handle(handle, [&](const step_handle_t& step) {
            
                        path = graph.get_path_handle_of_step(step);
                        ref_name = graph.get_path_name(path);
                        ref_name = reference_name_to_index.find(ref_name) == reference_name_to_index.end() ? "NA" : ref_name;
                        position = graph.get_position_of_step(step);

                        return ref_name == "NA"; // if this is not a reference path
                    });

                    allele_by_sample_t alleles_by_sample;
                    GenotypeTable empty_genotypes(std::unordered_map<std::string, size_t>(), 0);
                    node_info_t new_node_info(node, 
                                                ref_name,
                                                position,
                                                depth,
                                                empty_genotypes,
                                                all_sample_haplotypes,
                                                alleles_by_sample,
                                                node_sequence);
                    
                    // write node_info_t to file
                    if (out_filename != "") {
                        #pragma omp critical(node_collection)
                        {
                            write_node_data_collection_line(*temp_writer, new_node_info);
                        }
                    }
                });// end for each handle
            }// end omp task
        }// end omp single
    }// end omp shared

    stoat::LOG_INFO("Total number of node analyse: " + std::to_string(number_node_analyzed));
}

void NodeDataCollection::add_alleles_by_sample(
                const std::function<std::vector<size_t>(const node_info_t& node_data, 
                                                        const std::vector<stoat::sample_hap_t>& all_sample_haplotypes)>& find_alleles_by_sample,
                std::string chr) {

    #pragma omp parallel for schedule(dynamic)
    for (const node_info_internal_t& node_info : all_node_data) {

        // If we are limiting to a chromosome (by reference path), then skip anything not on this chromosome
        if (chr != "" && (node_info.reference_index == std::numeric_limits<size_t>::max() || chr != reference_names.at(node_info.reference_index))) {
            continue;
        }

        std::vector<size_t> empty_alleles_by_sample;
        
        // This might cause problems because it is a reference but it doesn't get used so I think its fine
        // I don't want to use the actual samples_to_index because then the empty genotype table with allocate memory for the vector
        GenotypeTable empty_genotypes(std::unordered_map<std::string, size_t>(), 0);
        std::string node_sequence;

        // Make the node_info_t from the information we have
        std::vector<PathTraversal> empty_walks (0); 
        std::vector<std::string> empty_sequences (0);
        node_info_t new_node_info(node_info.node,
                                node_info.reference_index == std::numeric_limits<size_t>::max() ? "NA" : reference_names.at(node_info.reference_index),
                                node_info.position,
                                node_info.depth,
                                empty_genotypes,
                                all_sample_haplotypes,
                                allele_by_sample_t(),
                                node_sequence);

        // Get the alleles from the node_info_t
        std::vector<size_t> new_alleles_by_sample = find_alleles_by_sample(new_node_info, all_sample_haplotypes);
        size_t max_allele = 0;
        for (size_t x : new_alleles_by_sample) {
            if (x != std::numeric_limits<size_t>::max()) {
                max_allele = std::max(x, max_allele);
            }
        }

        node_to_alleles_by_sample.emplace(node_info.node, allele_by_sample_t(max_allele+1, std::move(new_alleles_by_sample)));
        number_node_analyzed++;
    }
}

void NodeDataCollection::genotype_nodes_by_chr_from_vcf(std::vector<std::string>& sample_names, stoat_vcf::VCFParser& vcf_parser) {

    // we'll use this edge matrix object
    // TODO find the vector of sample names from the VCF header?
    stoat_vcf::AlleleBySampleMatrix<stoat::node_traversal_t> edge_matrix(sample_names, 0);

    // use the corresponding sample-haplotypes for this collection
    // remove any existing sample in the collection first
    all_sample_haplotypes.clear();

    // add two haplotypes per sample
    for (std::string sample_name: sample_names) {
        for (std::string hap_name: {"0", "1"}) {
            sample_hap_t samp_hap;
            samp_hap.sample = sample_name;
            samp_hap.haplotype = hap_name;
            all_sample_haplotypes.emplace_back(samp_hap);
        }
    }

    // Fill in sample_to_index
    size_t sample_index = 0;
    for (const sample_hap_t& sample_hap : all_sample_haplotypes) {
        if (!sample_to_index.count(sample_hap.sample)) {
            sample_to_index.emplace(sample_hap.sample, sample_index++);
        }
    }

    // now the index in the edge matrix should match the index in the collection sample-hap list
    // read the VCF by chunk, build the edge matrix and genotype each node
    // We assume that the vcf parser has read the header and is now ready to go through the nodes
    std::string chr = vcf_parser.get_next_chromosome_name();

    // Go through to the end of the VCF. Chunk by chromosome 
    while(chr != "") {

        // Skip chromosomes not in ref_chrs
        while (std::find(reference_names.begin(), reference_names.end(), chr) == reference_names.end()) {
            stoat::LOG_WARN("Chromosome " + chr + " not found in node paths file. Skipping.", 0);
            bool found_new_chr = false;

            // Just skip to the next one without doing anything
            vcf_parser.skip_to_next_chromosome(chr);
            
            chr = vcf_parser.get_next_chromosome_name();
            if (chr == "") {
                // If we've reached the end of the file, return
                return;
            }

            // chr is now the next chromosome we want to look at
        }

        // start analyzing this chromosome chr
        stoat::LOG_INFO("Analyzing chr : " + chr);
        auto timer_start_chr = std::chrono::high_resolution_clock::now();

        // prepare the edge matrix for this chromosome by reading the VCF
        // this will read to the end of this chr
        edge_matrix.load_vcf_chunk(vcf_parser, chr);

        auto timer_end_matrix = std::chrono::high_resolution_clock::now();
        stoat::LOG_INFO("Edge matrix construction for chr " + chr + " : " + std::to_string(std::chrono::duration<double>(timer_end_matrix - timer_start_chr).count()) + " s");

        add_alleles_by_sample([&] (const node_info_t& node_data, const std::vector<stoat::sample_hap_t>& all_sample_haplotypes) {
            // JEAN init with max of size_t which I believe means "absent"/"no allele"
            std::vector<size_t> allele_idx(all_sample_haplotypes.size(), std::numeric_limits<size_t>::max());

            for (size_t samp_hap_idx: edge_matrix.get_samples_on_node(node_data.node)) {
                allele_idx.at(samp_hap_idx) = 0;
            }

            return allele_idx;
        }, chr);

        stoat::LOG_INFO("Total number of node found in chr " + chr + " : " + std::to_string(number_node_analyzed));
        number_node_analyzed = 0; // reset for next chromosome

        auto timer_end_chr = std::chrono::high_resolution_clock::now();
        stoat::LOG_INFO("Node genotypes retrieved in chr " + chr + " : " + std::to_string(std::chrono::duration<double>(timer_end_chr - timer_end_matrix).count()) + " s");
        stoat::LOG_INFO("Total time for chr " + chr + " : " + std::to_string(std::chrono::duration<double>(timer_end_chr - timer_start_chr).count()) + " s");

        // The parser has now passed the current chromosome. Get the name of the next one
        chr = vcf_parser.get_next_chromosome_name();
    }
}

// Call interatee for all nodes
void NodeDataCollection::for_each_node(const std::function<void(node_info_t& node_info)>& iteratee) const {
    for (const node_info_internal_t& node_info : all_node_data) {
        run_iteratee_on_one_node(node_info, iteratee);
    }
}

void NodeDataCollection::for_each_node_in_file(stoat::Reader& in_reader, const std::function<void(node_info_t& node_info)>& iteratee) {

    load_node_data_collection_header(in_reader);

    std::string line;
    while (in_reader.getline(line)) {
        node_info_internal_t node_info = load_node_data_line(line);
        run_iteratee_on_one_node(node_info, iteratee);
    }
}

void NodeDataCollection::run_iteratee_on_one_node(const node_info_internal_t& internal_node_info, const std::function<void(node_info_t& node_info)>& iteratee) const {

    // GenotypeTable constructor takes a map from sample to index, and the number of alleles
    GenotypeTable genotypes(sample_to_index,
                        node_to_alleles_by_sample.count(internal_node_info.node) ? node_to_alleles_by_sample.at(internal_node_info.node).allele_count : 0);

    #ifdef DEBUG_NODE_DATA_COLLECTION
        std::cerr << " Make genotype table for " << sample_to_index.size() << " samples and " << (node_to_alleles_by_sample.count(internal_node_info.node) ? node_to_alleles_by_sample.at(internal_node_info.node).allele_count : 0) << " alleles" << std::endl; 
        for (const auto& pair : sample_to_index) {
            std::cerr << "\t" <<  pair.first << ": " << pair.second << std::endl;
        }
    #endif

    // Go through the alleles_by_sample vector for this node and add the counts to the genotype table
    // alleles_by_sample is a vector with the allele for each sample in all_sample_haplotypes
    if (node_to_alleles_by_sample.count(internal_node_info.node)) {
        const allele_by_sample_t alleles_by_sample = node_to_alleles_by_sample.at(internal_node_info.node);
        for (size_t sample_hap_i = 0; sample_hap_i < alleles_by_sample.alleles.size(); sample_hap_i++) {
            if (alleles_by_sample.alleles[sample_hap_i] != std::numeric_limits<size_t>::max()) {
                // JEAN ideally we would access the count in the collection by index but I'm not sure why so I'm using the map sample_to_index for now. Maybe Xian knows
                size_t sample_idx = this->sample_to_index.at(all_sample_haplotypes.at(sample_hap_i).sample);
                genotypes.increment_count(sample_idx, alleles_by_sample.alleles[sample_hap_i]);
            }
        }
    }

    allele_by_sample_t empty_alleles;
    std::vector<PathTraversal> empty_walks(0); 
    std::string empty_sequences;
    node_info_t new_node_info(internal_node_info.node, 
                          internal_node_info.reference_index == std::numeric_limits<size_t>::max() ? "NA" : reference_names.at(internal_node_info.reference_index),
                          internal_node_info.position,
                          internal_node_info.depth,
                          genotypes,
                          all_sample_haplotypes,
                          node_to_alleles_by_sample.count(internal_node_info.node) ? node_to_alleles_by_sample.at(internal_node_info.node) : empty_alleles,
                          node_to_sequences.count(internal_node_info.node) ? node_to_sequences.at(internal_node_info.node) : empty_sequences);
    iteratee(new_node_info);
}

void NodeDataCollection::write_node_data_collection_header(stoat::Writer& out_writer) const {
    std::stringstream outstream;
    
    // Write the header
    outstream << file_header << std::endl;

    // Next will be a list of reference path names.
    outstream << "#REFS" << std::endl;
    for (const auto& ref : reference_names) {
        outstream << "#" << ref << std::endl;
    }
    
    //Finally the nodes
    // Start with a header that will contain the names of all samples
    outstream << "#NODES" << std::endl;
    outstream << "#NODE\tREF_INDEX\tPOSITION\tDEPTH";

    // The header also includes a list of sample/haplotypes
    for (const auto& samp : all_sample_haplotypes) {
        outstream << "\t" << samp.sample << "#" + samp.haplotype;
    }
    outstream << std::endl;
    out_writer.write(outstream.str());
}

void NodeDataCollection::write_node_data_collection_line(stoat::Writer& out_writer, const node_info_t& node_data) const {

    std::stringstream outstream;
    
    // Start with just the contents of the node_info_internal_t
    outstream << node_data.node.to_string() << "\t"
              << node_data.ref_path << "\t"
              << node_data.position << "\t"
              << node_data.depth;

    outstream << std::endl;
    #pragma omp critical(write_node_collection) 
    {
        out_writer.write(outstream.str());
    }
}

void NodeDataCollection::write_node_data_collection(stoat::Writer& out_writer) const {
    write_node_data_collection_header(out_writer);

    // Now write the nodes, one per line
    for (const node_info_internal_t& node_data : all_node_data) {
        write_node_data_line(out_writer, node_data);
    }

    return;
}

void NodeDataCollection::write_node_data_line(stoat::Writer& out_writer, const node_info_internal_t& node_data) const {
    write_node_data_line(out_writer, node_data, node_to_sequences.empty() || node_to_sequences.count(node_data.node) == 0 ? nullptr : &node_to_sequences.at(node_data.node),
                          (node_to_alleles_by_sample.empty() || node_to_alleles_by_sample.count(node_data.node) == 0 || node_to_alleles_by_sample.at(node_data.node).alleles.size() == 0) ? nullptr : &node_to_alleles_by_sample.at(node_data.node)); 
}

// aucune instance de fonction surchargée "stoat::NodeDataCollection::write_node_data_line" ne correspond à la liste d'argumentsC/C++(304)
// node_data_collection.cpp(387, 5): les types d'arguments sont : (stoat::Writer, const stoat::NodeDataCollection::node_info_internal_t, const std::__1::string *, const stoat::allele_by_sample_t *)

void NodeDataCollection::write_node_data_line(stoat::Writer& out_writer, const node_info_internal_t& node_data, 
                                                const std::string* sequence, const allele_by_sample_t* alleles_by_sample) const {
    std::stringstream outstream;
    
    // Start with just the contents of the node_info_internal_t
    outstream << node_data.node.to_string() << "\t"
              << node_data.reference_index << "\t"
              << node_data.position << "\t"
              << node_data.depth << "\t";
    
    // Since writing can occur while building the collection, all the big maps could change so we need guards
    #pragma omp critical(node_collection)
    {

        // Add the sequences, if any, as comma separated strings
        if (sequence == nullptr) {
            outstream << ".";
        } else {
            outstream << sequence->c_str(); ;
        }
    
        // Next add the allele assignments, if there are any
        if (alleles_by_sample == nullptr) {
            // If we don't have alleles or we don't have alleles for this node or the alleles for this node are empty
            for (size_t i = 0 ; i <  all_sample_haplotypes.size() ; i++ ) {
                // . for a haplotype that doesn't traverse the node
                outstream << "\t.";
            }
        } else {
            for (size_t allele_num : alleles_by_sample->alleles) {
                if (allele_num == std::numeric_limits<size_t>::max()) {
                    // . for a haplotype that doesn't traverse the node
                    outstream << "\t.";
                } else {
                    outstream << "\t" << allele_num;
                }
            }
        }
    } // end critical
    
    outstream << std::endl;
    #pragma omp critical(write_node_collection) 
    {
        out_writer.write(outstream.str());
    }
}

void NodeDataCollection::load_node_data_collection_header(stoat::Reader& in_reader) {

    // Clear anything that has already been filled in, since we want it to match what was in the file
    all_node_data.clear();
    node_to_alleles_by_sample.clear();
    reference_names.clear();
    all_sample_haplotypes.clear();
    sample_to_index.clear();

    // Read the first line, which must match the header
    std::string line;
    in_reader.getline(line);
    if (line != file_header) {
        throw std::runtime_error("stoat: Node data file contains the wrong header: " + line);
    }

    // Next are the references
    in_reader.getline(line);
    if (line != "#REFS") {
        throw std::runtime_error("stoat: Node file is not formatted correctly");
    }

    in_reader.getline(line);
    
    // Get the reference path names from the file
    while (line != "#NODES") {
        std::string ref (line.begin()+1, line.end());
        reference_names.emplace_back(std::move(ref));
        in_reader.getline(line);
    }

    std::set<size_t> index_sample;
    for (auto idx_sample: sample_to_index) {
        index_sample.insert(idx_sample.second);
    }

    // The next header is  "#node\tEND_NODE\tREF\tSTART_OFFSET\tEND_OFFSET\tDEPTH\tALLELE_LENGTHS\tWALKS\tSEQUENCES", plus all of the sample/haplotypes
    in_reader.getline(line);
    {
        std::stringstream linestream(line);
        std::string sample_name;

        // First go through the first 9 things and ignore them
        for (size_t i = 0 ; i < 9 ; i++) {
            std::getline(linestream, sample_name, '\t');
        }

        while (std::getline(linestream, sample_name, '\t')) {
            all_sample_haplotypes.emplace_back(sample_name);
        }
    }

    // Fill in sample_to_index OR correct it if idx sample has been removed
    size_t sample_index = 0;
    for (const sample_hap_t& sample_hap : all_sample_haplotypes) {
        if (!sample_to_index.count(sample_hap.sample)) {
            sample_to_index.emplace(sample_hap.sample, sample_index++);
        }
    }
}

NodeDataCollection::node_info_internal_t NodeDataCollection::load_node_data_line(std::string& line) {

    node_info_internal_t node_info;
    std::stringstream linestream(line);
    std::string part;

    //Node start node traversal
    std::getline(linestream, part, '\t');
    node_info.node = stoat::node_traversal_t(part);

    // Index of reference
    std::getline(linestream, part, '\t');
    node_info.reference_index = std::stoull(part);
    
    // Start offset along reference
    std::getline(linestream, part, '\t');
    node_info.position = std::stoull(part);

    // Node depth
    std::getline(linestream, part, '\t');
    node_info.depth = std::stoull(part);

    // Fill this in with the contents of the line
    bool has_samples = false;
    std::vector<size_t> allele_assignments;
    allele_assignments.reserve(all_sample_haplotypes.size());
    size_t max_allele = 0;
    bool has_allele = false;
    while (std::getline(linestream, part, '\t')) {
        if (part == ".") {
            allele_assignments.emplace_back(std::numeric_limits<size_t>::max());
        } else {
            allele_assignments.emplace_back(std::stoull(part));
            if (allele_assignments.back() != std::numeric_limits<size_t>::max()) {
                has_allele = true;
                max_allele = std::max(max_allele, allele_assignments.back());
            }
            has_samples = true;
        }
    };
    
    #ifdef DEBUG_NODE_DATA_COLLECTION
    if (has_allele) {
        assert(max_allele+1 <= walk_count);
    }
    assert(allele_assignments.size() == all_sample_haplotypes.size()); 
    #endif

    if (has_samples) {
        node_to_alleles_by_sample[node_info.node] = allele_by_sample_t(has_allele ? max_allele+1 : 0, allele_assignments);
    }

    return node_info;
}

void NodeDataCollection::load_node_data_collection(stoat::Reader& in_reader, const bool header_only) {
    load_node_data_collection_header(in_reader);

    if (!header_only) {
        std::string line;

        // Get the nodes
        while (in_reader.getline(line)) {
            all_node_data.emplace_back(load_node_data_line(line));
        }
    }
}

std::unordered_map<std::string, size_t>& NodeDataCollection::get_sample_to_index_reference() {
    return sample_to_index;
}

std::unordered_map<std::string, size_t> NodeDataCollection::get_sample_to_index_copy() const {
    return sample_to_index;
}

void NodeDataCollection::get_all_walks_through_node(
        const handlegraph::PathPositionHandleGraph& graph, 
        const bdsg::SnarlDistanceIndex& distance_index,
        const net_handle_t& node, 
        const node_info_t& node_data, 
        std::vector<stoat::PathTraversal>& walks, 
        size_t walk_cycle_limit) {

    #ifdef DEBUG_NODE_DATA_COLLECTION
        std::cerr << "Get all possible walks through node " << distance_index.net_handle_as_string(node) << std::endl;
    #endif

    // Path exploration
    std::vector<std::vector<handlegraph::net_handle_t>> paths = {
        {distance_index.get_bound(node, false, true)}
    };

    std::vector<std::vector<handlegraph::net_handle_t>> walks_as_net_handles;

    // How many steps have we taken trying to enumerate paths? Includes all all paths
    size_t steps_taken = 0;
    bool break_node = false;

    // For each incomplete path in paths, walk out from the end and add a copy of the path plus each next step to paths
    // Do this until the path reaches the end
    while (!paths.empty()) {
        std::vector<handlegraph::net_handle_t> path = std::move(paths.back());
        paths.pop_back();

        std::unordered_map<handlegraph::net_handle_t, size_t> dict_path_occ;
        bool cycle = false;

        // TODO: Put this back
        for (const auto& net : path) {
            if (++dict_path_occ[net] > walk_cycle_limit + 1) {
                cycle = true;
                break;
            }
        }

        // Follow edges from the last element in path
        if (!path.empty()) {
            distance_index.follow_net_edges(path.back(), &graph, false, [&](const handlegraph::net_handle_t& next_child) {
                // If this is the bound of the node then we're done
                if (distance_index.is_sentinel(next_child)) {

                    size_t next_child_node_id = distance_index.node_id(distance_index.get_node_from_sentinel(next_child));
                    size_t first_element_path_node_id = distance_index.node_id(distance_index.get_node_from_sentinel(path[0]));

                    // Only keep the walk if it entered and exited the node at opposite sides
                    if (next_child_node_id != first_element_path_node_id) {
                        walks_as_net_handles.emplace_back(path);
                        walks_as_net_handles.back().push_back(next_child);
                    }

                } else {
                    //TODO: Look for the cycle sooner

                    if (cycle) { // Case where we find a loop
                        return false;
                    }
                    paths.emplace_back(path);
                    paths.back().push_back(next_child);
                }
                return true;
            });
        }
    }

    if (break_node) {
        walks_as_net_handles.clear();
    }

    #ifdef DEBUG_NODE_DATA_COLLECTION
        // Validate paths
        std::set<std::vector<handlegraph::net_handle_t>> found_walks;
        for (const auto& walk : walks_as_net_handles) {
            for (size_t i = 0 ; i < walk.size() - 1 ; i++) {
                assert(distance_index.distance_in_parent(node, walk[i], distance_index.flip(walk[i+1])) == 0);
            }
            assert(found_walks.count(walk) == 0);
            assert(walk.front() == distance_index.get_bound(node, false, true));
            assert(walk.back() == distance_index.get_bound(node, true, false));
            found_walks.insert(walk);
        }
    #endif

    walks = stoat::convert_path_traversals(distance_index, graph, walks_as_net_handles);  
    
    #ifdef DEBUG_NODE_DATA_COLLECTION
        // Validate paths
        std::cerr << "Found " << walks.size() << " paths through the node" << std::endl;
        for (const auto& walk : walks) {
            
            std::cerr << "\t" << walk.to_string() << std::endl;
        }
    #endif

    return;
}

}