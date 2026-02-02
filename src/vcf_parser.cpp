#include "vcf_parser.hpp"

#define DEBUG_VCF_PARSER

using namespace stoat;
namespace stoat_vcf{

std::vector<std::string> VCFParser::initialize_parser(const std::string& vcf_path) {

    // Open the VCF file
    ptr_vcf = bcf_open(vcf_path.c_str(), "r");
    
    // Read the VCF header
    hdr = bcf_hdr_read(ptr_vcf);
    if (!hdr) {
        bcf_close(ptr_vcf);
        throw std::invalid_argument("Could not read VCF header");
    }
    // Initialize a record
    rec = bcf_init();
    if (!rec) {
        bcf_hdr_destroy(hdr);
        bcf_close(ptr_vcf);
        throw std::invalid_argument("Failed to allocate memory for VCF record");
    }

    // If we want to untangle the snarls, then also open readers for the untangling steps
    if (untangle_snarls) {
        ptr_vcf_bounds = bcf_open(vcf_path.c_str(), "r");
        ptr_vcf_genotypes = bcf_open(vcf_path.c_str(), "r");
        hdr_bounds = bcf_hdr_read(ptr_vcf_bounds);
        hdr_genotypes = bcf_hdr_read(ptr_vcf_genotypes);
        rec_bounds = bcf_init();
        rec_genotypes = bcf_init();
    }

    std::vector<std::string> list_samples;
    // Get the samples names
    for (int i = 0; i < bcf_hdr_nsamples(hdr); i++) {
        list_samples.push_back(bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i));
    }
    //TOD: This assumes that the ploidy is 2 but idk if that is always true in a vcf
    hap_count = list_samples.size() * 2;

    // Read the current line
    read_status = bcf_read(ptr_vcf, hdr, rec);
    if (untangle_snarls) {
        bcf_read(ptr_vcf_bounds, hdr_bounds, rec_bounds);
        bcf_read(ptr_vcf_genotypes, hdr_genotypes, rec_genotypes);
    }
        
    return list_samples;
}

std::string VCFParser::get_next_chromosome_name() {
    if (read_status == -1) {
        return "";
    }
    return bcf_hdr_id2name(hdr, rec->rid);
}

void VCFParser::for_each_record_on_chromosome(const std::string& chr, const std::function<void(const vcf_info_t& vcf_info)>& iteratee) {

    if (untangle_snarls) {
        std::cerr << "Untangle snarls for chr " << chr << std::endl;
        // If we are going to untangle stuff, process the snarls first
        // Clear out all the data and get the next chromosome
        snarl_in_to_out.clear();
        genotypes.clear();
        snarl_count = 0;

        fill_in_nested_snarl_bounds(chr);
        fill_in_nested_genotypes(chr);
        std::cerr << "Found " << snarl_count << " snarls" << std::endl;
    }

    // Since we've already read the first line of this chunk, do a do-while loop and read the next at the end.
    // At the end of this loop, we'll be looking at the first line that is not this chromosome
    do {
        bcf_unpack(rec, BCF_UN_STR);

        int32_t *lv = nullptr;
        int n_lv = 0;
        size_t level;
        if (bcf_get_info_int32(hdr, rec, "LV", &lv, &n_lv) > 0)
        {
            level = lv[0];
        }
        free(lv);

        // Get the snarl bounds, which are saved in the VCF as the ID
        std::string snarl_bounds_string (rec->d.id);

        // This should be a vector of two node_traversal_t's of the snarl bounds, first one pointing in, second one pointing out
        std::vector<stoat::node_traversal_t> snarl_bounds = string_to_path_node_traversal(snarl_bounds_string);
        
        // extract genotypes GT
        int ngt = 0;
        int32_t *gt = nullptr;
        ngt = bcf_get_genotypes(hdr, rec, &gt, &ngt);
        
        if (ngt <= 0 || gt == nullptr)
        {
            throw std::invalid_argument("GT field is missing in VCF at position " + std::to_string(rec->pos + 1));
        }

        // Make the actual vector of genotypes
        // If we want to untangle the snarls, then check that the parent snarl actually was genotyped as having this child snarl
        std::vector<int> genotypes;
        genotypes.reserve(hap_count);
        for (size_t i = 0 ; i < hap_count ; i++) {
            if (!untangle_snarls || level == 0 || does_sample_have_snarl(i, snarl_bounds[0])) {
                // Always keep the genotype if we don't untangle snarls, or if this is a top-level snarl 
                genotypes.emplace_back(bcf_gt_allele(gt[i]));
            } else {
                genotypes.emplace_back((int)-1);
            }
        }
        free(gt);
        
        // extract AT field from INFO
        char *at = nullptr;
        int nat = 0;
        nat = bcf_get_info_string(hdr, rec, "AT", &at, &nat);
        if (nat <= 0 || !at)
        {
            // AT field is mandatory, throw an error
            throw std::invalid_argument("AT field is missing in VCF at position " + std::to_string(rec->pos + 1));
        }
        std::string at_str(at); // convert to C++ std::string
        free(at);

        // split by comma and save as a vector of edge lists [vector vector stoat::edge_t]
        // from: ">123>213<234", ">123<234", ">123<234<345"
        // to: [[edge_t(123, 213),stoat::edge_t(213, 234)], [...]]
        std::vector<std::vector<stoat::node_traversal_t>> paths;
        std::stringstream at_ss(at_str);
        std::string item;
        while (std::getline(at_ss, item, ','))
        {
            // If we are untangling snarls, then skip any nested snarls
            std::vector<stoat::node_traversal_t> path_as_nodes = string_to_path_node_traversal(item);
            if (untangle_snarls) {
                // If we want to resolve snarls, remove any nested snarls by copying everything except nested snarls into new path
                // Add a <0 node to represent the snarl
                std::vector<stoat::node_traversal_t> filtered_path;
                filtered_path.reserve(path_as_nodes.size());
                size_t path_i = 0;
                while (path_i < path_as_nodes.size()) {

                    // Add the current node
                    filtered_path.emplace_back(path_as_nodes[path_i]);

                    // Check if the current node is the start of a snarl
                    if (path_i != 0 && path_i != path_as_nodes.size()-1) {
                        stoat::node_traversal_t skip_to_node = get_opposite_snarl_bound(filtered_path.back());
                        if (skip_to_node != filtered_path.back()) {
                            // If this is the start of a snarl, a fake snarl node, skip to the end of the snarl, and add the end.
                            // The loop should continue on the node after the end node of the nested snarl
                            filtered_path.emplace_back(0, false);
                            while (path_as_nodes[path_i] != skip_to_node) {
                                path_i++;
                            }
                            assert(path_as_nodes[path_i] == skip_to_node);
                            filtered_path.emplace_back(path_as_nodes[path_i]);
                        }
                    }
                    path_i++;
                }
                paths.push_back(std::move(filtered_path));
            } else {
                paths.push_back(std::move(path_as_nodes));
            }
            
        }

        iteratee(vcf_info_t({level, genotypes, paths}));

        read_status = bcf_read(ptr_vcf, hdr, rec);
#ifdef DEBUG_VCF_PARSER
        std::cerr << read_status << " on chr " << chr << std::endl;
#endif

    } while ((read_status >= 0) && (chr == bcf_hdr_id2name(hdr, rec->rid)));

#ifdef DEBUG_VCF_PARSER
    std::cerr << " broke out of loop with " << read_status << " At chr " << bcf_hdr_id2name(hdr, rec->rid) << std::endl;
#endif
}

void VCFParser::skip_to_next_chromosome(const std::string& chr) {
#ifdef DEBUG_VCF_PARSER
        std::cerr << "Skip through chr " << chr << std::endl;
#endif

    // Since we've already read the first line of this chunk, do a do-while loop and read the next at the end.
    // At the end of this loop, we'll be looking at the first line that is not this chromosome
    do {
        //TODO: I think this is unecessary
        //bcf_unpack(rec, BCF_UN_STR);

        read_status = bcf_read(ptr_vcf, hdr, rec);
        if (untangle_snarls) {
            bcf_read(ptr_vcf_bounds, hdr_bounds, rec_bounds);
            bcf_read(ptr_vcf_genotypes, hdr_genotypes, rec_genotypes);
        }
#ifdef DEBUG_VCF_PARSER
        std::cerr << "\t" << read_status << " on chr " << chr << std::endl;
#endif

    } while ((read_status >= 0) && (chr == bcf_hdr_id2name(hdr, rec->rid)));

#ifdef DEBUG_VCF_PARSER
    std::cerr << " broke out of loop with " << read_status << " At chr " << bcf_hdr_id2name(hdr, rec->rid) << std::endl;
#endif
}


void VCFParser::fill_in_nested_snarl_bounds(const std::string& chr) {
    // This goes through the vcf for this chromosome and fills in snarl_in_to_out, to map each start bound of a snarl to its end bound (and end to start)
    // and gives an id to each snarl
    //TODO: Make sure all the reading matches that of the edge matrix
    
    // loop over the VCF file for each line and stop where chr is different
    do {

        // Unpack the vcf up to ALT field
        bcf_unpack(rec_bounds, BCF_UN_STR);
    
        // check the INFO field for LV (Level in the snarl tree) so we can skip LV=0
        int32_t *lv = nullptr;
        int n_lv = 0;
        if (bcf_get_info_int32(hdr_bounds, rec_bounds, "LV", &lv, &n_lv) > 0) {
            // Skip LV=0 snarls
            if (lv[0] == 0) {
                free(lv);
                continue;
            }
        }
        free(lv);

    
        // Get the snarl bounds, which are saved in the VCF as the ID
        std::string snarl_bounds_string (rec_bounds->d.id);

        // This should be a vector of two node_traversal_t's of the snarl bounds, first one pointing in, second one pointing out
        std::vector<stoat::node_traversal_t> snarl_bounds = string_to_path_node_traversal(snarl_bounds_string);
        #ifdef DEBUG_VCF_PARSER
        std::cerr << "Add snarl bounds " << snarl_bounds.at(0).to_string() << " and " << snarl_bounds.at(1).to_string() << std::endl;
        #endif

        size_t snarl_num = snarl_count;
        // Save start mapping to end
        snarl_in_to_out.emplace(std::make_pair(snarl_bounds.at(0), std::make_pair(snarl_bounds.at(1), snarl_count)));

        // Save end mapping to start, in the opposite direction
        snarl_in_to_out.emplace(std::make_pair(snarl_bounds.at(1).get_flipped(), std::make_pair(snarl_bounds.at(0).get_flipped(), snarl_count)));
        snarl_count++;
    
    
    } while ((bcf_read(ptr_vcf_bounds, hdr_bounds, rec_bounds) >= 0) && (chr == bcf_hdr_id2name(hdr_bounds, rec_bounds->rid)));
    
}

void VCFParser::fill_in_nested_genotypes(const std::string& chr) {
    // This goes through the vcf for this chromosome and fills in genotypes for each nested snarl
    // This assumes that fill_in_snarl_bounds has already been called

    // This is basically going to be a matrix of the presence/absence of each snarl for each sample/hap.
    // One row for each snarl, except flattened
    // Index into genotypes can be found with get_genotype_index()
    genotypes.resize(hap_count * snarl_count, 0);
    
    // loop over the VCF file for each line and stop where chr is different
     do {
        // Unpack the vcf up to ALT field
        bcf_unpack(rec_genotypes, BCF_UN_STR);

        // extract genotypes GT
        int ngt = 0;
        int32_t *gt = nullptr;
        ngt = bcf_get_genotypes(hdr_genotypes, rec_genotypes, &gt, &ngt);
        
        if (ngt <= 0 || gt == nullptr) {
            throw std::invalid_argument("GT field is missing in VCF at position " + std::to_string(rec_genotypes->pos + 1));
        }

        // extract AT field from INFO
        char *at = nullptr;
        int nat = 0;
        nat = bcf_get_info_string(hdr_genotypes, rec_genotypes, "AT", &at, &nat);
        if (nat <= 0 || !at) {
            // AT field is mandatory, throw an error
            throw std::invalid_argument("AT field is missing in VCF at position " + std::to_string(rec_genotypes->pos + 1));
        }
        std::string at_str(at); // convert to C++ std::string
        free(at);
        
        // split by comma and save as a vector of edge lists [vector vector stoat::edge_t]
        // from: ">123>213<234", ">123<234", ">123<234<345"
        // to: [[edge_t(123, 213),stoat::edge_t(213, 234)], [...]]
        std::vector<std::vector<stoat::node_traversal_t>> allele_paths;
        std::stringstream at_ss(at_str);
        std::string item;
        while (std::getline(at_ss, item, ',')) {
            std::vector<stoat::node_traversal_t> path_as_node_traversal = string_to_path_node_traversal(item);
            allele_paths.push_back(std::move(path_as_node_traversal));
        }
        
        // Now go through the paths and for each snarl in the path, remember how many copies of the snarl we see
        for (int sample_num = 0; sample_num < rec_genotypes->n_sample; ++sample_num){
            int ploidy = 2;
            for (int hap_num = 0; hap_num < ploidy; ++hap_num){
                // allele hap_num of that sample
                // JEAN here we are assuming diploid genotypes. check how to make sure we're really/always getting the genotype for sample sample_num with bcf_gt_allele
                size_t sample_hap_index = sample_num*2 + hap_num;

                int idx_path_allele = bcf_gt_allele(gt[sample_hap_index]);

                if (idx_path_allele != -1) { // If this has genotypes
                    #ifdef DEBUG_VCF_PARSER
                    std::cerr << "For sample number " << sample_hap_index << " Found path for allele number " << idx_path_allele << ": ";
                    for (const auto& x : allele_paths[idx_path_allele]) {
                        //std::cerr << x.to_string();
                    }
                    std::cerr << std::endl;
                    #endif
                    // Since this will be the path through the snarl, including the start of this snarl, skip the first and last node
                    for (size_t i = 1 ; i < allele_paths[idx_path_allele].size()-1 ; i++ ) {
                        node_traversal_t node = allele_paths[idx_path_allele][i];
                        if (snarl_in_to_out.count(node)) {

                            #ifdef DEBUG_VCF_PARSER
                            std::cerr << "\tadd snarl starting at " << node.to_string() << std::endl;
                            #endif

                            std::pair<stoat::node_traversal_t,size_t> snarl_end = snarl_in_to_out.at(node);

                            // Save the genotype
                            genotypes[genotype_index(sample_hap_index, node)] = 1;

                            // Skip to the end of the snarl
                            // The next iteration of the for loop should have node being the end node.
                            // TODO Double check that this is true
                            while (allele_paths[idx_path_allele][i+1] != snarl_end.first) {
                                i++;
                            }
                        }
                    }
                }
            }
        }
        free(gt);

    } while ((bcf_read(ptr_vcf_genotypes, hdr_genotypes, rec_genotypes) >= 0) && (chr == bcf_hdr_id2name(hdr_genotypes, rec_genotypes->rid)));
    
}

bool VCFParser::does_sample_have_snarl(size_t sample_hap_index, stoat::node_traversal_t snarl_bound) {
    if (snarl_in_to_out.count(snarl_bound)) {
        return genotypes.at(genotype_index(sample_hap_index, snarl_bound));
    } else {
        // If this snarl wasn't saved, then it must be a top-level snarl so it is always present
        return true;
    }
}

stoat::node_traversal_t VCFParser::get_opposite_snarl_bound(stoat::node_traversal_t snarl_bound) {
    if (snarl_in_to_out.count(snarl_bound)) {
        return snarl_in_to_out.at(snarl_bound).first;
    } else {
        return snarl_bound;
    }
}

size_t VCFParser::genotype_index(size_t sample_hap_index, const stoat::node_traversal_t& snarl_bound) {
    size_t snarl_index = snarl_in_to_out.at(snarl_bound).second;
    return snarl_index * hap_count + sample_hap_index;
}

void VCFParser::close_vcf(){
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(ptr_vcf);

    if (untangle_snarls) {
        bcf_destroy(rec_bounds);
        bcf_hdr_destroy(hdr_bounds);
        bcf_close(ptr_vcf_bounds);

        bcf_destroy(rec_genotypes);
        bcf_hdr_destroy(hdr_genotypes);
        bcf_close(ptr_vcf_genotypes);
    }
}


}//end namespace


