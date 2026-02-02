#include "vcf_untangler.hpp"
#include "arg_parser.hpp"

#define DEBUG_VCF_UNTANGLER

using namespace stoat;
namespace stoat_vcf{

VCFUntangler::VCFUntangler(const std::string& vcf_filename){
    // Read through the header for both sets of pointers
    std::vector<std::string> list_samples;
    
    std::tie(list_samples, ptr_vcf_bounds, hdr_bounds, rec_bounds) = stoat_vcf::parseHeader(vcf_filename);
    std::tie(list_samples, ptr_vcf_genotypes, hdr_genotypes, rec_genotypes) = stoat_vcf::parseHeader(vcf_filename);

    sample_hap_count = list_samples.size()*2;
    snarl_count = 0;
}

void VCFUntangler::fill_in_snarls_for_chromosome(const std::string& chr) {
    // Clear out all the data and get the next chromosome
    snarl_in_to_out.clear();
    genotypes.clear();
    snarl_count = 0;

    fill_in_nested_snarl_bounds(chr);
    fill_in_nested_genotypes(chr);
}

void VCFUntangler::fill_in_nested_snarl_bounds(const std::string& chr) {
    // This goes through the vcf for this chromosome and fills in snarl_in_to_out, to map each start bound of a snarl to its end bound (and end to start)
    // and gives an id to each snarl
    //TODO: Make sure all the reading matches that of the edge matrix
    
    // loop over the VCF file for each line and stop where chr is different
    while ((bcf_read(ptr_vcf_bounds, hdr_bounds, rec_bounds) >= 0) && (chr == bcf_hdr_id2name(hdr_bounds, rec_bounds->rid))){

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
        #ifdef DEBUG_VCF_UNTANGLER
        std::cerr << "Add snarl bounds " << snarl_bounds.at(0).to_string() << " and " << snarl_bounds.at(1).to_string() << std::endl;
        #endif

        size_t snarl_num = snarl_count;
        // Save start mapping to end
        snarl_in_to_out.emplace(std::make_pair(snarl_bounds.at(0), std::make_pair(snarl_bounds.at(1), snarl_count)));

        // Save end mapping to start, in the opposite direction
        snarl_in_to_out.emplace(std::make_pair(snarl_bounds.at(1).get_flipped(), std::make_pair(snarl_bounds.at(0).get_flipped(), snarl_count)));
        snarl_count++;
    
    
    } ;
    
}

void VCFUntangler::fill_in_nested_genotypes(const std::string& chr) {
    // This goes through the vcf for this chromosome and fills in genotypes for each nested snarl
    // This assumes that fill_in_snarl_bounds has already been called

    // This is basically going to be a matrix of the presence/absence of each snarl for each sample/hap.
    // One row for each snarl, except flattened
    // Index into genotypes can be found with get_genotype_index()
    genotypes.resize(sample_hap_count * snarl_count, 0);
    
    // loop over the VCF file for each line and stop where chr is different
     while ((bcf_read(ptr_vcf_genotypes, hdr_genotypes, rec_genotypes) >= 0) && (chr == bcf_hdr_id2name(hdr_genotypes, rec_genotypes->rid))) {
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
                #ifdef DEBUG_VCF_UNTANGLER
                std::cerr << "For sample number " << sample_hap_index << " Found path: ";
                for (const auto& x : allele_paths[idx_path_allele]) {
                    std::cerr << x.to_string();
                }
                std::cerr << std::endl;
                #endif

                if (idx_path_allele != -1) { // If this has genotypes
                    // Since this will be the path through the snarl, including the start of this snarl, skip the first and last node
                    for (size_t i = 1 ; i < allele_paths[idx_path_allele].size()-1 ; i++ ) {
                        node_traversal_t node = allele_paths[idx_path_allele][i];
                        if (snarl_in_to_out.count(node)) {

                            #ifdef DEBUG_VCF_UNTANGLER
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


       
    
    };
    
}

bool VCFUntangler::does_sample_have_snarl(size_t sample_hap_index, stoat::node_traversal_t snarl_bound) {
    return genotypes.at(genotype_index(sample_hap_index, snarl_bound));
}

stoat::node_traversal_t VCFUntangler::get_opposite_snarl_bound(stoat::node_traversal_t snarl_bound) {
    if (snarl_in_to_out.count(snarl_bound)) {
        return snarl_in_to_out.at(snarl_bound).first;
    } else {
        return snarl_bound;
    }
}

size_t VCFUntangler::genotype_index(size_t sample_hap_index, const stoat::node_traversal_t& snarl_bound) {
    size_t snarl_index = snarl_in_to_out.at(snarl_bound).second;
    return snarl_index * sample_hap_count + sample_hap_index;
}



}//end namespace
