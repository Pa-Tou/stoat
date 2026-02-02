#include "vcf_parser.hpp"

//#define DEBUG_VCF_PARSER

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

    std::vector<std::string> list_samples;
    // Get the samples names
    for (int i = 0; i < bcf_hdr_nsamples(hdr); i++) {
        list_samples.push_back(bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i));
    }
    //TOD: This assumes that the ploidy is 2 but idk if that is always true in a vcf
    hap_count = list_samples.size() * 2;

    // Read the current line
    read_status = bcf_read(ptr_vcf, hdr, rec);
        
    return list_samples;
}

std::string VCFParser::get_next_chromosome_name() {
    if (read_status == -1) {
        return "";
    }
    return bcf_hdr_id2name(hdr, rec->rid);
}

void VCFParser::for_each_record_on_chromosome(const std::string& chr, const std::function<void(const vcf_info_t& vcf_info)>& iteratee) {
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
        
        // extract genotypes GT
        int ngt = 0;
        int32_t *gt = nullptr;
        ngt = bcf_get_genotypes(hdr, rec, &gt, &ngt);
        
        if (ngt <= 0 || gt == nullptr)
        {
            throw std::invalid_argument("GT field is missing in VCF at position " + std::to_string(rec->pos + 1));
        }
        // Make the actual vector of genotypes
        std::vector<int> genotypes;
        genotypes.reserve(hap_count);
        for (size_t i = 0 ; i < hap_count ; i++) {
            genotypes.emplace_back(bcf_gt_allele(gt[i]));
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
            paths.push_back(string_to_path_node_traversal(item));
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
        bcf_unpack(rec, BCF_UN_STR);

        read_status = bcf_read(ptr_vcf, hdr, rec);
#ifdef DEBUG_VCF_PARSER
        std::cerr << "\t" << read_status << " on chr " << chr << std::endl;
#endif

    } while ((read_status >= 0) && (chr == bcf_hdr_id2name(hdr, rec->rid)));

#ifdef DEBUG_VCF_PARSER
    std::cerr << " broke out of loop with " << read_status << " At chr " << bcf_hdr_id2name(hdr, rec->rid) << std::endl;
#endif
}

}//end namespace
