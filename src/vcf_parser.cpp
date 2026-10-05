#include <algorithm>
#include <limits>
#include <atomic>
#include <omp.h>

#include "vcf_parser.hpp"

//#define DEBUG_VCF_PARSER

using namespace stoat;
namespace stoat_vcf{

void VCFParser::parse_haplotype_counts(const std::string& counts) {
    std::unordered_map<std::string, size_t> parsed_counts;
    size_t start = 0;

    while (start < counts.size()) {
        const size_t end = counts.find(',', start);
        const std::string entry = trim(counts.substr(
            start,
            end == std::string::npos ? std::string::npos : end - start
        ));

        const size_t colon = entry.find(':');
        if (colon == std::string::npos || colon == 0 || colon == entry.size() - 1 ||
            entry.find(':', colon + 1) != std::string::npos) {
            throw std::invalid_argument("Invalid entry; expected chromosome:count: " + entry);
        }

        const std::string chromosome = trim(entry.substr(0, colon));
        const std::string count_text = trim(entry.substr(colon + 1));
        if (chromosome.empty()) {
            throw std::invalid_argument("Invalid entry; expected chromosome:count: " + entry);
        }

        const size_t count = parse_count(count_text, ": " + count_text);
        if (!parsed_counts.emplace(chromosome, count).second) {
            throw std::invalid_argument("duplicate chromosome: " + chromosome);
        }

        if (end == std::string::npos) {
            break;
        }

        start = end + 1;

        if (start == counts.size()) {
            throw std::invalid_argument("Invalid entry; expected chromosome:count after trailing comma");
        }
    }

    if (parsed_counts.empty()) {
        throw std::invalid_argument("Expected at least one chromosome:count entry");
    }
    chr_haplotype_counts = std::move(parsed_counts);
}

void VCFParser::load_haplotype_counts_file(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open haplotype-count file: " + filename);
    }

    std::unordered_map<std::string, size_t> parsed_counts;
    std::string line;
    size_t line_number = 0;

    while (std::getline(file, line)) {
        ++line_number;

        if (trim(line).empty()) {
            continue;
        }

        if (trim(line) == "chromosome\thaplotype_count") {
            continue;
        }

        const size_t tab = line.find('\t');
        if (tab == std::string::npos || line.find('\t', tab + 1) != std::string::npos) {
            throw std::runtime_error("Expected exactly two tab-separated columns at line " +
                std::to_string(line_number));
        }

        const std::string chromosome = trim(line.substr(0, tab));
        const std::string count_text = trim(line.substr(tab + 1));
        if (chromosome.empty()) {
            throw std::runtime_error("Missing chromosome at line " + std::to_string(line_number));
        }

        size_t count;
        try {
            count = parse_count(count_text, " at line " + std::to_string(line_number));
        } catch (const std::invalid_argument& error) {
            throw std::runtime_error(error.what());
        }

        if (!parsed_counts.emplace(chromosome, count).second) {
            throw std::runtime_error("duplicate chromosome at line " + std::to_string(line_number));
        }
    }

    if (file.bad()) {
        throw std::runtime_error("Error reading ploidy count file: " + filename);
    }

    if (parsed_counts.empty()) {
        throw std::runtime_error("No haplotype counts found in file: " + filename);
    }
    chr_haplotype_counts = std::move(parsed_counts);
}

void VCFParser::initialize_parser(const std::string& vcf_path) {

    // Open the VCF file
    ptr_vcf = bcf_open(vcf_path.c_str(), "r");
    hts_set_threads(ptr_vcf, omp_get_max_threads());
    
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

    // Get the samples names
    for (int i = 0; i < bcf_hdr_nsamples(hdr); i++) {
        sample_names.push_back(bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i));
    }

    size_t sample_count = sample_names.size();
    if (sample_count == 0) {
        throw std::invalid_argument("No samples found in VCF file");
    }
    hap_count = sample_count * ploidy;

    // Read the current line
    read_status = bcf_read(ptr_vcf, hdr, rec);

    // If we want to untangle the snarls, then also open readers for the untangling steps
    if (resolve_nested_calls) {
        ptr_vcf_bounds = bcf_open(vcf_path.c_str(), "r");
        hts_set_threads(ptr_vcf_bounds, omp_get_max_threads());

        ptr_vcf_genotypes = bcf_open(vcf_path.c_str(), "r");
        hts_set_threads(ptr_vcf_genotypes, omp_get_max_threads());

        hdr_bounds = bcf_hdr_read(ptr_vcf_bounds);
        hdr_genotypes = bcf_hdr_read(ptr_vcf_genotypes);
        rec_bounds = bcf_init();
        rec_genotypes = bcf_init();

        bcf_read(ptr_vcf_bounds, hdr_bounds, rec_bounds);
        bcf_read(ptr_vcf_genotypes, hdr_genotypes, rec_genotypes);
    }
}

void VCFParser::set_chromosome_ploidy(const std::string& chromosome) {
    const auto user_ploidy_count = chr_haplotype_counts.find(chromosome);
    size_t inferred_ploidy = 0;

    if (read_status >= 0 && chromosome == bcf_hdr_id2name(hdr, rec->rid)) {
        int32_t* gt = nullptr;
        int genotype_count = 0;
        genotype_count = bcf_get_genotypes(hdr, rec, &gt, &genotype_count);

        if (genotype_count > 0 && gt != nullptr && !sample_names.empty()) {
            const size_t sample_count = sample_names.size();
            const size_t slots_per_sample = static_cast<size_t>(genotype_count) / sample_count;

            for (size_t sample = 0; sample < sample_count; ++sample) {
                size_t sample_ploidy = 0;

                while (sample_ploidy < slots_per_sample && gt[sample * slots_per_sample + sample_ploidy] != bcf_int32_vector_end) {
                    ++sample_ploidy;
                }

                inferred_ploidy = std::max(inferred_ploidy, sample_ploidy);
            }
        }

        free(gt);
    }

    // Chr ploidy was provided by the user
    if (user_ploidy_count != chr_haplotype_counts.end()) {
        size_t user_ploidy = user_ploidy_count->second;

        // Case where the user-provided haplotype count is higher than the inferred ploidy
        // report a warning but still use the user-provided count
        if (user_ploidy > inferred_ploidy) {
            stoat::LOG_WARN("User-provided haplotype count (" + std::to_string(user_ploidy) +
                ") for chromosome " + chromosome + " is greater than the inferred count (" +
                std::to_string(inferred_ploidy) + "); using the user-provided count.",
                "user_haplotype_count_higher_than_inferred");
        }

        // Case where the user-provided haplotype count is lower than the inferred ploidy
        // report a warning but still use the user-provided count
        if (user_ploidy < inferred_ploidy) {
            stoat::LOG_WARN("User-provided haplotype count (" + std::to_string(user_ploidy) +
                ") for chromosome " + chromosome + " is lower than the inferred count (" +
                std::to_string(inferred_ploidy) + "); using the user-provided count.",
                "user_haplotype_count_lower_than_inferred");
        }

    } else {
        chr_haplotype_counts.emplace(chromosome, inferred_ploidy);
    }

    ploidy = chr_haplotype_counts.at(chromosome);
    hap_count = sample_names.size() * ploidy;
}

std::string VCFParser::get_next_chromosome_name() {
    if (read_status == -1) {
        return "";
    }
    return bcf_hdr_id2name(hdr, rec->rid);
}

void VCFParser::for_each_record_on_chromosome(const std::string& chr, const std::function<void(const vcf_info_t& vcf_info)>& iteratee) {

    set_chromosome_ploidy(chr);

    if (resolve_nested_calls) {
        // If we are going to untangle stuff, process the snarls first
        // Clear out all the data and get the next chromosome
        snarl_in_to_out.clear();
        genotypes.clear();
        snarl_count = 0;

        fill_in_nested_snarl_bounds(chr);
        fill_in_nested_genotypes(chr);
    }

    // Process the chromosome chunk by chunk.
    while (read_status >= 0 && chr == bcf_hdr_id2name(hdr, rec->rid)) {

        // Buffer the next chunk of records in a vector of unique pointers so that they are freed when the vector is cleared.
        std::vector<Bcf1Ptr> raw_records;
        raw_records.reserve(CHUNK_SIZE);

        // ------------------------------------------------------------
        // Phase 1: serial VCF reading
        // ------------------------------------------------------------
        do {
            raw_records.emplace_back(bcf_dup(rec), &bcf_destroy);
            read_status = bcf_read(ptr_vcf, hdr, rec);

#ifdef DEBUG_VCF_PARSER
    std::cerr << read_status << " on chr " << chr << std::endl;
#endif
            // If there was a problem with the VCF, stop. The bcftools reader should have output its own more informative error message but doesn't seem to throw an error
            if (read_status < -1) {
                throw std::runtime_error("Unable to read VCF file");
            }
        } while (raw_records.size() < CHUNK_SIZE && read_status >= 0 && chr == bcf_hdr_id2name(hdr, rec->rid));

        // ------------------------------------------------------------
        // Phase 2: parallel processing of this chunk
        // ------------------------------------------------------------
        std::exception_ptr parse_exception = nullptr;
        bool has_error = false;

        #pragma omp parallel for schedule(static)
        for (size_t record_i = 0; record_i < raw_records.size(); ++record_i) {

            // If anything has failed, skip the rest of the chunk and throw the error later
            if (has_error) {
                continue;
            }

            try {
                // Get the raw record from the vector of unique pointers
                bcf1_t* record = raw_records.at(record_i).get();

                // CPU-intensive operation: do this in parallel.
                const vcf_info_t vcf_info = parse_record(record, chr);
                iteratee(vcf_info);

            // If anything has failed, skip the rest of the chunk and throw the error later
            // This is a bit of a hack to get around the fact that OpenMP doesn't support exceptions.
            // We just set a flag and throw the exception after the parallel region.
            } catch (...) {
                #pragma omp atomic write
                has_error = true;
                #pragma omp critical(vcf_parser_exception)
                {
                    parse_exception = std::current_exception();
                    has_error = true;
                }
            }
        }

        // Free the whole chunk before loading the next one.
        raw_records.clear();
        if (parse_exception) {
            std::rethrow_exception(parse_exception);
        }

#ifdef DEBUG_VCF_PARSER
    if (read_status >= 0) {
        std::cerr << "Finished chunk on chr " << chr << std::endl;
    }
#endif

    }

#ifdef DEBUG_VCF_PARSER
    if (read_status >= 0) {
        std::cerr << "Broke out of chromosome loop with " << read_status << " at chr " << bcf_hdr_id2name(hdr, rec->rid) << std::endl;
    }
#endif

}

vcf_info_t VCFParser::parse_record(bcf1_t* raw_record, const std::string& chr) {
    bcf_unpack(raw_record, BCF_UN_STR);

    int32_t *lv = nullptr;
    int n_lv = 0;
    // Default to LV=0 if it wasn't there
    size_t level = 0;
    if (bcf_get_info_int32(hdr, raw_record, "LV", &lv, &n_lv) > 0) {
        level = lv[0];
    }
    free(lv);

    // For a vg call vcf, the snarl id is the snarl bounds
    std::string snarl_id (raw_record->d.id);

    // Get the paths of the alleles. This is either from the AT and RT fields (vg call) or the ID field (pangenie)
    std::vector<std::vector<stoat::node_traversal_t>> paths;

    // extract AT field from INFO
    char *at = nullptr;
    int nat = 0;
    nat = bcf_get_info_string(hdr, raw_record, "AT", &at, &nat);

    char *id_field = nullptr;
    int nid = 0;
    nid = bcf_get_info_string(hdr, raw_record, "ID", &id_field, &nid);
    if ((nat > 0 && at) || (nid > 0 && id_field)) {
        std::string info_str;
        if (nat>0 && at) {
            // If there is an AT field, then this is a vg call vcf with the paths directly in the AT field
            std::string at_str(at); // convert to C++ std::string
            free(at);
            free(id_field);
            info_str= std::move(at_str);
        } else {
            std::string id_str(id_field); // convert to C++ std::string
            free(id_field);
            free(at);

            // If there is an ID field, then this is a pangenie vcf with the paths as part of the ID
            // Also add the RD field and put them together
            char *rd_field = nullptr;
            int nrd = 0;
            nrd = bcf_get_info_string(hdr, raw_record, "RD", &rd_field, &nrd);
            if (nrd <= 0 || !rd_field) {
                throw std::invalid_argument("VCF contains ID field but not RD field " + std::to_string(raw_record->pos + 1) + "\n\tThis pangenie VCF contains the alt paths in the ID field but not the reference path in the RD field");
            }
            std::string rd_str(rd_field);
            free(rd_field);
            rd_str+=",";
            rd_str.insert(rd_str.end(), std::make_move_iterator(id_str.begin()), std::make_move_iterator(id_str.end()));

            info_str= std::move(rd_str);
        }

        // split by comma and save as a vector of edge lists [vector vector stoat::edge_t]
        // from: ">123>213<234", ">123<234", ">123<234<345"
        // to: [[edge_t(123, 213),stoat::edge_t(213, 234)], [...]]
        std::stringstream info_ss(info_str);
        std::string item;
        while (std::getline(info_ss, item, ','))
        {
            // If we are untangling snarls, then skip any nested snarls
            std::vector<stoat::node_traversal_t> path_as_nodes;
            if (nat > 0 && at) {
                path_as_nodes  = string_to_path_node_traversal(item);
            } else {
                path_as_nodes = parse_pangenie_id(item);
            }

            if (resolve_nested_calls) {
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
                            // If this is the start of a snarl, add a fake snarl node and skip to the end of the snarl
                            // The loop should continue on the end node of the nested snarl
                            // TODO: This will include the boundary nodes between snarls which is wasteful but simpler
                            filtered_path.emplace_back(0, false);
                            while (path_as_nodes[path_i] != skip_to_node) {
                                path_i++;
                            }
                        } else {
                            path_i++;
                        }
                    } else {
                        path_i++;
                    }
                }
                paths.push_back(std::move(filtered_path));
            } else {
                paths.push_back(std::move(path_as_nodes));
            }

        }

    // End if ID field
    } else {
        // AT or ID field is mandatory, throw an error
        throw std::invalid_argument("AT and ID field is missing in VCF at position " + std::to_string(raw_record->pos + 1) + "\n\tstoat vcf requires the AT or ID fields containing graph walks from each allele. This can be obtained using vg call or pangenie");
    }

    // extract genotypes GT
    int ngt = 0;
    int32_t *gt = nullptr;
    ngt = bcf_get_genotypes(hdr, raw_record, &gt, &ngt);

    if (ngt <= 0 || gt == nullptr) {
        throw std::invalid_argument("GT field is missing in VCF at position " + std::to_string(raw_record->pos + 1));
    }

    const size_t sample_count = sample_names.size();
    if (sample_count == 0 || static_cast<size_t>(ngt) % sample_count != 0) {
        throw std::invalid_argument("GT field has an invalid number of genotype slots at position " +
            std::to_string(raw_record->pos + 1));
    }

    // Number of GT slots actually stored per sample in this record.
    const size_t gt_ploidy = ngt / sample_count;

    // Warn if the GT field has a different ploidy than expected for this chromosome.
    if (gt_ploidy != ploidy) {
        #pragma omp critical(vcf_parser_log)
        {
            stoat::LOG_WARN("GT field has " + std::to_string(gt_ploidy) +
                " slots per sample, but expected ploidy is " + std::to_string(ploidy) +
                (gt_ploidy > ploidy ? "; extra slots will be ignored" : "; remaining slots will be treated as missing") +
                " at position " + std::to_string(raw_record->pos + 1),
                "gt_field_ploidy_mismatch");
        }
    }

    // Make the actual vector of genotypes
    // If we want to untangle the snarls, then check that the parent snarl actually was genotyped as having this child snarl
    std::vector<int> record_genotypes;
    record_genotypes.reserve(hap_count);

    for (size_t i = 0; i < hap_count; ++i) {
        const size_t sample = i / ploidy;
        const size_t haplotype = i % ploidy;

        // This haplotype does not exist in the GT field for this sample.
        // Fill it with missing genotype.
        if (haplotype >= gt_ploidy) {
            record_genotypes.emplace_back(-1);
            continue;
        }

        const size_t gt_index = sample * gt_ploidy + haplotype;
        const int32_t encoded_gt = gt[gt_index];

        // Missing values
        if (encoded_gt == bcf_int32_vector_end || bcf_gt_is_missing(encoded_gt)) {
            record_genotypes.emplace_back(-1);
            continue;
        }

        int genotype = bcf_gt_allele(encoded_gt);

        // Invalid allele index.
        if (genotype < 0 || genotype >= static_cast<int>(paths.size())) {
            #pragma omp critical(vcf_parser_log)
            {
                stoat::LOG_WARN("VCF variant " + snarl_id + " at " + chr + ":" +
                    std::to_string(raw_record->pos + 1) + " has invalid genotype of " +
                    std::to_string(genotype), "bad_vcf_gt");
            }

            record_genotypes.emplace_back(-1);
            continue;
        }

        if (!resolve_nested_calls || level == 0 || does_sample_have_snarl(i, snarl_id)) {
            record_genotypes.emplace_back(genotype);
        } else {
            record_genotypes.emplace_back(-1);
        }
    }
    free(gt);

#ifdef DEBUG_VCF_PARSER
    std::cerr << " broke out of loop with " << read_status << " At chr " << bcf_hdr_id2name(hdr, rec->rid) << std::endl;
#endif

    return vcf_info_t{level, std::move(record_genotypes), std::move(paths)};
}

std::vector<stoat::node_traversal_t> VCFParser::parse_pangenie_id(std::string& allele_id_str) {
    // split by comma and save as a vector of edge lists [vector vector stoat::edge_t]
    // The id field is a comma separated list of alleles, and each allele is a :-separated list of ids, and each id is a --separated list of fields
    // Each id is formatted: [ref_path]-[offset]-[variant_type]-[path]-[size]:
    // from: ">123>213<234", ">123<234", ">123<234<345"
    // to: [[edge_t(123, 213),stoat::edge_t(213, 234)], [...]]
    std::stringstream id_str_ss(allele_id_str);
    std::string id_str;
    
    // Make a path for this allele. Since the allele may have multiple parts, separate each part by a fake >0 node
    std::vector<stoat::node_traversal_t> allele_path;
    bool first_part = true;
    while (std::getline(id_str_ss, id_str, ':')) {
        if (first_part) {
            first_part=false;
        } else {
            // If we've already seen an id for this allele, add a fake 
            allele_path.emplace_back(0, false);
        }
        // The path should be the fourth field in the id string
        std::stringstream id_field_str_ss(id_str);
        std::string path_str;
        std::getline(id_field_str_ss, path_str, '-');
        std::getline(id_field_str_ss, path_str, '-');
        std::getline(id_field_str_ss, path_str, '-');
        std::getline(id_field_str_ss, path_str, '-');
    
        // Get this part's traversal and add it to allele_path
        std::vector<stoat::node_traversal_t> path_as_nodes = string_to_path_node_traversal(path_str);
        allele_path.insert(allele_path.end(), std::make_move_iterator(path_as_nodes.begin()), std::make_move_iterator(path_as_nodes.end()));
        
    } //end looping through ids for one allele
    return allele_path;
}

void VCFParser::skip_to_next_chromosome(const std::string& chr) {
#ifdef DEBUG_VCF_PARSER
        std::cerr << "Skip through chr " << chr << std::endl;
#endif

    // Since we've already read the first line of this chunk, do a do-while loop and read the next at the end.
    // At the end of this loop, we'll be looking at the first line that is not this chromosome
    do {
        //TODO: I think this is unnecessary
        //bcf_unpack(rec, BCF_UN_STR);

        read_status = bcf_read(ptr_vcf, hdr, rec);
        if (resolve_nested_calls) {
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

        if (sample_names.empty() || static_cast<size_t>(ngt) % sample_names.size() != 0) {
            throw std::invalid_argument("GT field has an invalid number of genotype slots at position " +
                std::to_string(rec_genotypes->pos + 1));
        }

        const size_t gt_ploidy = static_cast<size_t>(ngt) / sample_names.size();
        std::vector<std::vector<stoat::node_traversal_t>> allele_paths;

        // extract AT or ID field from INFO
        char *at = nullptr;
        int nat = 0;
        nat = bcf_get_info_string(hdr_genotypes, rec_genotypes, "AT", &at, &nat);

        if (nat > 0 && at) {
            std::string at_str(at); // convert to C++ std::string
            free(at);
            
            // split by comma and save as a vector of edge lists [vector vector stoat::edge_t]
            // from: ">123>213<234", ">123<234", ">123<234<345"
            // to: [[edge_t(123, 213),stoat::edge_t(213, 234)], [...]]
            std::stringstream at_ss(at_str);
            std::string item;
            while (std::getline(at_ss, item, ',')) {
                std::vector<stoat::node_traversal_t> path_as_node_traversal = string_to_path_node_traversal(item);
                allele_paths.push_back(std::move(path_as_node_traversal));
            }
        } else {
            // AT field is mandatory, throw an error
            throw std::invalid_argument("AT fields are missing in VCF at position " + std::to_string(rec_genotypes->pos + 1) + "\n\tPangenie VCFs cannot be used with the --resolve-vcf option");
        }

        // Iterate only the configured haplotypes, ignoring any extra GT slots.
        for (int sample_num = 0; sample_num < rec_genotypes->n_sample; ++sample_num){
            for (int hap_num = 0; hap_num < ploidy; ++hap_num){
                // allele hap_num of that sample
                const size_t sample_hap_index = static_cast<size_t>(sample_num) * ploidy + hap_num;

                // HTSlib stores a rectangular GT array. Shorter genotypes are
                // padded with bcf_int32_vector_end.
                if (static_cast<size_t>(hap_num) >= gt_ploidy) {
                    continue;
                }

                const size_t gt_index = static_cast<size_t>(sample_num) * gt_ploidy + hap_num;
                const int32_t encoded_gt = gt[gt_index];

                if (encoded_gt == bcf_int32_vector_end || bcf_gt_is_missing(encoded_gt)) {
                    continue;
                }

                int idx_path_allele = bcf_gt_allele(encoded_gt);

                if (idx_path_allele > (int)-1 && idx_path_allele < (int)allele_paths.size() ) { // If this has acceptable genotypes
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
    
                } else if (sample_hap_index % ploidy == 1 &&
                           gt_index > 0 &&
                           gt[gt_index - 1] != bcf_int32_vector_end &&
                           !bcf_gt_is_missing(gt[gt_index - 1]) &&
                           bcf_gt_allele(gt[gt_index - 1]) >= 0) {
                    throw std::invalid_argument("VCF variant has undefined genotype of " + std::to_string(idx_path_allele));
                }
            }
        }
        free(gt);

    } while ((bcf_read(ptr_vcf_genotypes, hdr_genotypes, rec_genotypes) >= 0) && (chr == bcf_hdr_id2name(hdr_genotypes, rec_genotypes->rid)));
    
}

bool VCFParser::does_sample_have_snarl(size_t sample_hap_index, const std::string& snarl_id) {

    if (snarl_id == ".") {
        throw std::runtime_error("error: Trying to untangle a pangenie vcf");
        return true;
    }

    // This should be a vector of two node_traversal_t's of the snarl bounds, first one pointing in, second one pointing out
    std::vector<stoat::node_traversal_t> snarl_bounds = string_to_path_node_traversal(snarl_id);
    if (snarl_in_to_out.count(snarl_bounds[0])) {
        return genotypes.at(genotype_index(sample_hap_index, snarl_bounds[0]));
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

    if (resolve_nested_calls) {
        bcf_destroy(rec_bounds);
        bcf_hdr_destroy(hdr_bounds);
        bcf_close(ptr_vcf_bounds);

        bcf_destroy(rec_genotypes);
        bcf_hdr_destroy(hdr_genotypes);
        bcf_close(ptr_vcf_genotypes);
    }
}


}//end namespace
