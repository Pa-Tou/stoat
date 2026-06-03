#include <filesystem>
#include "log.hpp"
#include "arg_parser.hpp"

namespace fs = std::filesystem;

namespace stoat_vcf {

std::unordered_set<std::string> parse_chromosome_reference(const std::string& file_path) {
    std::unordered_set<std::string> reference;
    if (!std::filesystem::exists(file_path)) {
        stoat::LOG_WARN("given reference file " + file_path + " does not exist. Defaulting to using any reference- or generic-sense paths as references", 0);
        return reference;
    }
    std::ifstream file(file_path);
    std::string line;

    while (getline(file, line)) {
        reference.insert(line);
    }

    file.close();

    if (reference.size() == 0) {
        stoat::LOG_WARN("given reference file " + file_path + " is empty. Defaulting to using any reference- or generic-sense paths as references", 0);
    }
    return reference;
}

std::string methods_stats_prediction(const std::string& file_path, const bool& covariate) {
    std::ifstream file(file_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + file_path);
    }

    std::string line, sample_name, phenoStr;
    std::getline(file, line); // header

    std::istringstream header_stream(line);
    std::string col1, col2;
    header_stream >> col1 >> col2;

    // If not SAMPLE PHENO format, assume gene expression matrix
    if (col1 != "SAMPLE" || col2 != "PHENO") {
        return "linreg"; // assuming eqtl format
    }

    bool only_0_1 = true;

    while (std::getline(file, line)) {
        std::istringstream iss(line);

        if (!(iss >> sample_name >> phenoStr)) {
            throw std::invalid_argument("Malformed line: " + line);
        }

        double value;
        try {
            value = std::stod(phenoStr);
        } catch (...) {
            throw std::invalid_argument("Phenotype is not numeric: " + phenoStr);
        }

        if (!(value == 0.0 || value == 1.0)) {
            only_0_1 = false;
        }
    }

    file.close();

    if (only_0_1 && covariate) {
        return "logreg"; // binary + covariate
    } 
    else if (only_0_1) {
        return "chi2"; // binary
    } 
    else {
        return "linreg"; // quantitative
    }

}

std::unique_ptr<stoat::BinaryPhenotypeTable> parse_binary_pheno_table(const std::string& file_path, std::unordered_map<std::string, size_t>& sample_to_index) {

    // fill this map first
    std::unordered_map<std::string, bool> binary_pheno;
    // should we update the sample to index map? yes if it's empty at the start
    bool update_sample_to_index = sample_to_index.empty();

    // prepare to read the file
    std::ifstream file(file_path);
    std::string line;

    // read the header and check for expected column names (SAMPLE, then PHENO)
    std::getline(file, line);
    std::istringstream header_stream(line);
    std::string sample_name, phenoStr;
    header_stream >> sample_name >> phenoStr;
    
    if (sample_name != "SAMPLE" || phenoStr != "PHENO") {
        throw std::invalid_argument("Invalid header. Must be SAMPLE then PHENO, got: " + line);
    }
    
    // read each line and tally the number of cases and controls (for the log)
    int count_controls = 0;
    int count_cases = 0;
    size_t line_number = 1;
    size_t sample_idx = 0;
    while (std::getline(file, line)) {
        ++line_number;
        std::istringstream iss(line);
        if (!(iss >> sample_name >> phenoStr)) {
            throw std::invalid_argument("Malformed line " + std::to_string(line_number) + ": " + line);
        }

        // make sure the phenotype is an integer
        int pheno;
        try {
            pheno = std::stoi(phenoStr);
        } catch (...) {
            throw std::invalid_argument("Bad phenotype type at line " + std::to_string(line_number) + ": " + phenoStr);
        }

        // make sure phenotype is 0 or 1
        if (pheno == 0) {
            ++count_controls;
        } else if (pheno == 1) {
            ++count_cases;
        } else {
            throw std::invalid_argument("Binary phenotype must be 0 or 1 at line " + std::to_string(line_number) + ", got: " + phenoStr);
        }

        // insert and detect duplicate
        auto [it, inserted] = binary_pheno.emplace(sample_name, pheno == 1);
        if (!inserted) {
            throw std::invalid_argument("Duplicate sample in phenotype file: " + sample_name);
        }

        if (update_sample_to_index) {
            sample_to_index[sample_name] = sample_idx++;
        }
    }

    file.close();

    stoat::LOG_INFO("Binary phenotypes found: " + std::to_string(count_controls + count_cases)
        + " (Control: " + std::to_string(count_controls)
        + ", Case: " + std::to_string(count_cases) + ")");

    // Prepare the Table object to fill and output
    // ideally we could fill it when reading each line
   std::unique_ptr< stoat::BinaryPhenotypeTable> output_table ( new stoat::BinaryPhenotypeTable(sample_to_index));

    size_t unique_to_pheno = 0;
    size_t missing_in_pheno = 0;

    size_t used_cases = 0;
    size_t used_controls = 0;

    // phenotype -> genotype match
    // Prepare the Table object to fill and output
    // ideally we could fill it when reading each line
    for (const auto& [sample, pheno] : binary_pheno) {
        if (output_table->has_sample(sample)) {
            output_table->set_value_for_sample(sample, pheno);
            if (pheno) {
                ++used_cases;
            } else {
                ++used_controls;
            }
        } else {
            ++unique_to_pheno;
        }
    }

    // samples in the index but not in the phenotype file are marked as missing
    for (const auto& [sample, idx] : sample_to_index) {
        if (binary_pheno.find(sample) == binary_pheno.end()) {
            ++missing_in_pheno;
            output_table->add_missing_sample_index(idx);
        }
    }

    // print the warnings
    if (unique_to_pheno > 0) {
        stoat::LOG_WARN(std::to_string(unique_to_pheno) + " phenotype samples were not found in the sample set", 0);
    }

    if (missing_in_pheno > 0) {
        stoat::LOG_WARN(std::to_string(missing_in_pheno) + " samples have no phenotype", 0);
    }

    // print how many samples will be used
    stoat::LOG_INFO("Binary phenotypes used for GWAS: " + std::to_string(used_cases + used_controls) + 
        " (Control: " + std::to_string(used_controls) + ", Case: " + std::to_string(used_cases) + ")");

    return output_table;
}

std::unique_ptr<stoat::QuantitativePhenotypeTable> parse_quantitative_pheno_table(const std::string& file_path, std::unordered_map<std::string, size_t>& sample_to_index) {

    // fill this map first
    std::unordered_map<std::string, double> quantitative_pheno;

    // should we update the sample to index map? yes if it's empty at the start
    bool update_sample_to_index = sample_to_index.empty();

    // prepare to read the file
    std::ifstream file(file_path);
    std::string line;

    // read the header and check for expected column names (SAMPLE, then PHENO)
    std::getline(file, line);
    std::istringstream header_stream(line);
    std::string sample_name, phenoStr;
    header_stream >> sample_name >> phenoStr;

    if (sample_name != "SAMPLE" || phenoStr != "PHENO") {
        throw std::invalid_argument("Invalid header. Expected: SAMPLE PHENO, got: " + line);
    }
    
    // read each line and tally the number of cases and controls (for the log)
    size_t samp_with_pheno = 0;
    size_t line_number = 1;
    size_t sample_idx = 0;

    while (std::getline(file, line)) {
        ++line_number;

        std::istringstream iss(line);

        if (!(iss >> sample_name >> phenoStr)) {
            throw std::invalid_argument("Malformed line " + std::to_string(line_number) + ": " + line);
        }

        // make sure the phenotype is a double
        double pheno;
        try {
            // add the sample and phenotype to the temporary map
            pheno = std::stod(phenoStr);
        } catch (...) {
            throw std::invalid_argument("Bad phenotype type: " + phenoStr);
        }

        // insert and detect duplicate
        auto [it, inserted] = quantitative_pheno.emplace(sample_name, pheno);
        if (!inserted) {
            throw std::invalid_argument("Duplicate sample in phenotype file: " + sample_name);
        }

        samp_with_pheno++;
        if (update_sample_to_index) {
            sample_to_index.emplace(sample_name, sample_idx++);
        }
    }

    file.close();
    stoat::LOG_INFO("Quantitative phenotypes found in phenotype file: " + 
        std::to_string(samp_with_pheno) + " samples");

    // Prepare the Table object to fill and output
    // ideally we could fill it when reading each line
    std::unique_ptr<stoat::QuantitativePhenotypeTable> output_table ( new stoat::QuantitativePhenotypeTable(sample_to_index));
    
    size_t unique_to_pheno = 0;
    size_t missing_in_pheno = 0;

    for (const auto samp_pheno: quantitative_pheno){
        // add the sample and phenotype to the Table
        if (output_table->has_sample(samp_pheno.first)) {
            output_table->set_value_for_sample(samp_pheno.first, samp_pheno.second);
        } else {
            ++unique_to_pheno;
        }
    }

    // samples in the index but not in the phenotype file are marked as missing
    for (const auto& [sample, idx] : sample_to_index) {
        if (quantitative_pheno.find(sample) == quantitative_pheno.end()) {
            ++missing_in_pheno;
            output_table->add_missing_sample_index(idx);
        }
    }

    // Warnings
    if (unique_to_pheno > 0) {
        stoat::LOG_WARN(std::to_string(unique_to_pheno) + " phenotype samples were not found in the sample set", 0);
    }

    if (missing_in_pheno > 0) {
        stoat::LOG_WARN(std::to_string(missing_in_pheno) + " samples have no phenotype", 0);
    }

    // Final Pheno GWAS counts
    stoat::LOG_INFO("Quantitative phenotypes used for GWAS: " + 
        std::to_string(samp_with_pheno) + " samples");

    return output_table;
}

std::unique_ptr<stoat::GeneExpressionTable> parse_gene_expression_table(const std::string& gene_expression_path, const std::string& gene_position_path, std::unordered_map<std::string, size_t>& sample_to_index, std::unordered_map<std::string, size_t>& gene_to_index) {

    // should we update the sample to index map? yes if it's empty at the start
    bool update_sample_to_index = sample_to_index.empty();

    // currently we expect an empty gene_to_index map because we might have gotten samples from
    // another file but we won't have gotten genes from another file
    assert(gene_to_index.empty());

    // prepare to read and parse the gene expression file
    std::ifstream file(gene_expression_path);
    std::string line;
    std::string line_value;
    std::getline(file, line);

    // read the header and define the map sample-> index
    std::stringstream ss_header(line);
    std::getline(ss_header, line_value, '\t'); // Skip the first column (gene name)

    size_t samp_index = 0;
    std::vector<std::string> sample_names;
    std::unordered_set<std::string> seen_samples;

    while (std::getline(ss_header, line_value, '\t')) {

        if (line_value.empty()) {
            throw std::invalid_argument("Empty sample name in header");
        }

        // check duplicate sample
        if (!seen_samples.insert(line_value).second) {
            throw std::invalid_argument("Duplicate sample name in header: " + line_value);
        }

        // preserve column order
        sample_names.emplace_back(line_value);

        // optionally build sample->index map
        if (update_sample_to_index) {
            sample_to_index.emplace(line_value, samp_index++);
        }
    }

    // fill this map: sample name -> expression vector for all genes
    std::unordered_map<std::string, std::vector<double>> ge_map;
    size_t gene_idx = 0;
    size_t n_samples = sample_names.size();
    size_t line_number = 1;

    // read gene expressions for each gene
    while (std::getline(file, line)) {
        ++line_number;
        std::stringstream ss(line);

        // parse the gene name
        std::string gene_name;
        std::getline(ss, gene_name, '\t');

        if (gene_name.empty()) {
            throw std::invalid_argument("Missing gene name in expression file at line " + std::to_string(line_number));
        }

        size_t col = 0;
        while (std::getline(ss, line_value, '\t')) {
            if (col >= n_samples) {
                throw std::invalid_argument(
                    "Too many columns for gene " + gene_name + " at line " + std::to_string(line_number));
            }

            try {
                ge_map[sample_names[col]].push_back(std::stod(line_value));
            } catch (...) {
                throw std::invalid_argument(
                    "Invalid expression value for gene " + gene_name +
                    " at line " + std::to_string(line_number) + ": " + line_value);
            }

            ++col;
        }

        if (col != n_samples) {
            throw std::invalid_argument("Incorrect number of columns for gene " +
                                        gene_name + "at line " +
                                        std::to_string(line_number));
        }

        auto [it, inserted] = gene_to_index.emplace(gene_name, gene_idx++);
        if (!inserted) {
            throw std::invalid_argument("Duplicate gene name in expression file: " +
                                        gene_name + "at line " +
                                        std::to_string(line_number));
        }
    }

    file.close();

    stoat::LOG_INFO("Gene number found in expression file: " + std::to_string(gene_idx));

    // Build output table
    std::unique_ptr<stoat::GeneExpressionTable> output_table ( new stoat::GeneExpressionTable(sample_to_index, gene_to_index));

    size_t unique_to_expression = 0;
    size_t missing_in_expression = 0;
    size_t used_samples = 0;

    for (const auto& [sample, expr_vec] : ge_map) {
        if (output_table->has_sample(sample)) {
            ++used_samples;
            for (const auto& [gene, idx] : gene_to_index) {
                output_table->set_value_for_sample_and_feature(sample, gene, expr_vec[idx]);
            }
        } else {
            ++unique_to_expression;
        }
    }

    // samples in the index but not in the expression file are marked as missing
    for (const auto& [sample, idx] : sample_to_index) {
        if (ge_map.find(sample) == ge_map.end()) {
            ++missing_in_expression;
            output_table->add_missing_sample_index(idx);
        }
    }

    // Warnings
    if (unique_to_expression > 0) {
        stoat::LOG_WARN(std::to_string(unique_to_expression) +
            " samples with expression not found in current sample set", 0);
    }

    if (missing_in_expression > 0) {
        stoat::LOG_WARN(std::to_string(missing_in_expression) +
            " samples missing expression values", 0);
    }

    stoat::LOG_INFO("Samples with expression used in analysis: " +
        std::to_string(used_samples));

    // Load gene positions
    output_table->read_gene_positions_from_file(gene_position_path);
    return output_table;
}

std::unique_ptr<stoat::CovariateTable> parse_covariate_table(
    const std::string& file_path,
    std::unordered_map<std::string, size_t>& sample_to_index,
    std::unordered_map<std::string, size_t>& covar_to_index) {

    std::unordered_map<std::string, std::vector<double>> covariate_map;
    bool update_sample_to_index = sample_to_index.empty();

    std::ifstream file(file_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open covariate file: " + file_path);
    }

    std::string line;

    // Parse header
    if (!std::getline(file, line)) {
        throw std::invalid_argument("Empty covariate file: " + file_path);
    }

    std::vector<std::string> headers;
    std::istringstream header_stream(line);
    std::string head_val;

    while (header_stream >> head_val) {
        headers.push_back(head_val);
    }
    //TODO: Check that everything in covar_to_index is in headers

    // check for a SAMPLE column
    auto samp_head_it = std::find(headers.begin(), headers.end(), "SAMPLE");
    if (samp_head_it == headers.end()) {
        throw std::invalid_argument("Header must include 'SAMPLE' column");
    }

    // If we didn't already specify which covariates, use all of them
    size_t samp_head_i =  samp_head_it - headers.begin();
    size_t covar_count = 0;
    if (covar_to_index.empty()) {
        for (size_t i = 0 ; i < headers.size() ; i++) {
            if (i != samp_head_i) {
                covar_to_index[headers.at(i)] = covar_count++;
            }
        }
    }

    size_t samp_head_idx = std::distance(headers.begin(), samp_head_it);

    // save the position of each covariate column
    std::unordered_map<std::string, size_t> col_index;
    for (size_t i = 0; i < headers.size(); ++i) {
        col_index.emplace(headers[i], i);
    }

    // check specified covariate columns
    // also save the maximum covar index that will define the size of the vectors in covariate_map
    size_t max_covar_idx = 0;
    for (const auto& [covar_name, covar_idx] : covar_to_index) {
        if (col_index.find(covar_name) == col_index.end()) {
            throw std::invalid_argument(
                "Covariate column '" + covar_name + "' not found in file");
        }
        max_covar_idx = std::max(max_covar_idx, covar_idx);
    }

    // parse rows
    std::unordered_set<std::string> seen_samples;
    size_t sample_idx = 0;
    size_t line_number = 1;
    while (std::getline(file, line)) {
        ++line_number;
        std::istringstream line_stream(line);
        std::vector<std::string> line_vec;
        std::string line_val;
        while (line_stream >> line_val) {
            line_vec.push_back(line_val);
        }

        if (line_vec.size() <= samp_head_idx) {
            throw std::invalid_argument("Malformed line in covariate file at line " +
                                        std::to_string(line_number) + ": " + line);
        }

        // extract the specified covariables
        std::string samp_name = line_vec[samp_head_idx];

        if (!seen_samples.insert(samp_name).second) {
            throw std::invalid_argument("Duplicate sample in covariate file at line " +
                                        std::to_string(line_number) + ": " + samp_name);
        }

        // prepare a vector with enough elements to fill the value using covar_to_index
        std::vector<double> samp_covars(max_covar_idx + 1);
        try {
            for (const auto& [covar_name, covar_idx] : covar_to_index) {
                size_t col = col_index[covar_name];
                if (col >= line_vec.size()) {
                    throw std::invalid_argument(
                        "Missing value for covariate '" + covar_name +
                        "' in sample " + samp_name + "at line " +
                        std::to_string(line_number));
                }
                samp_covars[covar_idx] = std::stod(line_vec[col]);
            }
        } catch (...) {
            throw std::invalid_argument(
                "Sample " + samp_name + " has a non-numeric covariate value at line " + std::to_string(line_number));
        }

        covariate_map.emplace(samp_name, std::move(samp_covars));
        if (update_sample_to_index) {
            sample_to_index.emplace(samp_name, sample_idx++);
        }
    }

    // close the file connection
    file.close();

    // Build output table
    std::unique_ptr<stoat::CovariateTable> output_table (new stoat::CovariateTable(sample_to_index, covar_to_index));
    size_t unique_to_covar = 0;
    size_t genotype_missing_covars = 0;
    size_t used_samples = 0;

    for (const auto& [sample, covars] : covariate_map) {
        if (output_table->has_sample(sample)) {
            ++used_samples;
            for (const auto& [covar_name, covar_idx] : covar_to_index) {
                output_table->set_value_for_sample_and_feature(
                    sample,
                    covar_name,
                    covars[covar_idx]
                );
            }
        } else {
            ++unique_to_covar;
        }
    }

    // samples in the index but not in the covariable file are marked as missing
    for (const auto& [sample, idx] : sample_to_index) {
        if (covariate_map.find(sample) == covariate_map.end()) {
            ++genotype_missing_covars;
            output_table->add_missing_sample_index(idx);
        }
    }

    // Warnings
    if (unique_to_covar > 0) {
        stoat::LOG_WARN(std::to_string(unique_to_covar) +
            " samples with covariates not found in sample set", 0);
    }

    if (genotype_missing_covars > 0) {
        stoat::LOG_WARN(std::to_string(genotype_missing_covars) +
            " samples missing covariates", 0);
    }

    stoat::LOG_INFO("Samples with covariates used in analysis: " +
        std::to_string(used_samples));

    return output_table;
}

// Function to open a VCF file and return pointers to the file, header, and record
std::tuple<htsFile*, bcf_hdr_t*, bcf1_t*> parse_vcf(const std::string& vcf_path) {
    // Open the VCF file
    htsFile *ptr_vcf = bcf_open(vcf_path.c_str(), "r");

    // Read the VCF header
    bcf_hdr_t *hdr = bcf_hdr_read(ptr_vcf);
    if (!hdr) {
        bcf_close(ptr_vcf);
        throw std::invalid_argument("Could not read VCF header");
    }

    // Initialize a record
    bcf1_t *rec = bcf_init();
    if (!rec) {
        bcf_hdr_destroy(hdr);
        bcf_close(ptr_vcf);
        throw std::invalid_argument("Failed to allocate memory for VCF record");
    }

    // Return the three initialized pointers
    return std::make_tuple(ptr_vcf, hdr, rec);
}

std::tuple<std::vector<std::string>, htsFile*, bcf_hdr_t*, bcf1_t*> parseHeader(const std::string& vcf_path) {
    auto [ptr_vcf, hdr, rec] = parse_vcf(vcf_path);

    std::vector<std::string> list_samples;
    // Get the samples names
    for (int i = 0; i < bcf_hdr_nsamples(hdr); i++) {
        list_samples.push_back(bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i));
    }
        
    return std::make_tuple(list_samples, ptr_vcf, hdr, rec);
}

// Explicit instantiation for specific types
template void check_match_samples<bool>(const std::unordered_map<std::string, bool>&, const std::vector<std::string>&);
template void check_match_samples<double>(const std::unordered_map<std::string, double>&, const std::vector<std::string>&);
template void check_match_samples<std::vector<double>>(const std::unordered_map<std::string, std::vector<double>>&, const std::vector<std::string>&);
template void check_match_samples<std::tuple<std::string, int, int>>(const std::unordered_map<std::string, std::tuple<std::string, int, int>>&, const std::vector<std::string>&);

template <typename T>
void check_match_samples(const std::unordered_map<std::string, T>& map, const std::vector<std::string>& keys) {
    for (const auto& key : keys) {
        if (map.find(key) == map.end()) {
            throw std::invalid_argument("Sample '" + key + "' not found in the phenotype file");
        }
    }
    if (map.size() != keys.size()) {
        stoat::LOG_WARN("Number of samples found in VCF (" + std::to_string(keys.size()) + ") does not match the number of samples in the phenotype file (" + std::to_string(map.size()) + ").", 0);
    }
}

void check_methods(const std::string& method_name) {
    if (method_name != "chi2" && method_name != "linreg" && method_name != "logreg" && method_name != "exact") {
        throw std::invalid_argument("Invalid methods: " + method_name + ". Available methods are: chi2, linreg, logreg, exact.");
    }
}

void check_file(const std::string& file_path) {
    
    std::string line;

    // Check if file is a file
    if (!fs::is_regular_file(file_path)) {
        throw std::invalid_argument("File " + file_path + " does not exist.");
    }

    // Check if file can be open
    std::ifstream file(file_path);
    if (!file.is_open()) {
        throw std::invalid_argument("Unable to open the file " + file_path);
    }

    // Check if file can be read and not empty file
    if (!std::getline(file, line)) {
        throw std::invalid_argument("File " + file_path + "is empty or failed to read header.");
    }

    file.close();
}

} //end stoat namespace
