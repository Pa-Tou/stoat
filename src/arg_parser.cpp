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

std::vector<bool> parse_binary_pheno(
    const std::string& file_path,
    std::vector<std::string>& list_samples) {
    bool fill_in_samples = false;
    if (list_samples.size() == 0) {
        fill_in_samples = true;
    }
    
    std::unordered_map<std::string, bool> binary_pheno;
    std::ifstream file(file_path);
    std::string line;

    // Header is assumed already read and validated externally
    // JEAN but it's not, is it? we're doing this just below
    std::getline(file, line);
    std::istringstream header_stream(line);
    std::string sample_name, phenoStr;
    header_stream >> sample_name >> phenoStr;
    if (sample_name != "SAMPLE" || phenoStr != "PHENO") {
        throw std::invalid_argument("Invalid header. Must be SAMPLE then PHENO, got " + line);
    }
    
    // --- Read and process data ---
    int count_controls = 0;
    int count_cases = 0;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        if (!(iss >> sample_name >> phenoStr)) {
            throw std::invalid_argument("Malformed line: " + line);
        }

        int pheno;
        try {
            pheno = std::stoi(phenoStr);
        } catch (...) {
            throw std::invalid_argument("Bad phenotype type: " + phenoStr);
        }

        if (pheno == 0) {
            ++count_controls;
            binary_pheno[sample_name] = false;
        } else if (pheno == 1) {
            ++count_cases;
            binary_pheno[sample_name] = true;
        } else {
            throw std::invalid_argument("Binary phenotype must be 0 or 1, got: " + std::to_string(pheno));
        }
        if (fill_in_samples) {
            list_samples.emplace_back(sample_name);
        }
    }

    stoat::LOG_INFO("Binary phenotypes found: " + std::to_string(count_controls + count_cases)
        + " (Control: " + std::to_string(count_controls)
        + ", Case: " + std::to_string(count_cases) + ")");

    file.close();

    if (!fill_in_samples) {
        // If we were given samples, make sure that they check the phenotype file
        check_match_samples(binary_pheno, list_samples);
    }

    std::vector<bool> vector_binary_pheno;
    vector_binary_pheno.reserve(list_samples.size());

    for (const auto& sample : list_samples) {
        auto it = binary_pheno.find(sample);
        if (it != binary_pheno.end()) {
            vector_binary_pheno.push_back(it->second);
        }
    }

    return vector_binary_pheno;
}
    
stoat::BinaryPhenotypeTable* parse_binary_pheno_table(const std::string& file_path, std::unordered_map<std::string, size_t>& sample_to_index) {
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
        throw std::invalid_argument("Invalid header. Must be SAMPLE then PHENO, got " + line);
    }
    
    // read each line and tally the number of cases and controls (for the log)
    int count_controls = 0;
    int count_cases = 0;
    size_t sample_idx = 0;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        if (!(iss >> sample_name >> phenoStr)) {
            throw std::invalid_argument("Malformed line: " + line);
        }
        int pheno;
        // make sure the phenotype is an integer
        try {
            pheno = std::stoi(phenoStr);
        } catch (...) {
            throw std::invalid_argument("Bad phenotype type: " + phenoStr);
        }
        // make sure phenotype is 0 or 1
        if (pheno == 0) {
            ++count_controls;
        } else if (pheno == 1) {
            ++count_cases;
        } else {
            throw std::invalid_argument("Binary phenotype must be 0 or 1, got: " + std::to_string(pheno));
        }
        // add the sample and phenotype to the temporary map
        binary_pheno[sample_name] = pheno == 1;
        if (update_sample_to_index) {
            sample_to_index[sample_name] = sample_idx++;
        }
    }

    stoat::LOG_INFO("Binary phenotypes found: " + std::to_string(count_controls + count_cases)
        + " (Control: " + std::to_string(count_controls)
        + ", Case: " + std::to_string(count_cases) + ")");

    file.close();

    // Prepare the Table object to fill and output
    // ideally we could fill it when reading each line
    stoat::BinaryPhenotypeTable* output_table = new stoat::BinaryPhenotypeTable(sample_to_index);
    for (const auto samp_pheno: binary_pheno){
        // add the sample and phenotype to the Table
        if (output_table->has_sample(samp_pheno.first)) {
            output_table->set_value_for_sample(samp_pheno.first, samp_pheno.second);
        }
    }

    return (output_table);
}

std::string methods_stats_prediction(const std::string& file_path, const bool& covariate, const bool& eqtl) {
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

    std::set<double> unique_values;
    bool has_non_integer = false;

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
        unique_values.insert(value);
        if (value != static_cast<int>(value)) {
            has_non_integer = true;
        }
    }

    file.close();

    if (eqtl || has_non_integer) {
        return "linreg"; // quantitative
    } else if (unique_values.size() == 2 && covariate) {
        return "logreg"; // binary + covar
    } else {
        return "chi2"; // binary no covar
    }
}

stoat::QuantitativePhenotypeTable* parse_quantitative_pheno_table(const std::string& file_path, std::unordered_map<std::string, size_t>& sample_to_index) {
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
        throw std::invalid_argument("Invalid header. Must be SAMPLE then PHENO, got " + line);
    }
    
    // read each line and tally the number of cases and controls (for the log)
    int samp_with_pheno = 0;
    size_t sample_idx = 0;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        if (!(iss >> sample_name >> phenoStr)) {
            throw std::invalid_argument("Malformed line: " + line);
        }
        // make sure the phenotype is a double
        try {
            // add the sample and phenotype to the temporary map
            quantitative_pheno[sample_name] = std::stod(phenoStr);
        } catch (...) {
            throw std::invalid_argument("Bad phenotype type: " + phenoStr);
        }
        samp_with_pheno++;
        if (update_sample_to_index) {
            sample_to_index[sample_name] = sample_idx++;
        }
    }
    file.close();
    stoat::LOG_INFO("Quantitative phenotypes found for " + std::to_string(samp_with_pheno) + " samples");

    // Prepare the Table object to fill and output
    // ideally we could fill it when reading each line
    stoat::QuantitativePhenotypeTable* output_table = new stoat::QuantitativePhenotypeTable(sample_to_index);
    for (const auto samp_pheno: quantitative_pheno){
        // add the sample and phenotype to the Table
        if (output_table->has_sample(samp_pheno.first)) {
            output_table->set_value_for_sample(samp_pheno.first, samp_pheno.second);
        }
    }

    return (output_table);
}

stoat::GeneExpressionTable* parse_gene_expression_table(const std::string& gene_expression_path, const std::string& gene_position_path, std::unordered_map<std::string, size_t>& sample_to_index, std::unordered_map<std::string, size_t>& gene_to_index) {
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
    while (std::getline(ss_header, line_value, '\t')) {
        // always save the order of the samples to know which column corresponds to what
        sample_names.emplace_back(line_value);
        // eventually set the sample->index map
        if (update_sample_to_index) {
            sample_to_index[line_value] = samp_index++;
        }
    }

    // fill this map: sample name -> expression vector for all genes
    std::unordered_map<std::string, std::vector<double>> ge_map;
    // read gene expressions for each gene
    size_t gene_idx = 0;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        // parse the gene name
        std::string gene_name;
        std::getline(ss, gene_name, '\t');
        // append that gene's expression for each sample
        samp_index = 0;
        while (std::getline(ss, line_value, '\t')) {
            try {
                ge_map[sample_names[samp_index]].push_back(std::stod(line_value));
                samp_index++;
            } catch (...) {
                throw std::invalid_argument("Invalid expression value for gene " + gene_name + ": " + line_value);
            }
        }
        // save the gene index
        gene_to_index[gene_name] = gene_idx++;
    }

    // close the file connection
    file.close();

    // Prepare the Table object to fill and output
    // ideally we could fill it when reading each line
    stoat::GeneExpressionTable* output_table = new stoat::GeneExpressionTable(sample_to_index, gene_to_index);
    for (const auto samp_ge: ge_map){
        // add the sample and gene expression to the Table
        if (output_table->has_sample(samp_ge.first)) {
            for (const auto& gene_idx : gene_to_index) {
                output_table->set_value_for_sample_and_feature(samp_ge.first, gene_idx.first, samp_ge.second[gene_idx.second]);
            }
        }
    }

    // add gene positions from the file
    output_table->read_gene_positions_from_file(gene_position_path);
    return (output_table);

}
    
stoat::CovariateTable* parse_covariate_table(const std::string& file_path, std::unordered_map<std::string, size_t>& sample_to_index,
                                             std::unordered_map<std::string, size_t>& covar_to_index) {
    // fill this map first
    std::unordered_map<std::string, std::vector<double>> covariate_map;
    // should we update the sample to index map? yes if it's empty at the start
    bool update_sample_to_index = sample_to_index.empty();
    // currently we need a non empty sample_to_index map defining which covariables to use
    assert(!covar_to_index.empty());

    // prepare to read the file
    std::ifstream file(file_path);
    std::string line;

    // read the header and check for expected column names (at least SAMPLE)
    std::vector<std::string> headers;
    std::getline(file, line);
    std::istringstream header_stream(line);
    std::string head_val;
    while (header_stream >> head_val) {
        headers.push_back(head_val);
    }

    // check for a SAMPLE column
    auto samp_head_it = std::find(headers.begin(), headers.end(), "SAMPLE");
    if (samp_head_it == headers.end()) {
        throw std::invalid_argument("header must include 'SAMPLE' column.\n");
    }
    size_t samp_head_idx = std::distance(headers.begin(), samp_head_it);

    // save the position of each covariate column
    std::unordered_map<std::string, size_t> col_index;
    for (size_t i = 0; i < headers.size(); ++i) {
        col_index[headers[i]] = i;
    }

    // check specified covariate columns
    // also save the maximum covar index that will define the size of the vectors in covariate_map
    size_t max_covar_idx = 0;
    for (const auto& covar_idx : covar_to_index) {
        if (col_index.find(covar_idx.first) == col_index.end()) {
            throw std::invalid_argument("covariate column '" + covar_idx.first + "' not found in file.\n");
        }
        if (max_covar_idx < covar_idx.second) {
            max_covar_idx = covar_idx.second;
        }
    }
    
    // read covariates for each sample
    size_t sample_idx = 0;
    while (std::getline(file, line)) {
        // parse the line into a vector
        std::istringstream line_stream(line);
        std::vector<std::string> line_vec;
        std::string line_val;
        while (line_stream >> line_val) {
            line_vec.push_back(line_val);
        }
        
        if (line_vec.size() <= samp_head_idx) continue; // JEAN isn't that a sign that the file is wrong and we should raise an error?

        // extract the specified covariables
        std::string samp_name = line_vec[samp_head_idx];
        // prepare a vector with enough elements to fill the value using covar_to_index
        std::vector<double> samp_covars(max_covar_idx + 1);
        try {
            for (const auto& covar_idx : covar_to_index) {
                samp_covars[covar_idx.second] = std::stod(line_vec[col_index[covar_idx.first]]);
            }
        } catch (...) {
            throw std::invalid_argument("Individual " + samp_name + " got an non-numeric value\n");
        }
        covariate_map[samp_name] = samp_covars;
        if (update_sample_to_index) {
            sample_to_index[samp_name] = sample_idx++;
        }
    }

    // close the file connection
    file.close();

    // Prepare the Table object to fill and output
    // ideally we could fill it when reading each line
    stoat::CovariateTable* output_table = new stoat::CovariateTable(sample_to_index, covar_to_index);
    for (const auto samp_covar: covariate_map){
        // add the sample and phenotype to the Table
        if (output_table->has_sample(samp_covar.first)) {
            for (const auto& covar_idx : covar_to_index) {
                output_table->set_value_for_sample_and_feature(samp_covar.first, covar_idx.first, samp_covar.second[covar_idx.second]);
            }
        }
    }

    return (output_table);
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
