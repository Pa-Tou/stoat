#include "log.hpp"
#include "arg_parser.hpp"

namespace fs = std::filesystem;

namespace stoat_vcf {

std::unordered_set<std::string> parse_chromosome_reference(const std::string& file_path) {
    std::unordered_set<std::string> reference;
    ifstream file(file_path);
    std::string line;

    while (getline(file, line)) {
        reference.insert(line);
    }

    file.close();
    return reference;
}

std::vector<bool> parse_binary_pheno(
    const std::string& file_path,
    const std::vector<std::string>& list_samples) {
    
    std::unordered_map<std::string, bool> binary_pheno;
    
    std::ifstream file(file_path);
    std::string line;
    int count_controls = 0;
    int count_cases = 0;
    bool firstLine = true;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string fid, iid, phenoStr;

        if (!(iss >> fid >> iid >> phenoStr)) {
            stoat::LOG_FATAL("Malformed line: " + line);
        }

        if (firstLine) {
            firstLine = false;
            // Check that the header contains FID, IID, and PHENO
            if (fid != "FID" || iid != "IID" || phenoStr != "PHENO") {
                stoat::LOG_FATAL("Invalid header: " + line);
            }
            continue;
        }

        int pheno = -1;
        try {
            pheno = std::stoi(phenoStr);
        } catch (...) {
            stoat::LOG_FATAL("Bad phenotype type : " + phenoStr);
        }
        if (pheno == 1) {
            count_controls++;
            binary_pheno[iid] = static_cast<bool>(false);
        } else if (pheno == 2) {
            count_cases++;
            binary_pheno[iid] = static_cast<bool>(true);
        } else {
            stoat::LOG_FATAL("Binary phenotype must be 1 or 2");
        }
    }

    stoat::LOG_INFO("Binary phenotypes founds : " + std::to_string(count_controls+count_cases)
        + " (Control : " + std::to_string(count_controls) 
        + ", Case : " + std::to_string(count_cases) + ")");

    file.close();

    check_match_samples(binary_pheno, list_samples);
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

// Function to parse the phenotype file
std::vector<double> parse_quantitative_pheno(
    const std::string& file_path, 
    const std::vector<std::string>& list_samples) {

    std::unordered_map<std::string, double> quantitative_pheno;

    std::ifstream file(file_path);
    std::string line;
    int count_pheno = 0;
    bool firstLine = true;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string fid, iid, phenoStr;

        if (!(iss >> fid >> iid >> phenoStr)) {
            stoat::LOG_FATAL("In parsing phenotype, malformed line: " + line);
        }

        if (firstLine) {
            firstLine = false;
            // Check that the header contains FID, IID, and PHENO
            if (fid != "FID" || iid != "IID" || phenoStr != "PHENO") {
                stoat::LOG_FATAL("In parsing phenotype, invalid header: " + line);
            }
            continue;
        }

        try {
            quantitative_pheno[iid] = std::stod(phenoStr);
        } catch(...) {
            stoat::LOG_FATAL("Bad phenotype type : " + phenoStr);
        }
        count_pheno++;
    }

    stoat::LOG_INFO("Quantitative phenotypes founds : " + std::to_string(count_pheno));

    file.close();

    check_match_samples(quantitative_pheno, list_samples);
    std::vector<double> vector_quantitative_pheno;
    vector_quantitative_pheno.reserve(list_samples.size());

    for (const auto& sample : list_samples) {
        auto it = quantitative_pheno.find(sample);
        if (it != quantitative_pheno.end()) {
            vector_quantitative_pheno.push_back(it->second);
        }
    }

    return vector_quantitative_pheno;
}

// Function to open a VCF file and return pointers to the file, header, and record
std::tuple<htsFile*, bcf_hdr_t*, bcf1_t*> parse_vcf(const std::string& vcf_path) {
    // Open the VCF file
    htsFile *ptr_vcf = bcf_open(vcf_path.c_str(), "r");

    // Read the VCF header
    bcf_hdr_t *hdr = bcf_hdr_read(ptr_vcf);
    if (!hdr) {
        bcf_close(ptr_vcf);
        stoat::LOG_FATAL("Could not read VCF header");
    }

    // Initialize a record
    bcf1_t *rec = bcf_init();
    if (!rec) {
        bcf_hdr_destroy(hdr);
        bcf_close(ptr_vcf);
        stoat::LOG_FATAL("Failed to allocate memory for VCF record");
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
            stoat::LOG_FATAL("Key '" + key + "' not found in the phenotype file");
        }
    }
    if (map.size() != keys.size()) {
        stoat::LOG_WARN("Number of samples found in VCF (" + std::to_string(keys.size()) + ") does not match the number of samples in the phenotype file (" + std::to_string(map.size()) + ").");
    }
}

// dict chr:string : vector{(geneName:string, sample_expression:vector<double>, start_pos:size_t, end_pos:size_t)}
std::unordered_map<std::string, std::vector<Qtl_data>> parse_qtl_gene_file(
    const std::string& eqtl_path, 
    const std::string& gene_position_path, 
    const std::vector<std::string>& list_samples) {

    // dict sampleName:string : std::vector<double> sample_expression
    auto qtl = parse_qtl_file(eqtl_path, list_samples); // and check in the same time

    // dict geneName:string : tuple{chrom:string, start_pos:size_t, end_pos:size_t}
    auto gene_position = parse_gene_positions(gene_position_path);
    std::unordered_map<std::string, std::vector<Qtl_data>> qtl_map;

    for (const auto& [gene, expression_vector] : qtl) {
        auto it = gene_position.find(gene);
        if (it != gene_position.end()) {
            const auto& [chrom, start, end] = it->second;
            Qtl_data qtl_info(gene, expression_vector, start, end);
            qtl_map[chrom].emplace_back(qtl_info);
        } else {
            stoat::LOG_FATAL("Gene " + gene + " not found in gene positions.");
        }
    }
  
    // Warn if gene_position has more genes than qtl
    if (gene_position.size() > qtl.size()) {
        stoat::LOG_WARN("More genes present in the gene position file than in the QTL file.");
    }

    return qtl_map;
}

// Function to parse the gene positions file
// dict geneName:string : tuple{chrom:string, start_pos:size_t, end_pos:size_t}
std::unordered_map<std::string, std::tuple<std::string, size_t, size_t>> parse_gene_positions(
    const std::string& filename) {

    std::unordered_map<std::string, std::tuple<std::string, size_t, size_t>> geneMap;
    std::ifstream file(filename);
    std::string line;

    // parse and check header
    std::getline(file, line);
    std::stringstream ss_header(line);
    std::string gene, chrom, startStr, endStr;
    if (!std::getline(ss_header, gene, '\t') || !std::getline(ss_header, chrom, '\t') || !std::getline(ss_header, startStr, '\t') || !std::getline(ss_header, endStr, '\t') ||
         gene != "gene_name" || chrom != "chr" || startStr != "start" || endStr != "end") {
        stoat::LOG_FATAL("In parsing gene position file, invalid header format. Expected: gene_name\tchr\tstart\tend");
    }

    // Check for required columns
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::getline(ss, gene, '\t');
        std::getline(ss, chrom, '\t');
        std::getline(ss, startStr, '\t');
        std::getline(ss, endStr, '\t');

        try {
            int start = std::stoi(startStr);
            int end = std::stoi(endStr);
            geneMap[gene] = std::make_tuple(chrom, start, end);
        } catch (...) {
            stoat::LOG_FATAL("In parsing gene position file, invalid line " + line);
        }
    }

    file.close();
    return geneMap;
}

// Function to parse the qtl file
// dict sampleName:string : std::vector<double> sample_expression
std::unordered_map<std::string, std::vector<double>> parse_qtl_file(
    const std::string& filename, const std::vector<std::string>& list_samples) {

    std::ifstream file(filename);
    std::unordered_map<std::string, std::vector<double>> geneExpressions;

    std::string line;
    bool isHeader = true;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;

        if (isHeader) {
            std::getline(ss, token, '\t');  // Skip the first column (gene name)
            std::vector<std::string> sampleNames;
            while (std::getline(ss, token, '\t')) {
                sampleNames.push_back(token);
            }

            // Check if all sample names are present in the list_samples
            for (const auto& sample : sampleNames) {
                if (std::find(list_samples.begin(), list_samples.end(), sample) == list_samples.end()) {
                    stoat::LOG_FATAL("Sample " + sample + " not found in the list of samples.");
                }
            }

            // warning if the number of samples in the file does not match the number of samples in the list
            if (sampleNames.size() != list_samples.size()) {
                stoat::LOG_WARN("Number of samples in the qtl file is greater that the number of samples in the VCF.");
            }

            isHeader = false;  // Skip header
            continue;
        }

        std::string geneName;
        std::vector<double> expressions;

        std::getline(ss, geneName, '\t');
        while (std::getline(ss, token, '\t')) {
            try {
                expressions.push_back(std::stod(token));
            } catch (...) {
                stoat::LOG_FATAL("Invalid expression value for gene " + geneName + ": " + token);
            }
        }
        geneExpressions[geneName] = expressions;
    }

    file.close();
    return geneExpressions;
}

// Function to parse covariates into an unordered_map
std::vector<std::vector<double>> parse_covariates(
    const std::string& filename, 
    const std::vector<std::string>& covar_names,
    const std::vector<std::string>& list_samples) {

    std::ifstream file(filename);
    std::string line;
    std::vector<std::vector<double>> covariate;
    std::unordered_map<string, std::vector<double>>covariate_map;

    // Read header
    std::getline(file, line);
    std::istringstream headerStream(line);
    std::vector<std::string> headers;
    std::string col;
    while (headerStream >> col) {
        headers.push_back(col);
    }

    // Check for required columns
    auto it_iid = std::find(headers.begin(), headers.end(), "IID");
    if (it_iid == headers.end()) {
        stoat::LOG_FATAL("header must include 'IID' column.\n");
    }
    size_t iid_index = std::distance(headers.begin(), it_iid);

    std::unordered_map<std::string, size_t> col_index;
    for (size_t i = 0; i < headers.size(); ++i) {
        col_index[headers[i]] = i;
    }

    // Check header for covariate names
    for (const auto& name : covar_names) {
        if (col_index.find(name) == col_index.end()) {
            stoat::LOG_FATAL("covariate column '" + name + "' not found in file.\n");
        }
    }

    // Read data
    while (std::getline(file, line)) {
        std::istringstream lineStream(line);
        std::vector<std::string> tokens;
        std::string token;
        while (lineStream >> token) {
            tokens.push_back(token);
        }

        if (tokens.size() <= iid_index) continue;
        std::string iid = tokens[iid_index];

        std::vector<double> selected;
        try {
            for (const auto& name : covar_names) {
                double val = std::stod(tokens[col_index[name]]);
                selected.push_back(val);
            }
        } catch (...) {
            stoat::LOG_FATAL("Individual " + iid + " got an non-numeric value\n");
        }
        covariate_map[iid] = selected;
    }

    check_match_samples(covariate_map, list_samples);

    // Order covariate_map by list_samples
    for (const auto& sample : list_samples) {
        auto it = covariate_map.find(sample);
        if (it != covariate_map.end()) {
            covariate.push_back(it->second);
        } else {
            stoat::LOG_FATAL("Sample " + sample + " not found in the covariate file.");
        }
    }
    file.close();

    return covariate;
}

void check_file(const std::string& file_path) {
    
    std::string line;
    // Check if file is a file
    if (!fs::is_regular_file(file_path)) {
        stoat::LOG_FATAL("File " + file_path + " does not exist.");
    }

    // Check if file can be open
    std::ifstream file(file_path);
    if (!file.is_open()) {
        stoat::LOG_FATAL("Unable to open the file " + file_path);
    }

    // Check if file can be read and not empty file
    if (!std::getline(file, line)) {
        stoat::LOG_FATAL("File " + file_path + "is empty or failed to read header.");
    }

    file.close();
}

void KinshipMatrix::parseKinshipMatrix(const std::string& filename) {

    std::ifstream file(filename);
    std::string line;

    // Parse header line for IDs
    if (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        // Skip the empty top-left cell
        std::getline(ss, token, '\t');
        ids.clear();
        while (std::getline(ss, token, '\t')) {
            ids.push_back(token);
        }
    }

    matrix.clear();

    // Parse matrix rows
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string rowLabel;
        std::getline(ss, rowLabel, '\t'); // row label
        std::vector<double> row;
        std::string value;
        while (std::getline(ss, value, '\t')) {
            row.push_back(std::stod(value));
        }
        matrix.push_back(row);
    }
}

} //end stoat namespace
