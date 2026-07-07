#include "test.hpp"

#include <chrono>
#include <getopt.h>
#include <filesystem>

#include "../log.hpp"
#include "../banner.hpp"
#include "../snarl_analyzer.hpp"
#include "../arg_parser.hpp"
#include "../writer.hpp"

// #define USE_CALLGRIND

#ifdef USE_CALLGRIND
    #include <valgrind/callgrind.h>
#endif

namespace stoat_command {

// STOAT_VERSION are define in the CMAKELIST file
void print_help_test() {
    stoat::print_banner(std::string(STOAT_VERSION));
    std::cerr << "Usage: stoat test [options]\n\n"
              << "  -g, --genotype FILE             Path to the genotype file from stoat graph or stoat vcf\n"
              << "  -m, --method STR                Which test method to use: chi2 (Fisher/Chi-Squared), linreg (linear regression), \n"
              << "                                                            logreg (logistic regression), exact (exact phenotype/genotype match)\n"
              << "  -p, --phenotype FILE            Path to the phenotype file\n"
              << "  -P, --gene-position FILE        Path to the gene position file (activates the eQTL testing mode)\n"
              << "  -w, --max-gene-distance INT     Include snarls up to this distance from the gene when looking for eQTLs [1000000]\n"
              << "  -c, --covariate FILE            Path to the covariate file\n"
              << "  -C, --covar-name NAME           Covariate column name(s) used (default: use all covariates in file)\n"
              << "  -I, --min-individuals INT       Minimum number of individuals with at least one allele in the snarl [2]\n"
              << "  -M, --maf FLOAT                 Minimum allele frequency threshold [0.05]\n"
              << "  -t, --threads INT               Number of threads to use [1]\n"
              << "  -V, --verbose INT               Verbosity level (0=error, 1=warn, 2=info, 3=debug, 4=trace) [2]\n"
              << "  -o, --output FILE               Output directory name [stoat_output]\n"
              << "  -u, --no-bgzip                  Don't compress the output file with bgzip\n"
              << "  -a, --ascii                     Print the STOAT ascii art banner\n"
              << "  -h, --help                      Print this help message\n";
}

int main_stoat_test(int argc, char* argv[]) {
    
    // Declare variables to hold argument values
    std::string genotype_path, phenotype_path, covariate_path, gene_position_path;
    bool bgzip_output = true;
    bool ascii = false;

    // default value for the filter thresholds
    double maf_threshold = 0.05;
    size_t min_individuals = 2;
    size_t max_gene_dist = 1000000;

    // which method to use
    std::string method = "";
    
    // default output directory
    std::string output_dir = "stoat_output";

    // will store the names of the covariates to use
    std::vector<std::string> covar_names;

    size_t thread_count = 1;

    // Parse arguments
    int c;

    static struct option long_options[] = {
        {"genotype", required_argument, 0, 'g'},
        {"method", required_argument, 0, 'm'},
        {"phenotype", required_argument, 0, 'p'},
        {"gene-position", required_argument, 0, 'P'},
        {"max-gene-distance", required_argument, 0, 'w'},
        {"covariate", required_argument, 0, 'c'},
        {"covar-name", required_argument, 0, 'C'},
        {"min-individuals", required_argument, 0, 'I'},
        {"maf", required_argument, 0, 'M'},
        {"threads", required_argument, 0, 't'},
        {"verbose", required_argument, 0, 'V'},
        {"output", required_argument, 0, 'o'},
        {"ascii", no_argument, 0, 'a'},
        {"no-bgzip", no_argument, 0, 'u'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    while ((c = getopt_long(argc, argv, "g:m:p:P:w:c:C:I:M:t:V:o:uh", long_options, nullptr)) != -1) {
        switch (c) {
            case 'g': genotype_path = optarg; stoat_vcf::check_file(genotype_path); break;
            case 'a': ascii = true; break;
            case 'm': method = optarg; stoat_vcf::check_methods(method); break;
            case 'p': phenotype_path = optarg; stoat_vcf::check_file(phenotype_path); break;
            case 'c': covariate_path = optarg; stoat_vcf::check_file(covariate_path); break;
            case 'C': {
                std::stringstream ss(optarg);
                std::string token;
                while (std::getline(ss, token, ',')) covar_names.push_back(token);
                break;
            }
            case 'I':
                min_individuals = std::stoi(optarg);
                if (min_individuals < 2) {
                    throw std::runtime_error("Error: [stoat test] min_individuals threshold must be > 1");
                }
                break;
            case 'P': gene_position_path = optarg; stoat_vcf::check_file(gene_position_path); break;
            case 'w':
                max_gene_dist = std::stoi(optarg);
                if (max_gene_dist < 1) {
                    throw std::runtime_error("Error: [stoat test] Maximum gene distance (-w, --max-gene-distance) must be > 0");
                }
                break;
            case 'M':
                maf_threshold = std::stod(optarg);
                if (maf_threshold < 0 || maf_threshold > 1) {
                    throw std::runtime_error("Error: [stoat test] MAF must be in [0,1]");
                }
                break;
            case 't':
                if (std::stoi(optarg) < 1) {
                    stoat::LOG_ERROR("[stoat graph] Number of threads must be > 0");
                    return EXIT_FAILURE;
                }
                thread_count = std::stoi(optarg);
                omp_set_num_threads(thread_count);

                break;
            case 'V': 
                {
                int level = std::stoi(optarg);
                if (level < 0 || level > 4) {
                    throw std::runtime_error("Error: [stoat test] Invalid verbosity level. Use 0=Error, 1=Warn, 2=Info, 3=Debug, 4=Trace");
                }
                stoat::LogLevel logLevel = static_cast<stoat::LogLevel>(level);
                stoat::Logger::instance().setLevel(logLevel);                
                break;
                }
            case 'o': output_dir = optarg; break;
            case 'u': bgzip_output = false; break;
            case 'h': 
                print_help_test(); 
                return EXIT_SUCCESS; 
            default:
                stoat::LOG_ERROR("[stoat test] Unknown argument " + c);
                print_help_test();
                return EXIT_FAILURE;
        }
    }

    if (argc == 2) {
        print_help_test();
        return EXIT_FAILURE;
    }

    if (phenotype_path.empty()) {
        stoat::LOG_ERROR("[stoat test] a phenotype file must be provided with -p/--phenotype");
        print_help_test();
        return EXIT_FAILURE;
    }

    if (genotype_path.empty()) {
        stoat::LOG_ERROR("[stoat test] a genotype file, made with stoat graph or stoat vcf, must be provided with -g/--genotype");
        print_help_test();
        return EXIT_FAILURE;
    }


    std::filesystem::create_directory(output_dir);
    stoat::Logger::instance().setLogFile(output_dir + "/stoat.test.log");

    // add command launch in log file
    if (ascii) {
        print_ascii_banner(std::string(STOAT_VERSION));
    } else {
        print_banner(std::string(STOAT_VERSION));
    }

    std::stringstream ss;
    ss << "stoat ";
    for (int i = 0; i < argc; ++i) ss << argv[i] << " ";
    stoat::LOG_SILENTE(ss.str());

    auto start_total_timer = std::chrono::high_resolution_clock::now();

// Start tracking with callgrind
#ifdef USE_CALLGRIND
    CALLGRIND_START_INSTRUMENTATION;
#endif

    // Load the SnarlDataCollection
    stoat::SnarlDataCollection snarl_collection(0, 0, 0);
    // load the header from the snarl collection file. We'll use those sample indices
    std::shared_ptr<stoat::Reader> snarl_reader;
    if ((genotype_path.compare(genotype_path.length()-3, 3, ".gz") == 0) ||
        (genotype_path.compare(genotype_path.length()-4, 4, ".bgz") == 0)) {
        snarl_reader.reset(new BgzReader(genotype_path));
    } else {
        snarl_reader.reset(new StdReader(genotype_path));
    }
    snarl_collection.load_snarl_data_collection(*snarl_reader, true);
    snarl_reader->close();

    //////////////////////////////////////// Go through the genotypes and test against the phenotype    
    stoat::LOG_INFO("Starting GWAS analysis...");
    
    //////////////////// Load the phenotypes and covariate matrix from files
    std::unique_ptr<stoat::BinaryPhenotypeTable> binary_phenotype_table;
    std::unique_ptr<stoat::QuantitativePhenotypeTable> quantitative_phenotype_table;
    std::unique_ptr<stoat::GeneExpressionTable> gene_expression_table;

    // prepare the vector mapping samples to index here because other objects (phenotype or genotypes) will use it
    std::unordered_map<std::string, size_t> sample_to_index = snarl_collection.get_sample_to_index_copy();

    // prepare the vector mapping genes to index
    std::unordered_map<std::string, size_t> gene_to_index;

    // predict methods
    if (method.empty()) {
        method = stoat_vcf::methods_stats_prediction(phenotype_path, !covariate_path.empty());
        stoat::LOG_INFO("No statistical method provided; using " + method + " test.");
    }

    // read phenotype file
    if (!gene_position_path.empty()) {
        stoat::LOG_TRACE("Parsing eQTL phenotype file");
        gene_expression_table = stoat_vcf::parse_gene_expression_table(phenotype_path, gene_position_path, sample_to_index, gene_to_index);
    } else if (method == "chi2" || method == "logreg" || method == "exact") {
        stoat::LOG_TRACE("Parsing binary phenotype file");
        binary_phenotype_table = stoat_vcf::parse_binary_pheno_table(phenotype_path, sample_to_index);
    } else if (method == "linreg") {
        stoat::LOG_TRACE("Parsing quantitative phenotype file");
        quantitative_phenotype_table = stoat_vcf::parse_quantitative_pheno_table(phenotype_path, sample_to_index);
    } else {
        stoat::LOG_ERROR("Method : " + method + " not recognized");
    }

    // eventually parse the covariate file
    std::unique_ptr<stoat::CovariateTable> covariate_table = std::unique_ptr<stoat::CovariateTable>(new CovariateTable({}, {}));

    // prepare the vector mapping covariates to index
    // needs to be defined here to stay in memory because the CovariateTable don't store it
    std::unordered_map<std::string, size_t> covar_to_index;
    for (std::string covar: covar_names) {
        covar_to_index[covar] = covar_to_index.size();
    }

    if (!covariate_path.empty()) {
        stoat::LOG_TRACE("Parsing covariate file");
        covariate_table = stoat_vcf::parse_covariate_table(covariate_path, sample_to_index, covar_to_index);
    }
    
    // the object to orchestrate the testing of the snarls
    std::shared_ptr<stoat_vcf::SnarlAnalyzer> snarl_analyzer;

    if (method == "exact") {
        // Exact test between binary phenotype and genotypes
        snarl_analyzer.reset(new stoat_vcf::ExactBinarySnarlAnalyzer(snarl_collection, maf_threshold,
                                                                     *binary_phenotype_table, min_individuals));
    } else if (method == "chi2") {
        if (!covariate_path.empty()) {
            stoat::LOG_WARN("Covariate will not be used by the Chi2 test. To include them in the test use '-m logreg' instead.", 0);
        }
        // Binary using Chi2/Fisher (no covariate)
        snarl_analyzer.reset(new stoat_vcf::BinarySnarlAnalyzer(snarl_collection, maf_threshold,
                                                                *binary_phenotype_table, min_individuals));
    } else if (method == "logreg") {
        // Binary using logistic regression, possibly with covariates
        snarl_analyzer.reset(new stoat_vcf::BinaryCovarSnarlAnalyzer(snarl_collection, *covariate_table, maf_threshold,
                                                                     *binary_phenotype_table, min_individuals));
    } else if (!gene_position_path.empty()) {
        // expression QTL analysis
        snarl_analyzer.reset(new stoat_vcf::EQTLSnarlAnalyzer(snarl_collection, *covariate_table, maf_threshold, 
                                                              *gene_expression_table, max_gene_dist, min_individuals));
    } else if (method == "linreg") {
        // quantitative phenotype using linear regression
        snarl_analyzer.reset(new stoat_vcf::QuantitativeSnarlAnalyzer(snarl_collection, *covariate_table, maf_threshold, 
                                                                      *quantitative_phenotype_table, min_individuals));
    }

    // prepare to write the TSV
    // JEAN if we want to keep track of what was run, we might as well include a header in the file with the full info (all parameters, input files, etc)
    std::shared_ptr<stoat::Writer> out_writer;
    if (bgzip_output) {
        out_writer.reset(new BgzWriter(output_dir + "/stoat.assoc.pvalues.tsv.gz", thread_count));
    } else {
        out_writer.reset(new StdWriter(output_dir + "/stoat.assoc.pvalues.tsv", thread_count));
    }

    // guess if the input genotype file is bgzipped based on the suffix
    std::shared_ptr<stoat::Reader> gt_reader;
    if ((genotype_path.compare(genotype_path.length()-3, 3, ".gz") == 0) ||
        (genotype_path.compare(genotype_path.length()-4, 4, ".bgz") == 0)) {
        gt_reader.reset(new BgzReader(genotype_path));
    } else {
        gt_reader.reset(new StdReader(genotype_path));
    }

    // Test each snarl, line by line
    snarl_analyzer->test_snarls_from_file(*gt_reader, *out_writer);

    // close the file connections
    gt_reader->close();
    out_writer->close();

    auto end_total_timer = std::chrono::high_resolution_clock::now();
    stoat::LOG_INFO("stoat test took " + std::to_string(std::chrono::duration<double>(end_total_timer - start_total_timer).count()) + " s");
    return EXIT_SUCCESS;
}

} // end stoat
