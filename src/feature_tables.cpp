#include "feature_tables.hpp"
#include <limits>
#include <iostream>
#include <cassert>
#include <algorithm>

#define DEBUG_TABLES

namespace stoat {

template<class ValueType>
FeatureBySampleTable<ValueType>::FeatureBySampleTable(const std::unordered_map<std::string, size_t>& sample_to_index) :
    sample_to_index(sample_to_index),
    values_per_sample(sample_to_index.size()) {}

template<class ValueType>
CategoricalFeatureBySampleTable<ValueType>::CategoricalFeatureBySampleTable(const std::unordered_map<std::string, size_t>& sample_to_index, 
                                                                            const std::unordered_map<std::string, size_t>& feature_to_index) :
    FeatureBySampleTable<std::vector<ValueType>>::FeatureBySampleTable(sample_to_index),
    feature_to_index(feature_to_index) {

    this->values_per_sample.reserve(this->sample_to_index.size());
    for (size_t i = 0 ; i < this->sample_to_index.size() ; i++) {
        this->values_per_sample[i] = std::vector<ValueType>(feature_to_index.size());
    }
}

// Constructor for a categorical table of doubles fills everything in with a default value of std::numeric_limits<double>::max()
template<>
CategoricalFeatureBySampleTable<double>::CategoricalFeatureBySampleTable(const std::unordered_map<std::string, size_t>& sample_to_index,
                                                                         const std::unordered_map<std::string, size_t>& feature_to_index) :
    FeatureBySampleTable<std::vector<double>>::FeatureBySampleTable(sample_to_index),
    feature_to_index(feature_to_index) {

    this->values_per_sample.reserve(this->sample_to_index.size());
    for (size_t i = 0 ; i < this->sample_to_index.size() ; i++) {
        this->values_per_sample[i] = std::vector<double>(feature_to_index.size(), std::numeric_limits<double>::max());
    }
}
    

// Constructor for a genotype table fills everything in with a default value of 0 for the counts
GenotypeTable::GenotypeTable(const std::unordered_map<std::string, size_t>& sample_to_index, size_t allele_count) :
    FeatureBySampleTable<std::vector<size_t>>::FeatureBySampleTable(sample_to_index) {

    this->values_per_sample.reserve(this->sample_to_index.size());
    for (size_t i = 0 ; i < this->sample_to_index.size() ; i++) {
        this->values_per_sample[i] = std::vector<size_t>(allele_count, 0);
    }
}

// Getter for FeatureBySampleTable
template<class ValueType>
ValueType FeatureBySampleTable<ValueType>::get_value_for_sample(const std::string& sample) const {
    return values_per_sample.at(sample_to_index.at(sample));
}

// Setter for FeatureBySampleTable
template<class ValueType>
void FeatureBySampleTable<ValueType>::set_value_for_sample(const std::string& sample, ValueType value) {
    values_per_sample[sample_to_index.at(sample)] = value;
}

template<class ValueType>
bool FeatureBySampleTable<ValueType>::has_sample(const std::string& sample) const {
    return(sample_to_index.find(sample) != sample_to_index.end());
}
    
template<class ValueType>
std::vector<std::string> FeatureBySampleTable<ValueType>::get_sample_names() const {
    std::vector<std::string> output;
    // JEAN would that help?
    // output.reserve(sample_to_index.size());
    for(const auto& samp_name_idx: sample_to_index) {
        output.emplace_back(samp_name_idx.first);
    }
    return output;
}

template<class ValueType>
const std::unordered_map<std::string, size_t>& FeatureBySampleTable<ValueType>::get_sample_to_index() const {
    return sample_to_index;
}


// Getter for CategoricalFeatureBySampleTable
template<class ValueType>
ValueType CategoricalFeatureBySampleTable<ValueType>::get_value_for_sample_and_feature(const std::string& sample, const std::string& feature) const {
    return this->values_per_sample.at(this->sample_to_index.at(sample)).at(this->feature_to_index.at(feature));
}

// Setter for CategoricalFeatureBySampleTable
template<class ValueType>
void CategoricalFeatureBySampleTable<ValueType>::set_value_for_sample_and_feature(const std::string& sample, const std::string& feature, ValueType value) {
    this->values_per_sample[this->sample_to_index.at(sample)][this->feature_to_index.at(feature)] = value;
}

template<class ValueType>
std::vector<std::string> CategoricalFeatureBySampleTable<ValueType>::get_feature_names() const {
    std::vector<std::string> output;
    // JEAN would that help?
    // output.reserve(feature_to_index.size());
    for(const auto& feat_name_idx: feature_to_index) {
        output.emplace_back(feat_name_idx.first);
    }
    return output;
}

    
void GenotypeTable::increment_count(const std::string& sample, size_t allele_num) {
    #ifdef DEBUG_TABLES
    assert(this->sample_to_index.count(sample));
    assert(this->sample_to_index.at(sample) < this->values_per_sample.size());
    assert(allele_num < this->values_per_sample.at(this->sample_to_index.at(sample)).size());
    #endif
    this->values_per_sample[this->sample_to_index.at(sample)][allele_num]++;
}
size_t GenotypeTable::get_count_for_sample_and_allele(const std::string& sample, size_t allele_num) const {
    return this->values_per_sample[this->sample_to_index.at(sample)][allele_num];
}

std::string GenotypeTable::get_genotype_as_string(const std::string& sample) const {
    std::string genotype = "";
    for (size_t i = 0 ; i < this->values_per_sample[this->sample_to_index.at(sample)].size() ; i++) {
        size_t count = this->values_per_sample[this->sample_to_index.at(sample)][i];
        if (i != 0) {
            genotype += ",";
        }
        genotype += std::to_string(count);
    }
    return genotype;
}

std::string GenotypeTable::allele_paths_as_str() const {
    size_t n_alleles = get_allele_count();
    std::vector<size_t> allele_paths(n_alleles, 0);
    // tally the alleles across samples
    for (auto samp_idx: sample_to_index) {
        for (size_t al_idx = 0; al_idx < n_alleles; al_idx++) {
            if (values_per_sample[samp_idx.second][al_idx] > 0){
                allele_paths[al_idx] += values_per_sample[samp_idx.second][al_idx];
            }
        }
    }

    // join into a string with comma separators
    std::ostringstream allele_paths_str;
    for (size_t i = 0; i < allele_paths.size(); ++i) {
        if (i > 0) allele_paths_str << ",";
        allele_paths_str << allele_paths[i];
    }
    return allele_paths_str.str();
}
    
// Constructor for a gene_expression table fills everything except the gene positions
GeneExpressionTable::GeneExpressionTable(const std::unordered_map<std::string, size_t>& sample_to_index, const std::unordered_map<std::string, size_t>& gene_to_index) :
    CategoricalFeatureBySampleTable<double>::CategoricalFeatureBySampleTable(sample_to_index, gene_to_index) {}

void GeneExpressionTable::read_gene_positions_from_file(const std::string filename) {
    // empty gene positions, just in case
    gene_positions_by_chr.clear();
    // start reading the file
    std::ifstream file(filename);
    std::string line;

    // Read and validate header (assumes file and first line already checked)
    std::getline(file, line);
    std::stringstream ss_header(line);
    std::string gene, chrom, startStr, endStr;
    std::getline(ss_header, gene, '\t');
    std::getline(ss_header, chrom, '\t');
    std::getline(ss_header, startStr, '\t');
    std::getline(ss_header, endStr, '\t');

    if (gene != "gene_name" || chrom != "chr" || startStr != "start" || endStr != "end") {
        throw std::invalid_argument("Invalid gene position header. Expected: gene_name, chr, start, and end. Got " + line);
    }

    // read the gene positions and fill our map
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string gene_val, chrom_val, start_val, end_val;

        if (!std::getline(ss, gene_val, '\t') ||
            !std::getline(ss, chrom_val, '\t') ||
            !std::getline(ss, start_val, '\t') ||
            !std::getline(ss, end_val, '\t')) {
            throw std::invalid_argument("In parsing gene position file, malformed line: " + line);
        }

        try {
            // try to make the gene position object and save it
            size_t start = std::stoull(start_val);
            size_t end = std::stoull(end_val);
            gene_position_t gpos(gene_val, start, end);
            gene_positions_by_chr[chrom_val].emplace_back(gpos);
        } catch (...) {
            throw std::invalid_argument("In parsing gene position file, invalid numeric value in line: " + line);
        }
    }

    file.close();
}

std::vector<std::string> GeneExpressionTable::get_genes_around_pos(const std::string chrom, const size_t start_pos, const size_t end_pos, const size_t max_distance) const {
    // we'll save the nearby genes here
    std::vector<std::string> near_genes;
    
    // we look for genes in the specified range +- the maximum distance
    size_t min_start_pos = (start_pos > max_distance) ? start_pos - max_distance : 0;
    size_t max_end_pos = end_pos + max_distance;

    // loop over all genes and save the ones that are close enough
    for (const gene_position_t gene_pos: gene_positions_by_chr.at(chrom)) {
        // for (size_t gene_i = 0; gene_i < gene_positions_by_chr[chrom].size(); gene_i++) {
        //     gene_position_t gene_pos = gene_positions_by_chr[chrom][gene_i];
        if (!(gene_pos.end < min_start_pos || gene_pos.start > max_end_pos)) {
            near_genes.push_back(gene_pos.gene);
        }
    }
    
    return near_genes;
}
    
CombinedTable::CombinedTable(const GenotypeTable& genotypes){
    n_alleles = genotypes.get_allele_count();
    predictors.resize(n_alleles);
    sample_names = genotypes.get_sample_names();
    for(std::string samp_name: sample_names) {
        for (int al_i = 0; al_i < n_alleles; al_i++){
            size_t ac = genotypes.get_count_for_sample_and_allele(samp_name, al_i);
            predictors[al_i].push_back(static_cast<double>(ac));
        }
    }
}
    
int CombinedTable::get_n_alleles() const {
    return n_alleles;
}

int CombinedTable::get_n_covariates() const {
    return predictors.size() - n_alleles;
}

Eigen::MatrixXd CombinedTable::make_matrixXd_features() {
    // we'll prepare a matrix with an additional "intercept" first column, filled with ones
    Eigen::MatrixXd X = Eigen::MatrixXd::Constant(phenotype_vec.size(), predictors.size() + 1, 1);
    // fill other columns with the values all the predictors
    for (size_t pred_i = 0; pred_i < predictors.size(); ++pred_i) {
        for (size_t samp_i = 0; samp_i < X.rows(); ++samp_i) {
            X(samp_i, pred_i + 1) = predictors[pred_i][samp_i];
        }
    }
    return X;
}

Eigen::VectorXd CombinedTable::make_vectorxd_phenotype() {
    Eigen::VectorXd y(phenotype_vec.size());
    for (size_t samp_i = 0; samp_i < phenotype_vec.size(); ++samp_i) {
        y(samp_i) = phenotype_vec[samp_i];
    }
    return y;
}

std::vector<double> CombinedTable::get_phenotype() const {
    return phenotype_vec;
}
    
std::vector<std::vector<double>> CombinedTable::get_predictors() const {
    return predictors;
}
    
void CombinedTable::combine_binary_phenotype(const BinaryPhenotypeTable& phenotype) {
    assert(phenotype_vec.empty());
    // loop over the samples and fill up the phenotype vector
    for(std::string samp_name: sample_names) {
        if (phenotype.has_sample(samp_name)) {
            if (phenotype.get_value_for_sample(samp_name)){
                phenotype_vec.push_back(1);
            } else {
                phenotype_vec.push_back(0);
            }
        } else {
            // this sample is not in the matrix, mark for filtering and fill with 0s
            samples_to_exclude[samp_name] = true;
            phenotype_vec.push_back(0);
        }
    }
}

void CombinedTable::combine_quantitative_phenotype(const QuantitativePhenotypeTable& phenotype) {
    assert(phenotype_vec.empty());
    // loop over the samples and fill up the phenotype vector
    for(std::string samp_name: sample_names) {
        if (phenotype.has_sample(samp_name)) {
            phenotype_vec.push_back(phenotype.get_value_for_sample(samp_name));
        } else {
            // this sample is not in the matrix, mark for filtering and fill with 0s
            samples_to_exclude[samp_name] = true;
            phenotype_vec.push_back(0);
        }
    }
}

void CombinedTable::combine_covariates(const CovariateTable& covariates) {
    // don't do anything if there are no covariates
    if (covariates.get_feature_names().size() == 0) {
        return;
    }
    // there should only be alleles in the combined table at this point
    assert(predictors.size() == n_alleles);
    // loop over the samples and covariates and fill new columns in the table
    std::vector<std::string> covar_names = covariates.get_feature_names();
    size_t ncovar = covar_names.size();
    // extend the predictors
    for (int covar_i = 0; covar_i < ncovar; covar_i++){
        // init with 0?
        std::vector<double> covar_init(sample_names.size(), 0);
        predictors.emplace_back(covar_init);
    }
    // fill the table with the covariate values
    for (size_t samp_ii = 0; samp_ii < sample_names.size(); samp_ii++) {
        std::string samp_name = sample_names[samp_ii];
        if (covariates.has_sample(samp_name)) {
            for (int covar_i = 0; covar_i < ncovar; covar_i++){
                double covar = covariates.get_value_for_sample_and_feature(samp_name, covar_names[covar_i]);
                predictors[n_alleles + covar_i][samp_ii] = covar;
            }
        } else {
            // this sample is not in the matrix, mark for filtering and fill with 0s
            samples_to_exclude[samp_name] = true;
        }
    }
    // JEAN here we are marking samples to exclude because of missing values. Maybe we also want a check somewhere to warn the user when the different sample sets don't match? 
}

void CombinedTable::remove_missing_values() {
    // nothing to do if no samples to exclude
    if (samples_to_exclude.empty()) {
        return;
    }
    
    // remove marked samples in the phenotype vector
    auto pheno_it = phenotype_vec.begin();
    int samp_i = 0;
    while (pheno_it != phenotype_vec.end()) {
        if (samples_to_exclude.find(sample_names[samp_i]) != samples_to_exclude.end()) {
            pheno_it = phenotype_vec.erase(pheno_it);
        } else {
            pheno_it++;
        }
        samp_i++;
    }
    
    // remove marked samples in the predictor matrix
    for (size_t pred_i = 0; pred_i < predictors.size(); pred_i++) {
        auto pred_it = predictors[pred_i].begin();
        samp_i = 0;
        while (pred_it != predictors[pred_i].end()) {
            if (samples_to_exclude.find(sample_names[samp_i]) != samples_to_exclude.end()) {
                pred_it = predictors[pred_i].erase(pred_it);
            } else {
                pred_it++;
            }
            samp_i++;
        }
    }

    // remove sample from sample names
    auto sampnames_it = sample_names.begin();
    samp_i = 0;
    while (sampnames_it != sample_names.end()) {
        if (samples_to_exclude.find(sample_names[samp_i]) != samples_to_exclude.end()) {
            sampnames_it = sample_names.erase(sampnames_it);
        } else {
            sampnames_it++;
        }
        samp_i++;
    }
    
    // clear samples_to_exclude
    samples_to_exclude.clear();
}

void CombinedTable::remove_constant_predictors() {
    // just in case
    this->remove_missing_values();
    
    // check each predictor using an iterator to remove it efficiently(-ish) if needed
    auto preds_it = predictors.begin();
    int pred_idx = 0; // to figure out if we removed an allele or covariate
    while (preds_it != predictors.end()) {
        // we'll test if any value is different from the first value
        std::vector<double>& pred_it = *preds_it; // dereference the iterator to get a reference to predictor vector
        double first_value = pred_it[0];
        bool constant = true;
        for (double pred_value: pred_it) {
            if (pred_value != first_value) {
                // not constant, we can stop already
                constant = false;
                break;
            }
        }
        // remove if constant
        if (constant){
            preds_it = predictors.erase(preds_it);
            // if we've removed an allele, decrement the counter
            if (pred_idx < n_alleles){
                n_alleles--;
            }
        } else {
            preds_it++;
            pred_idx++;
        }
    }
}

void CombinedTable::remove_duplicated_predictors() {
    // just in case
    this->remove_missing_values();

    // first find which predictor should be removed, starting from the end so that we remove covariates over alleles
    std::vector<bool> to_remove(predictors.size(), false);
    size_t n_samples = predictors.front().size();
    bool any_removed = false;
    for (size_t pred_i = predictors.size() - 1; pred_i > 0; pred_i--) {
        // compare to all the other predictors before
        for (size_t pred_j = 0; pred_j < pred_i; pred_j++) {
            size_t samp_i = 0;
            while (samp_i < n_samples && predictors[pred_i][samp_i] == predictors[pred_j][samp_i]) {
                samp_i++;
            }
            if (samp_i == n_samples) {
                // duplicated predictors, mark and stop looking for more
                to_remove[pred_i] = true;
                any_removed = true;
                break;
            }
        }
    }

    // stop if none were marked for removal
    if (!any_removed) {
        return;
    }
    
    // remove the predictors efficiently(-ish) using an iterator
    auto preds_it = predictors.begin();
    int new_pred_idx = 0; // to figure out if we're removing an allele or covariate
    int old_pred_idx = 0; // to check if the original predictor was marked for removal
    while (preds_it != predictors.end()) {
        // remove if marked
        if (to_remove[old_pred_idx]){
            preds_it = predictors.erase(preds_it);
            // if we've removed an allele, decrement the counter
            if (new_pred_idx < n_alleles){
                n_alleles--;
            }
        } else {
            // move to the next element
            preds_it++;
            new_pred_idx++;
        }
        old_pred_idx++;
    }
}

void CombinedTable::remove_one_allele() {
    // if there is only one allele, do nothing
    if (n_alleles < 2) {
        return;
    }

    // remove the first predictor
    predictors.erase(predictors.begin());
    
    // decrement the number of alleles
    n_alleles--;
}
    
void CombinedTable::add_total_allele_count_covariable() {
    // just in case
    this->remove_missing_values();
    
    // to check if at least one sample has a different total allele count
    bool all_same_ac = true;
    // prepare the new "covariate" holding the total allele count per sample
    size_t n_samples = sample_names.size();
    std::vector<double> tot_ac(n_samples, 0);
    for (size_t samp_i = 0; samp_i < n_samples; samp_i++) {
        // compute the total allele count for this sample
        size_t samp_ac = 0;
        for (size_t pred_i = 0; pred_i < n_alleles; pred_i++) {
            if (predictors[pred_i][samp_i] > 0) {
                samp_ac += predictors[pred_i][samp_i];
            }
        }
        // append to new vector and check if different from the first one
        tot_ac[samp_i] = samp_ac;
        if (samp_i > 0 && samp_ac != tot_ac[0]) {
            all_same_ac = false;
        }
    }
    // add as a covariate if worth it
    if (!all_same_ac) {
        predictors.emplace_back(tot_ac);
    }
}
    
bool CombinedTable::passes_filters(const double maf, const size_t min_individuals) const {
    // JEAN this assumes we've removed samples with no genotype information (or that they've never were included in that combined table)

    // make sure there is/was at least two alleles
    if (n_alleles < 2) {
        stoat::LOG_DEBUG("Filtered: less than 2 independent alleles (" + std::to_string(n_alleles) + ")");
        return false;
    }

    // make sure there are enough individuals
    if(sample_names.size() < min_individuals) {
        stoat::LOG_DEBUG("Filtered: not enough individuals: " + std::to_string(sample_names.size()));
        return false;
    }
    
    // make sure the allele frequency is high enough
    if (maf == 0) {
        return true;
    }
    // first count the total number of observed alleles for each allele
    size_t total_allele_counts = 0;
    std::vector<size_t> allele_counts(n_alleles, 0);
    for (size_t al_ii = 0; al_ii < n_alleles; al_ii++) {
        size_t allele_count_ii = 0;
        for (double ac: predictors[al_ii]){
            allele_count_ii += ac;
        }
        total_allele_counts += allele_count_ii;
        allele_counts[al_ii] = allele_count_ii;
    }
    // compute the frequency of the second most frequent allele
    std::sort(allele_counts.begin(), allele_counts.end(), std::greater<size_t>());
    double af = static_cast<double>(allele_counts[1]) / total_allele_counts;
    // filter if too low
    if (af < maf) {
        stoat::LOG_DEBUG("Filtered: less than two alleles with frequency above " + std::to_string(maf));
        return false;
    }
    return true;
}

// Utils function to summarize the table in the output (binary phenotype)
// Write a std::string of: g0[0]:g1[1],g0[1]:g1[1],g0[2]:g1[2]...
std::string format_group_paths(const std::vector<size_t>& g0, const std::vector<size_t>& g1) {
    std::string result;
    size_t numb_col = g0.size();
    for (size_t index_col = 0; index_col < numb_col; ++index_col) {
        result += std::to_string(g0[index_col]) + ":" + std::to_string(g1[index_col]);
        if (index_col < numb_col - 1) {
            result += ","; // Separate row pairs with ','
        }
    }
    return result;
}

    
} //end stoat namespace

// Apparently these definitions are supposed to be done here, after all the members are defined
template class stoat::FeatureBySampleTable<bool>;
template class stoat::FeatureBySampleTable<double>;
template class stoat::FeatureBySampleTable<std::vector<size_t>>;
template class stoat::FeatureBySampleTable<std::vector<double>>;
template class stoat::CategoricalFeatureBySampleTable<double>;
