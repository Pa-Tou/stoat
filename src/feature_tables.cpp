#include "feature_tables.hpp"
#include <limits>
#include <iostream>
#include <cassert>
#include <algorithm>

// #define DEBUG_TABLES

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
    for (size_t i = 0; i < sample_to_index.size(); i++) {
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
    for (size_t i = 0; i < sample_to_index.size(); i++) {
        this->values_per_sample[i] = std::vector<double>(feature_to_index.size(), std::numeric_limits<double>::max());
    }
}

// Getter for FeatureBySampleTable
template<class ValueType>
ValueType FeatureBySampleTable<ValueType>::get_value_for_sample(const std::string& sample) const {
    return values_per_sample.at(sample_to_index.at(sample));
}

template<class ValueType>
ValueType FeatureBySampleTable<ValueType>::get_value_for_sample_id(size_t sample_idx) const {
    return values_per_sample.at(sample_idx);
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

template<class ValueType>
ValueType CategoricalFeatureBySampleTable<ValueType>::get_value_for_sample_and_feature_ids(size_t sample_idx, size_t feature_idx) const {
    return this->values_per_sample.at(sample_idx).at(feature_idx);
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

template<class ValueType>
size_t CategoricalFeatureBySampleTable<ValueType>::get_feature_number() const {
    return feature_to_index.size();
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

// Constructor for a genotype table fills everything in with a default value of 0 for the counts
GenotypeTable::GenotypeTable(const std::unordered_map<std::string, size_t>& sample_to_index, size_t allele_count) :
    FeatureBySampleTable<std::vector<double>>::FeatureBySampleTable(sample_to_index) {

    this->values_per_sample.reserve(this->sample_to_index.size());

    for (size_t i = 0; i < sample_to_index.size(); i++) {
        this->values_per_sample[i] = std::vector<double>(allele_count, 0);
    }

    // init the variable to keep track of the matrix dimensions
    n_alleles = allele_count;
    n_covariates = 0;
    n_samples = sample_to_index.size();

    // init the masks and the variables to keep track of active rows/columns
    col_mask = std::vector<bool>(n_alleles, false);
    row_mask = std::vector<bool>(n_samples, false);
    n_active_samples = n_samples;
    n_active_alleles = n_alleles;
    n_active_columns = n_alleles;

    // init the column with the total allele counts
    total_allele_counts_per_sample = std::vector<double>(n_samples, 0);
    use_total_ac = false;
}
    
void GenotypeTable::clear() {

    // we don't touch the alleles but reinit the covariates (for now)
    n_covariates = 0;
    n_active_samples = n_samples;
    n_active_alleles = n_alleles;
    n_active_columns = n_alleles;

    // reinit the masks
    col_mask = std::vector<bool>(n_alleles, false);
    row_mask = std::vector<bool>(n_samples, false);

    // reinit the total allele count column too
    use_total_ac = false;
}

void GenotypeTable::increment_count(size_t sample_idx, size_t allele_num) {
    // increment the appropriate allele column and total count
    this->values_per_sample[sample_idx][allele_num]++;
    total_allele_counts_per_sample[sample_idx]++;
}

double GenotypeTable::get_value(size_t row, size_t col) const {
    // depending on the column index, get the value from the genotype or covariates
    // the virtual matrix starts with the alleles and then the covariates
    if (col < n_alleles) {
        // we're looking for an allele column
        return values_per_sample.at(row).at(col);
    } else if (col < n_alleles + n_covariates) {
        // we're looking for a covariate
        return covariates->get_value_for_sample_and_feature_ids(row, col - n_alleles);
    } else {
        // we're overflowing
        stoat::LOG_ERROR("Column overflow of the table.");
    }
}

size_t GenotypeTable::get_count_for_sample_and_allele(const std::string& sample, size_t allele_num) const {
    return this->values_per_sample.at(this->sample_to_index.at(sample)).at(allele_num);
}

std::string GenotypeTable::allele_paths_as_str() const {
    std::vector<size_t> allele_paths(n_active_alleles, 0);
    // tally the alleles across samples
    for (size_t samp_i = 0; samp_i < n_samples; samp_i++) {
        // skip if sample was masked
        if (row_mask[samp_i]) {
            continue;
        }
        size_t al_c_i = 0;
        for (size_t al_i = 0; al_i < n_alleles; al_i++) {
            // skip if allele was masked
            if (col_mask[al_i]) {
                continue;
            }
            allele_paths[al_c_i] += values_per_sample.at(samp_i).at(al_i);
            al_c_i++;
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

std::string GenotypeTable::get_genotype_as_string(const std::string& sample) const {
    std::string genotype = "";
    for (size_t i = 0 ; i < this->values_per_sample.at(this->sample_to_index.at(sample)).size() ; i++) {
        size_t count = this->values_per_sample.at(this->sample_to_index.at(sample)).at(i);
        if (i != 0) {
            genotype += ",";
        }
        genotype += std::to_string(count);
    }
    return genotype;
}

size_t GenotypeTable::get_allele_count() const {
    return this->values_per_sample.size() == 0 ? 0 : this->values_per_sample.front().size(); 
}

void GenotypeTable::link_to_binary_phenotype(const BinaryPhenotypeTable& phenotype) {
    b_phenotype = &phenotype;
    is_binary_pheno = true;
}

void GenotypeTable::link_to_quantitative_phenotype(const QuantitativePhenotypeTable& phenotype) {
    q_phenotype = &phenotype;
    is_binary_pheno = false;
}

void GenotypeTable::link_to_covariates(const CovariateTable& in_covariates) {
    covariates = &in_covariates;
    n_covariates = in_covariates.get_feature_number();
    n_active_columns += n_covariates;
    col_mask.resize(n_alleles + n_covariates, false);
}

void GenotypeTable::remove_noncovered_samples() {

    // check the total allele count (filled previously) to decide if a sample should be masked
    for (size_t samp_i = 0; samp_i < n_samples; samp_i++) {
        // no allele counts present in sample or sample not have phenotype 
        if (total_allele_counts_per_sample[samp_i] == 0 && !row_mask[samp_i]) {
            row_mask[samp_i] = true;
            n_active_samples--;
        }
    }
    // JEAN TODO check that this is called before the test/regression
}

void GenotypeTable::remove_constant_predictors() {

    // don't do anything if there are no samples or predictors, it will be filtered out later
    if (n_alleles + n_covariates == 0 || n_samples == 0) {
        return;
    }

    // check each predictor
    for (size_t col_ii = 0; col_ii < n_alleles + n_covariates; col_ii++) {

        // skip if already masked
        if (col_mask[col_ii]) {
            continue;
        }

        // compare each value with the first value
        double first_value;
        bool constant = false;
        for (size_t row_ii = 0; row_ii < n_samples; row_ii++) {

            // skip if already masked
            if (row_mask[row_ii]) {
                continue;
            }

            // if first value to consider, save it
            if (!constant) {
                first_value = get_value(row_ii, col_ii);
                constant = true;
            }

            if (get_value(row_ii, col_ii) != first_value) {
                // not constant, we can stop already
                constant = false;
                break;
            }
        }

        // mask if constant
        if (constant){
            col_mask[col_ii] = true;

            // update the number of active columns
            n_active_columns--;
            if (col_ii < n_alleles) {
                n_active_alleles--;
            }
        }
    }
}

void GenotypeTable::remove_duplicated_predictors() {
    // first find which predictor should be removed, starting from the end so that we remove covariates over alleles
    for (size_t col_ii = 0; col_ii < n_alleles + n_covariates; col_ii++) {
        // skip if already masked
        if (col_mask[col_ii]) {
            continue;
        }
        // compare to all the other predictors before
        for (size_t col_jj = 0; col_jj < col_ii; col_jj++) {
            // skip if that one already masked
            if (col_mask[col_jj]) {
                continue;
            }
            // compare each (not masked) rows
            size_t samp_i = 0;
            while (samp_i < n_samples && (row_mask[samp_i] || get_value(samp_i, col_ii) == get_value(samp_i, col_jj))) {
                samp_i++;
            }
            // if one value is different, we will stop before reaching the end
            // if we reached the end, the columns are duplicated (on the active samples)
            if (samp_i == n_samples) {
                // duplicated predictors, mark and stop looking for more
                col_mask[col_ii] = true;
                n_active_columns--;
                if (col_ii < n_alleles) {
                    n_active_alleles--;
                }
                break;
            }
        }
    }
    // if the total allele count column is to be used, check if it's maybe duplicated (with a covariate most likely)
    if (use_total_ac) {
        // compare to all the other predictors
        for (size_t col_jj = 0; col_jj < n_alleles + n_covariates; col_jj++) {
            // skip if that one already masked
            if (col_mask[col_jj]) {
                continue;
            }
            // compare each (not masked) rows
            size_t samp_i = 0;
            while (samp_i < n_samples && (row_mask[samp_i] || total_allele_counts_per_sample[samp_i] == get_value(samp_i, col_jj))) {
                samp_i++;
            }
            // if one value is different, we will stop before reaching the end
            // if we reached the end, the columns are duplicated (on the active samples)
            if (samp_i == n_samples) {
                // duplicated, then don't use that total allele count column
                use_total_ac = false;
                n_active_columns--;
                break;
            }
        }
    }
}

void GenotypeTable::fill_contingency_table(std::vector<size_t>& g0, std::vector<size_t>& g1) const {
    g0.resize(n_active_alleles, 0);
    g1.resize(n_active_alleles, 0);
    size_t active_al_i = 0;
    for (size_t al_i = 0; al_i < n_alleles; al_i++) {
        if (col_mask[al_i]) {continue;}
        for (size_t samp_i = 0; samp_i < n_samples; samp_i++) {
            if (row_mask[samp_i]) {continue;}
            double value_sample = get_value(samp_i, al_i);
            if (value_sample > 0){
                // tally the allele counts in each group
                if (b_phenotype->get_value_for_sample_id(samp_i)){
                    g1[active_al_i] += value_sample;
                } else {
                    g0[active_al_i] += value_sample;
                }
            }
        }
        active_al_i++;
    }
}

// JEAN maybe these won't be used in the end. Remove at the end if not
size_t GenotypeTable::get_n_active_alleles() const {
    return n_active_alleles;
}

size_t GenotypeTable::get_n_active_columns() const {
    return n_active_columns;
}

size_t GenotypeTable::get_n_active_samples() const {
    return n_active_samples;
}

void GenotypeTable::remove_one_allele() {
    // if there is only one allele, do nothing
    if (n_active_alleles < 2) {
        return;
    }

    // remove the first active allele
    for (size_t ii = 0; ii < n_alleles; ii++) {
        if (!col_mask[ii]) {
            col_mask[ii] = true;
            n_active_columns--;
            n_active_alleles--;
            break;
        }
    }
}

// potentially add a new column with the total allele count
// return true if it did, otherwise false
void GenotypeTable::add_total_allele_count_covariable() {
    // to check if at least one sample has a different total allele count
    use_total_ac = true;
    double first_tot_ac;
    // we've filled the new "covariate" holding the total allele count per sample when creating the Table
    // now we just check if we should include it, i.e. if some of the samples that we want to use have different total allele counts
    for (size_t samp_i = 0; samp_i < n_samples; samp_i++) {
        // skip if that sample is masked
        if (row_mask[samp_i]) {
            continue;
        }
        if (use_total_ac) {
            // first active sample, save the first value
            first_tot_ac = total_allele_counts_per_sample[samp_i];
            // and start checking for any differences
            use_total_ac = false;
        } else {
            // check if different from the first one
            if (total_allele_counts_per_sample[samp_i] != first_tot_ac) {
                // "add" the new column
                use_total_ac = true;
                n_active_columns++;
                break;
            }
        }
    }
}

bool GenotypeTable::passes_filters(const double maf, const size_t min_individuals) const {
    // make sure there is/was at least two alleles
    if (n_active_alleles < 2) {
        stoat::LOG_DEBUG("Filtered: less than 2 independent alleles.");
        return false;
    }

    // make sure there are enough individuals
    if(n_active_samples < min_individuals) {
        stoat::LOG_DEBUG("Filtered: not enough individuals: " + std::to_string(n_active_samples));
        return false;
    }

    // make sure the allele frequency is high enough
    if (maf == 0) {
        return true;
    }
    // first count the total number of observed alleles for each allele
    size_t total_allele_counts = 0;
    std::vector<size_t> allele_counts(n_alleles, 0);
    for (size_t al_i = 0; al_i < n_alleles; al_i++) {
        // skip masked alleles
        if (col_mask[al_i]) {
            continue;
        }
        // count all active samples support
        size_t allele_count_ii = 0;
        for (size_t samp_i = 0; samp_i < n_samples; samp_i++) {
            if (!row_mask[samp_i]) {
                allele_count_ii += get_value(samp_i, al_i);
            }
        }
        total_allele_counts += allele_count_ii;
        allele_counts[al_i] = allele_count_ii;
    }
    // compute the frequency of the second most frequent allele
    std::sort(allele_counts.begin(), allele_counts.end(), std::greater<size_t>());
    // JEAN faster way to keep track of the second most frequent allele?
    double af = static_cast<double>(allele_counts[1]) / total_allele_counts;
    // filter if too low
    if (af < maf) {
        stoat::LOG_DEBUG("Filtered: less than two alleles with frequency above " + std::to_string(maf));
        return false;
    }
    return true;
}

Eigen::MatrixXd GenotypeTable::make_matrixXd_features() {
    // we'll prepare a matrix with an additional "intercept" first column filled with ones
    size_t ncols = n_active_columns + 1;
    size_t nrows = n_active_samples;
    Eigen::MatrixXd X = Eigen::MatrixXd::Constant(nrows, ncols, 1);
    // fill other columns with the values all the predictors
    size_t cur_row;
    size_t cur_col = 0;
    for (size_t col_i = 0; col_i < n_alleles + n_covariates; ++col_i) {
        if (col_mask[col_i]) {
            continue;
        }
        cur_col++;
        cur_row = -1;
        for (size_t samp_i = 0; samp_i < n_samples; ++samp_i) {
            if (row_mask[samp_i]) {
                continue;
            }
            cur_row++;
            X(cur_row, cur_col) = get_value(samp_i, col_i);
        }
    }
    // add the total allele count as a last column (if necessary)
    if (use_total_ac) {
        cur_row = -1;
        for (size_t samp_i = 0; samp_i < n_samples; ++samp_i) {
            if (row_mask[samp_i]) {
                continue;
            }
            cur_row++;
            X(cur_row, X.cols() - 1) = total_allele_counts_per_sample[samp_i];
        }           
    }
    return X;
}

Eigen::VectorXd GenotypeTable::make_vectorxd_phenotype() {
    size_t nrows = n_active_samples;
    Eigen::VectorXd y(nrows);
    size_t cur_row = -1;
    // would be nice to directly cast the phenotype no matter if it's binary or quantitative
    if (is_binary_pheno) {
        for (size_t samp_i = 0; samp_i < n_samples; ++samp_i) {
            if (row_mask[samp_i]) {
                continue;
            }
            cur_row++;
            y(cur_row) = static_cast<double>(b_phenotype->get_value_for_sample_id(samp_i));
        }           
    } else {
        for (size_t samp_i = 0; samp_i < n_samples; ++samp_i) {
            if (row_mask[samp_i]) {
                continue;
            }
            cur_row++;
            y(cur_row) = q_phenotype->get_value_for_sample_id(samp_i);
        }           
    }
    return y;
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
