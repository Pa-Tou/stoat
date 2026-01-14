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
    // // if sample not in the table, add it
    // // JEAN convenient but maybe not super efficient? Also had to remove the const state of sample_to_index
    // if (sample_to_index.count(sample) == 0) {
    //     sample_to_index[sample] = sample_to_index.size();
    //     values_per_sample.emplace_back(value);
    // } else {
    //     // otherwise set it
    //     values_per_sample[sample_to_index.at(sample)] = value;
    // }
    values_per_sample[sample_to_index.at(sample)] = value;
}

template<class ValueType>
bool FeatureBySampleTable<ValueType>::has_sample(const std::string& sample) const {
    return(sample_to_index.find(sample) != sample_to_index.end());
}
    
template<class ValueType>
std::vector<std::string> FeatureBySampleTable<ValueType>::get_sample_names() const {
    std::vector<std::string> output;
    for(const auto& samp_name_idx: sample_to_index) {
        output.emplace_back(samp_name_idx.first);
    }
    return output;
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

// template<class ValueType>
// double CategoricalFeatureBySampleTable<ValueType>::get_value_for_sample_as_double(const std::string& sample) const {
//     throw std::runtime_error("Not possible to call get_value_for_sample_as_double on a vector");
//     return 0;
// }

// template<class ValueType>
// double FeatureBySampleTable<ValueType>::get_value_for_sample_as_double(const std::string& sample) const {
//     return static_cast<double>(values_per_sample.at(sample_to_index.at(sample)));
// }

template<class ValueType>
double FeatureBySampleTable<ValueType>::get_value_for_sample_as_double(const std::string& sample) const {
    return 0;
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
            if (phenotype.get_value_for_sample_as_double(samp_name)){
                phenotype_vec.push_back(0);
            } else {
                phenotype_vec.push_back(1);
            }
        } else {
            // this sample is not in the matrix, mark for filtering
            samples_to_exclude[samp_name] = true;
        }
    }
}

void CombinedTable::remove_constant_predictors() {
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

bool CombinedTable::passes_filters(const size_t maf, const size_t min_individuals) const {
    // JEAN this assumes we've removed samples with no genotype information (or that they've never were included in that combined table)

    // make sure it has at least two alleles
    if (n_alleles < 2) {
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
    
} //end stoat namespace

// Apparently these definitions are supposed to be done here, after all the members are defined
template class stoat::FeatureBySampleTable<bool>;
template class stoat::FeatureBySampleTable<double>;
template class stoat::FeatureBySampleTable<std::vector<size_t>>;
template class stoat::FeatureBySampleTable<std::vector<double>>;
template class stoat::CategoricalFeatureBySampleTable<double>;
