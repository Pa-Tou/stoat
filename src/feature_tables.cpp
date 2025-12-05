#include "feature_tables.hpp"
#include <limits>
#include <iostream>

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
    

// Constructor for a genotype table fills everything in with a default value of 0 for the counts
GenotypeTable::GenotypeTable(const std::unordered_map<std::string, size_t>& sample_to_index, size_t allele_count) :
    FeatureBySampleTable<std::vector<size_t>>::FeatureBySampleTable(sample_to_index) {

    for (size_t i = 0 ; i < this->sample_to_index.size() ; i++) {
        this->values_per_sample[i] = std::vector<size_t>(allele_count);
        for(size_t j = 0 ; j < allele_count ; j++) {
            this->values_per_sample[i][j] = 0;
        }
    }
}

// Constructor for a gene expression table fills everything in with a default value of std::numeric_limits<double>::max()
GeneExpressionTable::GeneExpressionTable(const std::unordered_map<std::string, size_t>& sample_to_index, const std::unordered_map<std::string, size_t>& feature_to_index) :
    CategoricalFeatureBySampleTable<double>::CategoricalFeatureBySampleTable(sample_to_index, feature_to_index) {

    for (size_t i = 0 ; i < this->sample_to_index.size() ; i++) {
        for(size_t j = 0 ; j < this->values_per_sample[i].size() ; j++) {
            this->values_per_sample[i][j] = std::numeric_limits<double>::max();
        }
    }
}
// Constructor for a covariate table fills everything in with a default value of std::numeric_limits<double>::max()
CovariateTable::CovariateTable(const std::unordered_map<std::string, size_t>& sample_to_index, const std::unordered_map<std::string, size_t>& feature_to_index) :
    CategoricalFeatureBySampleTable<double>::CategoricalFeatureBySampleTable(sample_to_index, feature_to_index) {

    for (size_t i = 0 ; i < sample_to_index.size() ; i++) {
        for(size_t j = 0 ; j < values_per_sample[i].size() ; j++) {
            values_per_sample[i][j] = std::numeric_limits<double>::max();
        }
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

void GenotypeTable::increment_count(const std::string& sample, size_t allele_num) {
    this->values_per_sample[this->sample_to_index.at(sample)][allele_num]++;
}
size_t GenotypeTable::get_count_for_sample_and_allele(const std::string& sample, size_t allele_num) const {
    return this->values_per_sample[this->sample_to_index.at(sample)][allele_num];
}


} //end stoat namespace

// Apparently these definitions are supposed to be done here, after all the members are defined
template class stoat::FeatureBySampleTable<bool>;
template class stoat::FeatureBySampleTable<double>;
template class stoat::CategoricalFeatureBySampleTable<size_t>;
template class stoat::CategoricalFeatureBySampleTable<double>;

