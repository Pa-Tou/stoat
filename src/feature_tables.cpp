#include "feature_tables.hpp"
#include <limits>

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
        this->values_per_sample.emplace_back(feature_to_index.size());
    }
}
    

// Constructor for a genotype table fills everything in with a default value of 0 for the counts
GenotypeTable::GenotypeTable(const std::unordered_map<std::string, size_t>& sample_to_index, const std::unordered_map<std::string, size_t>& feature_to_index) :
    CategoricalFeatureBySampleTable<size_t>::CategoricalFeatureBySampleTable(sample_to_index, feature_to_index) {

    for (size_t i = 0 ; i < this->sample_to_index.size() ; i++) {
        for(size_t j = 0 ; j < this->values_per_sample[i].size() ; j++) {
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

    for (size_t i = 0 ; i < this->sample_to_index.size() ; i++) {
        for(size_t j = 0 ; j < this->values_per_sample[i].size() ; j++) {
            this->values_per_sample[i][j] = std::numeric_limits<double>::max();
        }
    }
}

template<class ValueType>
const ValueType& FeatureBySampleTable<ValueType>::get_value_for_sample(const std::string& sample) const {
    return this->values_per_sample.at(this->sample_to_index.at(sample));
}

template<class ValueType>
void FeatureBySampleTable<ValueType>::set_value_for_sample(const std::string& sample, ValueType value) {
    return this->values_per_sample.at(this->sample_to_index.at(sample)) = std::move(value);
}

template<class ValueType>
const ValueType& CategoricalFeatureBySampleTable<ValueType>::get_value_for_sample_and_feature(const std::string& sample, const std::string& feature) const {
    return this->values_per_sample.at(this->sample_to_index.at(sample)).at(this->feature_to_index.at(feature));
}

template<class ValueType>
void CategoricalFeatureBySampleTable<ValueType>::set_value_for_sample_and_feature(const std::string& sample, const std::string& feature, ValueType value) {
    return this->values_per_sample.at(this->sample_to_index.at(sample)).at(this->feature_to_index.at(feature)) = std::move(value);
}

}

// Apparently these definitions are supposed to be done here, after all the members are defined
//template class stoat::SampleFeatureTable<bool>;
//template class stoat::SampleFeatureTable<double>;


