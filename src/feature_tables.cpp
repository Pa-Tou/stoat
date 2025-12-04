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

    FeatureBySampleTable<std::vector<ValueType>>::values_per_sample.reserve(sample_to_index.size());
    for (size_t i = 0 ; i < sample_to_index.size() ; i++) {
        FeatureBySampleTable<std::vector<ValueType>>::values_per_sample.emplace_back(feature_to_index.size());
    }
}
    

// Constructor for a genotype table fills everything in with a default value of 0 for the counts
GenotypeTable::GenotypeTable(const std::unordered_map<std::string, size_t>& sample_to_index, const std::unordered_map<std::string, size_t>& feature_to_index) :
    CategoricalFeatureBySampleTable<size_t>::CategoricalFeatureBySampleTable(sample_to_index, feature_to_index) {

    for (size_t i = 0 ; i < sample_to_index.size() ; i++) {
        for(size_t j = 0 ; j < FeatureBySampleTable<std::vector<size_t>>::values_per_sample[i].size() ; j++) {
            FeatureBySampleTable<std::vector<size_t>>::values_per_sample[i][j] = 0;
        }
    }
}

GeneExpressionTable::GeneExpressionTable(const std::unordered_map<std::string, size_t>& sample_to_index, const std::unordered_map<std::string, size_t>& feature_to_index) :
    CategoricalFeatureBySampleTable<double>::CategoricalFeatureBySampleTable(sample_to_index, feature_to_index) {

    for (size_t i = 0 ; i < sample_to_index.size() ; i++) {
        for(size_t j = 0 ; j < FeatureBySampleTable<std::vector<double>>::values_per_sample[i].size() ; j++) {
            FeatureBySampleTable<std::vector<double>>::values_per_sample[i][j] = std::numeric_limits<double>::max();
        }
    }
}
CovariateTable::CovariateTable(const std::unordered_map<std::string, size_t>& sample_to_index, const std::unordered_map<std::string, size_t>& feature_to_index) :
    CategoricalFeatureBySampleTable<double>::CategoricalFeatureBySampleTable(sample_to_index, feature_to_index) {

    for (size_t i = 0 ; i < sample_to_index.size() ; i++) {
        for(size_t j = 0 ; j < FeatureBySampleTable<std::vector<double>>::values_per_sample[i].size() ; j++) {
            FeatureBySampleTable<std::vector<double>>::values_per_sample[i][j] = std::numeric_limits<double>::max();
        }
    }
}

}

// Apparently these definitions are supposed to be done here, after all the members are defined
//template class stoat::SampleFeatureTable<bool>;
//template class stoat::SampleFeatureTable<double>;


