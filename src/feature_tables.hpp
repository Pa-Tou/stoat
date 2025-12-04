#ifndef STOAT_FEATURE_TABLES_HPP
#define STOAT_FEATURE_TABLES_HPP

#include <vector>
#include <string>
#include <unordered_map>

namespace stoat {


// A generic table to store the value of a "feature" per sample. 
// The feature could be a genotype, phenotype
// This can be thought of as:
//                        sample 1  | sample 2 | sample 3 |
// values_per_sample:      value 1  |  value 2 |  value 3 | 
template<class ValueType>
class FeatureBySampleTable {

    public:

    // Constructor
    FeatureBySampleTable(const std::unordered_map<std::string, size_t>& sample_to_index);

    // Access the value saved for a sample
    const ValueType& get_value_for_sample(const std::string& sample) const;

    // Set the value for the sample
    void set_value_for_sample(const std::string& sample, ValueType value);

    protected:
    // Map from the samples that we have features for to their index in values_per_sample
    const std::unordered_map<std::string, size_t>& sample_to_index;

    // The values of the feature per sample in sample_to_index
    std::vector<ValueType> values_per_sample;
};

// Both phenotype tables assume that every sample has a phenotype
using BinaryPhenotypeTable = FeatureBySampleTable<bool>;
using QuantitativePhenotypeTable = FeatureBySampleTable<double>;


// A table to represent the value for each feature from a category of features. 
// values_per_sample is now a vector of vectors, with each entry in the inner vector being the value of a different feature.
// For example, for a gentoype, the feature is the allele and the value is the count. For gene expression, the feature is the gene and the value is the expression level 
// This can be thought of as:
//              |  feature 1        | feature 2         | feature 3 
//     -------------------------------------------------------------------
//     sample 1 | value samp1, cat1 | value samp1, cat3 | value samp1, cat3  
//     sample 2 | value samp2, cat1 | value samp2, cat2 | value samp2, cat3  
//     sample 3 | value samp3, cat1 | value samp3, cat3 | value samp3, cat3  
//
// Where each sample row is one vector in values_per_sample
template<class ValueType>
class CategoricalFeatureBySampleTable : public FeatureBySampleTable<std::vector<ValueType>> {

    public:

    // Constructor just calls base constructor
    CategoricalFeatureBySampleTable(const std::unordered_map<std::string, size_t>& sample_to_index, const std::unordered_map<std::string, size_t>& feature_to_index);

    // Access the value saved for a sample and feature
    const ValueType& get_value_for_sample_and_feature(const std::string& sample, const std::string& feature) const;

    // Set the value for the sample and feature
    void set_value_for_sample_and_feature(const std::string& sample, const std::string& feature, ValueType value);

    protected:

    // Map the name of the feather(what each value in the inner vector represents. e.g. gene name) to its index in the inner vector of values_per_sample 
    const std::unordered_map<std::string, size_t>& feature_to_index;

};


// Each of these derived classes just defines the default value

class GenotypeTable : public CategoricalFeatureBySampleTable<size_t> {
    public:

    GenotypeTable(const std::unordered_map<std::string, size_t>& sample_to_index, const std::unordered_map<std::string, size_t>& feature_to_index);

            CategoricalFeatureBySampleTable<size_t>::get_value_for_sample_and_feature;
};

class GeneExpressionTable : public CategoricalFeatureBySampleTable<double> {
    public:

    GeneExpressionTable(const std::unordered_map<std::string, size_t>& sample_to_index, const std::unordered_map<std::string, size_t>& feature_to_index);
};
class CovariateTable : public CategoricalFeatureBySampleTable<double> {
    public:

    CovariateTable(const std::unordered_map<std::string, size_t>& sample_to_index, const std::unordered_map<std::string, size_t>& feature_to_index);
};

} // namespace stoat


#endif 
