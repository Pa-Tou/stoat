#ifndef STOAT_FEATURE_TABLES_HPP
#define STOAT_FEATURE_TABLES_HPP

#include <vector>
#include <string>
#include <unordered_map>
#include "log.hpp"

namespace stoat {

/***

This file defines tables to represent per-sample values such as phenotypes, genotypes, gene expression, etc.

The base class is a FeatureBySampleTable, which is a class template to store a different type of value per sample.
This can be thought of a 1D matrix of values.
The classes of FeatureBySampleTable are:

-  BinaryPhenotypeTable: stores one bool per sample
-  QuantitativePhenotypeTable: stores one double per sample


The CategoricalFeatureBySampleTable inherits from the FeatureBySampleTable and provides a vector of values per sample.
This can be thought of as a 2D matrix. Each sample now has multiple values from a category of features. For example,
the category may be gene expression, and for each feature (gene), each sample has a gene expression value.
CategoricalFeatureBySampleTable is also a class template and its classes are:

- GeneExpressionTable: stores a vector of double's per sample
- CovariateTable: stores a vector of double's per sample


The GenotypeTable also extends FeatureBySampleTable, using a vector of size_t's like a CategoricalFeatureBySampleTable,
but unlike the Categorical table, it accesses alleles by index, rather than by a string name

- GenotypeTable: stores a vector of size_t's per sample, used as a count of alleles



Each of these tables has a const reference to an unordered map sample_to_index that maps the sample name to a unique index
used to place the sample in the vector.
CategoricalFeatureBySampleTable also has a const reference to feature_to_index that maps the feature name to a unique index.
Values are accessed by sample and feature names

***/


/// A generic table to store the value of a "feature" per sample. 
/// The feature could be a phenotype
/// This can be thought of as:
///                        sample 1  | sample 2 | sample 3 |
/// values_per_sample:      value 1  |  value 2 |  value 3 | 
template<class ValueType>
class FeatureBySampleTable {

    public:

    // Constructor
    // Note that using the generic constructor will fill in values_per_sample with undefined values. 
    FeatureBySampleTable(const std::unordered_map<std::string, size_t>& sample_to_index);

    // Access the value saved for a sample
    //TODO: I wanted this to be a reference in case the value is very big (eg, we want the whole vector) but it doesn't work with bools
    ValueType get_value_for_sample(const std::string& sample) const;

    double get_value_for_sample_as_double(const std::string& sample) const;
    
    // Set the value for the sample
    // JEAN if sample is not already present, add it?
    void set_value_for_sample(const std::string& sample, ValueType value);

    // is this sample in the table?
    bool has_sample(const std::string& sample) const;

    // return a vector with the names of all samples
    std::vector<std::string> get_sample_names() const;
    
    protected:
    // Map from the samples that we have features for to their index in values_per_sample
    const std::unordered_map<std::string, size_t>& sample_to_index;

    // The values of the feature per sample in sample_to_index
    std::vector<ValueType> values_per_sample;
};

// template<class ValueType>
// double FeatureBySampleTable<std::vector<ValueType>>::get_value_for_sample_as_double(const std::string& sample) const;

// Both phenotype tables assume that every sample has a phenotype
using BinaryPhenotypeTable = FeatureBySampleTable<bool>;
using QuantitativePhenotypeTable = FeatureBySampleTable<double>;

/// A CategoricalFeatureBySampleTable inherits from FeatureBySampleTable, but specifies that it is a 2D matrix with a vector of values per sample.
/// This is used to represent the value of each feature from a category of features. 
/// values_per_sample is now a vector of vectors, with each entry in the inner vector being the value of a different feature.
/// For example, for gene expression, the feature is the gene and the value is the expression level 
/// This can be thought of as:
///              |  feature 1        | feature 2         | feature 3 
///     -------------------------------------------------------------------
///     sample 1 | value samp1, feat1 | value samp1, feat3 | value samp1, feat3  
///     sample 2 | value samp2, feat1 | value samp2, feat2 | value samp2, feat3  
///     sample 3 | value samp3, feat1 | value samp3, feat3 | value samp3, feat3  
///
/// Where each row (corresponding to one sample) is one vector in values_per_sample
template<class ValueType>
class CategoricalFeatureBySampleTable : public FeatureBySampleTable<std::vector<ValueType>> {

    public:

    // Constructor
    // Note that using the generic constructor will fill in values_per_sample with undefined values
    // A derived class should specify the default value or make sure that all values are filled in.
    CategoricalFeatureBySampleTable(const std::unordered_map<std::string, size_t>& sample_to_index, const std::unordered_map<std::string, size_t>& feature_to_index);

    // Access the value saved for a sample and feature
    ValueType get_value_for_sample_and_feature(const std::string& sample, const std::string& feature) const;

    // Set the value for the sample and feature
    void set_value_for_sample_and_feature(const std::string& sample, const std::string& feature, ValueType value);

    // double get_value_for_sample_as_double(const std::string& sample) const;

    protected:

    // Map the name of the feature(what each value in the inner vector represents. e.g. gene name) to its index in the inner vector of values_per_sample 
    const std::unordered_map<std::string, size_t>& feature_to_index;

};
    
// Specialize a CategoricalFeatureBySampleTable constructor for doubles to fill in the matrix with std::numeric_limits<double>::max()
template<>
CategoricalFeatureBySampleTable<double>::CategoricalFeatureBySampleTable(const std::unordered_map<std::string, size_t>& sample_to_index, const std::unordered_map<std::string, size_t>& feature_to_index);
    
using GeneExpressionTable = CategoricalFeatureBySampleTable<double>;
using CovariateTable = CategoricalFeatureBySampleTable<double>;
    
/// A GenotypeTable is a 2D matrix of per-sample per-allele counts. The alleles are accessed by index, instead of by name
/// Technically this also has get_value_for_sample() which would return the vector but it shouldn't be used because it will return a copy of the vector 
class GenotypeTable : public FeatureBySampleTable<std::vector<size_t>> {
    public:

    GenotypeTable(const std::unordered_map<std::string, size_t>& sample_to_index, size_t allele_count);

    // Add 1 to the current value for sample and feature
    void increment_count(const std::string& sample, size_t allele_num);

    // Access the count value saved for a sample and allele
    size_t get_count_for_sample_and_allele(const std::string& sample, size_t allele_num) const;

    // Get the genotype of a sample as a string of counts.
    // This is only really useful to compare if two genotypes are the same
    std::string get_genotype_as_string(const std::string& sample) const;

    // How many alleles are there?
    size_t get_allele_count() const {
        return this->values_per_sample.size() == 0 ? 0 : this->values_per_sample.front().size(); 
    }
};

class CombinedTable {
public:
    CombinedTable(const GenotypeTable& genotypes);

    // combine a phenotype table with the current table (matching sample order etc)
    void combine_binary_phenotype(const BinaryPhenotypeTable& phenotype);
    // void combine_quantitative_phenotype(const QuantitativePhenotypeTable& phenotype);
    // void combine_gene_expression(const GeneExpressionTable& ge, std::string gene_name);

    // combine covariates
    // void combine_covariates(const CovariateTable& covariates);

    // remove predictors with the same values across all samples
    // (typically 0, i.e. alleles carried by no one)
    void remove_constant_predictors();

    // specific table operations used before a regression. Includes: adding a copy number
    // covariate if needed, merging duplicated predictors, removing the first allele
    // void prepare_for_regression();

    // should we test this table? does it pass the filters?
    bool passes_filters(const size_t maf, const size_t min_individuals) const;

    // get the (usually final) tables
    std::vector<double> get_phenotype() const;
    std::vector<std::vector<double>> get_predictors() const;

    // how many of the predictors are alleles? The others are covariates
    int get_n_alleles() const;
    
protected:
    // sample names in the current order
    std::vector<std::string> sample_names;
    // vector of phenotypes. For binary traits, 0 and 1 are used
    std::vector<double> phenotype_vec;
    // the allele counts + other covariates
    std::vector<std::vector<double>> predictors;
    // how many of the predictors are alleles (the rest will be covariates)
    int n_alleles;
    // which samples should be excluded from the matrix during filtering
    // (because some info is missing)
    std::unordered_map<std::string, bool> samples_to_exclude;
    
};

    
} // namespace stoat


#endif 
