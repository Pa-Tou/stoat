// This file is part of STOAT 0.0.1, copyright (C) 2024-2025 
// Authors : Matis Alias-Bagarre, Xian Hui Chang & Jean Monlong.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <string>
#include <unordered_map>
#include <chrono>
#include <Eigen/Dense>
#include <cstdlib>
#include <getopt.h>
#include <omp.h>

#include "subcommand/vcf.hpp"
#include "subcommand/graph.hpp"
#include "subcommand/test.hpp"
#include "subcommand/bh_correct.hpp"
#include "subcommand/change_reference.hpp"
#include "log.hpp"
#include "banner.hpp"

// STOAT_VERSION are define in the CMAKELIST file
void print_help() {
    stoat::print_banner(std::string(STOAT_VERSION));
    std::cerr   << "stoat: Snarl Tree Orchestrated Association Test, i.e GWAS on a pangenome's snarls.\n"
                << "usage: stoat <command> [options]\n\n"    
                << "  -- version       version information\n"
                << "\n"
                << "snarl genotyping:\n"
                << "  -- vcf           genotypes snarls using a VCF from a pangenome genotyping\n"
                << "  -- graph         genotypes snarls from the haplotypes in a pangenome graph\n"
                << "\n"
                << "association testing:\n"
                << "  -- test          performs the test between snarl genotypes and a phenotype\n"
                << "\n"
                << "post-processing:\n"
                << "  -- BHcorrect     apply the Benjamini-Hochberg procedure for multiple testing to a tsv file (DEPRECATED?)\n"
                << "  -- change-ref    change the reference coordinates of the output tsv to a new reference genome\n";     
}

int main(int argc, char* argv[]) {

    if (argc < 2) {
        print_help();
        return EXIT_FAILURE;
    }

    std::string subcommand = argv[1];
    stoat::LogLevel verbosity = stoat::LogLevel::Info;  // default level info

    // Shift argv to skip the subcommand itself
    argc -= 1;
    argv += 1;

    // Set the number of threads to 1 by default
    omp_set_num_threads(1);

    if (subcommand == "vcf") {
        stoat_command::main_stoat_vcf(argc, argv);

    } else if (subcommand == "graph") {
        stoat_command::main_stoat_graph(argc, argv);

    } else if (subcommand == "test") {
        stoat_command::main_stoat_test(argc, argv);

    } else if (subcommand == "BHcorrect") {
        stoat_command::main_stoat_bh_correct(argc, argv);

    } else if (subcommand == "change-ref") {
        stoat_command::main_stoat_change_reference(argc, argv);

    } else if (subcommand == "version") {
        std::cout << "stoat: GWAS analysis tool, version " << STOAT_VERSION;
        // stoat::LOG_INFO("Compiled with g++ (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0 on Linux)";
        // stoat::LOG_INFO("Linked against libstd++ 20230528)";

    } else {
        print_help();
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

// -------------------------------------------------------------- VCF --------------------------------------------------------------

// BINARY SNARL
// ./stoat vcf -g ../tests/test_data/input_data/binary/pg.full.pg -d ../tests/test_data/input_data/binary/pg.full.dist -R ../tests/test_data/input_data/binary/pg.chromosome --no-bgzip --output ../output_binary_snarl

// BINARY GENOTYPE
// ./stoat vcf -s ../output_binary_snarl/snarl_info.tsv -v ../tests/test_data/input_data/binary/merged_output.vcf.gz --no-bgzip --output ../output_binary_genotype

// BINARY TEST + Fisher/Chi-Squared
// ./stoat test -g ../output_binary_genotype/snarl_genotypes.tsv -p ../tests/test_data/input_data/binary/phenotype.tsv -m chi2 --no-bgzip --output ../output_binary_test

// BINARY + COVARIATE
// ./stoat vcf -g ../tests/test_data/input_data/binary/pg.full.pg -d ../tests/test_data/input_data/binary/pg.full.dist -r ../tests/test_data/input_data/binary/pg.chromosome -v ../tests/test_data/input_data/binary/merged_output.vcf.gz -b ../tests/test_data/input_data/binary/phenotype.tsv --covariate ../tests/test_data/input_data/binary/covariate.tsv --covar-name PC1,SEX,PC3 --output ../output_binary_covar

// BINARY + VCF
// ./stoat vcf -g ../tests/test_data/input_data/binary/pg.full.pg -d ../tests/test_data/input_data/binary/pg.full.dist -r ../tests/test_data/input_data/binary/pg.chromosome -v ../tests/test_data/input_data/binary/merged_output.vcf.gz -b ../tests/test_data/input_data/binary/phenotype.tsv --make-genotype vcf --output ../output_binary

// QUANTITATIVE
// ./stoat vcf -g ../tests/test_data/input_data/quantitative/pg.full.pg -d ../tests/test_data/input_data/quantitative/pg.full.dist -r ../tests/test_data/input_data/quantitative/pg.chromosome -v ../tests/test_data/input_data/quantitative/merged_output.vcf.gz -q ../tests/test_data/input_data/quantitative/phenotype.tsv --output ../output_quantitative

// QUANTITATIVE + COVARIATE
// ./stoat vcf -g ../tests/test_data/input_data/quantitative/pg.full.pg -d ../tests/test_data/input_data/quantitative/pg.full.dist -r ../tests/test_data/input_data/quantitative/pg.chromosome -v ../tests/test_data/input_data/quantitative/merged_output.vcf.gz -q ../tests/test_data/input_data/quantitative/phenotype.tsv  --covariate ../tests/test_data/input_data/quantitative/covariate.tsv --covar-name PC1,SEX,PC3 --output ../output_quantitative_covar

// EQTL
// ./stoat vcf -s ../tests/test_data/input_data/eqtl/snarl_info.tsv -v ../tests/test_data/input_data/eqtl/merged_output.vcf.gz -e ../tests/test_data/input_data/eqtl/qtl.tsv --gene-position ../tests/test_data/input_data/eqtl/gene_position.tsv --output ../output_eqtl

// EQTL + COVARIATE
// ./stoat vcf -s ../tests/test_data/input_data/eqtl/snarl_info.tsv -v ../tests/test_data/input_data/eqtl/merged_output.vcf.gz -e ../tests/test_data/input_data/eqtl/qtl.tsv --gene-position ../tests/test_data/input_data/eqtl/gene_position.tsv --covariate ../tests/test_data/input_data/eqtl/covariate.tsv --covar-name PC1,SEX,PC3 --output ../output_eqtl_covar

// SIMU TEST
// ./stoat vcf -g ../tests/graph_test/3th_snp.pg -d ../tests/graph_test/3th_snp.dist --output ../output

// BINARY-PLINK
// ./stoat vcf -g ../tests/test_data/input_data/binary/pg.full.pg -d ../tests/test_data/input_data/binary/pg.full.dist -v ../tests/test_data/input_data/binary/merged_output.vcf.gz --make-bed --output ../output

// QUANTITATIVE-PLINK
// ./stoat vcf -g ../tests/test_data/input_data/quantitative/pg.full.pg -d ../tests/test_data/input_data/quantitative/pg.full.dist -v ../tests/test_data/input_data/quantitative/merged_output.vcf.gz --make-bed --output ../output

// SIMULATION NEW
// ./stoat vcf -v ../tests/test_data/input_data/simu/merged_output.vcf.gz -s ../tests/test_data/input_data/simu/paths_snarl.tsv -b ../tests/test_data/input_data/simu/phenotypes.txt --covariate ../tests/test_data/input_data/simu/covar.tsv --covar-name AGE,SEX,PC1,PC2 --output ../output

// ./stoat vcf -v ../tests/test_data/input_data/simu/merged_output.vcf.gz -s ../tests/test_data/input_data/simu/paths_snarl.tsv -b ../tests/test_data/input_data/simu/phenotypes.txt --make-bed --output ../output
// plink --bfile ../output/output --pheno ../tests/test_data/input_data/simu/phenotypes.txt --pheno-name PHENO --assoc --allow-no-sex --allow-extra-chr --out ../output/stoat_plink

// -------------------------------------------------------------- GRAPH --------------------------------------------------------------

// BINARY
// ./stoat graph -g ../tests/test_data/input_data/binary/pg.full.pg -d ../tests/test_data/input_data/binary/pg.full.dist -r ref --output -u ../output_binary_graph
// ./stoat test -g ../output_binary_graph/snarl_genotypes.tsv -p ../tests/test_data/input_data/binary/phenotype_samples.tsv -m chi2 --no-bgzip --output ../output_binary_graph

// -------------------------------------------------------------- DECONSTRUCT --------------------------------------------------------------

// BINARY DECONSTRUCT
// ./stoat vcf -g ../tests/test_data/input_data/binary/pg.full.pg -d ../tests/test_data/input_data/binary/pg.full.dist -r ../tests/test_data/input_data/binary/pg.chromosome -v ../tests/test_data/input_data/binary/pg.full.deconstruct.change2.vcf -b ../tests/test_data/input_data/binary/phenotype.tsv --output ../output_binary_deconstruct

// -------------------------------------------------------------- OTHER --------------------------------------------------------------

// PLINK
// plink --vcf ../tests/test_data/input_data/simu/merged_output.vcf.gz --make-bed --allow-extra-chr --out ../output/genotype
// plink --bfile ../output/genotype --pheno ../tests/test_data/input_data/simu/phenotypes.txt --pheno-name PHENO --assoc --allow-no-sex --allow-extra-chr --out ../output/plink

// -------------------------------------------------------------- DEBUG --------------------------------------------------------------

// VALGRIND
// valgrind --tool=callgrind ./stoat -s ../tests/test_data/input_data/binary/snarl_paths.tsv -v ../tests/test_data/input_data/binary/merged_output.vcf.gz -b ../tests/test_data/input_data/binary/phenotype.tsv --output ../output
// kcachegrind callgrind.out.<id>

// GDB
// gdb --args ./stoat vcf -g ../tests/test_data/input_data/binary/pg.full.pg -d ../tests/test_data/input_data/binary/pg.full.dist -r ../tests/test_data/input_data/binary/pg.chromosome -v ../tests/test_data/input_data/binary/merged_output.vcf.gz -b ../tests/test_data/input_data/binary/phenotype.tsv --output ../output_binary
// (gdb) run

// BINARY SNARL
// gdb --args ./stoat vcf -g ../tests/test_data/input_data/binary/pg.full.pg -d ../tests/test_data/input_data/binary/pg.full.dist -r ../tests/test_data/input_data/binary/pg.chromosome --no-bgzip --output ../output_binary_snarl
// gdb --args ./stoat vcf -s ../output_binary_snarl/snarl_info.tsv -v ../tests/test_data/input_data/binary/merged_output.vcf.gz --no-bgzip --output ../output_binary_genotype
// gdb --args ./stoat test -g ../output_binary_genotype/snarl_genotypes.tsv -p ../tests/test_data/input_data/binary/sample_phenotype.tsv -m chi2 --no-bgzip --output ../output_binary_test
