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

// STOAT_VERSION are available in the CMAKELIST file
void print_help() {
    stoat::print_banner_v1(std::string(STOAT_VERSION));
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

// SNARL
// ./stoat vcf -g ../data/binary/pg.full.pg -d ../data/binary/pg.full.dist -r ../data/binary/pg.chromosome --output ../output_snarl

// BINARY
// ./stoat vcf -g ../data/binary/pg.full.pg -d ../data/binary/pg.full.dist -r ../data/binary/pg.chromosome -v ../data/binary/merged_output.vcf.gz -b ../data/binary/phenotype.tsv --output ../output_binary

// BINARY + COVARIATE
// ./stoat vcf -g ../data/binary/pg.full.pg -d ../data/binary/pg.full.dist -r ../data/binary/pg.chromosome -v ../data/binary/merged_output.vcf.gz -b ../data/binary/phenotype.tsv --covariate ../data/binary/covariate.tsv --covar-name PC1,SEX,PC3 --output ../output_binary_covar

// BINARY + VCF
// ./stoat vcf -g ../data/binary/pg.full.pg -d ../data/binary/pg.full.dist -r ../data/binary/pg.chromosome -v ../data/binary/merged_output.vcf.gz -b ../data/binary/phenotype.tsv --make-genotype vcf --output ../output_binary

// QUANTITATIVE
// ./stoat vcf -g ../data/quantitative/pg.full.pg -d ../data/quantitative/pg.full.dist -r ../data/quantitative/pg.chromosome -v ../data/quantitative/merged_output.vcf.gz -q ../data/quantitative/phenotype.tsv --output ../output_quantitative

// QUANTITATIVE + COVARIATE
// ./stoat vcf -g ../data/quantitative/pg.full.pg -d ../data/quantitative/pg.full.dist -r ../data/quantitative/pg.chromosome -v ../data/quantitative/merged_output.vcf.gz -q ../data/quantitative/phenotype.tsv  --covariate ../data/quantitative/covariate.tsv --covar-name PC1,SEX,PC3 --output ../output_quantitative_covar

// EQTL
// ./stoat vcf -s ../data/eqtl/snarl_analyse.tsv -v ../data/eqtl/merged_output.vcf.gz -e ../data/eqtl/qtl.tsv --gene-position ../data/eqtl/gene_position.tsv --output ../output_eqtl

// EQTL + COVARIATE
// ./stoat vcf -s ../data/eqtl/snarl_analyse.tsv -v ../data/eqtl/merged_output.vcf.gz -e ../data/eqtl/qtl.tsv --gene-position ../data/eqtl/gene_position.tsv --covariate ../data/eqtl/covariate.tsv --covar-name PC1,SEX,PC3 --output ../output_eqtl_covar

// SIMU TEST
// ./stoat vcf -g ../tests/graph_test/3th_snp.pg -d ../tests/graph_test/3th_snp.dist --output ../output

// BINARY-PLINK
// ./stoat vcf -g ../data/binary/pg.full.pg -d ../data/binary/pg.full.dist -v ../data/binary/merged_output.vcf.gz --make-bed --output ../output

// QUANTITATIVE-PLINK
// ./stoat vcf -g ../data/quantitative/pg.full.pg -d ../data/quantitative/pg.full.dist -v ../data/quantitative/merged_output.vcf.gz --make-bed --output ../output

// SIMULATION NEW
// ./stoat vcf -v ../data/simu/merged_output.vcf.gz -s ../data/simu/paths_snarl.tsv -b ../data/simu/phenotypes.txt --covariate ../data/simu/covar.tsv --covar-name AGE,SEX,PC1,PC2 --output ../output

// ./stoat vcf -v ../data/simu/merged_output.vcf.gz -s ../data/simu/paths_snarl.tsv -b ../data/simu/phenotypes.txt --make-bed --output ../output
// plink --bfile ../output/output --pheno ../data/simu/phenotypes.txt --pheno-name PHENO --assoc --allow-no-sex --allow-extra-chr --out ../output/stoat_plink

// -------------------------------------------------------------- GRAPH --------------------------------------------------------------

// BINARY
// ./stoat graph -g ../data/binary/pg.full.pg -d ../data/binary/pg.full.dist -b ../data/binary/phenotype_samples.tsv -T chi2 -r ref --output ../output_binary_graph
// echo "#CHR	START_POS	END_POS	SNARL	PATH_LENGTHS	P_FISHER	P_CHI2	GROUP_PATHS	DEPTH" > binary_table_graph.modify.tsv
// awk 'BEGIN{OFS=FS="\t"} !/^#/ {split($4, a, "_"); print a[1], $0}' binary_table_graph.tsv | sort -k1,1nr | cut -f2- >> binary_table_graph.modify.tsv
// mv binary_table_graph.modify.tsv binary_table_graph.tsv
//
// -------------------------------------------------------------- DECONSTRUCT --------------------------------------------------------------

// BINARY DECONSTRUCT
// ./stoat vcf -g ../data/binary/pg.full.pg -d ../data/binary/pg.full.dist -r ../data/binary/pg.chromosome -v ../data/binary/pg.full.deconstruct.change2.vcf -b ../data/binary/phenotype.tsv --output ../output_binary_deconstruct

// -------------------------------------------------------------- OTHER --------------------------------------------------------------

// PLINK
// plink --vcf ../data/simu/merged_output.vcf.gz --make-bed --allow-extra-chr --out ../output/genotype
// plink --bfile ../output/genotype --pheno ../data/simu/phenotypes.txt --pheno-name PHENO --assoc --allow-no-sex --allow-extra-chr --out ../output/plink

// DROSO
// ./stoat vcf -g ../../lab/droso/data/fly.pg -d ../../lab/droso/data/fly.dist -r ../../lab/droso/data/chromosome_ref.tsv --output ../output_droso
// sed -i 's/dm6#0#chr2L/1/g' ../output_droso/snarl_analyse.tsv
// sed -i 's/dm6#0#chr2R/2/g' ../output_droso/snarl_analyse.tsv
// sed -i 's/dm6#0#chr3L/3/g' ../output_droso/snarl_analyse.tsv
// sed -i 's/dm6#0#chr3R/4/g' ../output_droso/snarl_analyse.tsv
// sed -i 's/dm6#0#chr4/5/g' ../output_droso/snarl_analyse.tsv
// sed -i 's/dm6#0#chrX/6/g' ../output_droso/snarl_analyse.tsv
// sed -i 's/dm6#0#chrY/7/g' ../output_droso/snarl_analyse.tsv
// sed -i 's/dm6#0#chrM/8/g' ../output_droso/snarl_analyse.tsv
// ./stoat vcf -s ../output_droso/snarl_analyse.tsv -v ../../lab/droso/data/merging_stoat.vcf -q ../../lab/droso/data/dgrpool.dm6.male.phenotype.tsv --output ../output_droso

// -------------------------------------------------------------- DEBUG --------------------------------------------------------------

// VALGRIND
// valgrind --tool=callgrind ./stoat -s ../data/binary/snarl_paths.tsv -v ../data/binary/merged_output.vcf.gz -b ../data/binary/phenotype.tsv --output ../output
// kcachegrind callgrind.out.<id>

// GDB
// gdb --args ./stoat vcf -g ../data/binary/pg.full.pg -d ../data/binary/pg.full.dist -r ../data/binary/pg.chromosome -v ../data/binary/merged_output.vcf.gz -b ../data/binary/phenotype.tsv --output ../output_binary
// (gdb) run
