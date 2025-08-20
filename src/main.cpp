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
#include "subcommand/bh_correct.hpp"

// Global variable
const std::string VERSION = "v0.0.2";

void print_help() {
    std::cerr   << "stoat: gwas analysis tool, version " << VERSION << "\n"
                << "usage: stoat <command> [options]\n\n"    
                << "main usage:\n"
                << "  -- vcf           gwas analysis base on vcf pangenome calling\n"
                << "  -- graph         gwas analysis base on pangenome graph\n"
                << "  -- version       version information\n"
                << "\n"
                << "post-processing:\n"
                << "  -- BHcorrect    apply the Benjamini-Hochberg procedure for multiple testing to a tsv file\n"
                << "                   (this already done by `stoat vcf` and `stoat graph` by default)\n";     
}

int main(int argc, char* argv[]) {

    if (argc < 2) {
        print_help();
        return EXIT_FAILURE;
    }

   std::string subcommand = argv[1];

    // Shift argv to skip the subcommand itself
    argc -= 1;
    argv += 1;

    // Set the number of threads to 1 by default
    omp_set_num_threads(1);

    if (subcommand == "vcf") {
        stoat_command::main_stoat(argc, argv);

    } else if (subcommand == "graph") {
        stoat_command::main_stoat_graph(argc, argv);

    } else if (subcommand == "BHcorrect") {
        stoat_command::main_stoat_bh_correct(argc, argv);

    } else if (subcommand == "version") {
        std::cout << "stoat: GWAS analysis tool, version " << VERSION;
        // stoat::LOG_INFO("Compiled with g++ (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0 on Linux)";
        // stoat::LOG_INFO("Linked against libstd++ 20230528)";

    } else {
        print_help();
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

// -------------------------------------------------------------- VCF --------------------------------------------------------------

// BINARY
// ./stoat vcf -p ../data/binary/pg.full.pg -d ../data/binary/pg.full.dist -r ../data/binary/pg.chromosome -v ../data/binary/merged_output.vcf.gz -b ../data/binary/phenotype.tsv --output ../output

// BINARY + COVARIATE
// ./stoat vcf -p ../data/binary/pg.full.pg -d ../data/binary/pg.full.dist -r ../data/binary/pg.chromosome -v ../data/binary/merged_output.vcf.gz -b ../data/binary/phenotype.tsv --covariate ../data/binary/covariate.tsv --covar-name PC1,SEX,PC3 --output ../output

// QUANTITATIVE
// ./stoat vcf -p ../data/quantitative/pg.full.pg -d ../data/quantitative/pg.full.dist -r ../data/quantitative/pg.chromosome -v ../data/quantitative/merged_output.vcf.gz -q ../data/quantitative/phenotype.tsv --output ../output

// QUANTITATIVE + COVARIATE
// ./stoat vcf -p ../data/quantitative/pg.full.pg -d ../data/quantitative/pg.full.dist -r ../data/quantitative/pg.chromosome -v ../data/quantitative/merged_output.vcf.gz -q ../data/quantitative/phenotype.tsv  --covariate ../data/quantitative/covariate.tsv --covar-name PC1,SEX,PC3 --output ../output

// EQTL
// ./stoat vcf -s ../data/eqtl/snarl_analyse.tsv -v ../data/eqtl/merged_output.vcf.gz -e ../data/eqtl/qtl.tsv --gene-position ../data/eqtl/gene_position.tsv --output ../output

// EQTL + COVARIATE
// ./stoat vcf -s ../data/eqtl/snarl_analyse.tsv -v ../data/eqtl/merged_output.vcf.gz -e ../data/eqtl/qtl.tsv --gene-position ../data/eqtl/gene_position.tsv --covariate ../data/eqtl/covariate.tsv --covar-name SEX,AGE,PC1,PC2,PC3 --output ../output

// SIMU TEST
// ./stoat vcf -p ../tests/graph_test/3th_snp.pg -d ../tests/graph_test/3th_snp.dist --output ../output

// BINARY-PLINK
// ./stoat vcf -p ../data/binary/pg.full.pg -d ../data/binary/pg.full.dist -v ../data/binary/merged_output.vcf.gz --make-bed --output ../output

// QUANTITATIVE-PLINK
// ./stoat vcf -p ../data/quantitative/pg.full.pg -d ../data/quantitative/pg.full.dist -v ../data/quantitative/merged_output.vcf.gz --make-bed --output ../output

// SIMULATION NEW
// ./stoat vcf -v ../data/simu/merged_output.vcf.gz -s ../data/simu/paths_snarl.tsv -b ../data/simu/phenotypes.txt --covariate ../data/simu/covar.tsv --covar-name AGE,SEX,PC1,PC2 --output ../output

// ./stoat vcf -v ../data/simu/merged_output.vcf.gz -s ../data/simu/paths_snarl.tsv -b ../data/simu/phenotypes.txt --make-bed --output ../output
// plink --bfile ../output/output --pheno ../data/simu/phenotypes.txt --pheno-name PHENO --assoc --allow-no-sex --allow-extra-chr --out ../output/stoat_plink

// -------------------------------------------------------------- GRAPH --------------------------------------------------------------

// BINARY
// ./stoat graph -g ../data/binary/pg.full.pg -d ../data/binary/pg.full.dist -T chi2 -r ref -S ../data/binary/samples.g0.tsv -o ../output

// -------------------------------------------------------------- OTHER --------------------------------------------------------------

// PLINK
// plink --vcf ../data/simu/merged_output.vcf.gz --make-bed --allow-extra-chr --out ../output/genotype
// plink --bfile ../output/genotype --pheno ../data/simu/phenotypes.txt --pheno-name PHENO --assoc --allow-no-sex --allow-extra-chr --out ../output/plink

// DROSO
// ./stoat -p ../data_droso/fly.pg -d ../data_droso/fly.dist -v ../data_droso/merged.vcf -q ../data_droso/phenotype.tsv --output ../output_droso
   
// DROSO
// ./stoat -p ../data/droso/fly.pg -d ../data/droso/fly.dist -r ../data/droso/chromosome_ref.tsv --output ../output_droso
// sed -i 's/dm6#0#chr2L/1/g' ../output_droso/snarl_analyse.tsv
// sed -i 's/dm6#0#chr2R/2/g' ../output_droso/snarl_analyse.tsv
// sed -i 's/dm6#0#chr3L/3/g' ../output_droso/snarl_analyse.tsv
// sed -i 's/dm6#0#chr3R/4/g' ../output_droso/snarl_analyse.tsv
// sed -i 's/dm6#0#chr4/5/g' ../output_droso/snarl_analyse.tsv
// sed -i 's/dm6#0#chrX/6/g' ../output_droso/snarl_analyse.tsv
// sed -i 's/dm6#0#chrY/7/g' ../output_droso/snarl_analyse.tsv
// sed -i 's/dm6#0#chrM/8/g' ../output_droso/snarl_analyse.tsv
// ./stoat -s ../output_droso/snarl_analyse.tsv -v ../data/droso/merging_stoat.vcf -q ../data/droso/pangenome_pheno.tsv --output ../output_droso

// VALGRIND
// valgrind --tool=callgrind ./stoat -s ../data/binary/snarl_paths.tsv -v ../data/binary/merged_output.vcf.gz -b ../data/binary/phenotype.tsv --output ../output
// kcachegrind callgrind.out.<id>
