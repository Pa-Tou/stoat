#!/bin/bash
#SBATCH -J simu
#SBATCH -o simu.out
#SBATCH -e simu.err
#SBATCH -t 6:00:00
#SBATCH -p unlimitq
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=12G
#SBATCH --mail-user=matis.alias-bagarre@inserm.fr
#SBATCH --mail-type=END,FAIL

module load bioinfo/Snakemake/7.32.4 bioinfo/vg/1.65.0 bioinfo/Bcftools/1.21 bioinfo/bgzip/1.18 devel/python/Python-3.12.4 

## run pipeline
snakemake --cores 16 --snakefile Snakefile_simulation --configfile config/config_eqtl.yaml
