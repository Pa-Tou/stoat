import argparse
import pandas as pd
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate covariate file with optional PCs and SEX.")

    parser.add_argument("-i", "--input", required=True, help="Input phenotype file (FID IID PHENO).")
    parser.add_argument("-o", "--output", default="covariates.txt", help="Output covariate file.")
    parser.add_argument("--pcs", type=int, default=3, help="Number of principal components to generate.")
    parser.add_argument("--sex", action="store_true", default=True, help="Include SEX covariate (0/1 randomly).")

    # Mutually exclusive group for phenotype type
    pheno_group = parser.add_mutually_exclusive_group(required=True)
    pheno_group.add_argument("-b", "--binary", action="store_true", help="Specify that the phenotype is binary (0/1).")
    pheno_group.add_argument("-q", "--quantitative", action="store_true", help="Specify that the phenotype is quantitative (continuous).")

    args = parser.parse_args()

    # Load phenotype
    df = pd.read_csv(args.input, sep="\t")

    # Start covariate dataframe
    covar_df = pd.DataFrame()
    covar_df['IID'] = df['IID']

    # Add SEX if requested
    if args.sex:
        covar_df['SEX'] = np.random.randint(0, 2, size=len(df))

    # Add PCs
    for i in range(1, args.pcs + 1):
        covar_df[f'PC{i}'] = np.random.normal(loc=0, scale=1, size=len(df))

    # Save output
    covar_df.to_csv(args.output, sep="\t", index=False)
    print(f"Covariate file saved to {args.output}")

# python3 tests/scripts/simulate_covar.py -i data/binary/phenotype.tsv -o data/binary/covariate.tsv -b 
# python3 tests/scripts/simulate_covar.py -i data/quantitative/phenotype.tsv -o data/quantitative/covariate.tsv -q 
# python3 tests/scripts/simulate_covar.py -i data/eqtl/phenotype.tsv -o data/eqtl/covariate.tsv -q
