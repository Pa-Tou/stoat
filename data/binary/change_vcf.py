import sys
import gzip

def open_vcf(filename):
    return gzip.open(filename, 'rt') if filename.endswith('.gz') else open(filename, 'r')

def merge_vcf_samples(input_vcf, output_vcf):
    with open_vcf(input_vcf) as f_in, open(output_vcf, 'w') as f_out:
        for line in f_in:
            line = line.rstrip('\n')
            
            # Header line
            if line.startswith("#CHROM"):
                parts = line.split('\t')
                fixed = parts[:9]  # Keep standard VCF columns
                samples = parts[9:]

                # Group samples two by two (_h0 and _h1) and strip suffix
                merged_sample_names = []
                for i in range(0, len(samples), 2):
                    base_name = samples[i].replace('_h0', '').replace('_h1', '')
                    merged_sample_names.append(base_name)

                f_out.write('\t'.join(fixed + merged_sample_names) + '\n')

            # Meta lines
            elif line.startswith("#"):
                f_out.write(line + '\n')

            # Variant lines
            else:
                parts = line.split('\t')
                fixed = parts[:9]
                genotypes = parts[9:]

                merged_genotypes = []
                for i in range(0, len(genotypes), 2):
                    gt0 = genotypes[i].split(':')[0]  # Get GT only
                    gt1 = genotypes[i+1].split(':')[0]

                    # Combine alleles from two haplotypes
                    combined_gt = f"{gt0}/{gt1}"

                    # If other fields exist, reconstruct full format field (you can extend this)
                    merged_genotypes.append(combined_gt)

                f_out.write('\t'.join(fixed + merged_genotypes) + '\n')

# Example usage
merge_vcf_samples("pg.full.deconstruct.change.vcf", "pg.full.deconstruct.change2.vcf")

# python3 change_vcf.py