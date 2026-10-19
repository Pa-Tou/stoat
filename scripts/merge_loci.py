import argparse
import math
import sys

# Do the allele counts match exactly the genotype groups?
def is_perfectly_split(count_str):
    allele_counts=count_str.split(',')
    total_allele_count = 0
    for allele_count in allele_counts:
        counts = allele_count.split(':')
        nonzero_count = 0
        for count in counts:
            if int(count) > 0:
                nonzero_count += 1
                total_allele_count += int(count)
        if nonzero_count > 1:
            return False
    return total_allele_count >= 12



def print_loci(infile_name, gap_limit, only_perfect_split, p_cutoff):

    # Get a sorted list of variants as list of <chr, start, end>
    variants = []
    infile=sys.stdin
    if infile_name != "-":
        infile = open(infile_name, "r")
    p_index = 7
    for line in infile:
        if line[0] == "#":
            l = line.split()
            for i in range(len(l)):
                if line[i] == "P" or line[i] == "P_CHI2":
                    p_index = i
            continue
        l = line.split()
        if float(l[p_index]) < p_cutoff and (not only_perfect_split or is_perfectly_split(l[8])):
            variants.append((l[0], int(l[1]), int(l[2])))

    # Now sort the variants, by start pos then stably by chr
    variants.sort(key = lambda x : x[1])
    variants.sort(key = lambda x : x[0])

    # Run through the variants
    current_chr = ""
    current_start = math.inf
    current_end = 0
    for (next_chr, next_start, next_end) in variants:
        if next_chr == "NA":
            continue
        if current_chr == "":
            current_chr = next_chr
            current_start = next_start
            current_end = next_end
        is_close = abs(next_start - current_start) < gap_limit or abs(next_start - current_end) < gap_limit or \
                  abs(next_end - current_start) < gap_limit or abs(next_end - current_end) < gap_limit
        if next_chr != current_chr or not is_close:
            print(current_chr + ":" + str(current_start) + "-" + str(current_end))
            current_chr = next_chr
            current_start = next_start
            current_end = next_end
        else:
            current_end = max(current_end, next_end)
    print(current_chr + ":" + str(current_start) + "-" + str(current_end))
    infile.close()


def main():
    parser = argparse.ArgumentParser(description="given a stoat.assoc.pvalues.tsv file sorted by chr and position, merge nearby variants and print loci")
    parser.add_argument('input_file', help="stoat.assoc file name")
    parser.add_argument('--only-perfect-split', action='store_true', help = "Only include variants for which ALLELE_COUNT_PER_PHENO allele counts are divided perfectly among phenotype groups")
    parser.add_argument('--p-cutoff', type=float, default = math.inf, help = "Only include variants for which the p-value is less than this")
    parser.add_argument('--gap-limit', type=int, default = 100000, help = "Only include variants for which the p-value is less than this")
    args = parser.parse_args()

    print_loci(args.input_file, args.gap_limit, args.only_perfect_split, args.p_cutoff)




if __name__ == "__main__":
    main()

