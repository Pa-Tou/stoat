import random
import argparse
from scipy.special import expit  # sigmoid function
import numpy as np

# -*- coding: utf-8 -*-

# Written by Jean Monlong (jean.monlong@inserm.fr) modified by Matis Alias-Bagarre
# Simulate a simple pangenome graph with SNPs and indels. Some indels can have
# nested SNPs too. See below for the parameters used.

# Briefly, the pangenome is first build, one variant site at a time. Then, it
# is traversed by one (random) path that will be used as reference path. A
# first GFA file can be written at this point. At each beginning "fork" of a
# variant site, we had set random probability to follow one edge or the other.
# To create, difference between two groups of samples, we then select some
# variants and slightly shift those probabilities for the second group. Then
# we simulate paths for each haplotype of each sample and each group following
# those probabilities. A GFA file including all of those haplotype paths can
# also be outputted (to be used later to simulate reads from each sample)

# Limitation: simple pangenome with only "biallelic" snarls (two out-edges max)
# and only one level of nesting.

def randSeq(length):
    """Simulate a random nucleotide sequence"""
    return ''.join(random.choices(_nuc, k=length))

class Graph:
    def __init__(self, phenotype_type_, sex_arg=None):
        # phenotype type
        self.phenotype_type = phenotype_type_
        # map node name to sequence
        self.nodes = {}
        # map predecessor nodes to their successor nodes
        self.edges = {}
        # map path names to list of nodes
        self.paths = {}
        # map a snarl start to frequencies
        self.snarls_freq = {}
        # phenotypes
        self.phenotypes_value = []
        self.phenotypes_name = []
        # eqtl phenotypes
        self.eqtl_phenotypes = []
        # covariates
        self.covariates = []
        # keep track of a few things
        self.next_node_id = 1
        self.ngroups = 2
        self.sex = sex_arg

    def binaryPhenotype(self, nsamp):
        """Simulate a binary phenotype"""
        samp_t = nsamp//2
        for samp in range(0, nsamp):
            if samp < samp_t:
                self.phenotypes_value.append(1)
                self.phenotypes_name.append(f"samp_{samp}_g{0}")
            else:
                self.phenotypes_value.append(2)
                self.phenotypes_name.append(f"samp_{samp}_g{1}")

    def quantitativePhenotype(self, nsamp):
        """Simulate a quantitative phenotype"""
        samp_t = nsamp//2
        for samp in range(0, nsamp):
            if samp < samp_t:
                self.phenotypes_value.append(random.uniform(0.0, 1.0))
                self.phenotypes_name.append(f"samp_{samp}_g{0}")
            else:
                self.phenotypes_value.append(random.uniform(-1.0, 0.0))
                self.phenotypes_name.append(f"samp_{samp}_g{1}")

    def eqtlPhenotype(self, nsamp, gen_prob=0.1, ngene=100):
        """Simulate an eQTL phenotype"""
        for _ in range(ngene):
            gene_expr = random.uniform(0, 8.0)
            qtl_col = []
            gen_signi = random.random()
            bool_gene_signi = gen_signi <= gen_prob
            for samp in range(0, nsamp):
                if bool_gene_signi: # significant gene
                    qtl_col.append(gene_expr + (gene_expr * self.phenotypes_value[samp]))
                else: # not significant gene
                    qtl_col.append(gene_expr + (gene_expr * random.uniform(-1.0, 1.0)))
            self.eqtl_phenotypes.append(qtl_col)

    def covariate(self, nsamp, ncov) -> list[float]:
        """Simulate quantitative covariates associated with phenotype"""
        covar_effect = [random.uniform(0.1, 0.5) for _ in range(ncov)]  # strength of association
        for idx_cov in range(ncov):
            beta = covar_effect[idx_cov]
            cov = []
            for idx in range(nsamp):
                pheno = self.phenotypes_value[idx]
                noise = np.random.normal(0, 1)  # add noise to reduce perfect correlation
                cov_value = beta * pheno + noise
                cov.append(cov_value)
            self.covariates.append(cov)
        print("covar_effect : ", covar_effect)
    
    def covariate_sex(self, nsamp, sex_bias=0.0):
        """Simulate a sex covariate associated with phenotype"""
        cov = []
        beta_sex = np.random.uniform(-2, 2)  # random strength of association
        for idx in range(nsamp):
            pheno = self.phenotypes_value[idx]
            p = expit(beta_sex * pheno + sex_bias)  # logistic function
            sex = np.random.binomial(1, p)
            cov.append(sex)
        self.covariates.append(cov)
        print("beta_sex : ", beta_sex)

    def addNode(self, pred_nodes=[], min_size=50, max_size=300):
        # create a node and connect it to specified predecessors
        seq = randSeq(random.randint(min_size, max_size))
        # add node
        new_node = self.next_node_id
        self.nodes[new_node] = seq
        # add edges
        for pnod in pred_nodes:
            self.addEdge(pnod, new_node)
        # update next node ID
        self.next_node_id += 1
        return (new_node)

    def addEdge(self, pnode, snode):
        if pnode not in self.edges:
            self.edges[pnode] = {}
        self.edges[pnode][snode] = True

    def addSNP(self, pred_node):
        # pick two different alleles for the SNP
        als = random.sample(_nuc, 2)
        # create nodes and edges to predecessors
        node_a = self.next_node_id
        self.nodes[node_a] = als[0]
        self.next_node_id += 1
        node_b = self.next_node_id
        self.nodes[node_b] = als[1]
        self.next_node_id += 1
        # add edges
        self.addEdge(pred_node, node_a)
        self.addEdge(pred_node, node_b)
        # create a successor node
        suc_node = self.addNode([node_a, node_b])
        # init snarl frequency (same in all groups)
        freqs = {}
        freq_init = random.random()
        freqs[node_a] = [freq_init] * self.ngroups
        freqs[node_b] = [1 - freq_init] * self.ngroups
        self.snarls_freq[pred_node] = [freqs, 'N', 'N']
        return (suc_node)

    def addSNP_Npath(self, pred_node, npath):
        assert npath != 3 or npath != 4, "Error npath must be 3 or 4"
        # pick two different alleles for the SNP
        list_nucleotide = random.sample(_nuc, npath)
        # create nodes and edges to predecessors
        node_list = []
        for nuc in list_nucleotide :
            node_n = self.next_node_id
            self.nodes[node_n] = nuc
            self.next_node_id += 1
            node_list.append(node_n)

            # add edges
            self.addEdge(pred_node, node_n)

        # create a successor node
        suc_node = self.addNode(node_n)
        # init snarl frequency (same in all groups)
        freqs = {}
        freq_init = random.random()

        pathn_freq = 0.05 if npath == 3 else 0.025
        freqs[node_n[0]] = [freq_init-pathn_freq] * self.ngroups
        freqs[node_n[1]] = [1 - freq_init-pathn_freq] * self.ngroups
        for idex_node in range(len(node_list-2)) :
            freqs[node_list[idex_node]] = [pathn_freq*2] * self.ngroups

        self.snarls_freq[pred_node] = [freqs, 'N', 'N']
        return (suc_node)
    
    def addIndel(self, pred_node, snp_prob=.5):
        # create a middle node
        mid_node = self.addNode([pred_node])
        final_mid_node = mid_node
        # create a nested SNP sometimes
        while random.random() < snp_prob:
            final_mid_node = self.addSNP(final_mid_node)
        # create a successor node
        suc_node = self.addNode([final_mid_node])
        # add edge from predecessor directly to successor
        self.addEdge(pred_node, suc_node)
        # init snarl frequency (same in all groups)
        freqs = {}
        freq_init = random.random()
        freqs[mid_node] = [freq_init] * self.ngroups
        freqs[suc_node] = [1 - freq_init] * self.ngroups
        self.snarls_freq[pred_node] = [freqs, 'N', 'N']
        return (suc_node)

    def addPath(self, path_name='ref', group=0, samp=0):
        # start at node 1
        path = [1]
        while path[-1] in self.edges:
            # look for successor nodes
            snodes = self.edges[path[-1]]
            snodes = list(snodes.keys())
            # if there are some, pick one randomly
            if len(snodes) > 0:
                if path[-1] in self.snarls_freq:
                    # pick next node based on group frequencies
                    freqs = self.snarls_freq[path[-1]][0]
                    tot_freq = 0
                    rr = random.random()
                    # pick next node based on phenotype
                    if self.phenotype_type == 'binary':
                        # choose path base on group
                        for snode in freqs:
                            tot_freq += freqs[snode][group]
                            if rr < tot_freq:
                                break
                        path.append(snode)

                    else : # quantitative or eqtl
                        sample_pheno = (self.phenotypes_value[samp] / 2) + 0.5  # [0,1]
                        for snode in freqs: # loop only 2 times 
                            if rr < sample_pheno:
                                break
                        path.append(snode)
                else:
                    path.append(random.sample(snodes, 1)[0])
        self.paths[path_name] = path

    def createMarkerFreq(self, prop_markers, sex_markers=None, min_dev=.01, max_dev=.5):
        # loop over snarls and change a group's frequency sometimes
        for pnode in self.snarls_freq:
            bool_prop = False
            bool_sex = False
            if random.random() < prop_markers:
                bool_prop = True
            elif sex_markers and random.random() < sex_markers:
                bool_sex = True
            else : 
                continue
            # pick a group to deviate
            gp = random.randint(1, self.ngroups) - 1
            # pick a deviation size
            dev = random.uniform(min_dev, max_dev)
            # apply deviation
            dev_added = False
            for snode in self.snarls_freq[pnode][0]:
                self.snarls_freq[pnode][1] = "Y" if bool_sex else "N"
                self.snarls_freq[pnode][2] = "Y" if bool_prop else "N"
                if self.snarls_freq[pnode][0][snode][gp] < 1 - dev and not dev_added:
                    self.snarls_freq[pnode][0][snode][gp] += dev
                    dev_added = True
                else:
                    self.snarls_freq[pnode][0][snode][gp] -= dev

    def writeGfa(self, out_fn):
        outf = open(out_fn, 'wt')
        outf.write('H\tVN:Z:1.1\tRS:Z:ref\n')
        for nod in self.nodes:
            outf.write('S\t{}\t{}\n'.format(nod, self.nodes[nod]))
        for pnode in self.edges:
            for snode in self.edges[pnode]:
                outf.write('L\t{}\t+\t{}\t+\t0M\n'.format(pnode, snode))
        for pathn in self.paths:
            path = ','.join([str(nn) + '+' for nn in self.paths[pathn]])
            outf.write('P\t{}\t{}\t*\n'.format(pathn, path))

    def writeSnarlsFreq(self, out_fn, sex_cov):
        with open(out_fn, 'wt') as outf:
            if sex_cov:
                outf.write('start_node\tnext_node\tgroup\tfreq\tsex\tpheno\n')
                line_fmt = '{}\t{}\t{}\t{}\t{}\t{}\n'
            else:
                outf.write('start_node\tnext_node\tgroup\tfreq\n')
                line_fmt = '{}\t{}\t{}\t{}\n'

            for pnode, (snode_dict, *extra) in self.snarls_freq.items():
                sex = extra[0] if sex_cov else None
                pheno = extra[1] if sex_cov else None

                for snode, freqs in snode_dict.items():
                    for group, freq in enumerate(freqs):
                        freq_rounded = round(freq, 3)
                        if sex_cov:
                            outf.write(line_fmt.format(pnode, snode, group, freq_rounded, sex, pheno))
                        else:
                            outf.write(line_fmt.format(pnode, snode, group, freq_rounded))

    def writePhenotype(self, out_ph):
        outf = open(out_ph, 'wt')
        outf.write('FID\tIID\tPHENO\n')
        for samp_name, samp_val in zip(self.phenotypes_name, self.phenotypes_value):
            # write sample FID, ID and phenotype
            outf.write(f"{samp_name}\t{samp_name}\t{samp_val}\n")

    def writeCovariate(self, ncov:int, out_cov:str):
        outf = open(out_cov, 'wt')
        num_pc = '\t'.join([f"PC{i+1}" for i in range(0, ncov)])
        outf.write(f"FID\tIID\t{num_pc}\n")
        # write covariates
        for idx in range(len(self.covariates[0])):
            # write sample FID, ID and covariate
            cov_col = '\t'.join([str(cov[idx]) for cov in self.covariates])
            outf.write(f"{self.phenotypes_name[idx]}\t{self.phenotypes_name[idx]}\t{cov_col}\n")

    def writeEqtl(self, out_eqtl):
        outf = open(out_eqtl, 'wt')
        sample_names = '\t'.join([f'samp_{i}' for i in range(len(self.phenotypes_value))])
        outf.write(f"gene_name\t{sample_names}\n")   
        for idx, samp in enumerate(self.eqtl_phenotypes):
            # write sample FID, ID and eQTL phenotype
            eqtl_col = '\t'.join([str(gene) for gene in samp])
            outf.write(f'gene_{idx}\t{eqtl_col}\n')

    def writeGenePosition(self, out_gp):
        outf = open(out_gp, 'wt')
        outf.write('gene_name\tchr\tstart\tend\n')
        for idx in range(len(self.phenotypes_value)):
            # write gene name and position
            outf.write(f'gene_{idx}\tref\t{idx * 100}\t{(idx * 100)+10000}\n')

if "__main__" == __name__ :

    # parse command line arguments
    parser = argparse.ArgumentParser(description='Simulate a simple pangenome graph with SNPs and indels.')
    parser.add_argument('-v', '--nvar', type=int, default=1000, help='Number of top-level variants (default: 1000)')
    parser.add_argument('-n', '--nsamp', type=int, default=200, help='Number of samples in total (default: 100)')
    parser.add_argument('--snp_prop', type=float, default=0.7, help='Proportion of top-level variants that are SNPs (default: 0.7)')
    parser.add_argument('--nested_prop', type=float, default=0.8, help='Proportion of indels that we want to try to add nested SNPs into (default: 0.8)')
    parser.add_argument('--sv_prop', type=float, default=0.1, help='Proportion of SVs (default: 0.1)')
    parser.add_argument('--prop_markers', type=float, default=0.1, help='Proportion of markers (default: 0.1)')
    parser.add_argument('-c', '--ncov', type=int, default=3, help='Number of covariate (default: 3)')
    parser.add_argument('-s', '--sex_cov', action='store_true', default=True, help='Include sex as a covariate (default: False)')
    parser.add_argument('--sex_eq', type=float, default=0.1, help='Sexual group desiquilibre (default: 0.1)')
    parser.add_argument('--sex_markers', type=float, default=0.15, help='Proportion of sexual markers (default: 0.15)')
    parser.add_argument('-g', '--ngene', type=int, default=100, help='Number of gene (default: 100) for eQTL phenotype')
    parser.add_argument('--gene_prob', type=float, default=0.1, help='Probability of a gene being significatif (default: 0.1) for eQTL phenotype')
    parser.add_argument('-o', '--output', type=str, default="output", help='Output directory (default: output)')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-e', '--eqtl', action='store_true', help='eqtl phenotypre for each sample')
    group.add_argument('-b', '--binary', action='store_true', help='Binary phenotypre for each sample')
    group.add_argument('-q', '--quantitative', action='store_true', help='Quantitative phenotype for each sample')
    args = parser.parse_args()

    # output directory
    output_dir = args.output
    # number of top-level variants
    nvar = args.nvar
    # proportion of marker
    prop_markers = args.prop_markers
    # proportion of top-level variants that are SNPs also probability to add or keep adding SNPs when we want to add nested SNPs
    snp_prop = args.snp_prop
    # proportion of indels that we want to try to add nested SNPs into
    nested_prop = args.nested_prop
    # number of samples
    nsamp = args.nsamp
    # number of covariates
    ncov = args.ncov
    # sexual covariate
    sex_cov = None
    if args.sex_cov :
        sex_cov = args.sex_cov
        sex_eq = args.sex_eq
        sex_markers = args.sex_markers

    # number of genes for eQTL phenotype
    ngene = args.ngene
    # probability of a gene being significatif for eQTL phenotype
    gene_prob = args.gene_prob
    # random seed for reproducibility
    seed = 42
    random.seed(seed)

    if nsamp % 2 != 0:
        raise ValueError("Number of samples must be even")

    _nuc = ["A", "T", "C", "G"]

    # write phenotypes
    if args.binary:
        pg_gfa = f'{output_dir}/pg.binary.gfa'
        pg_gfa_full = f'{output_dir}/pg.binary.full.gfa'
        pg_snarl_freq = f'{output_dir}/pg.binary.snarls.freq.tsv'
        pg = Graph("binary", sex_cov)
        pg.binaryPhenotype(nsamp)
        pg.writePhenotype(f'{output_dir}/pg.binary.phenotypes.tsv')
        pg.covariate(nsamp, ncov)
        pg.covariate_sex(nsamp, sex_eq)
        pg.writeCovariate(ncov, f'{output_dir}/pg.binary.covariates.tsv')

    elif args.quantitative:
        pg_gfa = f'{output_dir}/pg.quantitative.gfa'
        pg_gfa_full = f'{output_dir}/pg.quantitative.full.gfa'
        pg_snarl_freq = f'{output_dir}/pg.quantitative.snarls.freq.tsv'
        pg = Graph("quantitative", sex_cov)
        pg.quantitativePhenotype(nsamp)
        pg.writePhenotype(f'{output_dir}/pg.quantitative.phenotypes.tsv')
        pg.covariate(nsamp, ncov)
        pg.covariate_sex(nsamp, sex_eq)
        pg.writeCovariate(ncov, f'{output_dir}/pg.quantitative.covariates.tsv')

    elif args.eqtl:
        pg_gfa = f'{output_dir}/pg.eqtl.gfa'
        pg_gfa_full = f'{output_dir}/pg.eqtl.full.gfa'
        pg_snarl_freq = f'{output_dir}/pg.eqtl.snarls.freq.tsv'
        pg = Graph("quantitative", sex_cov)
        pg.quantitativePhenotype(nsamp)
        pg.eqtlPhenotype(nsamp, ngene)
        pg.writeEqtl(f'{output_dir}/pg.eqtl.phenotypes.tsv')
        pg.writeGenePosition(f'{output_dir}/pg.eqtl.phenotypes.gene_position.tsv')
        pg.covariate(nsamp, ncov)
        pg.covariate_sex(nsamp, sex_eq)
        pg.writeCovariate(ncov, f'{output_dir}/pg.eqtl.covariates.tsv')

    # first node larger than read length
    pnod = pg.addNode(min_size=300, max_size=500)

    # add each variant
    for ii in range(nvar):
        if random.random() < snp_prop:
            # add a SNP
            pnod = pg.addSNP(pnod)
        else:
            # add an indel
            if random.random() < nested_prop:
                # with potentially some nested SNPs
                pnod = pg.addIndel(pnod, snp_prop)
            else:
                # no nested variants
                pnod = pg.addIndel(pnod, 0)

    # add last node larger than read length
    pnod = pg.addNode([pnod], min_size=300, max_size=500)
    # traverse the graph to create a reference path
    pg.addPath()
    # write simple GFA with just the pangenome and reference path
    pg.writeGfa(pg_gfa)
    # Shift the frequencies/probabilities at a subset of variant sites.
    pg.createMarkerFreq(prop_markers, sex_markers, min_dev=.01, max_dev=.5)
    # write the "truth" frequencies
    pg.writeSnarlsFreq(pg_snarl_freq, sex_cov)

    # traverse the graph to create samples for each group
    # two haplotypes for each sample
    # one sample per group
    for samp_ in range(nsamp//2):
        pg.addPath('samp_g0_' + str(samp_) + '_h0', group=0, samp=samp_)
        pg.addPath('samp_g0_' + str(samp_) + '_h1', group=0, samp=samp_)
        pg.addPath('samp_g1_' + str(samp_) + '_h0', group=1, samp=samp_)
        pg.addPath('samp_g1_' + str(samp_) + '_h1', group=1, samp=samp_)

    # write GFA containing all those haplotype paths
    pg.writeGfa(pg_gfa_full)

# python3 scripts/simulate_test_dataset.py -b
# python3 scripts/simulate_test_dataset.py -q
# python3 scripts/simulate_test_dataset.py -e
