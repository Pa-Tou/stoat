import pandas as pd # type: ignore
import matplotlib.pyplot as plt
import qmplot # type: ignore
import argparse

def qq_plot_quantitatif(file_path, output_qqplot="qq_plot.png") :
    
    data = pd.read_csv(file_path, sep="\t")
    data = data.dropna(subset=['P'])

    # Create a Q-Q plot
    _, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")
    qmplot.qqplot(data=data["P"],
           marker="o",
           xlabel=r"Expected $-log_{10}{(P)}$",
           ylabel=r"Observed $-log_{10}{(P)}$",
           ax=ax)

    plt.savefig(output_qqplot)

def plot_manhattan_quantitatif(file_path, output_manhattan="output_manhattan_plot.png") :
    
    data = pd.read_csv(file_path, delimiter="\t")

    cleaned_data = data.dropna(subset=['CHR', 'P', 'POS']).copy()  # Use .copy() to avoid potential slice issues
    cleaned_data.loc[:, 'POS'] = cleaned_data['POS'].apply(lambda x: int(str(x).split(',')[0]) if ',' in str(x) else int(x))
    cleaned_data.loc[:, 'P'] = cleaned_data['P'].apply(lambda x: max(x, 1e-300))  # Avoid log10(0)
    
    # Prepare data for plotting
    plot_data = cleaned_data[['CHR', 'POS', 'P']].sort_values(by=['CHR', 'POS'])

    _, ax = plt.subplots(figsize=(12, 4), facecolor='w', edgecolor='k')
    qmplot.manhattanplot(data=plot_data,
                chrom="CHR",
                pv="P",
                sign_marker_p=1e-6,  # Genome wide significant p-value
                sign_marker_color="r",
                snp="POS",
                xlabel="Chromosome",
                ylabel=r"$-log_{10}{(P)}$",
                text_kws={"fontsize": 12,  # The fontsize of annotate text
                            "arrowprops": dict(arrowstyle="-", color="k", alpha=0.6)},
                ax=ax)

    plt.savefig(output_manhattan)

def qq_plot_binary(file_path, output_qqplot="qq_plot.png") :
    
    data = pd.read_csv(file_path, sep="\t")
    data = data.dropna(subset=['P_FISHER'])

    # Create a Q-Q plot
    _, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")
    qmplot.qqplot(data=data["P_FISHER"],
           marker="o",
           xlabel=r"Expected $-log_{10}{(P_FISHER)}$",
           ylabel=r"Observed $-log_{10}{(P_FISHER)}$",
           ax=ax)

    plt.savefig(output_qqplot)

def plot_manhattan_binary(file_path, output_manhattan="output_manhattan_plot.png") :

    data = pd.read_csv(file_path, sep="\t")

    # Clean data
    cleaned_data = data.dropna(subset=['CHR', 'P_FISHER', 'POS']).copy()
    cleaned_data.loc[:, 'POS'] = cleaned_data['POS'].apply(lambda x: int(str(x).split(',')[0]) if ',' in str(x) else int(x))
    cleaned_data.loc[:, 'P_FISHER'] = cleaned_data['P_FISHER'].apply(lambda x: max(x, 1e-300))  # Avoid log10(0)
    plot_data = cleaned_data[['CHR', 'POS', 'P_FISHER']].sort_values(by=['CHR', 'POS'])

    _, ax = plt.subplots(figsize=(12, 4), facecolor='w', edgecolor='k')
    qmplot.manhattanplot(data=plot_data,
                chrom="CHR",
                pv = "P_FISHER",
                snp="POS",
                xlabel="Chromosome",
                ylabel=r"$-log_{10}{(P_FISHER)}$",
                text_kws={"fontsize": 12,
                            "arrowprops": dict(arrowstyle="-", color="k", alpha=0.6)},
                ax=ax)

    plt.tight_layout()
    plt.savefig(output_manhattan)

if __name__ == "__main__" :

    parser = argparse.ArgumentParser(description="Run the Pvalue Stoat GWAS analysis")
    parser.add_argument("--pvalue",type=str, required=True)
    parser.add_argument("--qq", type=str, default="qq.png", required=False)
    parser.add_argument("--manh", type=str, default="mahanttan.png",required=False)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-b", "--binary", action="store_true")
    group.add_argument("-q", "--quantitative", action="store_true")
    args = parser.parse_args()

    if args.binary :
        qq_plot_binary(args.pvalue, args.qq)
        plot_manhattan_binary(args.pvalue, args.manh)
    
    elif args.quantitative :
        qq_plot_quantitatif(args.pvalue, args.qq)
        plot_manhattan_quantitatif(args.pvalue, args.manh)

# python3 qq_