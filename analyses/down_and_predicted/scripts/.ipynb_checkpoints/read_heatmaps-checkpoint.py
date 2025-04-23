import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import zscore
from rodent_loss_src.rnaseq import reads_per_gene_in_condition


def read_and_filter_predicted_targets(path, organism, mirna, targetscan_cutoff=-0.2):
    """
    Returns significantly down-regulated genes specific to one miRNA and one organism. Only takes **TargetScanHuman** results into account!
    """
    df = pd.read_csv(path, sep='\t')
    df.mirna = df.mirna.apply(lambda x: x.replace('Mir-197-3p', 'mir197').replace('Mir-769-5p', 'mir769'))
    
    df = df[df['organism']  == organism]
    df = df[df['mirna'] == mirna]
    df = df.filter(['gene', 'TShuman_tc++s'])
    df = df[(df['TShuman_tc++s'] <= targetscan_cutoff)]
    if organism == 'mouse':
        df['gene'] = df['gene'].str.capitalize()

    # modify index
    df = df.set_index('gene')

    return df

def zscore_targetscan_heatmap(df, cmap='viridis'):
    sns.set(rc={'figure.figsize':(2.5,6)})
    sns.set(font_scale=0.8)
    tardf = df.filter(['TShuman_tc++s'])
    expdf = df.drop(columns=['TShuman_tc++s'])
    zdf = expdf.apply(zscore, axis=1)
    
    fig, (ax,ax2) = plt.subplots(ncols=2)
    sns.heatmap(tardf, cmap=cmap, ax=ax, vmin=-1, vmax=-0.2, linewidth=0.3)
    sns.heatmap(zdf, cmap="vlag", ax=ax2, vmin=-1, vmax=1.5, linewidth=0.3)

    ax2.yaxis.tick_right()
    ax2.tick_params(rotation=0)
    return fig


def filter_readsdf_mirna(df, mirna):
    if mirna == 'mir197':
        columns = [col for col in df.columns if not '-769_' in col]
    else:
        columns = [col for col in df.columns if not '-197_' in col]
    return df.filter(columns)


def main(args):
    reads_df = reads_per_gene_in_condition(args.in_path)
    mir_reads_df = filter_readsdf_mirna(reads_df, args.mirna)
    targetdf = read_and_filter_predicted_targets(tarpath, args.organism, args.mirna)
    heat_df = mir_reads_df.join(targetdf, how='inner')
    fig = zscore_targetscan_heatmap(heat_df)
    plt.savefig(args.out_path})


if __name__ == "__main__":    
    parser = argparse.ArgumentParser(description="Find differentially expressed genes from RNAseq data.")
    parser.add_argument("-in_path", help="the input file")
    parser.add_argument("-out_path", help="the output file")
    parser.add_argument("-organism", help="Organism in whose iPSCs the experiment was performed")
    parser.add_argument("-mirna", help="miRNA that was overexpressed")
    parser.add_argument("-targetscan_cutoff", required=False, default=-0.2, help="TargetScanHuman context++ score cutoff")
    args = parser.parse_args()

    main(args)