import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
pd.set_option('mode.chained_assignment', None)
import argparse
from rodent_loss_src.rnaseq import reads_per_gene_in_condition


def main(args):
    df = reads_per_gene_in_condition(args.in_path)
    tdf = df.transpose()

    # Separating out the features
    x = tdf.loc[:, tdf.columns].values
    y = tdf.index# Standardizing the features
    x = StandardScaler().fit_transform(x)
    pca = PCA(n_components=2)
    pca.fit(x)

    exp1, exp2 = pca.explained_variance_ratio_

    exP1 = round(exp1*100, 1)
    exP2 = round(exp2*100, 1)

    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents, columns = ['Dim 1', 'Dim 2'], index=tdf.index)
    pdf = principalDf.reset_index()
    pdf['probe'] = pdf['index'].apply(lambda x: x.split('_')[0])
    pdf['probe'][pdf['probe'] == 'Neg-Ctl'] = 'CTRL'
    #display(pdf)

    sns.set(rc={'figure.figsize':(6,4), 'ytick.left': True, 'xtick.bottom': True}, font_scale = 1.3, style='whitegrid')
    sns.scatterplot(data=pdf, x='Dim 1', y='Dim 2', hue='probe', s=90)
    plt.xlabel(f'Dim 1 ({exP1}%)')
    plt.ylabel(f'Dim 2 ({exP2}%)')
    # plt.ylim([-100, 125])
    # plt.xlim([-100, 125])
    plt.tight_layout()
    plt.savefig(args.out_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find differentially expressed genes from RNAseq data.")
    parser.add_argument("-in_path", help="the input file")
    parser.add_argument("-out_path", help="the output file")
    parser.add_argument("-organism", help="Organism in whose iPSCs the experiment was performed")
    args = parser.parse_args()

    main(args)
    