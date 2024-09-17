from pathlib import Path
import pandas as pd
from statsmodels.stats.multitest import multipletests
import argparse
import numpy as np


def parse_rnaseq(path, orthos, organism, mirna):
    cols = ['orthomap', 'organism', 'mirna', 'gene', 'log2FoldChange', 'pvalue', 'baseMean','bMctrl', 'bMtreat']
    
    df = pd.read_csv(path, delimiter='\t')
    df = df[df['Ensembl biotype'] == 'protein_coding']
    mapcol = find_orthostring(df['Ensembl gene id'], orthos)
    
    df['orthomap'] = mapcol
    
    logfold_col = [col for col in df.columns if col.startswith('log2Fold')][0]
    bmA = [col for col in df.columns if col.startswith('baseMeanA')][0]
    bmB = [col for col in df.columns if col.startswith('baseMeanB')][0]
    con1, con2 = logfold_col.split()[1].split('/')
    # print(f'# {con1}/{con2}')
    renamed = df.copy()
    renamed.rename(columns={
        bmA: 'bMctrl', 
        bmB: 'bMtreat', 
        logfold_col: 'log2FoldChange', 
        'Ensembl gene': 'gene', 
        'Ensembl gene id': 'Ensembl_id', 
                           }, inplace=True)
    renamed['organism'] = organism
    renamed['mirna'] = mirna
    renamed = renamed.filter(cols)
    #d, u = find_significant(renamed, minT, pT, fT)
    #return d, u
    return renamed


def calc_FDR(pvals):
    reject, padj, alphacSidak, alphacBonf = multipletests(pvals, alpha=0.05, method='fdr_bh')
    return padj


def find_significant(df, minT, pT, fT):
    df = df[df['bMctrl'] >= minT]
    df = df[(df['log2FoldChange'] > fT) | (df['log2FoldChange'] < -fT)]
    df['padj'] = calc_FDR(df['pvalue'].values)
    df = df[df['padj'] <= float(pT)]

    return df


def find_orthostring(ensidseries, orthodict):
    outcol = []
    for entry in ensidseries.values:
        try:
            mapper = orthodict[entry]
        except KeyError:
            mapper = np.nan
        outcol.append(mapper)
    return outcol


def read_oma(path):
    omap = {}
    hucount = set()
    mucount = set()
    
    with open(path, 'r') as fh:
        for line in fh:
            data = line.strip().split()
            human = data[0]
            mouse = data[1]
            mapstring = f'{human}|{mouse}'
            omap[human] = mapstring
            omap[mouse] = mapstring
            
            hucount.add(human)
            mucount.add(mouse)
    print(f'No. unique human orthologs: {len(hucount)}')
    print(f'No. unique mouse orthologs: {len(mucount)}')
                                           
    return omap


def main(args):
    # load OMA pairwise orthologs
    pairwise = read_oma(args.oma_path)
    assert len(pairwise) > 0

    # load RNAseq data
    rnadf = parse_rnaseq(args.in_path, pairwise, args.organism, args.mirna)
    differential = find_significant(rnadf, args.min_reads, args.min_adj_p, args.min_logfoldchange)

    # save
    differential.to_csv(args.out_path, sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find differentially expressed genes from RNAseq data.")
    parser.add_argument("-in_path", help="the input file")
    parser.add_argument("-out_path", help="the output file")
    parser.add_argument("-oma_path", help="Path to the pairwise ortholog assignments of human and murine genes by OMA")
    parser.add_argument("-organism", help="Organism in whose iPSCs the experiment was performed")
    parser.add_argument("-mirna", help="miRNA that was overexpressed")
    parser.add_argument("-min_reads", required=False, default=100, help="Minimum number of reads in the negative control for a gene to be considered.")
    parser.add_argument("-min_logfoldchange", required=False, default=0.5, help="Minimum log2fold-change (both negative and positive).")
    parser.add_argument("-min_adj_p", required=False, default=0.05, help="adjusted p-value threshold (after filtering for minT and fT).")
    args = parser.parse_args()

    main(args)