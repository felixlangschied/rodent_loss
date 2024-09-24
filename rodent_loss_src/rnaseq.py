import pandas as pd
from statsmodels.stats.multitest import multipletests
import numpy as np


def reads_per_gene_in_condition(path, index_key='name'):
    """
    Parses normalized count.matrix files and returns Pandas dataframe.

    index_key: 
        'name': Ensembl gene name
        'id': Ensembl Gene ID (ENSG...)
    """
    def check_index_key(column_names, index_key):
        if index_key not in column_names.values():
            possible_keys = ', '.join(column_names.values())
            raise ValueError(f'Unknown index_key "{index_key}". Use any of: {possible_keys}')
            
    column_names = {
        'Ensembl gene': 'name',
        'Ensembl gene id': 'id',
    }
    # check input
    check_index_key(column_names, index_key)

    # load
    df = pd.read_csv(path, sep='\t')
    df = df.rename(columns=column_names)

    # filter
    conditions = [name for name in df.columns if '_' in name]  # negative control + mir197 + mir769
    relevant_columns = conditions + [index_key]
    
    return df.filter(relevant_columns).set_index(index_key).astype(int)


def read_oma(path):
    """
    Reads OMA pairwise file with Ensembl GeneIDs and adds them to dictionary. Both human and mouse genes are used as keys.
    """
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


def find_significant(df, minT, pT, fT):
    """
    Filters for minimum number of reads in the negative Control and minimum amount of log2Fold-Change (positive as well as negative). 
    Then re-calculates the adjusted p-value using Benjamini-Hochberg and then filters for minimum p-adjusted
    """
    def calc_FDR(pvals):
        reject, padj, alphacSidak, alphacBonf = multipletests(pvals, alpha=0.05, method='fdr_bh')
        return padj
    
    df = df[df['bMctrl'] >= minT]
    df = df[(df['log2FoldChange'] > fT) | (df['log2FoldChange'] < -fT)]
    df['padj'] = calc_FDR(df['pvalue'].values)
    df = df[df['padj'] <= float(pT)]

    return df


def parse_rnaseq(path, orthos, organism, mirna):
    """
    Parses differential gene expression data from TSV file to pandas dataframe. Maps OMA pairwise orthologous groups to each Ensembl GeneID.
    """
    def find_orthostring(ensidseries, orthodict):
        outcol = []
        for entry in ensidseries.values:
            try:
                mapper = orthodict[entry]
            except KeyError:
                mapper = np.nan
            outcol.append(mapper)
        return outcol
        
    # load and filter for protein-coding genes
    df = pd.read_csv(path, delimiter='\t')
    df = df[df['Ensembl biotype'] == 'protein_coding']

    # add pairwise ortholog information
    mapcol = find_orthostring(df['Ensembl gene id'], orthos)
    df['orthomap'] = mapcol

    # clean column names
    cols = ['orthomap', 'organism', 'mirna', 'gene', 'log2FoldChange', 'pvalue', 'baseMean','bMctrl', 'bMtreat']
    logfold_col = [col for col in df.columns if col.startswith('log2Fold')][0]
    bmA = [col for col in df.columns if col.startswith('baseMeanA')][0]
    bmB = [col for col in df.columns if col.startswith('baseMeanB')][0]

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

    return renamed

