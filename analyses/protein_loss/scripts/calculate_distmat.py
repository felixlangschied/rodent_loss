import pandas as pd
import numpy as np
#from sklearn.feature_selection import mutual_info_regression
import os

###############################################################################################################################
orderpath = '../data/data/PP_taxidorder.txt'
specmappath = '../data/spec_map.txt'
mirfampath = '../data/output/phyloprofile/mirgenedb2_0/family_seedpp.long'
protpath = '../data/human_mammalia_orthologs_fdog.txt'
theria_taxa = '../data/theria_taxids.txt'
theria_only = True
outpath = '../results/theria_only_pearson.tsv'
###############################################################################################################################

def order_matrix(df, opath):
    """
    Will remove Carlito syrichta and Petromyzon from taxon sampling because they are not part of pre-miRNA species tree
    """
    smallorder = []
    with open(opath, 'r') as fh:
        order = [line for line in fh.read().split('\n') if line]
        #print(order)
    for taxon in order:
        if taxon in df.columns:
            smallorder.append(taxon)
    df = df[smallorder]
    return df


def parse_mirna_matrix(pppath, filterSeed=True, o='index'):
    col = {}
    with open(pppath, 'r') as fh:
        header = next(fh)
        for line in fh:
            mirfam, ncbi, mirid, seedcon = line.strip().split()
            if filterSeed and not bool(seedcon):
                continue
            if not mirfam in col:
                col[mirfam] = {}
            col[mirfam][ncbi] = 1
    df = pd.DataFrame.from_dict(col, orient=o).fillna(0).astype('int')
    return df


def parse_protein_matrix(path, o='index'):
    col = {}
    with open(path, 'r') as fh:
        for line in fh:
            rawgroup, fdogid, protid, score = line.strip().split('|')
            ncbi = f'ncbi{fdogid.split("@")[1]}'
            orthogroup = rawgroup.replace('>', '')
            if not orthogroup in col:
                col[orthogroup] = {}
            col[orthogroup][ncbi] = 1
    df = pd.DataFrame.from_dict(col, orient=o).fillna(0).astype('int')
    return df
    
    
def calc_distmat(df1, df2, opath, force=False):
    if not force and os.path.isfile(opath):
        print(f'Found file at {opath}. Nothing done')
        return None

    dist = np.corrcoef(df1.to_numpy(), df2.to_numpy())
    names = list(df1.index.values) + list(df2.index.values)
    fulldf = pd.DataFrame(dist, columns=names, index=names)
    distdf = fulldf.loc[list(df1.index.values), list(df2.index.values)]   
    
    distdf.to_csv(opath, sep='\t')
    
    

###############################################################################################################################

with open(theria_taxa) as fh:
    theria_taxids = [f'ncbi{line.strip()}' for line in fh if line]

# protein-coding
protdf = parse_protein_matrix(protpath)
protdf = order_matrix(protdf, orderpath)
if theria_only:
    protdf = protdf.filter(theria_taxids)

# miRNA, assume non-seed-conserved orthologs as no "functional" ortholog
mirdf = parse_mirna_matrix(mirfampath, orderpath)
mirdf = order_matrix(mirdf, orderpath)
# adjust to taxa in protein-coding dataframe
smirdf = mirdf.filter(list(protdf.columns), axis='columns')


calc_distmat(smirdf, protdf, 'pearson', outpath, force=True)


