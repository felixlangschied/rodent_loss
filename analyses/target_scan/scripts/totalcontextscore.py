import pandas as pd
import json
from collections import Counter
import numpy as np
import argparse
import sys
import os

def mmu_map_from_emblgff(path):
    mapdict = {}
    with open(path, 'r') as gfh:
        for line in gfh:
            if not line.startswith('#'):
                data = line.strip().split()
                if data[2] == 'mRNA':
                    infostr = data[8]
                    trans_id = infostr.split('transcript:')[1].split(';')[0]
                    try:
                        name = infostr.split('Name=')[1].split('-')[0]
                    except IndexError as e:
                        name = trans_id
                    mapdict[trans_id] = name
    return mapdict


def hsa_map_from_emblgff(path):
    mapdict = {}
    with open(path, 'r') as gfh:
        for line in gfh:
            if not line.startswith('#'):
                data = line.strip().split()
                if data[2] == 'mRNA':
                    infostr = data[-1]
                    trans_id = infostr.split('transcript:')[1].split(';')[0]
                    name = infostr.split('Name=')[1].split('-')[0]
                    mapdict[trans_id] = name
    return mapdict


def calc_cwcs(df):
    targetsite_count = Counter(df['target_gene'])
    multiple_targetsites = [key for key, value in targetsite_count.items() if value > 1]
    singlesites = [key for key, value in targetsite_count.items() if value == 1]
    
    singledf = df[df['target_gene'].isin(singlesites)]
    singledf = singledf.rename(columns={'context++ score': 'Total context++ score'})  # score is total score if only one target site is found

    singledf = singledf.filter(['target_gene', 'Species ID', 'Mirbase ID', 'Total context++ score'])
    # get some data for later use
    species = singledf.iloc[0, 1]
    mirna = singledf.iloc[0, 2]
    # calculate total score and fill dict
    multipledict = {'target_gene': [], 'Total context++ score': []}
    for mts in multiple_targetsites:
        mtsdf = df[df['target_gene'] == mts].sort_values(by='UTR start', ascending=False)
        clist = [0]
        for i in range(mtsdf.index.size):
            csi = mtsdf.iloc[i, 5]
            airi = mtsdf.iloc[i, 6]
            ci = clist[i] + (1 - 2**csi)*(airi - clist[i])
            clist.append(ci)
        cwcs = round(np.log2(1-clist[-1]), 2)
        multipledict['target_gene'].append(mts)
        multipledict['Total context++ score'].append(cwcs)
    multipledf = pd.DataFrame.from_dict(multipledict)
    multipledf['Species ID'] = species
    multipledf['Mirbase ID'] = mirna
    totaldf = pd.concat([singledf, multipledf]).reset_index(drop=True)
    return totaldf

    
def total_context(rawdf, taxid, weighted=True):
    specdf = rawdf[rawdf['Species ID'] == taxid]
    if weighted:
        specdf = specdf.filter(['target_gene', 'Species ID', 'Mirbase ID', 'UTR start', 'UTR end', 'weighted context++ score', 'AIR'])
    elif not weighted:
        specdf = specdf.filter(['target_gene', 'Species ID', 'Mirbase ID', 'UTR start', 'UTR end', 'context++ score', 'AIR'])
    mirnas = specdf['Mirbase ID'].unique()
    dflist = []
    for mirna in mirnas:
        mirdf = specdf[specdf['Mirbase ID'] == mirna]
        dflist.append(calc_cwcs(mirdf))
    finaldf = pd.concat(dflist).reset_index(drop=True)
    return finaldf


def map_transid(transid, transk2name):
    try:
        name = transk2name[transid]
    except KeyError:
        name = transid
    return name


##########################################################################################################


def main():
    
    input = sys.argv[1]
    
    parser = argparse.ArgumentParser(
                description='Calculate Total context++-score from TargetScan output.'
            )
    parser.add_argument(
        '-i', '--input', metavar='<path>', type=str, required=True,
        help='Path to context++-scores'
    )
    parser.add_argument(
        '-s', '--species', metavar='str', type=str, required=True,
        help='"mouse" or "human" or path to GFF-file'
    )
    parser.add_argument(
        '-o', '--outfile', metavar='<path>', type=str, required=True,
        help='Path to write'
    )
    args = parser.parse_args()
    if args.species == 'mouse':
        transk2name = mmu_map_from_emblgff('/share/project/felixl/mousetargetscan/Mus_musculus.GRCm39.106.gff3')
    elif args.species == 'human':
        transk2name = hsa_map_from_emblgff('/share/project/felixl/targetscan/Homo_sapiens.GRCh37.87.gff3')

        
##########################################################################################################

    df = pd.read_csv(args.input, sep='\t')
    df = df[df['Site Type'] != '6mer']
    df['Mirbase ID'] = df['Mirbase ID'].str.replace('Mir-197-5p', 'Mir-769-5p')
    df['Gene ID'] = df['Gene ID'].apply(lambda x: x.split('.')[0])
    df['target_gene'] = df['Gene ID'].apply(lambda x: map_transid(x, transk2name))

    allspecies = df['Species ID'].unique()
    dfcollect = []
    for species in allspecies:
        specdf = total_context(df, species, weighted=False)
        dfcollect.append(specdf)

    totaldf = pd.concat(dfcollect).reset_index(drop=True)
    totaldf.to_csv(args.outfile, sep='\t', index=False)
    
##########################################################################################################

if __name__ == "__main__":
    main()