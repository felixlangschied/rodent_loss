import pandas as pd


def enumerate_overlap(list1, list2):
    common_elements = set(list1) & set(list2)
    return len(common_elements)


def count_targets(df, organism, mirna, cutoff):
    """
    df: dataframe as outputted by 'read_and_filter_predicted_targets'
    """
    df = df[df['Total context++ score'] <= cutoff]
    fdf = df[(df['organism'] == organism) & (df['mirna'] == mirna)]
    # print(organism, mirna, fdf.index.size)
    return list(fdf.target_gene.values)


def read_and_filter_doapr(path, organism, mirna, targetscan_cutoff=-0.2):
    """
    Returns significantly down-regulated genes specific to one miRNA and one organism that are also predicted as targets.
    Only takes **TargetScanHuman** results into account!
    
    path typically: f'{PROJECTDIR}/milestones/data/targetscan/sigdown_and_predictedTargets_tsM_tsH.tsv'
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


