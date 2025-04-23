ORGANISMS = ["human", "mouse"]
MIRNAS = ["mir197", "mir769"]

rule all:
    input:
        expand(
            "milestones/data/rnaseq/{organism}_{mirna}_differential_genes.tsv",
            organism=ORGANISMS,
            mirna=MIRNAS,
        ),
        expand(
            "analyses/rnaseq/results/{organism}_reads_pca.png",
            organism=ORGANISMS,
            mirna=MIRNAS,
        ),


rule reads_pca:
    input:
        'external_data/counts_matrix/{organism}_counts.matrix',
    output:
        "analyses/rnaseq/results/{organism}_reads_pca.png",
    shell:
        "python analyses/rnaseq/scripts/reads_pca.py "
        "-in_path {input} "
        "-out_path {output} "
        "-organism {wildcards.organism} "


rule find_differential_genes:
    input:
        deseq="analyses/rnaseq/data/{organism}_results_Neg_vs_{mirna}.tsv",
        oma="analyses/rnaseq/data/human_mouse_omapairwise.txt"
    output:
        "milestones/data/rnaseq/{organism}_{mirna}_differential_genes.tsv",
    shell:
        "python analyses/rnaseq/scripts/differential_expressed_genes.py "
        "-in_path {input.deseq} "
        "-out_path {output} "
        "-oma_path {input.oma} "
        "-organism {wildcards.organism} "
        "-mirna {wildcards.mirna}"
