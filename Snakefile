ORGANISMS = ["human", "mouse"]
MIRNAS = ["mir197", "mir769"]


rule pca_figure:
    input:
        expand(
            "analyses/rnaseq/results/{organism}_reads_pca.png",
            organism=ORGANISMS,
            mirna=MIRNAS,
        ),


rule reads_pca:
    input:
        expand(
            "analyses/rnaseq/data/{organism}_results_Neg_vs_{mirna}.tsv",
            organism=lambda wildcards: wildcards.organism,
            mirna=MIRNAS,
        ),
    output:
        "analyses/rnaseq/results/{organism}_reads_pca.png",
    conda:
        "environment.yml"
    shell:
        "python analyses/rnaseq/scripts/reads_pca.py {input} {wildcards.organism} {output}"


rule all:
    input:
        expand(
            "milestones/data/rnaseq/{organism}_{mirna}_differential_genes.tsv",
            organism=ORGANISMS,
            mirna=MIRNAS,
        ),


rule find_differential_genes:
    input:
        "analyses/rnaseq/data/{organism}_results_Neg_vs_{mirna}.tsv",
    output:
        "milestones/data/rnaseq/{organism}_{mirna}_differential_genes.tsv",
    shell:
        "python analyses/rnaseq/scripts/differential_expressed_genes.py "
        "-in_path {input} "
        "-out_path {output} "
        "-oma_path analyses/rnaseq/data/human_mouse_omapairwise.txt "
        "-organism {wildcards.organism} "
        "-mirna {wildcards.mirna}"
