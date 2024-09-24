ORGANISMS = ["human", "mouse"]
MIRNAS = ["mir197", "mir769"]

rule all:
    input:
        expand(
            "analyses/down_and_predicted/results/read_heatmap/{organism}_{mirna}_read_heatmap.svg",
            organism=ORGANISMS,
            mirna=MIRNAS,
        ),


rule doapr_heatmaps:
    input:
        deseq="external_data/counts_matrix/{organism}_counts.matrix",
        targets="milestones/data/targetscan/sigdown_and_predictedTargets_tsM_tsH.tsv"
    output:
        "analyses/down_and_predicted/results/read_heatmap/{organism}_{mirna}_read_heatmap.svg",
    shell:
        "python analyses/down_and_predicted/scripts/read_heatmaps.py "
        "-deseq_path {input.deseq} "
        "-targets_path {input.targets} "
        "-out_path {output} "
        "-organism {wildcards.organism} "
        "-mirna {wildcards.mirna} "
        "-targetscan_cutoff -0.2"