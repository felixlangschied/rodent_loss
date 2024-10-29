# all paths should be relative paths starting from the project directory
# e.g. scripts/example.py
rule all:
    input:
        "results/context_scores_output.txt",


rule targetscan:
    input:
        mirna="data/mirna_targetscan_input.txt",
        utr="data/UTR_sequences_all.txt",  # download UTR sequences from https://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi
    output:
        "results/predicted_targets.txt",
    shell:
        "perl scripts/targetscan7/targetscan_70.pl "
        "{input.mirna} "
        "{input.utr} "
        "{output} "


rule binning:
    input:
        utr="data/UTR_sequences_all.txt",  # download UTR sequences from https://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi
    output:
        "results/UTRs_median_BLs_bins.txt",
    shell:
        "perl scripts/targetscan7/TargetScan7_BL_PCT/targetscan_70_BL_bins.pl "
        "{input.utr} "
        "{output} "


rule pct:
    input:
        targets="results/predicted_targets.txt",
        mirna="data/mirna_targetscan_input.txt",
        bins="results/UTRs_median_BLs_bins.txt",
    output:
        "results/UTRs_BL_PCT.txt",
    shell:
        "perl scripts/targetscan7/TargetScan7_BL_PCT/targetscan_70_BL_PCT.pl "
        "{input.mirna} "
        "{input.targets} "
        "{input.bins} "
        "> {output}"


rule orf:
    input:
        orfs="data/ORF_Sequences.txt",  # download ORF sequences from https://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi
        mirna="data/mirna_targetscan_input.txt",
    output:
        counts="results/ORF_8mer_counts.txt",
        lengths="data/ORF_Sequences.lengths.txt",
    shell:
        "perl scripts/targetscan7/TargetScan7_context_scores/targetscan_count_8mers.pl "
        "{input.mirna} "
        "{input.orfs} "
        ">| {output.counts}"


rule contextscore:
    # make sure that you have RNAplfold from the ViennaRNA Package 2 installed and in the Path
    input:
        matmirs="data/mirna_contextscore_input.txt",
        utr="data/UTR_sequences_all.txt",
        pct="results/UTRs_BL_PCT.txt",
        lengths="data/ORF_Sequences.lengths.txt",
        counts="results/ORF_8mer_counts.txt",
    output:
        "results/context_scores_output.txt",
    shell:
        "perl scripts/targetscan7/TargetScan7_context_scores/targetscan_70_context_scores.pl "
        "{input.matmirs} "
        "{input.utr}"
        "{input.pct} "
        "{input.lengths} "
        "{input.counts} "
        "{output} "
