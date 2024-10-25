# all paths should be relative paths starting from the project directory
# e.g. analyses/target_scan/scripts/example.py

rule targetscan:
    input:
        mirna='analyses/target_scan/data/mirna_targetscan_input.txt',
        utr='analyses/target_scan/data/UTR_sequences_all.txt'  # download UTR sequences from https://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi
    output:
        "analyses/target_scan/results/predicted_targets.txt",
    shell:
        "perl analyses/target_scan/scripts/targetscan7/targetscan_70.pl "
        "{input.mirna} "
        "{input.utr} "
        "{output} "


rule binning:
    input:
        utr='analyses/target_scan/data/UTR_sequences_all.txt'  # download UTR sequences from https://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi
    output:
        "analyses/target_scan/results/UTRs_median_BLs_bins.txt",
    shell:
        "perl analyses/target_scan/scripts/targetscan7/TargetScan7_BL_PCT/targetscan_70_BL_bins.pl "
        "{input.utr} "
        "{output} "


rule pct:
    input:
        targets='analyses/target_scan/results/predicted_targets.txt'
        mirna='analyses/target_scan/data/mirna_targetscan_input.txt',
        bins='"analyses/target_scan/results/UTRs_median_BLs_bins.txt"'
    output:
        "analyses/target_scan/results/UTRs_BL_PCT.txt",
    shell:
        "perl analyses/target_scan/scripts/targetscan7/TargetScan7_BL_PCT/targetscan_70_BL_PCT.pl "
        "{input.mirna} "
        "{input.targets} "
        "{input.bins} "
        "> {output}"

        