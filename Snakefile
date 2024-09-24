ORGANISMS = ["human", "mouse"]
MIRNAS = ["mir197", "mir769"]


module rnaseq:
    snakefile:
        "analyses/rnaseq/rnaseq.smk"
use rule * from rnaseq as rnaseq_*

module doapr:
    snakefile:
        "analyses/down_and_predicted/doapr.smk"
use rule * from doapr as doapr_*

rule all:
    input:
        # identify differentiall expressed genes from DESEQ2 results,
        # generate PCA plots of reads (each condition clusters in plot)
        rules.rnaseq_all.input,
        # find genes that are down-regulated and predicted as targets
        rules.doapr_all.input,
    default_target: True
