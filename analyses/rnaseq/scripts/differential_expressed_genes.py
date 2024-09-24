import argparse
from rodent_loss_src.rnaseq import read_oma, parse_rnaseq, find_significant


def main(args):
    # load OMA pairwise orthologs
    pairwise = read_oma(args.oma_path)
    assert len(pairwise) > 0

    # load RNAseq data
    rnadf = parse_rnaseq(args.in_path, pairwise, args.organism, args.mirna)
    differential = find_significant(rnadf, args.min_reads, args.min_adj_p, args.min_logfoldchange)

    # save
    differential.to_csv(args.out_path, sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find differentially expressed genes from RNAseq data.")
    parser.add_argument("-in_path", help="the input file")
    parser.add_argument("-out_path", help="the output file")
    parser.add_argument("-oma_path", help="Path to the pairwise ortholog assignments of human and murine genes by OMA")
    parser.add_argument("-organism", help="Organism in whose iPSCs the experiment was performed")
    parser.add_argument("-mirna", help="miRNA that was overexpressed")
    parser.add_argument("-min_reads", required=False, default=100, help="Minimum number of reads in the negative control for a gene to be considered.")
    parser.add_argument("-min_logfoldchange", required=False, default=0.5, help="Minimum log2fold-change (both negative and positive).")
    parser.add_argument("-min_adj_p", required=False, default=0.05, help="adjusted p-value threshold (after filtering for minT and fT).")
    args = parser.parse_args()

    main(args)