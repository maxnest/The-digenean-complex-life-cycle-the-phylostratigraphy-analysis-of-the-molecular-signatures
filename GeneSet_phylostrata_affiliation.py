import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--fasta', type=argparse.FileType('r'), required=True,
                    help="Fasta file with genes under consideration")
parser.add_argument('--phylostratr', type=argparse.FileType('r'), required=True,
                    help="Table with results of phylostratigraphic analysis carried out with Phylostratr."
                         "Gene IDs should be in table!")
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def phylostratr_parsing(phylostratr, phylostratr_dict):
    header = phylostratr.readline()
    for line in phylostratr:
        description = line.strip().split("\t")
        geneID, mrca, ps, mrca_name = description[0], description[1], description[2], description[3]
        phylostratum = "{ps}:{mrca_name}".format(ps=ps, mrca_name=mrca_name)
        if phylostratum not in phylostratr_dict:
            phylostratr_dict[phylostratum] = []
        phylostratr_dict[phylostratum].append(geneID)


def fasta_parsing(fasta):
    fasta_list = []

    fasta_seqs = SeqIO.parse(fasta, "fasta")
    for record in fasta_seqs:
        fasta_list.append(record.id)

    return fasta_list


def output_writing(out, phylostratr_dict, fasta_list):
    summary_dict = {phylostratum: [] for phylostratum in phylostratr_dict}
    for phylostratum, genes in phylostratr_dict.items():
        for gene in genes:
            if gene in fasta_list:
                summary_dict[phylostratum].append(gene)

    with open("{out}.phylostrata_summary_table.tsv".format(out=out), 'a') as output:
        output.write("Phylostratum\tGeneCount\tPercent\n")
        for phylostratum, genes in summary_dict.items():
            output.write("{phylo}\t{count}\t{percent}\n".format(phylo=phylostratum, count=len(genes),
                                                                percent=round((len(genes)/len(fasta_list)) * 100, 2)))

    phylostrata = [phylostratum for phylostratum in phylostratr_dict]
    with open("{out}.phylostrata_pres_abs_matrix.tsv".format(out=out), 'a') as output_matrix:
        output_matrix.write("GeneIDs\t{phylostrata}\n".format(phylostrata="\t".join(phylostrata)))
        for gene in fasta_list:
            gene_pres_abs = ["1" if gene in phylostratr_dict[phylostratum] else "0" for phylostratum in phylostrata]
            output_matrix.write("{id}\t{gene_pres_abs}\n".format(id=gene, gene_pres_abs="\t".join(gene_pres_abs)))


if __name__ == "__main__":
    phylostratr_dict = {}
    phylostratr_parsing(args.phylostratr, phylostratr_dict)
    wanted_seq = fasta_parsing(args.fasta)
    output_writing(args.out, phylostratr_dict, wanted_seq)