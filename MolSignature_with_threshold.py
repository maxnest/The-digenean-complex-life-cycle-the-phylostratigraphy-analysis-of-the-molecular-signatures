import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True,
                    help="A table with expression values in TPM averaged between replicates")
parser.add_argument('--threshold', type=int, required=True,
                    help="Expression threshold for inclusion of a sequence in a molecular signature")
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()


def table_parsing(table, header, gene_dict):
    header.extend(table.readline().strip().split("\t"))
    print("Samples: {header}".format(header=" ".join(header[1:])))
    for line in table:
        description = line.strip().split("\t")
        geneID, values = description[0], description[1:]
        gene_dict[geneID] = {sample: 0 for sample in header[1:]}
        for sample in header[1:]:
            gene_dict[geneID][sample] += float(values[header[1:].index(sample)])


def get_molsign(header, gene_dict, out, threshold):
    for sample in header[1:]:
        with open("{out}.{sample}_molsign.threshold_{threshold}TPM.txt".format(
                out=out, sample=sample, threshold=threshold), 'a') as output:
            for geneID, values in gene_dict.items():
                if values[sample] >= threshold:
                    output.write("{geneID}\n".format(geneID=geneID))


def molsign_summary(header, gene_dict, out, threshold):
    with open("{out}.molsign_summary.threshold_{threshold}TPM.tsv".format(
            out=out, threshold=threshold), 'a') as output:
        output.write("GeneIDs\t{samples}\n".format(samples="\t".join(header[1:])))
        for geneID, values in gene_dict.items():
            sample_string = ["1" if float(values[sample]) >= threshold else "0" for sample in header[1:]]
            output.write("{geneID}\t{sample_string}\n".format(geneID=geneID, sample_string="\t".join(sample_string)))


if __name__ == "__main__":
    header, gene_dict = [], {}
    print("***** Input table parsing *****")
    table_parsing(args.tab, header, gene_dict)
    print("***** Identification of molecular signatures *****")
    get_molsign(header, gene_dict, args.out, args.threshold)
    print("***** Summary *****")
    molsign_summary(header, gene_dict, args.out, args.threshold)