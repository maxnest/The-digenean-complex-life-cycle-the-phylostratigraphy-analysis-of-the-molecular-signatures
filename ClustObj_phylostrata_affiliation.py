import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--clust_obj', type=argparse.FileType('r'), required=True,
                    help="Modified Clusters_Objects.tsv (Clust output) table:"
                         "blank cells were replaced to '-' and 1 row was removed")
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


def clust_obj_parsing(clust_obj, clusters, clust_dict):
    clusters.extend(clust_obj.readline().strip().split("\t"))
    print("The following {num} clusters were identified: {clusters}".format(num=len(clusters),
                                                                            clusters="\t".join(clusters)))
    for cluster in clusters:
        clust_dict[cluster] = []

    for line in clust_obj:
        genes = line.strip().split("\t")
        for gene in genes:
            if gene != '-':
                clust_dict[gene] = clusters[genes.index(gene)]


def output_writing(out, phylostratr_dict, clusters, clust_dict):
    summary_dict = {phylostratum: {cluster: [] for cluster in clusters} for phylostratum in phylostratr_dict}

    for phylostratum, genes in phylostratr_dict.items():
        for gene in genes:
            if gene in clust_dict:
                summary_dict[phylostratum][clust_dict[gene]].append(gene)

    with open("{out}.phylostrata_summary_table.gene_count.tsv".format(out=out), 'a') as output_count:
        output_count.write("Phylostratum\t{clusters}\n".format(clusters="\t".join(clusters)))
        for phylostratum, values in summary_dict.items():
            count_num = [str(len(values[cluster])) for cluster in clusters]
            output_count.write("{phylo}\t{count}\n".format(phylo=phylostratum, count="\t".join(count_num)))

    with open("{out}.phylostrata_summary_table.percents.tsv".format(out=out), 'a') as output_percent:
        output_percent.write("Phylostratum\t{clusters}\n".format(clusters="\t".join(clusters)))
        for phylostratum, values in summary_dict.items():
            percent = [
                str(round((len(values[cluster])/
                           len([gene for gene in clust_dict if clust_dict[gene] == cluster])) * 100, 2))
                for cluster in clusters]
            output_percent.write("{phylo}\t{percent}\n".format(phylo=phylostratum, percent="\t".join(percent)))


if __name__ == "__main__":
    phylostratr_dict, clusters, cluster_dict = {}, [], {}
    phylostratr_parsing(args.phylostratr, phylostratr_dict)
    clust_obj_parsing(args.clust_obj, clusters, cluster_dict)
    output_writing(args.out, phylostratr_dict, clusters, cluster_dict)