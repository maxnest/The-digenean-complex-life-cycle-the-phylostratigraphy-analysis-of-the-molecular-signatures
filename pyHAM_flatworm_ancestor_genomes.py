# --------------------------------------------------------------------------------------------------------
# Python code used to reconstruct Platyhelminthes and Digenea ancestor genome models in a publication
# 'The digenean complex life cycle: the phylostratigraphy analysis of the molecular signatures' written by
# Maksim A. Nesterenko, Sergei V. Shchenkov, Sofia A. Denisova and Viktor V. Starunov
#
# E-mail: maxnest.research@gmail.com
# --------------------------------------------------------------------------------------------------------

# NB!:
# PyHAM deals with ancestral genomes as well as extant genomes.
# It basically uses the extant genomes and orthology information to reconstruct the ancestral genomes.
# An ancestral genome can be considered as all the hogs present at a given taxonomic level.

# Load in all necessary scientific libraries
import pandas as pd
import pyham
from Bio import SeqIO

# Input files:
working_dir = 'Ancestral_pyHAM/'
nwk_file = working_dir + "ManualSpeciesTree_with_internal_nodes.nwk"  # Phylogenetic tree with internal nodes names
orthoxml_file = working_dir + "HierarchicalGroups.orthoxml"  # OMA standalone output file

fgig_fasta = working_dir + "DB/Fasciola_gigantica.fa"
fhep_fasta = working_dir + "DB/Fasciola_hepatica.fa"
psim_fasta = working_dir + "DB/Psilotrema_simillimum.fa"
sman_fasta = working_dir + "DB/Schistosoma_mansoni.fa"
treg_fasta = working_dir + "DB/Trichobilharzia_regenti.fa"
tszi_fasta = working_dir + "DB/Trichobilharzia_szidati.fa"


def ancestor_genome_with_sp_threshold(ham_analysis, ancestral_genome_name, sp_threshold, ancestral_genes_filtered):
    """
    Application of a threshold (sp_threshold) to the constructed by pyHAM model of the ancestor genome.
    The threshold is a value corresponding to the number of modern species that must have a gene
    in order for it to be included in the model.

    Parameters:
    :param ham_analysis: pyham.Ham function result
    :param ancestral_genome_name: python string corresponding to the name of one of the internal node of the nwk-tree
    :param sp_threshold: value corresponding to the number of modern species that must have a gene
    :param ancestral_genes_filtered: python list with selected genes
    """

    # get ancestral genome using pyham
    ancestral_genome = ham_analysis.get_ancestral_genome_by_name(ancestral_genome_name)

    # Get its ancestral genes/HOGs
    ancestral_genes = ancestral_genome.genes

    # filtering:
    for ancestral_gene in ancestral_genes:
        if len(ancestral_gene.get_all_descendant_genes_clustered_by_species()) >= sp_threshold:
            ancestral_genes_filtered.append(ancestral_gene)

    print("After applying the filter in >= {sp_threshold} species, "
          "{ancestor} ancestor genome model includes {num} genes".format(
            sp_threshold=str(sp_threshold), ancestor=ancestral_genome_name, num=len(ancestral_genes_filtered)))


def ancestor_genome_df(ancestral_genes):
    """
    Pandas dataframe formation

    :param ancestral_genes: python list with ancestor genes
    :return: pandas dataframe
    """

    ancestral_genome_df = pd.DataFrame(ancestral_genes, columns={"ancestral_gene"})

    # get descendant genes of the filtered hog
    ancestral_genome_df['descendant_genes'] = ancestral_genome_df.apply(
        lambda x: x["ancestral_gene"].get_all_descendant_genes(), axis=1)

    # modify df so that each descendant gene has its own row
    ancestral_genome_df = \
        ancestral_genome_df.set_index(['ancestral_gene'])['descendant_genes'].apply(pd.Series).stack().reset_index(
            level=1, drop=True).reset_index()

    # rename the column
    ancestral_genome_df = ancestral_genome_df.rename({0: "descendant_gene"}, axis=1)

    # make a string id of ancestral gene
    ancestral_genome_df['ancestral_gene_id'] = ancestral_genome_df.apply(
        lambda x: str(x['ancestral_gene'].get_all_descendant_genes()), axis=1)

    # convert pyham gene name to cross-reference gene name
    ancestral_genome_df['descendant_gene_xref'] = ancestral_genome_df.apply(
        lambda x: x['descendant_gene'].get_dict_xref()['protId'].split(" ")[0], axis=1)

    return ancestral_genome_df


def get_retained_ancestral_gene_ids(retained, first_ancestor_sp_threshold, second_ancestor_sp_threshold):
    """
    Applying the filters of the minimum number of modern species with the considered gene
    to the list of genes that have retained

    :param retained: python dictionary with retained genes (vertical comparison result)
    :param first_ancestor_sp_threshold: value corresponding to the number of modern species that must have a gene
    :param second_ancestor_sp_threshold: value corresponding to the number of modern species that must have a gene
    :return: python list with IDs of genes that retained
    """

    descendant_retained_gene_ids = []
    for first_ancestor_gene, second_ancestor_gene in retained.items():
        if len(first_ancestor_gene.get_all_descendant_genes_clustered_by_species()) >= first_ancestor_sp_threshold and \
                len(second_ancestor_gene.get_all_descendant_genes_clustered_by_species()) \
                >= second_ancestor_sp_threshold:
            # get all descendant extant gene ids
            descendant_retained_gene_ids.append(second_ancestor_gene.get_all_descendant_genes())
    return descendant_retained_gene_ids


def get_duplicated_ancestral_gene_ids(duplicated, copy_threshold, first_ancestor_sp_threshold,
                                      second_ancestor_sp_threshold):
    """
    Applying the filters to the list of genes that arose through duplication

    :param duplicated: python dictionary with duplicated genes (vertical comparison result)
    :param copy_threshold: minimum number of duplicated gene copies
    :param first_ancestor_sp_threshold: value corresponding to the number of modern species that must have a gene
    :param second_ancestor_sp_threshold: value corresponding to the number of modern species that must have a gene
    :return: python list with IDs of genes that duplicated
    """

    descendant_duplicated_gene_ids = []
    for first_ancestor_gene, second_ancestor_genes in duplicated.items():
        if len(first_ancestor_gene.get_all_descendant_genes_clustered_by_species()) >= first_ancestor_sp_threshold:
            second_ancestor_genes_filtered = [second_ancestor_gene for second_ancestor_gene in second_ancestor_genes if
                                              len(second_ancestor_gene.get_all_descendant_genes_clustered_by_species())
                                              >= second_ancestor_sp_threshold]
            if len(second_ancestor_genes_filtered) >= copy_threshold:
                descendant_duplicated_gene_ids.extend([second_ancestor_gene.get_all_descendant_genes()
                                                       for second_ancestor_gene in second_ancestor_genes_filtered])
    return descendant_duplicated_gene_ids


def get_gained_ancestral_gene_ids(gained, second_ancestor_sp_threshold):
    """
    Applying the filter to the list of genes that appeared (gained)

    :param gained: python list with appeared (gained) genes (vertical comparison results)
    :param second_ancestor_sp_threshold: value corresponding to the number of modern species that must have a gene
    :return: python list with IDs of genes that appeared (gained)
    """

    descandant_gained_gene_ids = []
    ancestral_gained_genes_filtered = [second_ancestor_gene for second_ancestor_gene in gained if
                                       len(second_ancestor_gene.get_all_descendant_genes_clustered_by_species())
                                       >= second_ancestor_sp_threshold]
    descandant_gained_gene_ids.extend([gene.get_all_descendant_genes() for gene in ancestral_gained_genes_filtered])
    return descandant_gained_gene_ids


def vertical_comparison_with_filters(ham_analysis, first_ancestor, second_ancestor,
                                    retained_gene_ids, duplicated_gene_ids, gained_gene_ids,
                                    copy_threshold, first_ancestor_sp_threshold, second_ancestor_sp_threshold):
    """
    Vertical comparison between reconstucred genome models of ancestors
    NB!:
    pyHam allows tracing of HOGs/genes provided in the orthoxml file and through the protein files along each branch of
    a phylogenetic tree (provided in the newick format) and report the genes from these HOGs that arose
    through duplication on each branch, got lost on each branch, each appeared on that branch, or were simply retained.

    :param ham_analysis: pyham.Ham function result
    :param first_ancestor: python string corresponding to the name of one of the internal node of the nwk-tree
    :param second_ancestor: python string corresponding to the name of other one of the internal node of the nwk-tree
    :param retained_gene_ids: empty python list for genes that retained
    :param duplicated_gene_ids: empty python list for genes that arose through duplication
    :param gained_gene_ids: empty python list for genes that appeared (gained)
    :param copy_threshold: minimum number of duplicated gene copies
    :param first_ancestor_sp_threshold: value corresponding to the number of modern species that must have a gene
    :param second_ancestor_sp_threshold: value corresponding to the number of modern species that must have a gene
    """

    first_ancestor_genome = ham_analysis.get_ancestral_genome_by_name(first_ancestor)
    second_ancestor_genome = ham_analysis.get_ancestral_genome_by_name(second_ancestor)

    vert_comp = ham_analysis.compare_genomes_vertically(first_ancestor_genome, second_ancestor_genome)

    # The identical genes (that stay single copies)
    # one HOG at first_ancestor -> one descendant gene in second_ancestor
    retained_genes = vert_comp.get_retained()   # {<HOG(17901)>: <HOG()>, <HOG(17902)>: <HOG()>}

    # The duplicated genes (that have duplicated)
    # one HOG at first_ancestor -> list of its descendants gene in second_ancestor
    duplicated_genes = vert_comp.get_duplicated()   # {<HOG(17912)>: [<HOG()>, <HOG()>, <HOG()>]}

    # The gained genes (that emerged in between)
    # list of gene that appeared after second_ancestor taxon
    gained_genes = vert_comp.get_gained()   # [<HOG(14230)>, <HOG(14231)>, <HOG(14232)>]

    # Applying filters and extract descendant_gene_ids
    retained_gene_ids.extend(get_retained_ancestral_gene_ids(
        retained_genes, first_ancestor_sp_threshold, second_ancestor_sp_threshold))
    duplicated_gene_ids.extend(get_duplicated_ancestral_gene_ids(
        duplicated_genes, copy_threshold, first_ancestor_sp_threshold, second_ancestor_sp_threshold))
    gained_gene_ids.extend(get_gained_ancestral_gene_ids(gained_genes, second_ancestor_sp_threshold))

    print("{first} VS. {second} vertical comparison results:\t"
          "{ret_num} retained, {dup_num} duplicated, and {gain_num} gained genes".format(
            first=first_ancestor, second=second_ancestor, ret_num=len(retained_gene_ids),
            dup_num=len(duplicated_gene_ids), gain_num=len(gained_gene_ids)))


def fasta_parsing(fasta):
    fasta_dict = {}

    fasta_seqs = SeqIO.parse(fasta, "fasta")
    for record in fasta_seqs:
        fasta_dict[record.id] = record.seq

    return fasta_dict


def ancestor_genome_fasta(fasta_dict, ancestor_df, working_dir, ref_tag, ancestor_tag):
    """
    Parsing a fasta file to extract proteins from an ancestor genome model

    :param fasta_dict: python dictionary for fasta file: protein ID is a key, sequence is the value
    :param ancestor_df: pandas dataframe
    :param working_dir: python string with path to working directory
    :param ref_tag: python string with the name of the species whose fasta file is being parsing
    :param ancestor_tag: python string with ancestor name
    """

    with open("{wd}/{ref_tag}.{ancestor_tag}_ancestor.fasta".format(
                wd=working_dir, ref_tag=ref_tag, ancestor_tag=ancestor_tag), 'a') as output_fasta:

        for id, seq in fasta_dict.items():
            if id in list(ancestor_df['descendant_gene_xref']):
                output_fasta.write(">{name}\n{seq}\n".format(name=id, seq=seq))


def comparison_results_fasta(fasta_dict, results, working_dir, ref_tag, first_ancestor, second_ancestor, results_tag):
    """
    Parsing a fasta file to extract proteins from an vertical comparison results

    :param fasta_dict: python dictionary for fasta file: protein ID is a key, sequence is the value
    :param results: python list with IDs of genes that retained/duplicated/gained
    :param working_dir: python string with path to working directory
    :param ref_tag: python string with the name of the species whose fasta file is being parsing
    :param first_ancestor: python string with ancestor name
    :param second_ancestor: python string with ancestor name
    :param results_tag: python string indicating whether these genes that retained/duplicated/gained
    """

    results_xref = []
    for gene_list in results:
        for gene in gene_list:
            results_xref.append(gene.get_dict_xref()['protId'].split(" ")[0])

    with open("{wd}/{first}_vs_{second}_vert_comp.{results_tag}.{ref_tag}.fasta".format(
            wd=working_dir, ref_tag=ref_tag, first=first_ancestor, second=second_ancestor,
            results_tag=results_tag), 'a') as output_fasta:
        for id, seq in fasta_dict.items():
            if id in results_xref:
                output_fasta.write(">{name}\n{seq}\n".format(name=id, seq=seq))


def comparison_results_txt(results, working_dir, first_ancestor, second_ancestor, results_tag):
    """
    Сreating a text file with a list of genes (one ID per line) that retained/duplicated/gained

    :param results: python list with IDs of genes that retained/duplicated/gained
    :param working_dir: python string with path to working directory
    :param first_ancestor: python string with ancestor name
    :param second_ancestor: python string with ancestor name
    :param results_tag: python string indicating whether these genes that retained/duplicated/gained
    """

    results_xref = []
    for gene_list in results:
        for gene in gene_list:
            results_xref.append(gene.get_dict_xref()['protId'].split(" ")[0])

    with open("{wd}/{first}_vs_{second}_vert_comp.{results_tag}.txt".format(
            wd=working_dir, first=first_ancestor, second=second_ancestor, results_tag=results_tag), 'a') as output_txt:
        for xref in results_xref:
            output_txt.write("{xref}\n".format(xref=xref))


if __name__ == "__main__":
    # All genes in ancestor genomes according to pyHAM:
    digenea_ancestor_all, platy_ancestor_all,  = [], []
    # Filtered sets of genes from ancestor genomes:
    digenea_ancestor_filtered, platy_ancestor_filtered = [], []
    # For vertical comparison:
    retained_gene_ids, duplicated_gene_ids, gained_gene_ids = [], [], []

    # pyHAM analysis:
    ham_analysis = pyham.Ham(nwk_file, orthoxml_file, use_internal_name=True)

    # Ancestor genome reconstruction:
    digenea_ancestor_genome = ham_analysis.get_ancestral_genome_by_name("Digenea")
    platy_ancestor_genome = ham_analysis.get_ancestral_genome_by_name("Platyhelminthes")

    digenea_ancestor_all.extend(digenea_ancestor_genome.genes)
    print("According to pyHAM, Digenea ancestor genome model includes {num} genes".format(
        num=len(digenea_ancestor_all)))

    platy_ancestor_all.extend(platy_ancestor_genome.genes)
    print("According to pyHAM, Platyhelminthes ancestor genome model includes {num} genes".format(
        num=len(platy_ancestor_all)))

    # Species threshold applying:
    # Given the possible incompleteness of available protein sets, we chose 75% of species as an appropriate threshold
    # In general, 11 Digenea gene set are analyzed, and 75% equal to 8.25 (8) species:
    ancestor_genome_with_sp_threshold(ham_analysis, "Digenea", 8, digenea_ancestor_filtered)
    # In general, 14 flatworm gene set are analyzed, and 75% equal to 10.5 (11) species:
    ancestor_genome_with_sp_threshold(ham_analysis, "Platyhelminthes", 11, platy_ancestor_filtered)

    # Make pandas dataframes:
    digenea_ancestor_all_df = ancestor_genome_df(digenea_ancestor_all)
    digenea_ancestor_filtered_df = ancestor_genome_df(digenea_ancestor_filtered)
    platy_ancestor_all_df = ancestor_genome_df(platy_ancestor_all)
    platy_ancestor_filtered_df = ancestor_genome_df(platy_ancestor_filtered)

    # Write pandas dataframes as TSV-files:
    digenea_ancestor_all_df.to_csv(working_dir + "Digenea_ancestor.all_genes.pyHAM.pd_df.tsv", sep="\t")
    digenea_ancestor_filtered_df.to_csv(working_dir + "Digenea_ancestor.8_sp_threshold.pyHAM.pd_df.tsv", sep="\t")
    platy_ancestor_all_df.to_csv(working_dir + "Platyhelminthes_ancestor.all_genes.pyHAM.pd_df.tsv", sep="\t")
    platy_ancestor_filtered_df.to_csv(working_dir + "Platyhelminthes_ancestor.11_sp_threshold.pyHAM.pd_df.tsv", sep="\t")

    # Vertical comparison between 2 ancestor genomes:
    vertical_comparison_with_filters(ham_analysis, "Platyhelminthes", "Digenea",
                                     retained_gene_ids, duplicated_gene_ids, gained_gene_ids,
                                     copy_threshold=2, first_ancestor_sp_threshold=11, second_ancestor_sp_threshold=8)

    # Fasta parsing:
    fgig_fasta_dict = fasta_parsing(fgig_fasta)
    fhep_fasta_dict = fasta_parsing(fhep_fasta)
    psim_fasta_dict = fasta_parsing(psim_fasta)
    sman_fasta_dict = fasta_parsing(sman_fasta)
    treg_fasta_dict = fasta_parsing(treg_fasta)
    tszi_fasta_dict = fasta_parsing(tszi_fasta)

    # Digenea ancestor genome model, fasta output:
    ancestor_genome_fasta(fgig_fasta_dict, digenea_ancestor_filtered_df, working_dir, "Fgigantica_100aa", "Digenea")
    ancestor_genome_fasta(fhep_fasta_dict, digenea_ancestor_filtered_df, working_dir, "Fhepatica_100aa", "Digenea")
    ancestor_genome_fasta(psim_fasta_dict, digenea_ancestor_filtered_df, working_dir, "Psimillimum_100aa", "Digenea")
    ancestor_genome_fasta(treg_fasta_dict, digenea_ancestor_filtered_df, working_dir, "Tregenti_100aa", "Digenea")
    ancestor_genome_fasta(tszi_fasta_dict, digenea_ancestor_filtered_df, working_dir, "Tszidati_100aa", "Digenea")
    # Digenea and Platyhelminthes ancestor genome models, Schistosoma mansoni fasta output:
    ancestor_genome_fasta(sman_fasta_dict, digenea_ancestor_filtered_df, working_dir, "Smansoni_100aa", "Digenea")
    ancestor_genome_fasta(sman_fasta_dict, platy_ancestor_filtered_df, working_dir, "Smansoni_100aa", "Platyhelminthes")

    # Vertical comparison output:
    comparison_results_fasta(sman_fasta_dict, retained_gene_ids, working_dir,
                             "Smansoni_100aa", "Platyhelminthes", "Digenea", "retained_genes")
    comparison_results_fasta(sman_fasta_dict, duplicated_gene_ids, working_dir,
                             "Smansoni_100aa", "Platyhelminthes", "Digenea", "duplicated_genes")
    comparison_results_fasta(sman_fasta_dict, gained_gene_ids, working_dir,
                             "Smansoni_100aa", "Platyhelminthes", "Digenea", "gained_genes")
    comparison_results_txt(retained_gene_ids, working_dir, "Platyhelminthes", "Digenea", "retained_genes")
    comparison_results_txt(duplicated_gene_ids, working_dir, "Platyhelminthes", "Digenea", "duplicated_genes")
    comparison_results_txt(gained_gene_ids, working_dir, "Platyhelminthes", "Digenea", "gained_genes")

    # TreeProfile:
    ham_analysis.create_tree_profile(outfile="{wd}/Flatworms_pyHAM_TreeProfile.all_genes_without_filters.html".format(
        wd=working_dir), as_html=True, export_with_histogram=True)