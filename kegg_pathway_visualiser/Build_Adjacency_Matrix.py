from Bio.KEGG.KGML.KGML_parser import read
import json
import numpy as np
import pandas as pd


def build_adjacency_matrix(pathway_file_name, symbols_path):
    '''
    This function is used to build an adjacency matrix for a given pathway. The adjacency matrix is built using a
    custom conversion dictionary that asigns a value to each relation type.

    param pathway_file_name: the name of the pathway file
    :param symbols_path: the path to the symbols json file
    :return:
    '''

    # loading geneIDs2Names.json
    def load_geneIDs2Names(path):
        with open(path, 'r') as fp:
            gene_IDs = json.load(fp)
        return gene_IDs


    def load_pathway(pathway_name):
        pathway = read(open(pathway_name, 'r'))
        return pathway


    # reading the pathway using BioPython
    pathway = load_pathway(f'media/{pathway_file_name}')


    def extract_pathway_genes(pathway):
        pathway_entries = []
        for entry_id in pathway.entries:
            li = pathway.entries[entry_id].name.split(" ")
            for item in li:
                if item != "undefined":
                    pathway_entries.extend(li)
        return pathway_entries

    pathway_genes = extract_pathway_genes(pathway)


    def count_relations(pathway):
        relations = {}
        for relation in pathway.relations:
            relation_type = relation.subtypes[0][0]
            if relation_type in relations:
                relations[relation_type] += 1
            else:
                relations[relation_type] = 1

        return relations

    relation_count = count_relations(pathway)
    print(relation_count)

    # print(len(pathway_genes))

    # visualising the pathway
    #from Bio.KEGG.KGML.KGML_pathway import Pathway
    #from Bio.Graphics.KGML_vis import KGMLCanvas
    #canvas = KGMLCanvas(pathway, import_imagemap=True, show_maps=True)
    #canvas.draw('hsa05200.pdf')

    # we get an input from the patient which are the genes that are mutated
    #mutated = ["caspase 9", "erb-b2 receptor tyrosine kinase 2"]

    # find id of mutated genes
    def get_gene_IDs(gene_names):
        gene_IDs = load_geneIDs2Names(symbols_path)
        selected_gene_IDs = []
        for gene_name in gene_names:
            if gene_name in gene_IDs.values():
                gene_ID = list(gene_IDs.keys())[list(gene_IDs.values()).index(gene_name)] # gene_IDs.key() that corresponds to gene_name
                selected_gene_IDs.append(gene_ID)
        return selected_gene_IDs


    #mutated_gene_IDs = get_gene_IDs(mutated)
    #print(mutated_gene_IDs)


    def relation_value_conversion(relation):
        conversion_dict = {"activation": 1, "inhibition": -1, "expression": 2, "repression": -2, "indirect effect": 3,
                           "state change": 4, "binding/association": 5, "dissociation": -5, "missing interaction": 6,
                           "phosphorylation": 7, "dephosphorylation": -7, "glycosylation": 8, "ubiquitination": 9,
                           "methylation": 10}

        relation_type = relation.subtypes[0][0]
        if relation_type in conversion_dict:
            v = conversion_dict[relation_type]
            return int(v)
        else:
            return 0

    def extract_genes_from_relations(pathway):
        rel_genes = []
        for relation in pathway.relations:
            entry1 = relation.entry1.name
            item1 = entry1.split(" ")
            entry2 = relation.entry2.name
            item2 = entry2.split(" ")
            rel_genes.append(item1)
            rel_genes.append(item2)

        # splitting the list of lists into a list
        rel_genes = [item for sublist in rel_genes for item in sublist]

        return rel_genes


    rel_genes = extract_genes_from_relations(pathway)


    # comparing the rel genes with the pathway genes
    def compare_genes(rel_genes, pathway_genes):
        different_genes = []
        for gene in rel_genes:
            if gene not in pathway_genes:
                different_genes.append(gene)

        return different_genes


    extra_items = compare_genes(pathway_genes, rel_genes)
    print(len(extra_items), "items from relations are not present in the pathway genes")
    print(extra_items[:-10])


    def pathway_adjacency_matrix(pathway, pathway_genes, rel_items):
        # creating a dataframe with the genes as columns and rows
        extra_items = compare_genes(pathway_genes, rel_items)
        all_names = pathway_genes + extra_items


        all_names = np.unique(all_names) # removing duplicates
        # removing all items that does not star with hsa
        # all_names = [x for x in all_names if x[:3] == "hsa"]

        df = pd.DataFrame(0, index=all_names, columns=all_names) # initialize the dataframe with zeros

        counter = 0
        # filling the dataframe with the relations
        for relation in pathway.relations:
            counter += 1
            rel1 = relation.entry1.name
            entry1 = rel1.split(" ")
            rel2 = relation.entry2.name
            entry2 = rel2.split(" ")
            for id1 in entry1:
                for id2 in entry2:
                    if id1 in df.index and id2 in df.columns:
                        df.loc[id1, id2] = relation_value_conversion(relation)

        return df

    matrix = pathway_adjacency_matrix(pathway, pathway_genes, rel_genes)

    # changing NaN to 0
    matrix = matrix.fillna(0)

    # changing all numbers to integers
    matrix = matrix.astype(int)

    print(matrix)
    name1 = pathway_file_name.split(".xml")[0]
    matrix_path = f"media/{name1}_adjacency_matrix.csv"
    matrix.to_csv(matrix_path)

    return matrix_path

