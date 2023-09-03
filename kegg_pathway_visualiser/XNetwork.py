import json

import networkx as nx
import matplotlib.pyplot as plt
import csv
import numpy as np
import pandas as pd
import os
import networkx.algorithms.community
import gravis as gv

matrix = "adjacency_matrix.csv"
with open(matrix, newline='') as f:
    reader = csv.reader(f)
    data = list(reader)


def get_gene_IDs(gene_names):
    # loading geneIDs2Names.json
    with open('genesIDs2Names.json', 'r') as fp:
        gene_IDs = json.load(fp)
    selected_gene_IDs = []
    for gene_name in gene_names:
        if gene_name in gene_IDs.values():
            gene_ID = list(gene_IDs.keys())[list(gene_IDs.values()).index(gene_name)] # gene_IDs.key() that corresponds to gene_name
            selected_gene_IDs.append(gene_ID)
    return selected_gene_IDs


def get_gene_name(gene_ID):
    # loading geneIDs2Names.json
    with open('genesIDs2Names.json', 'r') as fp:
        gene_IDs = json.load(fp)
    if gene_ID in gene_IDs.keys():
        gene_name = gene_IDs[gene_ID]
        return gene_name
    else:
        return gene_ID


def count_relations(pathway, target_genes=None):
    relations = {}
    for relation in pathway.relations:
        relation_type = relation.subtypes[0][0]
        if target_genes is None or relation.entry1.name in target_genes or relation.entry2.name in target_genes:
            if relation_type in relations:
                relations[relation_type] += 1
            else:
                relations[relation_type] = 1

    return relations


def relation_name(int):
    conversion_dict = {"activation": 1, "inhibition": -1, "expression": 2, "repression": -2, "indirect effect": 3,
                           "state change": 4, "binding/association": 5, "dissociation": -5, "missing interaction": 6,
                           "phosphorylation": 7, "dephosphorylation": -7, "glycosylation": 8, "ubiquitination": 9,
                           "methylation": 10}

    return list(conversion_dict.keys())[list(conversion_dict.values()).index(int)]


def Graph(subset, targets=None, target_rel=None, target_generations=2):
    if target_rel is None:
        target_rel = []
    if targets is None:
        targets = []

    # creating a graph from the matrix using DiGraph
    G = nx.DiGraph()

    relations = {}
    generation = 0
    while generation < target_generations:
        generation += 1
        children = []

        # ignoring the first row
        for row in range(1, len(subset)):
            # ignoring the first column
            for col in range(1, len(subset)):
                row_name = subset[row][0]
                col_name = subset[col][0]
                if row_name in targets or col_name in targets or targets == []:
                    rel_value = int(subset[row][col])
                    if rel_value != 0:
                        rel_name = relation_name(rel_value)
                        if target_rel == [] or rel_name in target_rel:

                            if row_name not in targets:
                                children.append(row_name)
                            if col_name not in targets:
                                children.append(col_name)

                            G.add_edge(get_gene_name(row_name), get_gene_name(col_name), weight=rel_value, label=rel_name)

                            if rel_name in relations:
                                relations[rel_name] += 1
                            else:
                                relations[rel_name] = 1

        targets = children + targets
        print("Generation: " + str(generation), "Targets: " + str(len(targets)), "Children: " + str(len(children)))

    return G, relations


# target genes
target_genes = ["TP53", "MDM2"]

target_IDs = get_gene_IDs(target_genes)

subset = data

target_relations = ["activation", "inhibition", "expression", "repression", "indirect effect", "state change",
                    "binding/association", "dissociation", "missing interaction", "phosphorylation", "dephosphorylation",
                    "glycosylation", "ubiquitination", "methylation"]

G, relations = Graph(subset, targets=target_IDs, target_rel=target_relations, target_generations=6)

print(relations)

# increasing space between nodes
pos = nx.spring_layout(G, k=2, seed=8)

color = {1: "red", -1: "blue", 2: "green", -2: "yellow", 3: "orange", 4: "purple", 5: "pink", -5: "brown", 6: "black",
         7: "gray", -7: "cyan", 8: "magenta", 9: "lime", 10: "olive"}


conversion_dict = {"activation": 1, "inhibition": -1, "expression": 2, "repression": -2, "indirect effect": 3,
                   "state change": 4, "binding/association": 5, "dissociation": -5, "missing interaction": 6,
                   "phosphorylation": 7, "dephosphorylation": -7, "glycosylation": 8, "ubiquitination": 9, "methylation": 10}

reversed_conversion_dict = {v: k for k, v in conversion_dict.items()}

options = {
    'node_color': 'white',
    'node_size': 1000,
    'width': 1,
    'arrowstyle': '-|>',
    'arrowsize': 12,
    'edge_color': [color[G[u][v]['weight']] for u, v in G.edges],
    'font_size': 10,
    'font_color': 'black',
    'font_weight': 'bold',
    'style': 'dashed',
    'labels': {e: e for e in G.nodes}
}

nx.draw(G, pos, **options)

# legend is target relations
# the color is in colors dictionary and values are in reversed_conversion_dict
plt.legend(handles=[plt.Line2D([0], [0], color=list(color.values())[list(reversed_conversion_dict.values()).index(i)], lw=4, label=i) for i in target_relations], loc='upper right')


plt.show()



# change edge colors according to the color dictionary
color_dict = {1: "red", -1: "blue", 2: "green", -2: "yellow", 3: "orange", 4: "purple", 5: "pink", -5: "brown", 6: "black",
            7: "gray", -7: "cyan", 8: "magenta", 9: "lime", 10: "olive"}
for u,v in G.edges:
    G.edges[u,v].update({'color': color_dict[G.edges[u,v]['weight']]})

gravis_fig = gv.d3(G, graph_height=200)

html = gravis_fig.to_html()
with open("gravis.html", "w", encoding="utf-8") as file:
    file.write(html)



