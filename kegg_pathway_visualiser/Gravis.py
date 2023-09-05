import json
import networkx as nx
import matplotlib.pyplot as plt
import csv
import gravis as gv


def generate_gravis_html(matrix_path, symbols_path, target_genes=None, target_relations=None, target_generations=2):
    """
    This function takes a matrix and generates a graph using gravis. The graph is then converted to html and returned.

    :param matrix_path: adjacency matrix path
    :param symbols_path: symbols json file path
    :param target_genes: Eg. TP53,EGFR
    :param target_relations: Eg. activation,inhibition
    :param target_generations: Eg. 2
    :return:
    """


    target_genes = target_genes[0]
    print("target_entries: ", target_genes)

    if target_genes is None:
        target_genes = []
        with open(symbols_path, 'r') as fp:
            gene_IDs = json.load(fp)
        for gene_ID in gene_IDs.keys():
            target_genes.append(gene_IDs[gene_ID])

    if target_relations is None:
        target_relations = ["activation", "inhibition", "expression", "repression", "indirect effect", "state change",
                            "binding/association", "dissociation", "missing interaction", "phosphorylation",
                            "dephosphorylation", "glycosylation", "ubiquitination", "methylation"]

    target_relations = target_relations[0]
    print(target_relations, type(target_relations))

    if target_generations is None:
        target_generations = 10


    with open(matrix_path, newline='') as f:
        reader = csv.reader(f)
        data = list(reader)

    def get_gene_IDs(gene_names):
        # loading geneIDs2Names.json
        with open(symbols_path, 'r') as fp:
            gene_IDs = json.load(fp)
        selected_gene_IDs = []
        for gene_name in gene_names:
            if gene_name in gene_IDs.values():
                gene_ID = list(gene_IDs.keys())[
                    list(gene_IDs.values()).index(gene_name)]  # gene_IDs.key() that corresponds to gene_name
                selected_gene_IDs.append(gene_ID)
        return selected_gene_IDs

    def get_gene_name(gene_ID):
        # loading geneIDs2Names.json
        with open(symbols_path, 'r') as fp:
            gene_IDs = json.load(fp)
        if gene_ID in gene_IDs.keys():
            gene_name = gene_IDs[gene_ID]
            return gene_name
        else:
            return gene_ID

    def count_relations(pathway, target_entries=None):
        relations = {}
        for relation in pathway.relations:
            relation_type = relation.subtypes[0][0]
            if target_entries is None or relation.entry1.name in target_entries or relation.entry2.name in target_entries:
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

    def Graph(subset, targets=None, target_rel=None, target_generations=None):
        if target_rel is None:
            target_rel = []
        if targets is None:
            targets = []
        if target_generations is None:
            target_generations = 10


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

                                G.add_edge(get_gene_name(row_name), get_gene_name(col_name), weight=rel_value,
                                           label=rel_name)

                                if rel_name in relations:
                                    relations[rel_name] += 1
                                else:
                                    relations[rel_name] = 1

            targets = children + targets
            print("Generation: " + str(generation), "Targets: " + str(len(targets)), "Children: " + str(len(children)))

        return G, relations

    # target genes
    target_genes = target_genes

    target_IDs = get_gene_IDs(target_genes)

    subset = data

    target_rels = target_relations

    G, relations = Graph(subset, targets=target_IDs, target_rel=target_rels, target_generations=target_generations)

    print(relations)

    # increasing space between nodes
    pos = nx.spring_layout(G, k=2, seed=8)

    color = {1: "red", -1: "blue", 2: "green", -2: "yellow", 3: "orange", 4: "purple", 5: "pink", -5: "brown",
             6: "black",
             7: "gray", -7: "cyan", 8: "magenta", 9: "lime", 10: "olive"}

    conversion_dict = {"activation": 1, "inhibition": -1, "expression": 2, "repression": -2, "indirect effect": 3,
                       "state change": 4, "binding/association": 5, "dissociation": -5, "missing interaction": 6,
                       "phosphorylation": 7, "dephosphorylation": -7, "glycosylation": 8, "ubiquitination": 9,
                       "methylation": 10}

    reversed_conversion_dict = {v: k for k, v in conversion_dict.items()}

    options = {
        'node_color': 'white',
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
    # node['shape'] = 'rectangle' if node in target_entries else 'circle'
    for node in G.nodes:
        if node in target_genes:
            G.nodes[node]['shape'] = 'rectangle'
            G.nodes[node]['color'] = 'red'
        else:
            G.nodes[node]['shape'] = 'circle'

    nxfig = nx.draw(G, pos, **options)

    # legend is target relations
    # the color is in colors dictionary and values are in reversed_conversion_dict
    plt.legend(handles=[
        plt.Line2D([0], [0], color=list(color.values())[list(reversed_conversion_dict.values()).index(i)], lw=4,
                   label=i) for i in target_relations], loc='upper right')

    # nxfig to svg
    #svg_path = f"{symbols_path.split('_symbols.json')[0]}.svg"
    #nxsvg = plt.savefig(svg_path, format="svg")

    # nxfig to png
    # matrix name: hsa05200_yEZHf3m_adjacency_matrix.csv
    #file_name = matrix_path.split('_adjacency_matrix.csv')[0]
    #png_path = f"{file_name}.png"
    print("symptoms_path :", symbols_path)
    print("matrix_path :", matrix_path)

    file_name = matrix_path.split("/")[-1].split("_adjacency_matrix.csv")[0]
    print("file_name :", file_name)
    if file_name.split("_") is not None:
        file_name = file_name.replace("/", "")

    png_path = f"media/{file_name}.png"
    plt.savefig(png_path, format="png")
    print("saved png file")

    # change edge colors according to the color dictionary
    color_dict = {1: "red", -1: "blue", 2: "green", -2: "yellow", 3: "orange", 4: "purple", 5: "pink", -5: "brown",
                  6: "black",
                  7: "gray", -7: "cyan", 8: "magenta", 9: "lime", 10: "olive"}
    for u, v in G.edges:
        G.edges[u, v].update({'color': color_dict[G.edges[u, v]['weight']]})

    gravis_fig = gv.d3(G, graph_height=400)

    html = gravis_fig.to_html()
    html.encode('utf-8')

    return html, png_path
