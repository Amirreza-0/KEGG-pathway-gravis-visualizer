import plotly.graph_objects as go
import json
import networkx as nx
import csv


def generate_plotly_html(matrix_path, symbols_path, target_genes=None, target_relations=None, target_generations=2):

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

    def get_gene_symbols(gene_ID):
        # loading geneIDs2Names.json
        with open(symbols_path, 'r') as fp:
            gene_IDs = json.load(fp)
        if gene_ID in gene_IDs.keys():
            gene_name = gene_IDs[gene_ID]
            return gene_name
        else:
            return gene_ID



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

                                G.add_edge(get_gene_symbols(row_name), get_gene_symbols(col_name), weight=rel_value,
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

    G2, relations = Graph(subset, targets=target_IDs, target_rel=target_rels, target_generations=target_generations)
    G = nx.random_geometric_graph(200, 0.125)

    # print(relations)

    def generate_plotly_Network_Graph(Graph):
        # extracting the positions of the nodes
        pos = nx.spring_layout(Graph)

        # adding positions to the nodes
        for node in Graph.nodes:
            Graph.nodes[node]['pos'] = list(pos[node])

        # extracting the positions of the nodes
        edge_x = []
        edge_y = []
        for edge in Graph.edges():
            x0, y0 = Graph.nodes[edge[0]]['pos']
            x1, y1 = Graph.nodes[edge[1]]['pos']
            edge_x.append(x0)
            edge_x.append(x1)
            edge_x.append(None)
            edge_y.append(y0)
            edge_y.append(y1)
            edge_y.append(None)

        edge_trace = go.Scatter(
            x=edge_x, y=edge_y,
            line=dict(width=0.5, color='#888'),
            hoverinfo='none',
            mode='lines')

        node_x = []
        node_y = []
        for node in Graph.nodes():
            x, y = Graph.nodes[node]['pos']
            node_x.append(x)
            node_y.append(y)

        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode='markers',
            hoverinfo='text',
            marker=dict(
                showscale=True,
                # colorscale options
                # 'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu'
                # 'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet'
                # 'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
                colorscale='YlGnBu',
                reversescale=True,
                color=[],
                size=10,
                colorbar=dict(
                    thickness=15,
                    title='Node Connections',
                    xanchor='left',
                    titleside='right'
                ),
                line_width=2))

        node_adjacencies = []
        node_text = []
        for node, adjacencies in enumerate(Graph.adjacency()):
            node_adjacencies.append(len(adjacencies[1]))
            node_text.append(get_gene_symbols(adjacencies[0]) + ': ' + str(len(adjacencies[1])))

        node_trace.marker.color = node_adjacencies
        node_trace.text = node_text

        fig = go.Figure(data=[edge_trace, node_trace],
                        layout=go.Layout(
                            title='<br>Network graph made with Python',
                            showlegend=False,
                            hovermode='closest',
                            margin=dict(b=20, l=5, r=5, t=40),
                            annotations=[dict(
                                text="Python code: <a href='https://plotly.com/ipython-notebooks/network-graphs/'> https://plotly.com/ipython-notebooks/network-graphs/</a>",
                                showarrow=False,
                                xref="paper", yref="paper",
                                x=0.005, y=-0.002)],
                            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))

        return fig

    figure = generate_plotly_Network_Graph(G2)
    figure.write_html("undirected-graph.html")

    def generate_plotly_directed_Network_Graph(Graph):
        # extracting the positions of the nodes
        pos = nx.spring_layout(Graph)

        # adding positions to the nodes
        for node in Graph.nodes:
            Graph.nodes[node]['pos'] = list(pos[node])

        # including all edge types
        egde_types = ["activation", "inhibition", "expression", "repression", "indirect effect", "state change",
                      "binding/association", "dissociation", "missing interaction", "phosphorylation",
                      "dephosphorylation", "glycosylation", "ubiquitination", "methylation"]
        edge_colors = ["#00FF00", "#FF0000", "#0000FF", "#FF00FF", "#00FFFF", "#FFFF00", "#000000", "#FFFFFF",
                       "#808080", "#800000", "#008000", "#000080", "#800080", "#008080"]

        def edge_trace(edge_type, edge_color):
            edge_x = []
            edge_y = []
            for edge in Graph.edges():
                if Graph.edges[edge]['label'] == edge_type:
                    x0, y0 = Graph.nodes[edge[0]]['pos']
                    x1, y1 = Graph.nodes[edge[1]]['pos']
                    edge_x.append(x0)
                    edge_x.append(x1)
                    edge_x.append(None)
                    edge_y.append(y0)
                    edge_y.append(y1)
                    edge_y.append(None)

            edge_trace = go.Scatter(
                x=edge_x, y=edge_y,
                line=dict(width=0.5, color=edge_color),
                hoverinfo='none',
                mode='lines')

            return edge_trace

        edge_traces = []
        for i in range(len(egde_types)):
            edge_traces.append(edge_trace(egde_types[i], edge_colors[i]))

        node_x = []
        node_y = []
        for node in Graph.nodes():
            x, y = Graph.nodes[node]['pos']
            node_x.append(x)
            node_y.append(y)

        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode='markers',
            hoverinfo='text',
            marker=dict(
                showscale=True,
                # colorscale options
                # 'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu'
                # 'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet'
                # 'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
                colorscale='YlGnBu',
                reversescale=True,
                color=[],
                size=10,
                colorbar=dict(
                    thickness=15,
                    title='Node Connections',
                    xanchor='left',
                    titleside='right'
                ),
                line_width=2))

        node_adjacencies = []
        node_text = []
        for node, adjacencies in enumerate(Graph.adjacency()):
            node_adjacencies.append(len(adjacencies[1]))
            node_text.append(get_gene_symbols(adjacencies[0]) + ': ' + str(len(adjacencies[1])))

        node_trace.marker.color = node_adjacencies
        node_trace.text = node_text

        fig = go.Figure(data=edge_traces + [node_trace],
                        layout=go.Layout(
                            title='<br>Network graph made with Python',
                            showlegend=False,
                            hovermode='closest',
                            margin=dict(b=20, l=5, r=5, t=40),
                            annotations=[dict(
                                text="Python code: <a href='https://plotly.com/ipython-notebooks/network-graphs/'> https://plotly.com/ipython-notebooks/network-graphs/</a>",
                                showarrow=False,
                                xref="paper", yref="paper",
                                x=0.005, y=-0.002)],
                            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))

        return fig

    figure = generate_plotly_directed_Network_Graph(G2)
    # figure.write_html("directed-graph.html")
    html = figure.to_html()
    html.encode("utf-8")

    return html


"""
G = nx.random_geometric_graph(200, 0.125)



edge_x = []
edge_y = []
for edge in G.edges():
    x0, y0 = G.nodes[edge[0]]['pos']
    x1, y1 = G.nodes[edge[1]]['pos']
    edge_x.append(x0)
    edge_x.append(x1)
    edge_x.append(None)
    edge_y.append(y0)
    edge_y.append(y1)
    edge_y.append(None)

edge_trace = go.Scatter(
    x=edge_x, y=edge_y,
    line=dict(width=0.5, color='#888'),
    hoverinfo='none',
    mode='lines')

node_x = []
node_y = []
for node in G.nodes():
    x, y = G.nodes[node]['pos']
    node_x.append(x)
    node_y.append(y)

node_trace = go.Scatter(
    x=node_x, y=node_y,
    mode='markers',
    hoverinfo='text',
    marker=dict(
        showscale=True,
        # colorscale options
        #'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
        #'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
        #'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
        colorscale='YlGnBu',
        reversescale=True,
        color=[],
        size=10,
        colorbar=dict(
            thickness=15,
            title='Node Connections',
            xanchor='left',
            titleside='right'
        ),
        line_width=2))

node_adjacencies = []
node_text = []
for node, adjacencies in enumerate(G.adjacency()):
    node_adjacencies.append(len(adjacencies[1]))
    node_text.append('# of connections: '+str(len(adjacencies[1])))

node_trace.marker.color = node_adjacencies
node_trace.text = node_text

fig = go.Figure(data=[edge_trace, node_trace],
             layout=go.Layout(
                title='<br>Network graph made with Python',
                titlefont_size=16,
                showlegend=False,
                hovermode='closest',
                margin=dict(b=20,l=5,r=5,t=40),
                annotations=[ dict(
                    text="Python code: <a href='https://plotly.com/ipython-notebooks/network-graphs/'> https://plotly.com/ipython-notebooks/network-graphs/</a>",
                    showarrow=False,
                    xref="paper", yref="paper",
                    x=0.005, y=-0.002 ) ],
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                )


# exporting the plot in a HTML file
fig.write_html("network_graph.html")
"""
