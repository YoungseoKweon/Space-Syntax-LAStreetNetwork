import networkx as nx
from itertools import combinations

### connectivity
def calculate_connectivity(G_d):
    """
    Calculate the connectivity value for each edge in the graph.

    Parameters:
    G_d (networkx.Graph): The graph for which to calculate edge connectivity.

    Returns:
    dict: A dictionary where keys are edges (u, v, key) and values are the connectivity values.
    """
    # Create an empty dictionary to store connectivity values
    connectivity_dict = {}

    # Iterate through each edge in G_d
    for u, v, key, data in G_d.edges(keys=True, data=True):
        # Calculate the number of out-edges from the starting and ending nodes of the edge
        connectivity_u = len(list(G_d.out_edges(u)))
        connectivity_v = len(list(G_d.out_edges(v)))
        # Sum the out-edge counts for the starting and ending nodes to get the edge's connectivity
        connectivity_sum = connectivity_u + connectivity_v
        # Store the connectivity value in the dictionary
        connectivity_dict[(u, v, key)] = connectivity_sum

    return connectivity_dict




### choice
def calculate_edge_choice(G_d):
    """
    Calculate the choice value for each edge in the graph.

    Parameters:
    G_d (networkx.DiGraph): The directed graph for which to calculate edge choice values.

    Returns:
    dict: A dictionary where keys are edges (u, v, k) and values are the choice values.
    """
    # Initialize edge counts
    edge_counts = {(u, v, k): 0 for u, v, k in G_d.edges(keys=True)}

    nodes = list(G_d.nodes())
    for source, target in combinations(nodes, 2):
        if source != target:
            try:
                # Calculate all shortest paths from source to target
                shortest_paths = list(nx.all_shortest_paths(G_d, source=source, target=target))
                for path in shortest_paths:
                    for u, v in zip(path[:-1], path[1:]):
                        for key in G_d[u][v]:
                            edge_counts[(u, v, key)] += 1
            except nx.NetworkXNoPath:
                pass

    # Calculate the edge choice values
    edge_choice = {}
    total_paths = sum(edge_counts.values())
    for edge, count in edge_counts.items():
        edge_choice[edge] = count / total_paths if total_paths > 0 else 0

    return edge_choice




### centrality
def calculate_centrality(G_d, weight='length'):
    """
    Calculate the closeness centrality for each node in the graph.

    Parameters:
    G_d (networkx.Graph): The graph for which to calculate node closeness centrality.
    weight (str): The edge attribute to use as weight. Defaults to 'length'.

    Returns:
    dict: A dictionary where keys are nodes and values are the closeness centrality.
    """
    # Calculate shortest path lengths
    shortest_path_lengths = dict(nx.shortest_path_length(G_d, weight=weight))

    # Calculate closeness centrality
    cc = {}
    for n in G_d.nodes():
        if n in shortest_path_lengths:
            if len(shortest_path_lengths[n]) > 1:
                cc[n] = sum(shortest_path_lengths[n].values()) / (len(shortest_path_lengths[n]) - 1)
            else:
                cc[n] = 0
        else:
            cc[n] = 0

    return cc
