import networkx as nx
from shapely.geometry import Point, LineString
from geopy.distance import great_circle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap as lsc

def gdf_to_multidigraph(gdf, source_col='INT_ID_FRO', target_col='INT_ID_TO', crs=None):
    """
    Convert a GeoDataFrame to a NetworkX MultiDiGraph.

    Parameters:
    gdf (GeoDataFrame): Input GeoDataFrame containing street data.
    source_col (str): Column name for source nodes.
    target_col (str): Column name for target nodes.
    crs (str): Coordinate Reference System to set for the graph. Defaults to the CRS of the GeoDataFrame.

    Returns:
    MultiDiGraph: A NetworkX MultiDiGraph representing the street network.
    """
    # Convert the input GeoDataFrame 'gdf' to a NetworkX MultiDiGraph object
    graph = nx.from_pandas_edgelist(gdf, source=source_col, target=target_col, create_using=nx.MultiDiGraph())

    # Iterate through each row of the GeoDataFrame and add the attributes to the corresponding edge in the graph
    for _, row in gdf.iterrows():
        source = row[source_col]
        target = row[target_col]
        attr = row.to_dict()
        del attr[source_col]
        del attr[target_col]
        graph.add_edge(source, target, **attr)

    # Set the CRS of the graph, either from input or from the GeoDataFrame
    if crs is not None:
        graph.graph["crs"] = crs
    else:
        graph.graph["crs"] = gdf.crs.to_string() if gdf.crs is not None else None

    return graph




def add_node_coordinates(G, gdf_edges, source_col='INT_ID_FRO', target_col='INT_ID_TO', geom_col='geometry'):
    """
    Add node coordinates to a NetworkX graph from a GeoDataFrame.

    Parameters:
    G (networkx.Graph): The graph to which node coordinates will be added.
    gdf_edges (GeoDataFrame): The GeoDataFrame containing edge data with geometry.
    source_col (str): Column name for source nodes in gdf_edges.
    target_col (str): Column name for target nodes in gdf_edges.
    geom_col (str): Column name for geometry in gdf_edges.

    Returns:
    None
    """
    # Create a dictionary to store the coordinates of each node in the graph
    coords = {}

    # Iterate through each row of the GeoDataFrame containing the edges
    for _, row in gdf_edges.iterrows():
        source = row[source_col]
        target = row[target_col]
        source_geom = row[geom_col].coords[0]
        target_geom = row[geom_col].coords[-1]

        # Store the coordinates of the source and target nodes in the dictionary
        coords[source] = source_geom
        coords[target] = target_geom

    # Add the coordinates from the dictionary as attributes to the corresponding nodes in the graph
    for node, (x, y) in coords.items():
        G.nodes[node]['x'] = x
        G.nodes[node]['y'] = y




def collapse_multidigraph_to_digraph(multi_digraph):
    """
    Collapse a MultiDiGraph to a DiGraph, preserving edge attributes.

    Parameters:
    multi_digraph (networkx.MultiDiGraph): The input MultiDiGraph to be collapsed.

    Returns:
    networkx.DiGraph: A DiGraph with collapsed edges and preserved attributes.
    """
    # Create a new directed graph 'digraph' to store the collapsed edges
    digraph = nx.DiGraph()

    # Iterate through the edges of the input MultiDiGraph 'multi_digraph'
    for u, v, data in multi_digraph.edges(data=True):
        # If an edge already exists between nodes u and v in 'digraph', update its attributes
        if digraph.has_edge(u, v):
            digraph[u][v].update(data)
        # Otherwise, add the edge with the corresponding attributes to 'digraph'
        else:
            digraph.add_edge(u, v, **data)

    # Add nodes from the input 'multi_digraph' to the new 'digraph', preserving their attributes
    digraph.add_nodes_from(multi_digraph.nodes(data=True))

    # If the input 'multi_digraph' has a "crs" attribute, copy it to the new 'digraph'
    if "crs" in multi_digraph.graph:
        digraph.graph["crs"] = multi_digraph.graph["crs"]

    return digraph




def nearest_edge(point, G):
    """
    Find the nearest edge in a graph to a given point.

    Parameters:
    point (tuple): A tuple representing the coordinates of the point (latitude, longitude).
    G (networkx.Graph): The graph containing the edges.

    Returns:
    tuple: The nodes and key of the nearest edge (u, v, k).
    """
    min_distance = float('inf')
    nearest_u, nearest_v, nearest_k = None, None, None
    
    # Iterate through the edges of the graph
    for u, v, k, data in G.edges(data=True, keys=True):
        # Get the coordinates of the nodes u and v
        u_coord = (G.nodes[u]['y'], G.nodes[u]['x'])
        v_coord = (G.nodes[v]['y'], G.nodes[v]['x'])
        
        # Create a LineString representing the edge between nodes u and v
        line = LineString([u_coord, v_coord])
        
        # Calculate the distance between the point and the nearest point on the line
        distance = great_circle(point, line.interpolate(line.project(Point(point), normalized=True)).coords[0]).meters
        
        # Update the minimum distance and the nearest edge if the current distance is smaller than the previously found minimum distance
        if distance < min_distance:
            min_distance = distance
            nearest_u, nearest_v, nearest_k = u, v, k
            
    return nearest_u, nearest_v, nearest_k




# Helper function to generate colormaps for GWR plots
def get_cmap(values, cmap_name='coolwarm', n=256):
    """
    Generate a colormap for GWR plots based on the range of values.

    Parameters:
    values (array-like): Array of values to determine the colormap range.
    cmap_name (str): Name of the colormap to use.
    n (int): Number of color bins in the colormap.

    Returns:
    LinearSegmentedColormap: The generated colormap.
    """
    name = f'{cmap_name}_new'
    cmap = plt.cm.get_cmap(cmap_name)
    vmin = values.min()
    vmax = values.max()

    if vmax < 0:
        # if all values are negative, use the negative half of the colormap
        return lsc.from_list(name, cmap(np.linspace(0, 0.5, n)))
    elif vmin > 0:
        # if all values are positive use the positive half of the colormap
        return lsc.from_list(name, cmap(np.linspace(0.5, 1, n)))
    else:
        # otherwise there are positive and negative values so use zero as midpoint
        # and truncate the colormap such that if the original spans ± the greatest
        # absolute value, use colors from it spanning vmin to vmax
        abs_max = max(values.abs())
        start = (vmin + abs_max) / (abs_max * 2)
        stop = (vmax + abs_max) / (abs_max * 2)
        return lsc.from_list(name, cmap(np.linspace(start, stop, n)))