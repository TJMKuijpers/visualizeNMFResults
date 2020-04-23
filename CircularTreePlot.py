import networkx as nx
import matplotlib.pyplot as plt

try:
    from networkx import graphviz_layout
except ImportError:
    raise ImportError("This example needs Graphviz and either PyGraphviz or Pydot")

G.add_edges_from(edges)

nx.draw_networkx(G, with_label = True)