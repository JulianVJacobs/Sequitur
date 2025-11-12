"""Simple example showing how to call sequitur_rs.analyse_reads from the project venv.

This prints adjacency, confidences, and detected cycles as JSON. If networkx is installed, it will
also draw the adjacency graph with edge widths proportional to confidence.
"""
import json

try:
    import sequitur_rs as s
except Exception as e:
    print("Failed to import sequitur_rs:", e)
    raise

reads = ["ACGT", "GTAA", "TACG"]
info = s.analyse_reads(reads)
print(json.dumps(info, indent=2))

# optional visualization
try:
    import networkx as nx
    import matplotlib.pyplot as plt
    G = nx.DiGraph()
    for i, row in enumerate(info['adjacency']):
        for entry in row:
            dst, weight, overlap = entry
            conf_row = info['confidences'][i]
            # find conf for dst
            conf = next((c for (d, c) in conf_row if d == dst), 0.0)
            G.add_edge(i, dst, weight=weight, conf=conf, overlap=overlap)
    pos = nx.spring_layout(G)
    widths = [G[u][v]['conf']*5 for u,v in G.edges()]
    nx.draw(G, pos, with_labels=True, width=widths)
    edge_labels = {(u,v): f"{G[u][v]['overlap']}" for u,v in G.edges()}
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
    plt.show()
except ImportError:
    pass
