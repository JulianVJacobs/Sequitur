def construct_debruijn_graph(reads,k=3):
    import networkx as nx
    
    G = nx.MultiDiGraph()
    for read in reads:
        for i in range(len(read)-k+1):
            G.add_edge(read[i:i+k-1],read[i+1:i+k])
    return G