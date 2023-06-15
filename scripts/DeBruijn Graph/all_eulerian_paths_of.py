def all_eulerian_paths_of(G):
    paths = []
    g = nx.DiGraph()
    while len(list(nx.selfloop_edges(G))):
        g.add_edge(list(nx.selfloop_edges(G))[0][0],list(nx.selfloop_edges(G))[0][1])
        G.remove_edges_from(g.edges)
        paths += [g.copy()]
        g.clear()
    n = min(i[1] for i in G.in_degree())
    while n <= max(i[1] for i in G.in_degree()):
        if not len(g): 
            if len(nx.subgraph_view(G,filter_node=lambda node: G.in_degree(node)==n)): 
                if len(G.in_edges(nx.subgraph_view(G,filter_node=lambda node: G.in_degree(node)==n).nodes())): 
                    edge = list(G.in_edges(nx.subgraph_view(G,filter_node=lambda node: G.in_degree(node)==n).nodes()))[0]
                    if len(set(G.in_edges(nx.subgraph_view(G,filter_node=lambda node: G.in_degree(node)==n).nodes())).intersection(G.in_edges(edge[0]))): 
                        edge = list(set(G.in_edges(nx.subgraph_view(G,filter_node=lambda node: G.in_degree(node)==n).nodes())).intersection(G.in_edges(edge[0])))[0]
                else: edge = list(G.out_edges(nx.subgraph_view(G,filter_node=lambda node: G.in_degree(node)==n).nodes()))[0]
            else: 
                n+=1
                continue
        else: edge = list(G.out_edges([edge[1]]))[0]
        g.add_edge(edge[0],edge[1])
        if G.out_degree(edge[1]) != 1:
            G.remove_edges_from(g.edges)
            paths += [g.copy()]
            g.clear()
            if not len(G.edges): break
            n = min(i[1] for i in nx.subgraph_view(G,filter_node=lambda node: G.out_degree(node) > 0).in_degree())
    return paths