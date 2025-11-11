def assemble(G):
    contigs = []
    for g in G:
        seq = ''
        init = True
        for n in nx.eulerian_path(g):
            if init: 
                seq = n[0] + n[1][-1]
                init = False
                continue
            seq += n[1][-1]
        contigs += [seq]
    return contigs