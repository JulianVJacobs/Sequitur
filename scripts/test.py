#!/usr/bin/env python3
def generate_reads(sequence,min_subseq_len,max_subseq_len,min_overlap,max_overlap,min_coverage=None,circularise=False,seed=None,shuffle=True):
    """
    DESCRIPTION 
        Utility function that chops a sequence into several reads with bounded random lengths that 
        have a bounded random overlap
    INPUT
        sequence       | a sequence of characters that will be divided into overlapping subsequences
        min_subseq_len | the shortest length a subsequence can have
        max_subseq_len | the longest length a subsequence can have
        min_overlap    | the shortest overlap two subsequences can share
        max_overlap    | the longest overlap two subsequences can share
        circularize    | boolean indicating whether to add a random amount of the end of the sequence
                    | to the beginning and vice versa
        seed           | random seed for the random function for reproducibility
    OUTPUT
        A list of overlapping reads of random bounded size which share a bounded random amount of
        overlap
    """
    import random

    random.seed(seed)
    if circularise: sequence = sequence[-random.randint(min_overlap,max_overlap):] + sequence + sequence[:random.randint(min_overlap,max_overlap)]
    reads = []
    while 1: 
        start = 0
        end = random.randint(min_subseq_len,max_subseq_len)
        reads += [sequence[start:end]]
        while end < len(sequence):
            start = random.randint(end-max_overlap,end-min_overlap)
            if (len(sequence) - start)/max_subseq_len < 2:
                if (len(sequence) - start)/max_subseq_len < 1:
                    end = len(sequence)
                else:
                    a = 0
                    while (len(sequence) - start)/(min_subseq_len+a) > 2: a+=1
                    end = random.randint(start+min_subseq_len+a,start+max_subseq_len) 
            else: end = random.randint(start+min_subseq_len,start+max_subseq_len) 
            reads += [sequence[start:end]]
        if min_coverage is None or len(set(reads))*(sum(len(read) for read in set(reads))/len(set(reads)))/len(sequence) >= min_coverage:
            if not shuffle: return reads
            reads_ = reads[:]
            random.shuffle(reads_)
            return reads_, list(reads_.index(read) for read in reads)

def generate_genome_sequence(n,palindrome=False,seed=None):
    """
    DESCRIPTION 
        Utility function that creates a random sequence containing only the letters A, T, G, and C
    INPUT
        n          | the length of the sequence
        palindrome | a boolean indicating whether the sequence must be a palidrome or not
        seed       | random seed for the random function for reproducibility
    OUTPUT
        A random sequence of length n
    """
    import random
    from math import fmod, ceil
    
    random.seed(seed)
    nucleotides = {1:'A',2:'C',3:'G',4:'T'}
    seq = ''
    if palindrome: n = ceil(n/2)
    for _ in range(n):
        seq += nucleotides[random.randint(1,4)]
    if palindrome: seq += ''.join(reversed(seq[:int(n-fmod(n,2))]))
    return seq

def find_longest_overlap(reads):
    import pylcs

    overlaps = []
    for read in reads:
        overlaps += pylcs.lcs2_of_list(read, list(set(reads).symmetric_difference([read])))
    return min(overlaps),max(overlaps)

def create_de_bruijn_graph(k, sequences,do_time = False):
    """
    Create a de Bruijn graph from a set of DNA sequences.
    
    Parameters:
    - k (int): k-mer size
    - sequences (list): List of DNA sequences
    - do_time (bool): should the function be time of not (default False)
    
    Returns:
    - Union[nx.Graph,(nx.Graph,float)]: De Bruijn graph or a tuple of the De Bruijn graph and the time
                                        to create it.
    """
    import networkx as nx
    import time

    if do_time: start_time = time.time()
    graph = nx.DiGraph()

    for sequence in sequences:
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            prefix = kmer[:-1]
            suffix = kmer[1:]
            
            if not graph.has_edge(prefix, suffix):
                graph.add_edge(prefix, suffix, weight=1)
            else:
                graph[prefix][suffix]['weight'] += 1

    if do_time: return graph, time.time() - start_time
    return graph

def eulerian_path(graph,do_time = False):
    """
    Find an Eulerian path in the given graph.
    
    Parameters:
    - graph (nx.Graph): De Bruijn graph
    - do_time (bool): should the function be time of not (default False)
    
    Returns:
    - Union[str,(str,float)]: string describing the eulerian path or a tuple of this 
                              string and the execution time.
    """
    import networkx as nx
    import time

    if do_time: start_time = time.time()
    path = []
    
    if nx.has_eulerian_path(graph):
        for node in nx.eulerian_path(graph):
            if len(path): path.append(node[0][-1])
            else: path.append(node[0])
        
        if len(node[1]): path.append(node[1][-1])
        else: path.append(node[1])
        
        if do_time: return ''.join(path), time.time() - start_time
        return ''.join(path)
    else: 
        if time: return '', time.time() - start_time

def normalised_damerau_levenshtein_distance(read,overlap):
    """
    Find the Damerau-Levenshtein edit distance of two strings normalised to the length
    of the shorter string. This normalisation is because we want to path prefixes to
    suffixes and this means that in general we will be comparing a full string to a
    portion of another string.
    
    Parameters:
    - read (str): string for comparison, usually the longer string 
    - overlap (str): string for comparison, usually the shorter string
    
    Returns:
    - float: the normalised Demarau-Levenshtein edit distance of the input strings
    """
    from jellyfish import damerau_levenshtein_distance
    return damerau_levenshtein_distance(read.__str__()[:min(len(overlap),len(read))],overlap.__str__()[:min(len(overlap),len(read))])/min(len(overlap),len(read))

def build_suffix_array(reads, min_suf_len = 3,do_time = False):
    import time

    if do_time: start = time.time()
    suf_arr = []
    for read in reads:
        read += '$' + str(reads.index(read))
        for i in range(len(read)-min_suf_len-1):
            # if len(read[i:]) < min_suf_len + 2: continue 
            suf_arr += [read[i:]]
    suf_arr.sort()
    suf_arr_ind = []
    for s in range(len(suf_arr)):
        suf_arr_ind += [int(suf_arr[s].split('$')[-1].__str__())]
        suf_arr[s] = suf_arr[s][:suf_arr[s].find('$')+1]
    if do_time: return suf_arr, suf_arr_ind, time.time() - start
    return suf_arr,suf_arr_ind

def create_bipartite_adjacency_matrix(reads, suf_arr = None, suf_arr_ind = None, do_time = False,max_diff = 0.25, min_suf_len = 3):
    import time
    
    if do_time: start = time.time()
    elif suf_arr is None or suf_arr_ind is None: suf_arr,suf_arr_ind = build_suffix_array(reads,min_suf_len=min_suf_len)
    reads_map = dict(zip(reads,list(range(len(reads)))))
    B = {}
    for read in reads:
        for j in range(min_suf_len + 1):
            i = suf_arr.index(read[j:]+'$') - 1
            while normalised_damerau_levenshtein_distance(read,suf_arr[i][:-1]) <= 0.5:
                if not reads[suf_arr_ind[i]] == read and \
                   normalised_damerau_levenshtein_distance(read,suf_arr[i][:-1]) < max_diff and \
                   read.startswith(suf_arr[i][:-1]):
                    if (reads_map[reads[suf_arr_ind[i]]],reads_map[read]) not in B: B[(reads_map[reads[suf_arr_ind[i]]],reads_map[read])] = len(suf_arr[i][:-1])
                    else: B[(reads_map[reads[suf_arr_ind[i]]],reads_map[read])] = max(len(suf_arr[i][:-1]),B[(reads_map[reads[suf_arr_ind[i]]],reads_map[read])])
                i -= 1
    if do_time: return B, time.time() - start
    return B

def move_col(B, cols):
    for c in range(len(B.col)):
        B.col[c] = cols[B.col[c]]
            
def move_row(B,rows):
    for r in range(len(B.row)):
        B.row[r] = rows[B.row[r]]

def find_lower_diagonal_path(B,reads_map,cols,rows,do_time = False):
    import time
    import numpy as np

    if do_time: start = time.time()
    argpen = lambda l: np.argpartition(l,-2)[-2]

    new_cols = cols[:]
    if B.sum(axis=0).min() == 0: new_cols = list(c for c in new_cols if c not in [new_cols[B.sum(axis=0).argmin()]]) + [new_cols[B.sum(axis=0).argmin()]]
    if B.sum(axis=1).min() == 0: 
        if B.sum(axis=1).argmin() == B.sum(axis=0).argmin():
            new_cols = [new_cols[-1]] + list(c for c in new_cols[:-1] if c not in [cols[B.getrow(rows.index(new_cols[-1])).argmax()]]) + [cols[B.getrow(rows.index(new_cols[-1])).argmax()]]
        else: new_cols = [rows[B.sum(axis=1).argmin()]] + list(c for c in new_cols if c not in [rows[B.sum(axis=1).argmin()]])

    cols_map = dict((cols.index(c),new_cols.index(c)) for c in range(len(cols)))
    move_col(B,cols_map)
    cols = new_cols

    new_rows = cols[:]
    rows_map = dict((rows.index(r),new_rows.index(r)) for r in range(len(rows)))
    move_row(B,rows_map)
    rows = new_rows

    i,j,k = len(rows), len(cols) - 1, B.sum(axis=1).argmin() if B.sum(axis=1).min() == 0 else None

    while j > (k if B.sum(axis=1).min() == 0 else 0):
        if k is not None and B.getrow(rows.index(cols[j])).argmax() == k: 
            cols_,c_ = [], 0

            while j + c_ + 1 < len(rows):
                c_ += 1
                if len(B.getrow(j+c_).nonzero()[1]) > 1:
                    cols_ = np.argpartition(B.getrow(j+c_).toarray().flatten(),-2)[::-1][:2]
                    if cols[cols_[1]] in cols[:j] and B.getcol(cols_[1]).argmax() == j+c_: break
            
            if j + c_ + 1 == len(cols): new_cols = cols[:k+1] + cols[j:] + cols[k+1:j]
            else: new_cols = cols[:k+1] + cols[j:j+c_] + list(c for c in cols[k+1:j] if c not in [cols[min(cols_)]]) + [cols[min(cols_)]] + cols[j+c_:]
            cols_map = dict((cols.index(c),new_cols.index(c)) for c in range(len(cols)))
            move_col(B,cols_map)
            cols = new_cols

            if j + c_ + 1 == len(rows): new_rows = cols[:]
            else: new_rows = cols[:k+c_+1] + list(r for r in rows[k:j+c_] if r not in cols[:k+c_+1] + cols[j+c_:]) + cols[j+c_:]
            rows_map = dict((rows.index(r),new_rows.index(r)) for r in range(len(rows)))
            move_row(B,rows_map)
            rows = new_rows

            i,j,k = j + c_ + 1, j + c_, k + c_
        else:
            cmax = B.getrow(rows.index(cols[j])).argmax()
            if len(B.getrow(rows.index(cols[j])).nonzero()[1]) > 1:
                cpen = argpen(B.getrow(rows.index(cols[j])).toarray().flatten()) 
                if cmax > j: 
                    if len(B.getrow(cmax+1).nonzero()[1]) > 1 and \
                    B.getrow(cmax+1).getcol(argpen(B.getrow(cmax+1).toarray().flatten())).data[0] >=  B.getrow(rows.index(cols[j])).getcol(cpen).data[0]: 
                        crange = [argpen(B.getrow(cmax).toarray().flatten()),cmax]
                    else: crange = [cpen]
                else: crange = [cmax]
            else: crange = [cmax]
            while crange[0] > j:
                if len(B.getrow(crange[0]).nonzero()[1]) > 1:
                    crange = [argpen(B.getrow(crange[0]).toarray().flatten())] + crange
                else:
                    crange = [B.getrow(crange[0]).argmax()] + crange
                if crange[0] == j: crange = [B.getrow(crange[1]).argmax()] + crange[1:]

            new_cols = list(c for c in cols[:j] if c not in list(cols[cr] for cr in crange)) + list(cols[cr] for cr in crange) + list(c for c in cols[j:] if c not in list(cols[cr] for cr in crange))
            cols_map = dict((cols.index(c),new_cols.index(c)) for c in range(len(cols)))
            move_col(B,cols_map)
            cols = new_cols

            new_rows = list(r for r in rows[:i] if r not in cols[j:]) + cols[j:]
            rows_map = dict((rows.index(r),new_rows.index(r)) for r in range(len(rows)))
            move_row(B,rows_map)
            rows = new_rows
        j -= 1
        i -= 1

    seq = ''
    for s,d in zip(list(reads_map[k] for k in rows)[:-1],B.diagonal(-1)):
        seq += s[:-d]
    seq += list(reads_map[k] for k in rows)[-1]
    if do_time: return seq, time.time() - start
    return seq

if __name__ == "__main__":
    import time, os
    from scipy.sparse import coo_matrix, vstack, hstack, save_npz # pip install scipy
    from jellyfish import damerau_levenshtein_distance # pip install jellyfish
    from Bio import SeqIO # pip install Bio

    # generated sequences
    n = 50
    m = 1000
    for seed in range(n):  
        with open('data/output/generated_sequence_seed_'+str(seed)+'_sequitur.csv','a') as f:
            f.write('edit_distance,target_sequence_length,output_sequence_length,suffix_array_construction_time,adjacency_matrix_construction_time,sequence_reconstruction_time\n')
            seq = generate_genome_sequence(10000,seed=seed)
            reads,_ = generate_reads(seq,250,500,50,100,seed=seed)
            reads_map = dict(zip(list(range(len(reads))),reads))
            rows = list(range(len(reads)))
            cols = list(range(len(reads)))
            for _ in range(m):
                suf_arr,suf_arr_ind,t1 = build_suffix_array(reads,do_time=True)
                B,t2 = create_bipartite_adjacency_matrix(reads,suf_arr=suf_arr,suf_arr_ind=suf_arr_ind,do_time=True)
                start = time.time()
                B = coo_matrix((list(B.values()),list(zip(*B.keys())))).T
                if B.shape[0] < len(rows): B = vstack([B,coo_matrix((1,B.shape[1]),dtype=B.dtype)])
                if B.shape[1] < len(cols): B = hstack([B,coo_matrix((B.shape[0], 1),dtype=B.dtype)])
                t2 += time.time() - start
                seq_,t3 = find_lower_diagonal_path(B,reads_map,cols,rows,do_time=True)
                f.write('{},{},{},{},{},{}\n'.format(damerau_levenshtein_distance(seq_,seq),len(seq),len(seq_),t1,t2,t3))
    
    # euler
    for seed in range(n): 
        with open('data/output/generated_sequence_seed_'+str(seed)+'_euler.csv','a') as f:
            f.write('k,edit_distance,target_sequence_length,output_sequence_length,de_bruijn_graph_construction_time,euler_path_reconstruction_time\n')
            seq = generate_genome_sequence(10000,seed=seed)
            reads,_ = generate_reads(seq,250,500,50,100,seed=seed)
            outputs_seq = {}
            outputs_dbg_time = {}
            outputs_euler_time = {}
            a,b = find_longest_overlap(reads)
            for _ in range(m):
                for k in range(a,b):
                    G,outputs_dbg_time[k] = create_de_bruijn_graph(k,reads,do_time=True)
                    outputs_seq[k],outputs_euler_time[k] = eulerian_path(G,do_time=True)
                for k in outputs_seq:
                    f.write('{},{},{},{},{},{}\n'.format(k,damerau_levenshtein_distance(outputs_seq[k],seq),len(seq),len(outputs_seq[k]),outputs_dbg_time[k],outputs_euler_time[k]))
    
    # real genomic sequence, generated reads
    # sequitur
    n = 50
    m = 1000
    for seed in range(n):
        with open('data/output/Raphanus sativus_NC_018551.1_seed_'+str(seed)+'_sequitur.csv','a') as f:
            f.write('edit_distance,target_sequence_length,output_sequence_length,suffix_array_construction_time,adjacency_matrix_construction_time,sequence_reconstruction_time\n')
            for record in SeqIO.parse("data/input/Raphanus sativus_NC_018551.1.fasta",'fasta'): seq = record.seq
            reads,_ = generate_reads(seq,250,250,50,50,seed=seed,min_coverage=None)
            reads_map = dict(zip(list(range(len(reads))),reads))
            rows = list(range(len(reads)))
            cols = list(range(len(reads)))
            for _ in range(m):
                suf_arr,suf_arr_ind,t1 = build_suffix_array(reads,do_time=True)
                B,t2 = create_bipartite_adjacency_matrix(reads,suf_arr=suf_arr,suf_arr_ind=suf_arr_ind,do_time=True)
                start = time.time()
                B = coo_matrix((list(B.values()),list(zip(*B.keys())))).T
                if B.shape[0] < len(rows): B = vstack([B,coo_matrix((1,B.shape[1]),dtype=B.dtype)])
                if B.shape[1] < len(cols): B = hstack([B,coo_matrix((B.shape[0], 1),dtype=B.dtype)])
                t2 += time.time() - start
                if not os.path.exists('data/input/matrices/seed_'+str(seed)+'.npz'):
                    save_npz('data/input/matrices/seed_'+str(seed)+'.npz', B)
                seq_,t3 = find_lower_diagonal_path(B,reads_map,cols,rows,do_time=True)
                f.write('{},{},{},{},{},{}\n'.format(damerau_levenshtein_distance(seq_,seq),len(seq),len(seq_),t1,t2,t3))
    
    # euler
    for seed in range(n):
        with open('data/output/Raphanus sativus_NC_018551.1_seed_'+str(seed)+'_euler.csv','a') as f:
            f.write('k,edit_distance,target_sequence_length,output_sequence_length,de_bruijn_graph_construction_time,euler_path_reconstruction_time\n')
            for record in SeqIO.parse("data/input/Raphanus sativus_NC_018551.1.fasta",'fasta'): seq = record.seq
            reads,_ = generate_reads(seq,250,250,50,50,seed=seed,min_coverage=None)
            outputs_seq = {}
            outputs_dbg_time = {}
            outputs_euler_time = {}
            a,b = find_longest_overlap(reads)
            for _ in range(m):
                for k in range(a,b):
                    G,outputs_dbg_time[k] = create_de_bruijn_graph(k,reads,do_time=True)
                    outputs_seq[k],outputs_euler_time[k] = eulerian_path(G,do_time=True)
                for k in outputs_seq:
                    f.write('{},{},{},{},{},{}\n'.format(k,damerau_levenshtein_distance(outputs_seq[k],seq),len(seq),len(outputs_seq[k]),outputs_dbg_time[k],outputs_euler_time[k]))