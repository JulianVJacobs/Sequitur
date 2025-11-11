import time, os
import argparse
import pylcs # pip install pylcs
from typing import Union, Tuple
import networkx as nx # pip install networkx
from itertools import chain
from Bio import SeqIO # pip install Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from fastDamerauLevenshtein import damerauLevenshtein as damerau_levenshtein_distance # pip install fastDamerauLevenshtein

def find_longest_overlap(reads: list) -> int:
    overlaps = []
    for read in reads:
        overlaps += pylcs.lcs2_of_list(read, list(set(reads).symmetric_difference([read])))
    return min(overlaps),max(overlaps)

def create_de_bruijn_graph(
        k: int, 
        sequences: list,
        do_time: bool = False
    ) -> Union[nx.Graph,Tuple[nx.Graph,float]]:
    """
    Create a de Bruijn graph from a set of DNA sequences.
    
    Parameters:
    - k (int): k-mer size
    - sequences (list): List of DNA sequences
    - do_time (bool): should the function be time of not (default False)
    
    Returns:
    - Union[nx.Graph,Tuple[nx.Graph,float]]: De Bruijn graph or a tuple of the De Bruijn graph and the time
                                        to create it.
    """
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

def eulerian_path(
        graph: nx.DiGraph, 
         do_time: bool = False,
         file: str = ""
    ) -> Union[None, float]:
    """
    Find an Eulerian path in the given graph or an approximate Eulerian path if none exists.
    
    Parameters:
    - graph (nx.Graph): De Bruijn graph
    - do_time (bool): should the function be timed or not (default False)
    
    Returns:
    - Union[str, Tuple[str, float]]: string describing the Eulerian path or a tuple of this 
                                string and the execution time.
    """
    if do_time: start_time = time.time()
    records = []
    seq = []
    
    i=0
    for subgraph in [graph.subgraph(c).copy() for c in nx.connected_components(graph.to_undirected())]:
        for node in nx.eulerian_path(graph.subgraph(nx.eulerize(subgraph))):
            if len(seq): seq.append(node[0][-1])
            else: seq.append(node[0])
        records.append(SeqRecord(
            Seq(''.join(seq)),
            id=f"contig_{i}",
            description=""
            ))
        i+=1
        seq = []
        
    if do_time: 
        t = time.time() - start_time
        SeqIO.write(records, f"/workspace/sequitur/data/output/Raphanus sativus_NC_018551.1/{file}.dbg.out.fasta", "fasta")
        return t

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str)
    args = parser.parse_args()

    for record in SeqIO.parse("/workspace/sequitur/data/input/Raphanus sativus_NC_018551.1/Raphanus sativus_NC_018551.1.fasta",'fasta'): seq = record.seq
    if not os.path.exists(f"/workspace/sequitur/data/output/{args.input}.euler.csv") or os.path.getsize(f"/workspace/sequitur/data/output/{args.input}.euler.csv") == 0:
        with open(f"/workspace/sequitur/data/output/{args.input}.euler.csv",'a') as f:
            f.write('edit_distance,target_sequence_length,output_sequence_length,suffix_array_construction_time,adjacency_matrix_construction_time,sequence_reconstruction_time\n')
    with open(f"/workspace/sequitur/data/output/{args.input}.euler.csv",'a') as f:
        f.write('k,de_bruijn_graph_construction_time,euler_path_reconstruction_time\n')
        reads= chain(
            (str(read.seq.upper()) for read in SeqIO.parse(f"/workspace/sequitur/data/input/Raphanus sativus_NC_018551.1/{args.input}.1.fastq",'fastq')),
            (str(read.seq.upper().reverse_complement()) for read in SeqIO.parse(f"/workspace/sequitur/data/input/Raphanus sativus_NC_018551.1/{args.input}.2.fastq",'fastq'))
            )
        outputs_seq = {}
        outputs_dbg_time = {}
        outputs_euler_time = {}
        a,b = find_longest_overlap(list(reads))
        for _ in range(100):
            for k in range(a,b):
                G,outputs_dbg_time[k] = create_de_bruijn_graph(k,list(str(read) for read in reads),do_time=True)
                outputs_euler_time[k] = eulerian_path(G,do_time=True,file=args.input)
            for k in outputs_seq:
                f.write('{},{},{}\n'.format(k,outputs_dbg_time[k],outputs_euler_time[k]))