import time, os
import pylcs # pip install pylcs
from typing import Union, Tuple
import networkx as nx # pip install networkx
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
        for node in nx.eulerian_path(nx.eulerize(subgraph)):
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
        SeqIO.write(records, f"/workspace/sequitur/data/output/natural_language_{file}.out.fasta", "fasta")
        return t

if __name__ == "__main__":
    nat_lang_seq = [
		('betty_bought_butter_the_butter_was_bitter_betty_bought_better_butter_to_make_the_bitter_butter_better',
		['betty_bought_butter_th',
							'tter_the_butter_was_',
								'he_butter_was_bitter_',
										'as_bitter_betty_bought',
													'tty_bought_better_butter_t',
														'ught_better_butter_to',
																	'r_butter_to_make_the_',
																				'ke_the_bitter_butter_better']),
		('you say hello world, i bellow go to hell',
		['you say hel',
					' say hello wo',
							'lo world, i be',
								'ld, i bellow go t',
											'ow go to hell']),
		('she_sells_sea_shells_on_the_sea_shore',
		['she_sells_s',
					'lls_sea_shel',
							'ea_shells_o',
							'shells_on_the_s',
										'he_sea_s',
											'ea_shore'])
    ]
    for i, (seq, reads) in enumerate(nat_lang_seq):
        if not os.path.exists('/workspace/sequitur/data/output/natural_language_sequences.euler.csv') or os.path.getsize('data/output/natural_language_sequences.euler.csv') == 0:
            with open('/workspace/sequitur/data/output/natural_language_sequences.euler.csv', 'a') as f:
                f.write('natural language_sequence,k,de_bruijn_graph_construction_time,euler_path_reconstruction_time\n')
            with open('/workspace/sequitur/data/output/natural_language_sequences.euler.csv','a') as f:
                for _ in range(100):
                    outputs_seq = {}
                    outputs_dbg_time = {}
                    outputs_euler_time = {}
                    a,b = find_longest_overlap(reads)
                    for k in range(a,b):
                        G,outputs_dbg_time[k] = create_de_bruijn_graph(k,reads,do_time=True)
                        outputs_euler_time[k] = eulerian_path(G,do_time=True,file=i)
                    for k in outputs_seq:
                        f.write('{},{},{},{}\n'.format(i,k,outputs_dbg_time[k],outputs_euler_time[k]))