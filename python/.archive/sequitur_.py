import time, os, argparse
from memory_profiler import memory_usage
from typing import Union, Tuple, Generator
from scipy.sparse import coo_matrix, csr_matrix, vstack, hstack, save_npz, load_npz # pip install scipy
from scipy.optimize import linear_sum_assignment
import numpy as np # pip install numpy
from itertools import chain
from Bio import SeqIO # pip install Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from math import ceil
import pickle
from hungarian_algorithm import algorithm
import copy
from collections import namedtuple
from itertools import product
from fastDamerauLevenshtein import damerauLevenshtein as damerau_levenshtein_distance # pip install fastDamerauLevenshtein

def normalised_damerau_levenshtein_distance(
        suf: str,
        pre: str
    ) -> Tuple[float,int]:
    """
    Find the Damerau-Levenshtein edit distance of two strings normalised to the length
    of the shorter string. This normalisation is because we want to path prefixes to
    suffixes and this means that in general we will be comparing a full string to a
    portion of another string.
    
    Parameters:
    - suf (str): suffix string
    - pre (str): prefix string
    
    Returns:
    - float: the normalised Demarau-Levenshtein edit distance of the input strings as a ratio
    - int  : the normalised Demarau-Levenshtein edit distance of the input strings as an integer
    """
    m = min(len(suf),len(pre))
    d = damerau_levenshtein_distance(str(suf)[-m:],str(pre)[:m],similarity=False)
    return (d/m,m-d)

def build_suffix_array(
        reads: list, 
        min_suf_len: int = 3,
        test: bool = False
    ) -> Tuple[list, list]:
    if test: start = time.time()
    suf_arr = []
    reads = list(reads)
    for index,read in enumerate(reads):
        read += '$' + str(index)
        for i in range(read.index('$')-min_suf_len+1):
            suf_arr += [read[i:]]
        read = read.replace('$','^')
        for i in range(read.index('^')-min_suf_len+1):
            suf_arr += [read[:read.index('^')-i]+read[read.index('^'):]]
    suf_arr.sort()
    if test: return suf_arr, time.time() - start
    return suf_arr

def create_bipartite_adjacency_matrix(
        reads: Generator, 
        suf_arr: list | None = None, 
        test: bool = False,
        max_diff: float = 0.25, 
        min_suf_len: int = 3
    ) -> dict:
    
    if test: start = time.time()
    if suf_arr is None: suf_arr = build_suffix_array(reads,min_suf_len=min_suf_len)
    reads = list(reads)
    B, O = {}, {}
    for suffix_i,suffix in enumerate(reads):
        _u, _d = False, False
        if suffix_i not in B: B[suffix_i], O[suffix_i] = {}, {}
        for j in range(min_suf_len, len(suffix) - 1):
            # get index of the value after the suffix
            i = suf_arr.index(f"{suffix[-j:]}${str(suffix_i)}") + 1
            while '$' in suf_arr[i] or suf_arr[i].endswith(str(suffix_i)): i += 1
            if i < len(suf_arr):
                prefix, prefix_i = suf_arr[i].split('^')
                prefix_i = int(str(prefix_i))
                edf, edi = normalised_damerau_levenshtein_distance(suffix[-j:], prefix)
                _d = edf > max_diff
                while edf <= max_diff and i < len(suf_arr):
                    if prefix_i != suffix_i and edf < max_diff:
                        if prefix_i not in B[suffix_i]: 
                            B[suffix_i][prefix_i] = edi
                            O[suffix_i][prefix_i] = min(len(suffix[-j:]), len(prefix))
                        else: 
                            if edi == B[suffix_i][prefix_i] and min(len(suffix[-j:]), len(prefix)) > O[suffix_i][prefix_i]:
                                O[suffix_i][prefix_i] = min(len(suffix[-j:]), len(prefix))
                            elif edi > B[suffix_i][prefix_i]:
                                B[suffix_i][prefix_i] = edi
                                O[suffix_i][prefix_i]  = min(len(suffix[-j:]), len(prefix))
                    i += 1
                    while i < len(suf_arr) and ('$' in suf_arr[i] or suf_arr[i].endswith(str(suffix_i))): i += 1
                    if i >= len(suf_arr): break
                    prefix, prefix_i = suf_arr[i].split('^')
                    prefix_i = int(str(prefix_i))
                    edf, edi = normalised_damerau_levenshtein_distance(suffix[-j:], prefix)
            k = suf_arr.index(f"{suffix[-j:]}${str(suffix_i)}") - 1
            while '$' in suf_arr[k] or suf_arr[k].endswith(str(suffix_i)): k -= 1
            if k >= 0:
                prefix, prefix_k = suf_arr[k].split('^')
                prefix_k = int(str(prefix_k))
                edf, edi = normalised_damerau_levenshtein_distance(suffix[-j:], prefix)
                _u = edf > max_diff
                while edf <= max_diff and k >= 0:
                    if prefix_k != suffix_i and edf < max_diff:
                        if prefix_k not in B[suffix_i]: 
                            B[suffix_i][prefix_k] = edi
                            O[suffix_i][prefix_k] = min(len(suffix[-j:]), len(prefix))
                        else: 
                            if edi == B[suffix_i][prefix_k] and min(len(suffix[-j:]), len(prefix)) > O[suffix_i][prefix_k]:
                                O[suffix_i][prefix_k] = min(len(suffix[-j:]), len(prefix))
                            elif edi > B[suffix_i][prefix_k]:
                                B[suffix_i][prefix_k] = edi
                                O[suffix_i][prefix_k] = min(len(suffix[-j:]), len(prefix))
                    k -= 1
                    while k >= 0 and ('$' in suf_arr[k] or suf_arr[k].endswith(str(suffix_i))): k -= 1
                    if k < 0: break
                    prefix, prefix_k = suf_arr[k].split('^')
                    prefix_k = int(str(prefix_k))
                    edf, edi = normalised_damerau_levenshtein_distance(suffix[-j:], prefix)
            if _u and _d: break 
    if test: return B, O, time.time() - start
    return B, O

def find_lower_diagonal_path(
        B: csr_matrix,
        O: dict,
        reads_map: dict,
        cols: list,
        rows: list,
        file: str,
        test: bool = False,
        quality_map: dict | None = None,
        target_seq: Seq = None
    ) -> Union[str, Tuple[str, float]]:
    if test: start = time.time()
    
    B_ = {i: {j: B.todok()[i, j] for j in range(B.shape[1]) if (i, j) in B.todok()} for i in range(B.shape[0])}
    
    l = algorithm.find_matching(B_, matching_type = 'max', return_type = 'list')
    
    # ~ 1. Subtract the smallest entry in each row from all the other entries in the row. 
    for row_index in range(B.shape[0]):
        B[row_index, :] -= B[row_index, :].min()
        
    # ~ 2. Subtract the smallest entry in each column from all the other entries in the 
    for col_index in range(B.shape[1]):
        B[:, col_index] -= B[:, col_index].min()
    
    # ~ 3. Draw lines through the row and columns that have the 0 entries such that the 
    # ~     fewest lines possible are drawn.
    # ~ 3.1. For each column, find the row indices for every 0.
    marked_cols, marked_rows = [], []
    for col_index in range(B.shape[1]):
        row_indices = np.where(B[:, col_index].toarray().flatten() == 0)[0]
        for row_index in row_indices:
            if (
                sum(B[row_index, :].toarray().flatten() == 0) > sum(B[:, col_index].toarray().flatten() == 0) and
                row_index not in marked_rows
            ):
                marked_rows += [row_index]
            else:
                marked_cols += [col_index]
    marked_cols = list(set(marked_cols))
    marked_rows = list(set(marked_rows))

    # ~ 4. If there are nn lines drawn, an optimal assignment of zeros is possible and the 
    # ~     algorithm is finished. If the number of lines is less than nn, then the optimal 
    # ~     number of zeroes is not yet reached. Go to the next step.
    
    # ~ 5. Find the smallest entry not covered by any line. Subtract this entry from each 
    # ~     row that isn’t crossed out, and then add it to each column that is crossed out. 
    # ~     Then, go back to Step 3.
    
    # row_ind, col_ind = min_weight_full_bipartite_matching(B,True)

    seq = ''
    reads_map = dict(reads_map)
    quality_map = dict(quality_map)
    b = True
    for c,r in zip(row_ind, col_ind):
        if b: 
            seq += reads_map[c][:-O[r][c]]
            b = False
        prefix,pre_q = reads_map[c][-O[r][c]:],quality_map[c][-O[r][c]:]
        suffix,suf_q = reads_map[r][:-O[r][c]],quality_map[r][:-O[r][c]]
        for base in range(O[r][c]):
            if pre_q[base] >= suf_q[base]: seq += prefix[base]
            else: seq += suffix[base]
    seq += reads_map[rows[-1]][-O[rows[-1]][cols[-1]]:]
    
    
    if test: 
        t = time.time() - start
        SeqIO.write([
            SeqRecord(
                Seq(seq),
                id="sequitur",
                description=""
                )], f"/workspaces/Sequitur/data/output/Raphanus sativus_NC_018551.1/{file}.sequitur.out.fasta", "fasta")
        return seq, t
    return seq

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str)
    args = parser.parse_args()
    test_env = False
 
    for record in SeqIO.parse("/workspaces/Sequitur/data/input/Raphanus sativus_NC_018551.1/Raphanus sativus_NC_018551.1.fasta",'fasta'): seq = record.seq
    def get_reads():
        return chain(
            (read.seq.upper().reverse_complement() for read in SeqIO.parse(f"/workspaces/Sequitur/data/input/Raphanus sativus_NC_018551.1/{args.input}.1.fastq",'fastq')),
            (read.seq.upper() for read in SeqIO.parse(f"/workspaces/Sequitur/data/input/Raphanus sativus_NC_018551.1/{args.input}.2.fastq",'fastq'))
            )
    def get_qualities():
        return chain(
            (read.letter_annotations["phred_quality"] for read in SeqIO.parse(f"/workspaces/Sequitur/data/input/Raphanus sativus_NC_018551.1/{args.input}.1.fastq",'fastq')),
            (read.letter_annotations["phred_quality"] for read in SeqIO.parse(f"/workspaces/Sequitur/data/input/Raphanus sativus_NC_018551.1/{args.input}.2.fastq",'fastq'))
            )
    def dict_to_mat(B,test=False):
        if test: start = time.time()
        rows = []
        cols = []
        data = []

        for row, col_dict in B.items():
            for col, value in col_dict.items():
                rows.append(row)
                cols.append(col)
                data.append(value)
        B = coo_matrix((np.array(data), (np.array(rows), np.array(cols)))).T
        if B.shape[0] < len(list(index for index, _ in enumerate(get_reads()))): B = vstack([B,coo_matrix((1, B.shape[1]),dtype=B.dtype)])
        if B.shape[1] < len(list(index for index, _ in enumerate(get_reads()))): B = hstack([B,coo_matrix((B.shape[0], 1),dtype=B.dtype)])
        if test: return B.tocsr(), time.time() - start
        return B.tocsr()
    if not os.path.exists(f"/workspaces/Sequitur/data/output/{args.input}.sequitur.csv") or os.path.getsize(f"/workspaces/Sequitur/data/output/{args.input}.sequitur.csv") == 0:
        with open(f"/workspaces/Sequitur/data/output/{args.input}.sequitur.csv",'a') as f:
            f.write('edit_distance,target_sequence_length,output_sequence_length,suffix_array_construction_time,suffix_array_construction_mem_peak,suffix_array_construction_mem_avg,adjacency_matrix_construction_time,adjacency_matrix_construction_mem_peak,adjacency_matrix_construction_mem_avg,sequence_reconstruction_time,sequence_reconstruction_mem_peak,sequence_reconstruction_mem_avg\n')
    with open(f"/workspaces/Sequitur/data/output/{args.input}.sequitur.csv",'a') as f:
        for _ in range(100):
            if not test_env or not os.path.exists(f"/workspaces/Sequitur/data/output/test/{args.input}.B.sequitur.npz"):
                suf_arr,t1 = build_suffix_array(get_reads(),test=True)
                B, O,t2 = create_bipartite_adjacency_matrix(get_reads(),suf_arr=suf_arr,test=True)
                l = algorithm.find_matching(B, matching_type = 'max', return_type = 'list')
                raise # ! set test_env = True
                with open(f"/workspaces/Sequitur/data/output/test/{args.input}.O.sequitur.pkl", 'wb') as f:
                    pickle.dump(O, f)
                B, t2_ = dict_to_mat(B,test=True)
                t2 += t2_
                save_npz(f"/workspaces/Sequitur/data/output/test/{args.input}.B.sequitur.npz",B)
            else:
                B = load_npz(f"/workspaces/Sequitur/data/output/test/{args.input}.B.sequitur.npz")
                with open(f"/workspaces/Sequitur/data/output/test/{args.input}.O.sequitur.pkl", 'rb') as f:
                    O = pickle.load(f)
            seq_,t3 = find_lower_diagonal_path(
                    B,O,
                    zip((index for index,_ in enumerate(get_reads())),get_reads()),
                    (index for index,_ in enumerate(get_reads())),
                    (index for index,_ in enumerate(get_reads())),
                    file=args.input,test=True,
                    quality_map=zip((index for index,_ in enumerate(get_reads())),get_qualities()),
                    target_seq=seq
                )
            f.write(f"{damerau_levenshtein_distance(seq_,seq)},{len(seq)},{len(seq_)},{t1},0,0,{t2},0,0,{t3},0,0\n")
            # m1, (suf_arr,suf_arr_ind,t1) = memory_usage(build_suffix_array(get_reads(),test=True),retval=True)
            # m2, (B, O,t2) = memory_usage(create_bipartite_adjacency_matrix(get_reads(),zip(get_reads(),(index for index,_ in enumerate(get_reads()))),suf_arr=suf_arr,suf_arr_ind=suf_arr_ind,test=True),retval=True)
            # m2_, (B, t2_) = memory_usage(dict_to_mat(B,test=True),retval=True)
            # m2 += m2_
            # t2 += t2_
            # m3,(seq_,t3) = memory_usage(find_lower_diagonal_path(B,O,zip(get_ids(),get_reads()),get_ids(),get_ids(),file=args.input,test=True,quality_map=zip(get_ids(),get_reads())),retval=True)
            # f.write(f"{damerau_levenshtein_distance(seq_,seq)},{len(seq)},{len(seq_)},{t1},{max(m1)},{sum(m1)/len(m1)},{t2},{max(m2)},{sum(m2)/len(m2)},{t3},{max(m3)},{sum(m3)/len(m3)}\n")