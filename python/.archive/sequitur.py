import time, os, argparse, pickle, concurrent.futures, copy
from memory_profiler import memory_usage
from typing import Union, Tuple, Generator
from scipy.sparse import coo_matrix, csr_matrix, vstack, hstack, save_npz, load_npz # pip install scipy
from scipy.sparse.csgraph import min_weight_full_bipartite_matching, maximum_bipartite_matching
import numpy as np # pip install numpy
from itertools import chain
from Bio import SeqIO # pip install Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from math import ceil
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
    # for suffix_i,suffix in enumerate(reads):
    def process_suffix(suffix_i, suffix):
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
            
    with concurrent.futures.ThreadPoolExecutor(max_workers=16) as executor:
        futures = [executor.submit(process_suffix, suffix_i, suffix) for suffix_i, suffix in enumerate(reads)]
        concurrent.futures.wait(futures) 
        
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
    
    def move_col(
        B: coo_matrix, 
        cols: dict
        ) -> csr_matrix:
        for c in range(len(B.col)):
            B.col[c] = cols[B.col[c]]
        return B.tocsr()
                
    def move_row(
            B: coo_matrix,
            rows: dict
        ) -> csr_matrix:
        for r in range(len(B.row)):
            B.row[r] = rows[B.row[r]]
        return B.tocsr()    
    
    class ArgIterator:
        def __init__(self,B,row):
            self.B = B
            self.row_index = row
            self.row = self.B.getrow(row)
            self.argsorted_row = list(np.lexsort((np.arange(self.row.shape[1])[::-1], self.row.toarray().flatten())))[-self.row.count_nonzero() - 2:]
            # self.argsorted_row = list(np.argsort(self.row.toarray().flatten(), stable=True))
            # self.argsorted_row = list(np.argsort(self.row.toarray().flatten()))
            self.memory = []
            self.next_best = None
            self.curr_best = None
            self.generator = self.get_next_best()
            self.rotareneg = self.get_prev_best()
            next(self)
            next(self)
        
        def __iter__(self):
            return self.generator
        
        def __reti__(self):
            return self.rotareneg
        
        def __str__(self):
            return f"ArgIterator(row_index={self.row_index}, curr_best={self.curr_best}, next_best={self.next_best})"
        
        def __repr__(self):
            return f"ArgIterator(row_index={self.row_index}, curr_best={self.curr_best}, next_best={self.next_best})"
        
        def __next__(self):
            try:
                return next(self.generator)
            except: 
                raise StopIteration
        
        def __prev__(self):
            try:
                return next(self.rotareneg)
            except: 
                raise StopIteration
            
        def __eq__(self,other):
            return (
                type(self) == type(other) and 
                self.curr_best == other.curr_best and 
                self.next_best == other.next_best and 
                self.row_index == other.row_index and 
                self.memory == other.memory
            )
            
        def __hash__(self):
            return hash((self.curr_best, self.next_best, self.row_index, tuple(self.memory)))
        
        def __deepcopy__(self, memo):
            new_copy = type(self)(copy.deepcopy(self.B, memo), copy.deepcopy(self.row_index, memo))
            new_copy.row = copy.deepcopy(self.row, memo)
            new_copy.argsorted_row = copy.deepcopy(self.argsorted_row, memo)
            new_copy.memory = copy.deepcopy(self.memory, memo)
            new_copy.curr_best = copy.deepcopy(self.curr_best, memo)
            new_copy.next_best = copy.deepcopy(self.next_best, memo)
            new_copy.generator = self.get_next_best(start_from = new_copy.argsorted_row)
            new_copy.rotareneg = self.get_prev_best(start_from = new_copy.memory)
        
        def get_next_best(self, start_from: list = None):
            if start_from is not None: self.argsorted_row = start_from
            while len(self.argsorted_row):
                arg = self.argsorted_row.pop()
                if self.row[:,arg] > 0: 
                    self.memory += [arg]
                    self.rotareneg = self.get_prev_best(start_from=self.memory)
                    self.curr_best = self.next_best
                    self.next_best = arg
                    yield 
                elif self.next_best is not None:
                    self.curr_best = self.next_best
                    self.next_best = None
                    yield
                elif self.curr_best is not None:
                    self.curr_best = None
                    yield
                else: 
                    break
                    
        def get_prev_best(self, start_from: list = None):
            if start_from is not None: self.memory = start_from
            while len(self.memory):
                arg = self.memory.pop()
                if self.curr_best is None:
                    self.curr_best = arg
                    yield
                elif self.next_best is None:
                    self.next_best = self.curr_best
                    self.curr_best = arg
                    yield
                else:
                    self.argsorted_row += [self.next_best]
                    self.generator = self.get_next_best(start_from=self.argsorted_row)
                    self.next_best = self.curr_best
                    self.curr_best = arg
                    yield
        
    def prev(o: ArgIterator):
        return o.__prev__()
    
    def reset(o: ArgIterator):
        o.__init__(o.B, o.row_index)
            
    # rows = list(rows)
    # cols = list(cols)
    
    # new_cols = cols[:]
    # cols_map = dict((c,new_cols.index(c)) for c in cols)
    # B = move_col(B.tocoo(),cols_map) 
    # cols = new_cols

    # new_rows = rows[:]
    # rows_map = dict((r,new_rows.index(r)) for r in rows)
    # B = move_row(B.tocoo(),rows_map)
    # rows = new_rows
    
    # last_col_cand = list(np.argsort(np.array(B.sum(axis=0)).flatten()))
    # first_row_cand = list(np.argsort(np.array(B.sum(axis=1)).flatten()))
    # first_row = first_row_cand.pop(0)
    # last_col = last_col_cand.pop(0)
    # while len(last_col_cand):
    # B_ = hstack([B[:, :last_col], B[:, last_col + 1:]])
    reads_map = dict(reads_map)
    matching = maximum_bipartite_matching(B)
    B_ = vstack([B[first_row, :], hstack([B[:, :last_col], B[:, last_col + 1:], B[:, last_col]])[:first_row, :], hstack([B[:, :last_col], B[:, last_col + 1:], B[:, last_col]])[first_row + 1:, :]])
    for i in range(0,B.shape[0],10):
        continue
        # if not any(r + 1 == c for r,c in matching): break
    b = vstack([B[first_row, :], hstack([B[:, :last_col], B[:, last_col + 1:], B[:, last_col]])[:first_row, :], hstack([B[:, :last_col], B[:, last_col + 1:], B[:, last_col]])[first_row + 1:, :]])
    b
    b
    B_ = dict((c, ArgIterator(vstack([B[first_row, :], hstack([B[:, :last_col], B[:, last_col + 1:], B[:, last_col]])[:first_row, :], hstack([B[:, :last_col], B[:, last_col + 1:], B[:, last_col]])[first_row + 1:, :]]),r + 1)) for r,c in matching)
    B_[len(B_)] = ArgIterator(vstack([B[first_row, :], hstack([B[:, :last_col], B[:, last_col + 1:], B[:, last_col]])[:first_row, :], hstack([B[:, :last_col], B[:, last_col + 1:], B[:, last_col]])[first_row + 1:, :]]),0)
    for c in range(len(B_)):
        if B_[c].row_index == 0:
            while B_[c].curr_best is not None:
                next(B_[c])
            B_ = dict(zip(range(len(B_)),[B_[c]] + list(v for v in B_.values() if v != B_[c])))
            continue
        while B_[c].curr_best != c:
            try:
                next(B_[c])
            except:
                break
            
    row_indices = list(v.row_index for v in B_.values())
    B_ = dict(list(zip(range(len(B_)),list(v for v in B_.values() if v != B_[row_indices.index(last_col)]) + [B_[row_indices.index(last_col)]])))
    
    i,j = len(rows), len(cols) - 1
    p = j
    k = []
    b_copy = None
    start_from = 0
    
    while j > -1:
        row_indices = list(v.row_index for v in B_.values())
        if B_[j].curr_best == B_[j].row_index:
            next(B_[j])
            continue
        if B_[j].curr_best in k:
            k_ = k.index(B_[j].curr_best) + 1
            if k_ < len(k):
                b_copy = copy.deepcopy(B_[k_])
                reset(B_[k_])
            else:
                reset(B_[k_ - 1])
                k = k[:-1]
        j_ = row_indices.index(B_[j].curr_best) if B_[j].curr_best is not None else 0
        if (
                B_[j].curr_best in row_indices[j:] or b_copy is not None or
                (j_ == 0 and B_[j_].curr_best is None) or j_ is None
        ):
            j_ += 1
            if j_ < B.shape[0]:
                if (
                        B_[j_].next_best is not None and B_[j].next_best is not None and
                        B[B_[j_].row_index,B_[j_].next_best] > B[B_[j].row_index,B_[j].next_best]
                ):
                    if b_copy is not None:
                        for ind in range(k_ - 1, len(k)):
                            reset(B_[ind])
                        k = k[:k_ - 1]
                        b_copy = None
                        j_ -= 1
                    else:
                        try:
                            next(B_[j_])
                        except:
                            B_[j_]
                        j = j_
                        continue
                else:
                    try:
                        if b_copy is not None:
                            B_[k_] = b_copy
                            b_copy = None
                        next(B_[j])
                        while B_[j].curr_best is None:
                            if j == 0:
                                B_[0]
                            if j == len(B_) - 1:
                                while True:
                                    if not len(last_col_cand):
                                        last_col_cand = list(np.argsort(np.array(B.sum(axis=0)).flatten()))
                                        start_from += 1
                                        for _ in range(start_from):
                                            last_col = last_col_cand.pop(0)
                                        first_row = first_row_cand.pop(0)
                                    else: last_col = last_col_cand.pop(0)
                                    B_ = hstack([B[:, :last_col], B[:, last_col + 1:]])
                                    B_ = vstack([B_[:first_row, :], B_[first_row + 1:, :]])
                                    matching = sorted(list(zip(*min_weight_full_bipartite_matching(B_,True))), key= lambda x: x[1])
                                    if not any(r + 1 == c for r,c in matching): break
                                B_ = dict((c, ArgIterator(vstack([B[first_row, :], hstack([B[:, :last_col], B[:, last_col + 1:], B[:, last_col]])[:first_row, :], hstack([B[:, :last_col], B[:, last_col + 1:], B[:, last_col]])[first_row + 1:, :]]),r + 1)) for r,c in matching)
                                B_[len(B_)] = ArgIterator(vstack([B[first_row, :], hstack([B[:, :last_col], B[:, last_col + 1:], B[:, last_col]])[:first_row, :], hstack([B[:, :last_col], B[:, last_col + 1:], B[:, last_col]])[first_row + 1:, :]]),0)
                                for c in range(len(B_)):
                                    if B_[c].row_index == 0:
                                        while B_[c].curr_best is not None:
                                            next(B_[c])
                                        B_ = dict(zip(range(len(B_)),[B_[c]] + list(v for v in B_.values() if v != B_[c])))
                                        continue
                                    while B_[c].curr_best != c:
                                        try:
                                            next(B_[c])
                                        except:
                                            break
                                row_indices = list(v.row_index for v in B_.values())
                                B_ = dict(list(zip(range(len(B_)),list(v for v in B_.values() if v != B_[row_indices.index(last_col)]) + [B_[row_indices.index(last_col)]])))
                                break
                            reset(B_[j])
                            j += 1
                            next(B_[j])
                        continue
                    except:
                        if B_[0].curr_best is None:
                            for ind in range(len(k)):
                                reset(B_[ind])
                            B_ = dict(list(zip(range(len(B_) - j),list(B_.values())[j:])) + list(zip(range(len(B_) - j,len(B_)),list(B_.values())[:j])))
                            k = list(c.row_index for c in list(B_.values())[:len(B_) - j])
                            j = len(cols) - 1
                            continue
                        else:
                            j_ -= 1
            else: j_ -= 1
        B_ = dict(list(zip(range(len(B_)),list(v for v in list(B_.values())[:j] if v not in [B_[j_], B_[j]]) + [B_[j_], B_[j]] + list(v for v in list(B_.values())[j:] if v not in [B_[j_], B_[j]]))))
        j -= 1
        if j < p: 
            p = j

    i,j = len(rows), len(cols) - 1
    while j > - 1:
        endpoints = np.argsort(B.sum(axis=0).A1,stable=True)
        crange = []   
        crange_mem = []
        cand_col = ArgIterator(B,rows.index(cols[j]))
        p = 0
        while cand_col.curr_best is None or cand_col.curr_best > j:
            p += 1
            b = False
            crange_ = list(c.curr_best for c in crange)
            if cand_col.curr_best in crange_: 
                ind = crange_.index(cand_col.curr_best) + 1
            while cand_col.curr_best in crange_:
                if (
                        ind < len(crange) and 
                        # crange[ind].next_best is not None and cand_col.next_best is not None and
                        B[crange[ind].row_index,crange[ind].next_best] > B[cand_col.row_index,cand_col.next_best]
                ):
                    if crange[ind].next_best in crange_[ind:]:
                        ind = crange_.index(crange[ind].next_best) + 1
                        continue
                    if e > 1 and endpoints[e - 1] in crange_[:ind]: 
                        e -= 1
                    crange = [crange[ind]] + crange[ind + 1:]
                    next(crange[0])
                    cand_col = ArgIterator(B,rows.index(cols[crange[0].curr_best]))
                    crange_ = list(c.curr_best for c in crange)
                    if cand_col.curr_best in crange_: 
                        ind = crange_.index(cand_col.curr_best) + 1
                    if hash(tuple(crange)) not in crange_mem: crange_mem += [hash(tuple(crange))]
                    else:
                        crange_mem
                    continue
                elif (
                        # cand_col.next_best is not None and crange[-1].next_best is not None and
                        B[crange[-1].row_index,crange[-1].next_best] > B[cand_col.row_index,cand_col.next_best]
                ):
                    try:
                        next(crange[-1])
                        if e > 1 and endpoints[e - 1] in crange[:-1]: 
                            e -= 1
                        crange = [crange[-1]]
                        cand_col = ArgIterator(B,rows.index(cols[crange[0].curr_best]))
                        if hash(tuple(crange)) not in crange_mem: crange_mem += [hash(tuple(crange))]
                        else:
                            crange_mem
                        break
                    except:
                        prev(crange[-1])
                try:
                    next(cand_col) 
                except:
                    cand_col = crange[0]
                    next(cand_col)
                    if e > 1 and endpoints[e - 1] == crange[0]: 
                        e -= 1
                    crange = crange[1:]
                    if hash(tuple(crange)) not in crange_mem: crange_mem += [hash(tuple(crange))]
                    else:
                        crange_mem
                crange_ = list(c.curr_best for c in crange)
                if cand_col.curr_best in crange_: 
                    ind = crange_.index(cand_col.curr_best) + 1
                    
            if cand_col.curr_best + 1 < B.shape[0]: 
                cand_col_succ = ArgIterator(B,rows.index(cols[cand_col.curr_best + 1]))
                if rows.index(cols[cand_col_succ.curr_best]) != rows.index(cols[cand_col_succ.curr_best]):
                    b = True
            elif cand_col.curr_best == endpoints[e]:
                e += 1
                cand_col_succ = ArgIterator(B,rows.index(cols[endpoints[e] + 1]))
            else:
                b = True
                p = 0
                
            if (
                    (
                        b and cand_col.curr_best is not None
                    ) or (
                        cand_col_succ.next_best is not None and cand_col.next_best is not None and
                        B[cand_col_succ.row_index,cand_col_succ.next_best] > B[cand_col.row_index,cand_col.next_best]
                    ) 
            ):
                    crange = [cand_col] + crange
                    if hash(tuple(crange)) not in crange_mem: crange_mem += [hash(tuple(crange))]
                    else:
                        crange_mem
                    cand_col = ArgIterator(B,rows.index(cols[crange[0].curr_best]))
            else: 
                if e > 1 and cand_col.curr_best == endpoints[e - 1]: 
                    e -= 1
                try:
                    next(cand_col)
                    if cand_col.curr_best is None: raise
                except:
                    cand_col = crange[0]
                    next(cand_col)
                    crange = crange[1:]
                    if hash(tuple(crange)) not in crange_mem: crange_mem += [hash(tuple(crange))]
                    else:
                        crange_mem
        if cand_col.curr_best < j and cand_col.curr_best not in list(c.curr_best for c in crange): 
            crange = [cand_col] + crange
        crange = list(c.curr_best for c in crange)
        
        new_cols = list(c for c in cols[:j] if c not in list(cols[cr] for cr in crange)) + list(cols[cr] for cr in crange) + list(c for c in cols[j:] if c not in list(cols[cr] for cr in crange))
        cols_map = dict((cols.index(c),new_cols.index(c)) for c in range(len(cols)))
        B = move_col(B.tocoo(),cols_map)
        cols = new_cols

        new_rows = list(r for r in rows if r not in cols[j:]) + cols[j:]
        rows_map = dict((rows.index(r),new_rows.index(r)) for r in range(len(rows)))
        B = move_row(B.tocoo(),rows_map)
        rows = new_rows
        j -= 1
        i -= 1

    seq = ''
    reads_map = dict(reads_map)
    if quality_map is None:
        for r,c in zip(rows[:-1],cols[1]):
            seq += reads_map[c][:-O[r][c]]
        seq += reads_map[cols[-1]]
        # for s,d in zip(list(reads_map[k] for k in rows)[:-1],B.diagonal(-1)):
        #     seq += s[:-d]
        # seq += list(reads_map[k] for k in rows)[-1]
    else:
        quality_map = dict(quality_map)
        b = True
        for r,c in zip(rows[:-1],cols[1]):
        # for c,r in zip(cols,rows):
            # if c == r: continue
            if b: 
                seq += reads_map[c][:-O[r][c]]
                b = False
            prefix,pre_q = reads_map[c][-O[r][c]:],quality_map[c][-O[r][c]:]
            suffix,suf_q = reads_map[r][:-O[r][c]],quality_map[r][:-O[r][c]]
            for base in range(O[r][c]):
                if pre_q[base] >= suf_q[base]: seq += prefix[base]
                else: seq += suffix[base]
        seq += reads_map[rows[-1]][-O[rows[-1]][cols[-1]]:]
        # for i,d in enumerate(B.diagonal(-1)):
        #     if i == 0: seq += reads_map[rows[i]][:-d]
        #     r1,q1 = reads_map[rows[i]][-d:],quality_map[rows[i]][-d:]
        #     r2,q2 = reads_map[rows[i+1]][:-d],quality_map[rows[i+1]][:-d]
        #     for c in range(d):
        #         if q1[c] >= q2[c]: seq+=r1[c]
        #         else: seq+=r2[c]
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
    test_env = True
 
    for record in SeqIO.parse("/workspaces/Sequitur/data/input/Raphanus sativus_NC_018551.1/Raphanus sativus_NC_018551.1.fasta",'fasta'): seq = record.seq
    def get_reads():
        return chain(
            (read.seq.upper() if read.seq.upper() in seq else read.seq.upper().reverse_complement() for read in SeqIO.parse(f"/workspaces/Sequitur/data/input/Raphanus sativus_NC_018551.1/{args.input}.1.fastq",'fastq')),
            (read.seq.upper() if read.seq.upper() in seq else read.seq.upper().reverse_complement() for read in SeqIO.parse(f"/workspaces/Sequitur/data/input/Raphanus sativus_NC_018551.1/{args.input}.2.fastq",'fastq'))
            )
    def get_qualities():
        return chain(
            (read.letter_annotations["phred_quality"] if read.seq.upper() in seq else read.letter_annotations["phred_quality"][::-1] for read in SeqIO.parse(f"/workspaces/Sequitur/data/input/Raphanus sativus_NC_018551.1/{args.input}.1.fastq",'fastq')),
            (read.letter_annotations["phred_quality"] if read.seq.upper() in seq else read.letter_annotations["phred_quality"][::-1] for read in SeqIO.parse(f"/workspaces/Sequitur/data/input/Raphanus sativus_NC_018551.1/{args.input}.2.fastq",'fastq'))
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