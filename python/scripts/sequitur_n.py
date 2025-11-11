import time, os
from typing import Union, Tuple
from scipy.sparse import coo_matrix, vstack, hstack # pip install scipy
import numpy as np # pip install numpy
from fastDamerauLevenshtein import damerauLevenshtein as damerau_levenshtein_distance # pip install fastDamerauLevenshtein

def normalised_damerau_levenshtein_distance(
        read: str,
        overlap: str
    ) -> Tuple[float,int]:
    """
    Find the Damerau-Levenshtein edit distance of two strings normalised to the length
    of the shorter string. This normalisation is because we want to path prefixes to
    suffixes and this means that in general we will be comparing a full string to a
    portion of another string.
    
    Parameters:
    - read (str)   : string for comparison, usually the longer string 
    - overlap (str): string for comparison, usually the shorter string
    
    Returns:
    - float: the normalised Demarau-Levenshtein edit distance of the input strings as a ratio
    - int  : the normalised Demarau-Levenshtein edit distance of the input strings as an integer
    """
    m = min(len(overlap),len(read))
    d = damerau_levenshtein_distance(str(read)[:m],str(overlap)[:m],similarity=False)
    return (d/m,m-d)

def build_suffix_array(
        reads: list, 
        min_suf_len: int = 3,
        do_time: bool = False
    ) -> Tuple[list, list]:
    if do_time: start = time.time()
    suf_arr = []
    for read in reads:
        read += '$' + str(reads.index(read))
        for i in range(len(read)-min_suf_len-1):
            suf_arr += [read[i:]]
    suf_arr.sort()
    suf_arr_ind = []
    for s in range(len(suf_arr)):
        suf_arr_ind += [int(suf_arr[s].split('$')[-1].__str__())]
        suf_arr[s] = suf_arr[s][:suf_arr[s].find('$')+1]
    if do_time: return suf_arr, suf_arr_ind, time.time() - start
    return suf_arr,suf_arr_ind

def create_bipartite_adjacency_matrix(
        reads: list, 
        suf_arr: list | None = None, 
        suf_arr_ind: list | None = None, 
        do_time: bool = False,
        max_diff: float = 0.25, 
        min_suf_len: int = 3
    ) -> dict:
    
    if do_time: start = time.time()
    if suf_arr is None or suf_arr_ind is None: suf_arr,suf_arr_ind = build_suffix_array(reads,min_suf_len=min_suf_len)
    reads_map = dict(zip(reads,list(range(len(reads)))))
    B, O = {}, {}
    for read in reads:
        for j in range(min_suf_len + 1):
            i = suf_arr.index(read[j:]+'$') - 1
            edf, edi = normalised_damerau_levenshtein_distance(read,suf_arr[i][:-1])
            while edf <= max_diff and i >= 0:
                if not reads[suf_arr_ind[i]] == read and edf < max_diff and \
                   read.startswith(suf_arr[i][:-1]):
                    if (reads_map[reads[suf_arr_ind[i]]],reads_map[read]) not in B: 
                        B[(reads_map[reads[suf_arr_ind[i]]],reads_map[read])] = edi
                        O[(reads_map[reads[suf_arr_ind[i]]],reads_map[read])] = len(suf_arr[i][:-1])
                    else: 
                        if edi == B[(reads_map[reads[suf_arr_ind[i]]],reads_map[read])] and len(suf_arr[i][:-1]) > O[(reads_map[reads[suf_arr_ind[i]]],reads_map[read])]:
                            O[(reads_map[reads[suf_arr_ind[i]]],reads_map[read])] = max(len(suf_arr[i][:-1]),O[(reads_map[reads[suf_arr_ind[i]]],reads_map[read])])
                        elif edi > B[(reads_map[reads[suf_arr_ind[i]]],reads_map[read])]:
                            B[(reads_map[reads[suf_arr_ind[i]]],reads_map[read])] = edi
                            O[(reads_map[reads[suf_arr_ind[i]]],reads_map[read])] = len(suf_arr[i][:-1])
                i -= 1
    if do_time: return B, O, time.time() - start
    return B, O

def move_col(
        B: coo_matrix, 
        cols: dict
    ) -> None:
    for c in range(len(B.col)):
        B.col[c] = cols[B.col[c]]
            
def move_row(
        B: coo_matrix,
        rows: dict
    ) -> None:
    for r in range(len(B.row)):
        B.row[r] = rows[B.row[r]]

def find_lower_diagonal_path(
        B: coo_matrix,
        reads_map: dict,
        cols: list,
        rows: list,
        do_time: bool = False
    ) -> Union[str, Tuple[str, float]]:
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
    # path = '/workspace/sequitur/data/output'
    path = './data/output'
    for i, (seq, reads) in enumerate(nat_lang_seq):
        if not os.path.exists(f"{path}/natural_language_sequences.sequitur.csv") or os.path.getsize(f"{path}/natural_language_sequences.sequitur.csv") == 0:
            with open(f"{path}/natural_language_sequences.sequitur.csv", 'a') as f:
                f.write('natural_language_sequence,edit_distance,target_sequence_length,output_sequence_length,suffix_array_construction_time,adjacency_matrix_construction_time,sequence_reconstruction_time\n')
        with open(f"{path}/natural_language_sequences.sequitur.csv",'a') as f:
            reads_map = dict(zip(list(range(len(reads))),reads))
            rows = list(range(len(reads)))
            cols = list(range(len(reads)))
            for _ in range(100):
                suf_arr,suf_arr_ind,t1 = build_suffix_array(reads,min_suf_len=2,do_time=True)
                B,t2 = create_bipartite_adjacency_matrix(reads,suf_arr=suf_arr,suf_arr_ind=suf_arr_ind,do_time=True)
                start = time.time()
                B = coo_matrix((list(B.values()),list(zip(*B.keys())))).T
                if B.shape[0] < len(rows): B = vstack([B,coo_matrix((1, B.shape[1]),dtype=B.dtype)])
                if B.shape[1] < len(cols): B = hstack([B,coo_matrix((B.shape[0], 1),dtype=B.dtype)])
                t2 += time.time() - start
                seq_,t3 = find_lower_diagonal_path(B,reads_map,cols,rows,do_time=True)
                f.write('{},{},{},{},{},{},{}\n'.format(i,damerau_levenshtein_distance(seq_,seq),len(seq),len(seq_),t1,t2,t3))