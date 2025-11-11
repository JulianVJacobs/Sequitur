import argparse
from Bio import SeqIO
import glob
import os
import sqlite3
from scipy.sparse import lil_matrix, save_npz
import pickle

from fastDamerauLevenshtein import damerauLevenshtein as damerau_levenshtein_distance

def normalised_damerau_levenshtein_distance(read: str,overlap: str) -> (float,int):
    m = min(len(overlap),len(read))
    d = damerau_levenshtein_distance(read[:m],overlap[:m],similarity=False)
    return (d/m,m-d)

def next_suffix(cursor: sqlite3.Cursor, offset: int, batch_size: int = 100):
    cursor.execute(f'''
        SELECT seq_id, suffix FROM ERS218597 ORDER BY suffix LIMIT -1 OFFSET ?
    ''', (offset,))
    while True:
        rows = cursor.fetchmany(batch_size)
        if not rows: break
        for row in rows: yield row

def main(min_suf_len: int = 3, max_diff: float = 0.25):
    conn = sqlite3.connect('/workspace/sequitur/data/input/ERS218597/suffix.db')
    cursor = conn.cursor()
    
    cursor.execute(f'''
        CREATE INDEX IF NOT EXISTS idx_ERS218597_suffix ON ERS218597(suffix)
    ''')
    
    conn.commit()
    
    pattern = os.path.join('/workspace/sequitur/data/input/ERS218597/', 'ERR234359*.fastq')
    B = {}
    for file_path in glob.glob(pattern):
        for record in SeqIO.parse(file_path, 'fastq'):
            s_seq_id = record.id
            if s_seq_id not in B: B[s_seq_id] = {}
            for i in range(min_suf_len + 1):
                cursor.execute(f'''
                    SELECT COUNT(*) FROM ERS218597 WHERE suffix < ?
                ''', (str(record.seq[i:]),))
                offset = cursor.fetchone()[0] + 1
            
                for seq_id, suffix in next_suffix(cursor, offset):
                    edit_distance_f, edit_distance_i = normalised_damerau_levenshtein_distance(str(record.seq),suffix)
                    if edit_distance_f > max_diff: break
                    if seq_id.startswith(s_seq_id): continue
                    p_seq_id = ''.join(seq_id.split(':')[:-2])
                    B[s_seq_id][p_seq_id] = max(B[s_seq_id][p_seq_id], edit_distance_i) if p_seq_id in B[s_seq_id] else edit_distance_i
    
    row_indices = []
    col_indices = []
    data = []
    seq_id_to_index = {seq_id: idx for idx, seq_id in enumerate(B.keys())}
    
    for s_seq_id, overlaps in B.items():
        for p_seq_id, edit_distance_i in overlaps.items():
            row_indices.append(seq_id_to_index[s_seq_id])
            col_indices.append(seq_id_to_index[p_seq_id])
            data.append(edit_distance_i)
    
    sparse_matrix = lil_matrix((len(seq_id_to_index), len(seq_id_to_index)))
    sparse_matrix[row_indices, col_indices] = data
    
    sparse_matrix_path = os.path.join('/workspace/sequitur/data/input/ERS218597/', 'ERS218597.matrix.npz')
    save_npz(sparse_matrix_path, sparse_matrix)
    
    seq_id_to_index_path = os.path.join('/workspace/sequitur/data/input/ERS218597/', 'ERS218597.index_map.pkl')
    with open(seq_id_to_index_path, 'wb') as f:
        pickle.dump(seq_id_to_index, f)
    
    conn.close()
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--min_suf_len", help="Minimum length any suffix can be (default 3)", default=3, type=int)
    parser.add_argument("--max_diff", help="""
                        Maximum edit distance by which two strings can differ and still make a connection 
                        in the bipartite adjacency matrix (default 0.25)
                        """, default=0.25, type=float)

    args = parser.parse_args()
    main(args.min_suf_len,args.max_diff)