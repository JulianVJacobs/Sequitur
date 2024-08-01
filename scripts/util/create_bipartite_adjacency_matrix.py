import argparse
from Bio import SeqIO
import glob
import os
import sqlite3
from fastDamerauLevenshtein import damerauLevenshtein as damerau_levenshtein_distance

def normalised_damerau_levenshtein_distance(read: str,overlap: str) -> float:
    m = min(len(overlap),len(read))
    d = damerau_levenshtein_distance(read[:m],overlap[:m],similarity=False)
    return (d/m,m-d)

def main(db: str, in_dir:str, acc: str, min_suf_len: int = 3, max_diff: float = 0.25):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    
    cursor.execute(f'''
        CREATE INDEX IF NOT EXISTS idx_{acc}_suffix ON {acc}(suffix)
    ''')
    
    conn.commit()
    
    pattern = os.path.join(in_dir, f'{acc}*.fastq')
    for file_path in glob.glob(pattern):
        B = {}
        for record in SeqIO.parse(file_path, 'fastq'):
            s_seq_id = f'{record.id}:{'1' if '1' == file_path[file_path.find('.fastq')-1] else '2'}'
            if s_seq_id not in B: B[s_seq_id] = {}
            for i in range(min_suf_len + 1):
                cursor.execute(f'''
                    SELECT COUNT(*) FROM {acc} WHERE suffix < ?
                ''', (str(record.seq[i:]),))
                offset = cursor.fetchone()[0] + 1
            
                cursor.execute(f'''
                    SELECT seq_id, suffix FROM {acc} ORDER BY suffix LIMIT -1 OFFSET ?
                ''', (offset,))
            
                for seq_id, suffix in cursor:
                    edit_distance_f, edit_distance_i = normalised_damerau_levenshtein_distance(str(record.seq[i:]),suffix)
                    if edit_distance_f > max_diff: break
                    if seq_id.startswith(s_seq_id): continue
                    p_seq_id = ''.join(seq_id.split(':')[:-1])
                    B[s_seq_id][p_seq_id] = max(B[s_seq_id][p_seq_id], edit_distance_i) if p_seq_id in B[s_seq_id] else edit_distance_i
    
    conn.close()
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
                                     Takes a SQLite database of suffixes, sorts it lexicographically
                                     and builds a bipartite adjacency matrix according to the
                                     Demerau-Levenshtein distance normalised to the length of the 
                                     shorter suffix. Connections are only made to the suffixes whose
                                     normalised edit distance is less than max_diff. Once max_diff is
                                     read, the program stops looking for connections to make.
                                                 """)
    parser.add_argument("in_dir", help="Directory where the inputs can be found")
    parser.add_argument("acc", help="Accession of FASTQ")
    parser.add_argument("db", help="SQLite DB which stores the suffixes")
    parser.add_argument("--max_diff", help="""
                        Maximum edit distance by which two strings can differ and still make a connection 
                        in the bipartite adjacency matrix (default 0.25)
                        """, default=0.25, type=float)
    parser.add_argument("--min_suf_len", help="Minimum length any suffix can be (default 3)", default=3, type=int)

    args = parser.parse_args()
    main(args.db,args.in_dir,args.acc,args.min_suf_len,args.max_diff)