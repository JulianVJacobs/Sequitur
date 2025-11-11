import os
import sqlite3
from Bio import SeqIO
import argparse
import glob

def main():
    min_suf_len: int = 3

    os.makedirs(os.path.dirname('/workspace/sequitur/data/input/ERS218597/suffix.db'), exist_ok=True)
    
    conn = sqlite3.connect('/workspace/sequitur/data/input/ERS218597/suffix.db')
    cursor = conn.cursor()
    
    cursor.execute(f'''
        CREATE TABLE IF NOT EXISTS ERS218597 (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            seq_id TEXT,
            suffix TEXT,
            UNIQUE(seq_id)
        )
    ''')
    
    conn.commit()
    i = 0
    for record in SeqIO.parse('/workspace/sequitur/data/input/ERS218597/ERR234359.1.fastq', 'fastq'):
        cursor.executemany(f'''
                                INSERT OR IGNORE INTO ERS218597 (seq_id,suffix) VALUES (?,?)
                            ''', [(f'{record.id}:2:{i}', 
                                str(record.seq.upper()[i:])) for i in range(len(record.seq) - (min_suf_len - 1))])
        i += 1
        if i == 10000:
            conn.commit()
            i = 0
    conn.commit()
    i = 0    
    for record in SeqIO.parse('/workspace/sequitur/data/input/ERS218597/ERR234359.2.fastq', 'fastq'):
        cursor.executemany(f'''
                                INSERT OR IGNORE INTO ERS218597 (seq_id,suffix) VALUES (?,?)
                            ''', [(f'{record.id}:2:{i}',  
                                    str(record.seq.upper().reverse_complement()[i:])) for i in range(len(record.seq) - (min_suf_len - 1))])
        i += 1
        if i == 10000:
            conn.commit()
            i = 0              
    
    conn.commit()
    conn.close()

if __name__ == "__main__":
    main()