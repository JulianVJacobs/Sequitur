import os
import sqlite3
from Bio import SeqIO
import argparse
import glob
from typing import Generator

def create_database(db_path: str, acc: str):
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute(f'''
        CREATE TABLE IF NOT EXISTS {acc} (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            seq_id TEXT,
            suffix TEXT,
            UNIQUE(seq_id)
        )
    ''')
    conn.commit()
    conn.close()

def insert_suffixes(db_path: str, acc:str, suffixes: Generator):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.executemany(f'''
            INSERT OR IGNORE INTO {acc} (seq_id,suffix) VALUES (?,?)
        ''', [(seq_id,suffix) for seq_id, suffix in suffixes])
    conn.commit()
    conn.close()

def main(acc, in_dir, db_dir, min_suf_len: int = 3):
    db_path = os.path.join(db_dir, 'suffix.db')

    create_database(db_path, acc)

    pattern = os.path.join(in_dir, f'{acc}*.fastq')
    for file_path in glob.glob(pattern):
        for record in SeqIO.parse(file_path, 'fastq'):
            if '1' == file_path[file_path.find('.fastq')-1]:
                insert_suffixes(db_path, acc, ((f'{record.id}:1:{i}',str(record.seq.upper()[i:])) for i in range(len(record.seq) - (min_suf_len - 1))))
            else:
                insert_suffixes(db_path, acc, ((f'{record.id}:2:{i}',str(record.seq.upper().reverse_complement()[i:])) for i in range(len(record.seq) - (min_suf_len - 1))))

    # retrieve_sorted_suffixes(db_path, output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
                                     Takes two sets of paired read FASTQ data, orients them to face the same direction
                                     and then writes all of the reads and their suffixes to an SQLite Database. The
                                     databse is then sorted to produce a suffix array.
                                                 """)
    parser.add_argument("in_dir", help="Directory where the inputs can be found")
    parser.add_argument("acc", help="Accession of FASTQ")
    parser.add_argument("db_dir", help="Directory of SQLite DB to store suffixes")
    parser.add_argument("--min_suf_len", help="Minimum length any suffix can be (default 3)", default=3, type=int)

    args = parser.parse_args()
    main(args.acc, args.in_dir,args.db_dir,args.min_suf_len)