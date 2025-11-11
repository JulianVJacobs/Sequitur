import sqlite3
from Bio import SeqIO

def count_sequences_in_fastq(fastq_file):
    records = SeqIO.parse(fastq_file, 'fastq')
    count = sum(1 for _ in records)
    return count

def count_in_db(db_path):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    cursor.execute("SELECT COUNT(*) FROM suffix WHERE seq_id LIKE '%0'")
    count = cursor.fetchone()[0]
    
    conn.close()
    return count

if __name__ == "__main__":
    db_path = '/workspace/sequitur/data/input/ERS218597/suffix.db'
    fastq_file1 = '/workspace/sequitur/data/input/ERS218597/ERR234359.1.fastq'
    fastq_file2 = '/workspace/sequitur/data/input/ERS218597/ERR234359.2.fastq'
    
    count_db = count_in_db(db_path)
    print(f"Number of reads in db: {count_db}")
    
    count_seq = count_sequences_in_fastq(fastq_file1) + count_sequences_in_fastq(fastq_file2)
    print(f"Number of reads in both fastq files: {count_seq}")
    