# generate/import sequence
# split sequence
# create Sequitur object and assemble
# output to file

if __name__ == "__main__":
    import generate_reads
    import Sequitur
    from Bio import SeqIO

    for record in SeqIO.parse("data/Raphanus sativus_NC_018551.1.fasta",'fasta'):
        sequence = record.seq
    reads = generate_reads(sequence,250,500,50,100,seed=0)
    sequitur = Sequitur(reads,assemble=True,biphasic=False)
    f = open("output.txt", "x")
    f.write(sequitur.sequence)
    f.write(sequence)
