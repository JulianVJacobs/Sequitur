'''
DESCRIPTION 
    Utility function that chops a sequence into several reads with bounded random lengths that 
    have a bounded random overlap
INPUT
    sequence       | a sequence of characters that will be divided into overlapping subsequences
    min_subseq_len | the shortest length a subsequence can have
    max_subseq_len | the longest length a subsequence can have
    min_overlap    | the shortest overlap two subsequences can share
    max_overlap    | the longest overlap two subsequences can share
    circularize    | boolean indicating whether to add a random amount of the end of the sequence
                   | to the beginning and vice versa
    seed           | random seed for the random function for reproducibility
OUTPUT
    A list of overlapping reads of random bounded size which share a bounded random amount of
    overlap
'''
def generate_reads(sequence,min_subseq_len,max_subseq_len,min_overlap,max_overlap,circularise=False,seed=None):
    import random

    random.seed(seed)
    if circularise: 
        sequence = sequence[-random.randint(min_overlap,max_overlap):] + sequence + sequence[:random.randint(min_overlap,max_overlap)]
        print(len(sequence))
    reads = []
    start = 0
    end = random.randint(min_subseq_len,max_subseq_len)
    reads += [sequence[start:end]]
    while end < len(sequence):
        start = random.randint(end-max_overlap,end-min_overlap)
        if (len(sequence) - start)/max_subseq_len < 2:
            if (len(sequence) - start)/max_subseq_len < 1:
                end = len(sequence)
            else:
                a = 0
                while (len(sequence) - start)/(min_subseq_len+a) > 2: a+=1
                end = random.randint(start+min_subseq_len+a,start+max_subseq_len) 
        else: end = random.randint(start+min_subseq_len,start+max_subseq_len) 
        reads += [sequence[start:end]]
    return reads