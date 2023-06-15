'''
DESCRIPTION 
    Utility function that creates a random sequence containing only the letters A, T, G, and C
INPUT
    n          | the length of the sequence
    palindrome | a boolean indicating whether the sequence must be a palidrome or not
    seed       | random seed for the random function for reproducibility
OUTPUT
    A random sequence of length n
'''
def generate_genome_sequence(n,palindrome=False,seed=None):
    import random,math
    
    random.seed(seed)
    nucleotides = {1:'A',2:'C',3:'G',4:'T'}
    seq = ''
    if palindrome: n = math.ceil(n/2)
    for _ in range(n):
        seq += nucleotides[random.randint(1,4)]
    if palindrome: seq += ''.join(reversed(seq[:int(n-math.fmod(n,2))]))
    return seq