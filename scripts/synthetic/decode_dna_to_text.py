#!/usr/bin/env python3
"""
Decode DNA sequences back to text using the reverse mapping from text_to_fastq.py.

Usage:
    python decode_dna_to_text.py assembled.fasta output.txt
"""

import argparse
import sys
from pathlib import Path


def char_to_dna(char):
    """Same fixed 4-base encoding as text_to_fastq.py."""
    code = ord(char)
    
    if code < 32 or code > 126:
        code = ord('?')
    
    idx = code - 32
    bases = ['A', 'C', 'G', 'T']
    
    # Fixed 4-base encoding (256 possible values, we use 95 for ASCII 32-126)
    b1 = bases[idx // 64]
    b2 = bases[(idx // 16) % 4]
    b3 = bases[(idx // 4) % 4]
    b4 = bases[idx % 4]
    return b1 + b2 + b3 + b4


def build_reverse_map():
    """Build DNA codon -> character reverse mapping."""
    reverse = {}
    
    # Map all printable ASCII
    for code in range(32, 127):
        char = chr(code)
        dna = char_to_dna(char)
        if dna not in reverse:  # First mapping wins
            reverse[dna] = char
    
    return reverse


def decode_sequence(dna_seq, reverse_map):
    """Decode DNA sequence to text with fixed 4-base codon support."""
    result = []
    i = 0
    
    while i < len(dna_seq):
        # Fixed 4-base codon
        if i + 4 <= len(dna_seq):
            codon = dna_seq[i:i+4]
            if codon in reverse_map:
                result.append(reverse_map[codon])
            else:
                result.append('?')
            i += 4
        else:
            # Handle incomplete codon at end
            result.append('?')
            i += 1
    
    return ''.join(result)


def read_fasta(path):
    """Read FASTA file, return list of (header, sequence) tuples."""
    sequences = []
    current_header = None
    current_seq = []
    
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    sequences.append((current_header, ''.join(current_seq)))
                current_header = line[1:]
                current_seq = []
            elif line:
                current_seq.append(line)
    
    if current_header is not None:
        sequences.append((current_header, ''.join(current_seq)))
    
    return sequences


def main():
    parser = argparse.ArgumentParser(
        description='Decode DNA sequences back to text',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('input', type=Path, help='Input FASTA file with DNA sequences')
    parser.add_argument('output', type=Path, help='Output text file')
    parser.add_argument('--show-unknown', action='store_true',
                       help='Show unknown codons as ? instead of omitting')
    
    args = parser.parse_args()
    
    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        return 1
    
    # Build reverse mapping
    print("Building DNA -> text reverse mapping...")
    reverse_map = build_reverse_map()
    print(f"  Mapped {len(reverse_map)} codons")
    
    # Read sequences
    print(f"Reading DNA sequences from {args.input}")
    sequences = read_fasta(args.input)
    print(f"  Found {len(sequences)} sequence(s)")
    
    # Decode and write
    print(f"Decoding to text and writing to {args.output}")
    with open(args.output, 'w') as out:
        for i, (header, dna_seq) in enumerate(sequences):
            print(f"\n  Sequence {i+1}: {header}")
            print(f"    DNA length: {len(dna_seq)} bases ({len(dna_seq)//3} codons)")
            
            text = decode_sequence(dna_seq, reverse_map)
            
            unknown_count = text.count('?')
            if unknown_count > 0:
                print(f"    Warning: {unknown_count} unknown codons decoded as '?'")
            
            # Write separator if multiple sequences
            if i > 0:
                out.write("\n" + "="*60 + "\n")
                out.write(f"SEQUENCE {i+1}: {header}\n")
                out.write("="*60 + "\n\n")
            
            out.write(text)
            
            if not args.show_unknown:
                # Show stats about decoded text
                visible_chars = len([c for c in text if c != '?'])
                print(f"    Decoded {visible_chars} visible characters")
    
    print(f"\nDone! Decoded text written to {args.output}")
    print(f"Open the file to see what the assembler produced.")


if __name__ == '__main__':
    sys.exit(main())
