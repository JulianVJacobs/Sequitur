#!/usr/bin/env python3
"""
Validate assembly against reference using edlib alignment.

Works with both DNA sequences (FASTA) and text converted to DNA.
Reports edit distance, accuracy, coverage, and alignment visualization.

Usage:
    python validate_assembly.py reference.fasta assembled.fasta [--detailed]
"""

import argparse
import sys
from pathlib import Path

try:
    import edlib
except ImportError:
    print("Error: edlib not installed. Run: pip install edlib", file=sys.stderr)
    sys.exit(1)


def read_fasta(path):
    """Read FASTA file, return (header, sequence)."""
    with open(path, 'r') as f:
        lines = f.readlines()
    
    header = None
    sequence = []
    
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if header is None:
                header = line[1:]
            else:
                # Multiple sequences, concatenate
                pass
        elif line:
            sequence.append(line)
    
    return header, ''.join(sequence)


def format_alignment(ref, asm, cigar, max_width=80):
    """Format alignment with CIGAR for visualization."""
    lines = []
    ref_pos = 0
    asm_pos = 0
    ref_line = []
    match_line = []
    asm_line = []
    
    # Parse CIGAR
    i = 0
    while i < len(cigar):
        # Read number
        num_str = ''
        while i < len(cigar) and cigar[i].isdigit():
            num_str += cigar[i]
            i += 1
        
        count = int(num_str) if num_str else 1
        op = cigar[i] if i < len(cigar) else 'M'
        i += 1
        
        for _ in range(count):
            if op == '=' or op == 'M':  # Match
                ref_line.append(ref[ref_pos])
                match_line.append('|')
                asm_line.append(asm[asm_pos])
                ref_pos += 1
                asm_pos += 1
            elif op == 'X':  # Mismatch
                ref_line.append(ref[ref_pos])
                match_line.append('X')
                asm_line.append(asm[asm_pos])
                ref_pos += 1
                asm_pos += 1
            elif op == 'I':  # Insertion in assembly
                ref_line.append('-')
                match_line.append(' ')
                asm_line.append(asm[asm_pos])
                asm_pos += 1
            elif op == 'D':  # Deletion in assembly
                ref_line.append(ref[ref_pos])
                match_line.append(' ')
                asm_line.append('-')
                ref_pos += 1
            
            # Wrap lines at max_width
            if len(ref_line) >= max_width:
                lines.append('Ref: ' + ''.join(ref_line))
                lines.append('     ' + ''.join(match_line))
                lines.append('Asm: ' + ''.join(asm_line))
                lines.append('')
                ref_line = []
                match_line = []
                asm_line = []
    
    # Add remaining
    if ref_line:
        lines.append('Ref: ' + ''.join(ref_line))
        lines.append('     ' + ''.join(match_line))
        lines.append('Asm: ' + ''.join(asm_line))
    
    return '\n'.join(lines)


def validate_assembly(ref_path, asm_path, detailed=False):
    """Validate assembly against reference."""
    print(f"Reading reference: {ref_path}")
    ref_header, ref_seq = read_fasta(ref_path)
    print(f"  Reference: {ref_header}")
    print(f"  Length: {len(ref_seq)} bases")
    
    print(f"\nReading assembly: {asm_path}")
    asm_header, asm_seq = read_fasta(asm_path)
    print(f"  Assembly: {asm_header}")
    print(f"  Length: {len(asm_seq)} bases")
    
    # Align
    print("\nAligning with edlib...")
    result = edlib.align(ref_seq, asm_seq, mode="NW", task="path")
    
    edit_dist = result['editDistance']
    cigar = result['cigar']
    
    # Calculate metrics
    max_len = max(len(ref_seq), len(asm_seq))
    accuracy = 1 - (edit_dist / max_len) if max_len > 0 else 0
    len_ratio = len(asm_seq) / len(ref_seq) if len(ref_seq) > 0 else 0
    
    # Count operations from CIGAR
    matches = cigar.count('=') if cigar else 0
    mismatches = cigar.count('X') if cigar else 0
    insertions = cigar.count('I') if cigar else 0
    deletions = cigar.count('D') if cigar else 0
    
    # Report
    print("\n" + "="*60)
    print("VALIDATION RESULTS")
    print("="*60)
    print(f"Edit distance:     {edit_dist:,}")
    print(f"Accuracy:          {accuracy*100:.2f}%")
    print(f"Length ratio:      {len_ratio:.3f} (asm/ref)")
    print(f"\nAlignment operations:")
    print(f"  Matches:         {matches:,}")
    print(f"  Mismatches:      {mismatches:,}")
    print(f"  Insertions:      {insertions:,}")
    print(f"  Deletions:       {deletions:,}")
    
    if detailed and cigar:
        print("\n" + "="*60)
        print("ALIGNMENT VISUALIZATION (first 500 bases)")
        print("="*60)
        # Show first 500 bases of alignment
        snippet_cigar = cigar  # Could truncate if needed
        vis = format_alignment(ref_seq[:500], asm_seq[:500], snippet_cigar)
        print(vis)
    
    # Success criteria
    print("\n" + "="*60)
    if accuracy >= 0.95:
        print("✓ PASS: Assembly accuracy >= 95%")
        return 0
    elif accuracy >= 0.90:
        print("⚠ WARN: Assembly accuracy 90-95% (acceptable but not ideal)")
        return 0
    else:
        print("✗ FAIL: Assembly accuracy < 90%")
        return 1


def main():
    parser = argparse.ArgumentParser(
        description='Validate assembly against reference using edlib',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('reference', type=Path, help='Reference FASTA file')
    parser.add_argument('assembly', type=Path, help='Assembled FASTA file')
    parser.add_argument('--detailed', action='store_true',
                       help='Show detailed alignment visualization')
    
    args = parser.parse_args()
    
    if not args.reference.exists():
        print(f"Error: Reference file not found: {args.reference}", file=sys.stderr)
        return 1
    
    if not args.assembly.exists():
        print(f"Error: Assembly file not found: {args.assembly}", file=sys.stderr)
        return 1
    
    return validate_assembly(args.reference, args.assembly, args.detailed)


if __name__ == '__main__':
    sys.exit(main())
