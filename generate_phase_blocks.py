#!/usr/bin/env python3

"""
Generate unique phase blocks from a file with 3 columns (SNP chromosome,
position, and phase set) and write the phase block coordinates to a BED file.
"""

import sys

def main():
    # Check if the input file is provided as a command-line argument
    if len(sys.argv) != 2:
        print("Usage: python3 generate_phase_blocks.py <input_file>")
        sys.exit(1)

    # Read in the input file and generate phase blocks
    phase_blocks = {}
    with open(sys.argv[1], 'r') as f:
        # Initialize the current chromosome and position
        curr_chr = None
        curr_pos = None

        for line in f:
            snp_chr, pos, phase_set = line.strip().split()

            # Check if we have moved on to a new chromosome
            if snp_chr != curr_chr:
                # Write the phase blocks for the previous chromosome
                if curr_chr is not None:
                    write_phase_blocks(phase_blocks, curr_chr)

                # Reset the phase blocks dictionary for the new chromosome
                phase_blocks = {}
                curr_chr = snp_chr
                curr_pos = None

            # Check if we have moved on to a new position
            if curr_pos is None or int(pos) > curr_pos:
                curr_pos = int(pos)

            if phase_set not in phase_blocks:
                phase_blocks[phase_set] = [snp_chr, curr_pos, curr_pos]
            else:
                # Update the phase block end coordinate
                end_pos = int(pos)
                if end_pos > int(phase_blocks[phase_set][2]):
                    phase_blocks[phase_set][2] = end_pos

    # Write the phase blocks for the final chromosome
    write_phase_blocks(phase_blocks, curr_chr)

def write_phase_blocks(phase_blocks, chr_name):
    # Write the phase blocks for a single chromosome to a BED file
    with open(f"phase_blocks_{chr_name}.bed", 'w') as f:
        for phase_set, coords in phase_blocks.items():
            chrom, start, end = coords
            f.write(f"{chrom}\t{start}\t{end}\t{phase_set}\n")

if __name__ == '__main__':
    main()