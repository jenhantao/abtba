#!/usr/bin/env python
import sys


def homerPeak_to_bed(input_path, output_path):
    """
    reads a homer peak file and then writes a bed file
    """
    with open(input_path) as f:
        data = f.readlines()
    out_file = open(output_path, 'w')
    for line in data:
        if not line[0] == '#':
            tokens = line.strip().split()
            name = tokens[0]
            chrom = tokens[1]
            start = tokens[2]
            end = tokens[3]
            
            out_file.write('\t'.join([chrom, start, end, name, '\n'])
    out_file.close()
            


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Usage')
        print('convert_homer_peaks.py input_path output_path')
    input_path = sys.argv[1]
    output_path = sys.argv[2]

    homerPeak_to_bed(input_path, output_path)

