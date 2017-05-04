#!/usr/bin/env python
"""
Given a set of genomic coordinates in BED format:
chr start end
...

Extracts the genomic sequences
"""

### imports ###
import sys


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Usage')
        print('fimo2homer.py <fimo_motif> <output_file_path>')
        sys.exit(0)
    else:
        fimo_path= sys.argv[1]
        output_path = sys.argv[2]

    with open(fimo_path) as f:
        data = f.readlines()
    
    nameline = data[8].strip()
    name_tokens = nameline.split()
    motif_name = name_tokens[1]
    alternate_motif_name = name_tokens[2]

    out_file = open(output_path, 'w')
    out_file.write('\t'.join(['>' + alternate_motif_name, motif_name, '\n']))
    for line in data[10:]:
        tokens = line.strip().split()
        if len(tokens) > 1:
            out_line = '\t'.join(tokens) + '\n'
            out_file.write(out_line)
    out_file.close()
     
