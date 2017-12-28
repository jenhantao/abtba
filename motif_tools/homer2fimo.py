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
        print('fimo2homer.py <homer_motif> <output_file_path>')
        sys.exit(0)
    else:
        homer_path= sys.argv[1]
        output_path = sys.argv[2]

    out_file = open(output_path, 'w')
    with open(homer_path) as f:
        data = f.readlines()
    motif_name = data[0].split()[1]
    motif_id = data[0].split()[0][1:]
    out_file.write('MEME version 4.4\n\nALPHABET= ACGT\n\nstrands: + -\n\n' +
                   'Background letter frequencies (from uniform background):\nA 0.25000 C 0.25000 G 0.25000 T 0.25000\n' +
                   'MOTIF '+ motif_name + ' ' + motif_id  + '\n')
    out_file.write('letter-probability matrix: nsites= 20 alength= 4 w= '+str(len(data)-1)+' E= 0 \n')
    for line in data[1:]:
        out_file.write('  ' + line)
    out_file.write('\n')
    out_file.close()

 
