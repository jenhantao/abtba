#!/usr/bin/env python
"""
converts a homer motif to a jaspar motif file
"""

### imports ###
import sys
import numpy as np


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Usage')
        print('homer2jaspar.py <homer_motif> <output_file_path>')
        sys.exit(0)
    else:
        homer_path= sys.argv[1]
        output_path = sys.argv[2]

    with open(homer_path) as f:
        data = f.readlines()
    
    name_line = data[0]
    name_tokens = name_line.strip().split('\t')
    motif_id = name_tokens[0][1:]
    motif_name = name_tokens[1]

    A_freqs = []
    C_freqs = []
    G_freqs = []
    T_freqs = []
    
    for line in data[1:]:
        tokens = line.strip().split('\t')
        A_freqs.append(str(int(float(tokens[0])*10000)))
        C_freqs.append(str(int(float(tokens[1])*10000)))
        G_freqs.append(str(int(float(tokens[2])*10000)))
        T_freqs.append(str(int(float(tokens[3])*10000)))
    
    out_file = open(output_path, 'w')

    out_file.write('>'+motif_id + '\t' + motif_name + '\n')
    out_file.write('A [ ' + ' '.join(A_freqs) + ' ]\n') 
    out_file.write('C [ ' + ' '.join(C_freqs) + ' ]\n') 
    out_file.write('G [ ' + ' '.join(G_freqs) + ' ]\n') 
    out_file.write('T [ ' + ' '.join(T_freqs) + ' ]\n') 
    out_file.close()
    
