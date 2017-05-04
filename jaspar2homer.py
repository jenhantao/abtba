#!/usr/bin/env python
"""
converts a jaspar motif to a homer motif file
"""

### imports ###
import sys
import numpy as np


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Usage')
        print('jaspar2homer.py <jaspar_motif> <output_file_path>')
        sys.exit(0)
    else:
        jaspar_path= sys.argv[1]
        output_path = sys.argv[2]

    with open(jaspar_path) as f:
        data = f.readlines()
    
    name_line = data[0]
    A_line = data[1]
    C_line = data[2]
    G_line = data[3]
    T_line = data[4]
    
    A_freqs = np.array([float(x) for x in A_line[4:-3].split()])
    C_freqs = np.array([float(x) for x in C_line[4:-3].split()])
    G_freqs = np.array([float(x) for x in G_line[4:-3].split()])
    T_freqs = np.array([float(x) for x in T_line[4:-3].split()])
    
    freqs = np.array([A_freqs, C_freqs, G_freqs, T_freqs])
    normed_freqs = freqs/freqs.sum(axis=0)
    normed_freqs = normed_freqs.T
    motif_name = name_line.strip().split()[-1]
    
    out_file = open(output_path, 'w')
    out_file.write(name_line)
    for i in range(normed_freqs.shape[0]):
        out_file.write('\t'.join([str(x) for x in normed_freqs[i]]) + '\n')
    out_file.close()
