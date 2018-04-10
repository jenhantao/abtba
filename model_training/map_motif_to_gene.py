#!/usr/bin/env python
"""
Given a TBA coefficients file or a TBA significance file, maps the 
motif names to gene names
"""

### imports ###
import argparse
import numpy as np
import os    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Given a TBA coefficients file or \
    a TBA significance file, maps the motif names to gene names')
    parser.add_argument("result_path",
        help="path to a TBA coefficients or significance file",
        type = str)
    parser.add_argument("output_path",
        help="file path where output should be written",
        type=str)
    parser.add_argument("motif_files",
        help="list of motif files",
        type=str,
        nargs="+")
    
    # parse arguments
    args = parser.parse_args()

    result_path = args.result_path
    motif_files = args.motif_files
    output_path = args.output_path
    
    motif_genes_dict = {}
    for mf in motif_files:
        with open(mf) as f:
            data = f.readlines()
        gene_names = data[0].replace('_var.2','').split()[0][1:].split('|')
        motif_name = data[0].strip().split()[1]
        motif_genes_dict[motif_name] = gene_names
    
    with open(result_path) as f:
        results = f.readlines()
    
    out_file = open(output_path, 'w')
    out_file.write(results[0])
    
    for line in results[1:]:
        motif_name = line.split()[0]
        motif_genes = motif_genes_dict[motif_name]
        for g in motif_genes:
            out_file.write(line.replace(motif_name,g))
    
    out_file.close()
