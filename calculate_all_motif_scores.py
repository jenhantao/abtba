#!/usr/bin/env python
"""
Given a fasta file, and a set of PWMs in Homer format calculates the top motif
match for each motif for each sequence:
"""

### imports ###
import argparse
import os
import numpy as np
import time
import multiprocessing
import pickle
import pandas as pd

### functions ###
def read_motif_file(motifPath, pseudocount):
    '''
    reads all motif files in a directory 
    inputs: path to a directory containing homer motif files
    outputs: an array of tuples representing each motif
    '''
    name_metadata_dict = {}
    with open(motifPath) as f:
        data = f.readlines()
    name = '.'.join(motifPath.split('/')[-1].split('.')[:-1])
    matrix = []
    metadata = data[0].strip()
    for line in data[1:]:
        tokens = line.strip().split("\t")
        if len(tokens) > 1:
            scores = np.array([float(x) for x in tokens])
            scores = scores + pseudocount
            scores= scores/np.sum(scores)
            matrix.append(scores)
    return (name,np.array(matrix))

def read_fasta(file_path):
    '''
    reads in a fasta file and returns a list of sequence ids and a list of sequences
    inputs: filepath - path to a fasta file
    outputs: sequence_list - a list of sequences
             id_list - a list of ids
    '''
    with open(file_path) as f:
        data = f.readlines()
    id_list = []
    sequence_list = []
    # loop through each sequence
    for line in data:
        if '>' in line:
            seq_id = line.strip()[1:]
            id_list.append(seq_id)
        else:
            seq = line.strip()
            sequence_list.append(seq)
    return sequence_list, id_list

def convert_sequences_to_array(sequences):
    '''
    inputs: sequence of nucleotides represented as a string composed of A, C, G, T
    outputs: a list of numpy array representations of a sequence with:
             A = [1, 0, 0, 0]
             C = [0, 1, 0, 0]
             G = [0, 0, 1, 0]
             T = [0, 0, 0, 1]
             
    '''
    nucleotide_array_dict = {'A': [1, 0, 0, 0],
                             'C': [0, 1, 0, 0],
                             'G': [0, 0, 1, 0],
                             'T': [0, 0, 0, 1],
                             'N': [0, 0, 0, 0]}
    sequence_array_list = []
    for seq in sequences:
        seq_array = []
        for nuc in seq:
            seq_array.append(nucleotide_array_dict[nuc])
        seq_array = np.array(seq_array)
        sequence_array_list.append(seq_array)
    return sequence_array_list

def calculate_all_motif_scores_async(sequence_array_list,
                                 pwm,
                                 motif_name,
                                 output_path
                                 ):
    '''
    identifies the highest scoring match to pwm for each sequence
    inputs: pwm - a numpy array representing a pwm
            sequence_array_list - a list of numpy array representations of a sequence with:
                                 A = [1, 0, 0, 0]
                                 C = [0, 1, 0, 0]
                                 G = [0, 0, 1, 0]
                                 T = [0, 0, 0, 1]
    outputs: top_scores - a list of the best motif scores in each sequence
    '''
    start = time.time()

    background_frequency = 0.25

    all_scores = [] # motif score of best match for each sequence
    all_rc_scores = []
    pwm_length = pwm.shape[0]
    # calculate scores for each motif at each position
    for seq_array in sequence_array_list:
        seq_length = seq_array.shape[0]
        scores = []
        for i in range(seq_length - pwm_length + 1):
            # get substring represented as matrix
            subseq_array = seq_array[i: i + pwm_length]
            # get corresponding pwm frequencies
            frequencies = (pwm * subseq_array).sum(axis=1)
            # 0.25 background freq
            ratios = (frequencies)/(background_frequency)
            # calculate log likelihood ratios
            llr = np.log2(ratios)
            # sum to calculate motif score
            score = np.sum(llr)
            scores.append(score)
        
        # calculate reverse complement and scores for reverse complement
        rc_seq_array = seq_array[::-1, ::-1]
        rc_scores = []
        for i in range(seq_length - pwm_length + 1):
            # get substring represented as matrix
            subseq_array = rc_seq_array[i: i + pwm_length]
            # get corresponding pwm frequencies
            frequencies = (pwm  * subseq_array).sum(axis=1)
            # 0.25 background freq
            ratios = (frequencies)/(background_frequency)
            # calculate log likelihood ratios
            llr = np.log2(ratios)
            # sum to calculate motif score
            score = np.sum(llr)
            rc_scores.append(score)
            
        rc_scores=rc_scores[::-1]
        all_scores.append(scores)
        all_rc_scores.append(rc_scores)

    pickle.dump(all_scores,
                open(output_path+'/'+motif_name+'.pickle','wb'))
    pickle.dump(all_rc_scores,
                open(output_path+'/'+motif_name+'.rc.pickle','wb'))

    end = time.time()
    print(motif_name, 'calculation time:', end-start)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='given a fasta file and a list \
                                     and a list of motif files, calculates the \
                                     the best match to each motif for each \
                                     sequence' )
    parser.add_argument("fasta_path",
        help="path to a fasta_file containing sequences to score",
        type = str)
    parser.add_argument("output_path",
        help="directory where output file should be written",
        default="./", type=str)
    parser.add_argument("motif_files",
        help="list of motif files",
        type=str,
        nargs="+")
    parser.add_argument("-num_procs", 
        help="number of processor cores to use",
        type=int,
        default=4)
    parser.add_argument("-pseudocount", 
        help="pseudocount for calculating motif scores",
        type=float,
        default=0.001)

    # parse arguments
    args = parser.parse_args()

    fasta_path = args.fasta_path
    output_path = args.output_path
    motif_files = args.motif_files
    num_processors = args.num_procs
    pseudocount = args.pseudocount
    
    if not os.path.isdir(output_path):
        os.mkdir(output_path)

    # read in motif files
    all_motifs = []
    for m in motif_files:
        motif = read_motif_file(m, pseudocount)
        all_motifs.append(motif)
    # sort motifs by name
    all_motifs.sort(key=lambda x:x[0])

    start = time.time()
    
    sequence_list, id_list = read_fasta(fasta_path)

    # convert strings to arrays
    sequence_array_list = convert_sequences_to_array(sequence_list)
    
    # calculate motif scores
    pool = multiprocessing.Pool(processes=num_processors)
    manager = multiprocessing.Manager()
       
    for motif in all_motifs:
        pwm = motif[1]
        motif_name = motif[0]
        
        pool.apply_async(
        calculate_all_motif_scores_async ,args=(sequence_array_list, 
                                                pwm, 
                                                motif_name, 
                                                output_path
                                           )
                        )
    pool.close()
    pool.join()
    
    end = time.time()
    print('total time', end-start)
