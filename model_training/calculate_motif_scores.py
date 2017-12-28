#!/usr/bin/env python
"""
Given a fasta file, and a set of PWMs in Homer format calculates the top motif
match for each motif for each sequence:
"""

### imports ###
import argparse
import numpy as np
import time
import multiprocessing
import pandas as pd
import os

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
    current_seq_tokens = []
    for line in data:
        if '>' in line:
            if len(current_seq_tokens) > 0:
                seq = ''.join(current_seq_tokens)
                id_list.append(seq_id)
                sequence_list.append(seq)
                current_seq_tokens = [] 
            seq_id = line.strip()[1:]
        else:
            current_seq_tokens.append(line.strip())
    sequence_list.append(seq)
    id_list.append(seq_id)
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

def calculate_top_motif_matches_async(sequence_array_list, 
                                      pwm, 
                                      motif_name, 
                                      motif_score_dict, 
                                      motif_start_dict,
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
             top_starts - a list of the start position of the best motif match in each sequence
    '''
    start = time.time()
    
    background_frequency = 0.25
    
    top_scores = [] # motif score of best match for each sequence
    top_starts = [] # start position of best match for each sequence
    
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
            scores.append((score, i, '+'))
            
        # calculate reverse complement and scores for reverse complement
        rc_seq_array = seq_array[::-1, ::-1]
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
            scores.append((score, seq_length - i, '-')) 
            
        scores.sort(key = lambda x:x[0], reverse = True)
        top_hit = scores[0]
        if top_hit[0] > 0:
            top_scores.append(top_hit[0])
            top_starts.append(str(top_hit[1]) + ' ' + top_hit[2])
        else:
            top_scores.append(0)
            top_starts.append('-1 ?')
    motif_score_dict[motif_name] = top_scores
    motif_start_dict[motif_name] = top_starts
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
        default=0.1)

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
    motif_score_dict= manager.dict() # {motif_name:[scores]}
    motif_start_dict = manager.dict() # {motif_name:[motif_start_pos]}
       
    for motif in all_motifs:
        pwm = motif[1]
        motif_name = motif[0]
        
        pool.apply_async(
        calculate_top_motif_matches_async,args=(sequence_array_list, 
                                                pwm, 
                                                motif_name, 
                                                motif_score_dict, 
                                                motif_start_dict
                                           )
                        )
    pool.close()
    pool.join()
    motif_score_dict = dict(motif_score_dict)
    motif_start_dict = dict(motif_start_dict)
    motif_score_frame = pd.DataFrame(motif_score_dict, 
                                     index = id_list)
    motif_start_frame = pd.DataFrame(motif_start_dict, 
                                     index = id_list)
    name_root = fasta_path.split('/')[-1].split('.')[0]
    score_tsv_path = output_path + '/' + name_root + '_motif_scores.tsv'
    start_tsv_path = output_path + '/' + name_root + '_motif_starts.tsv'
    motif_score_frame.to_csv(score_tsv_path, sep='\t')
    motif_start_frame.to_csv(start_tsv_path, sep='\t')
    
    end = time.time()
    print('total time', end-start)
