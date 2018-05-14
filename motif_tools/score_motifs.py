#!/usr/bin/env python 
# given paths to two motif files, aligns them using either needleman-wunsch 
# and/or smith-waterman. Scores are assigned in the score matrix using PCC

### imports ###
import sys
import numpy as np
import os
import matplotlib 
import pandas as pd
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import argparse
from motif_utilities import *
import multiprocessing

def score_motif_against_others(motifs, 
                               index, 
                               other_indices, 
                               result_dict
                               ):
    motif1 = motifs[index]
    cleaned_motif1 = (motif1[0], cleanMatrix(motif1[1]))

    for j in other_indices:
        motif2 = motifs[j]
        cleaned_motif2 = (motif2[0], cleanMatrix(motif2[1]))
        # calc scores for original orientation 
        alignment_fwd, alignScore_fwd = local_align_motifs(cleaned_motif1, cleaned_motif2)
        r_fwd = calcCorrelation(alignment_fwd[0], alignment_fwd[1])
        # calc scores for one motif reversed
        alignment_rev, alignScore_rev = local_align_motifs(cleaned_motif1, revCompMotif(cleaned_motif2))
        r_rev = calcCorrelation(alignment_rev[0], alignment_rev[1])
        # select largest score
        r = np.max([r_fwd, r_rev])
        result_dict[str(index) + '_' + str(j)] = r

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='calculates pairwise \
                                      scores between pairs of motifs')
    parser.add_argument('outputPath',
        help='path to output directory',
        type = str)
    parser.add_argument('motifFiles',
        help='space separated list of motif files',
        type=str,
        nargs='+')

    parser.add_argument('-num_procs', 
        help='number of cores to use',
        type=int,
        default=8)

    # parse arguments
    args = parser.parse_args()

    outputPath = args.outputPath
    motifFiles = args.motifFiles
    num_procs = args.num_procs

    # make output directory if it doesn't est
    if not os.path.isdir(outputPath):
        os.mkdir(outputPath)

    # read in motifs
    # find all motifs in input directory
    allMotifs = []
    motifNames = []
    motifIndexDict = {} # key: motif name, value: motif index
    print('Reading motif files...')
    counter = 0
    for mf in sorted(motifFiles):
        motif = readMotifFile(mf) # (name, PWM)
        allMotifs.append(motif)
        motif_name = motif[0]
        motifNames.append(motif_name)
        motifIndexDict[motif_name] = counter
        counter += 1


    # align motifs and calculate scores (pearson correlation)
    result_matrix = np.zeros((len(allMotifs), len(allMotifs)))
    print('Calculating alignments between motifs and scoring motifs')
    pool = multiprocessing.Pool(processes=num_procs)
    manager = multiprocessing.Manager()
    result_dict = manager.dict()
    for i in range(len(allMotifs) - 1 ):
        motif1 = allMotifs[i]
        other_indices = list(range(i+1, len(allMotifs)))
        pool.apply_async(score_motif_against_others,
                         args = (allMotifs, i, other_indices, result_dict)
                        )
    pool.close()
    pool.join()
    
    # fill in matrix with result dict
    for key in result_dict.keys():
        tokens = key.split('_')
        i = int(tokens[0])
        j = int(tokens[1])

        result_matrix[i][j] = result_dict[key]
        result_matrix[j][i] = result_dict[key]
 
    # fill in diagonal scores
    for i in range(len(allMotifs)):
        result_matrix[i][i] = 1.0

    result_matrix = np.array(result_matrix)

    pearson_array = []
    print('Creating visualizations...')
    # plot the score distribution
    for i in range(result_matrix.shape[0] - 1):
        for j in range(i+1, result_matrix.shape[0]):
            pearson_array.append(result_matrix[i][j])
    plt.hist(pearson_array, bins = 20)
    plt.ylabel('Frequency (KDE)')
    plt.xlabel('PCC')
    plt.title('PCC Distribution')
    plt.savefig(outputPath + '/correlation_distribution.pdf')
    plt.close()
            
    # save matrix as spreadsheet
    frame = pd.DataFrame(data = result_matrix,
        index=motifNames,
        columns=motifNames)
    frame.to_csv(outputPath+'/correlation.tsv', sep='\t')

