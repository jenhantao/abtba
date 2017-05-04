#!/usr/bin/env python 
# given paths to two motif files, aligns them using either needleman-wunsch and/or smith-waterman. Scores are assigned in the score matrix using a linear combination of the nucleotide frequenceies

### imports ###
import sys
import numpy as np
import os
import matplotlib 
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import argparse
from motif_utilities import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='calculates pairwise \
                                      scores between pairs of motifs')
    parser.add_argument("outputPath",
        help="path to output directory",
        type = str)
    parser.add_argument("motifFiles",
        help="space separated list of motif files",
        type=str,
        nargs="+")

    # parse arguments
    args = parser.parse_args()

    outputPath = args.outputPath
    motifFiles = args.motifFiles

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

    # align motifs and calculate scores (euclidian, r-value)
    alignScores = np.zeros((len(allMotifs), len(allMotifs)))
    distances = np.zeros((len(allMotifs), len(allMotifs)))
    rs = np.zeros((len(allMotifs), len(allMotifs)))
    print('Calculating alignments between motifs and scoring motifs')
    for i in range(len(allMotifs) - 1 ):
        motif1 = allMotifs[i]
        print('scoring against ' +  str(i+1) + '/' +str(len(allMotifs)))
        for j in range(i + 1, len(allMotifs)):
            # calc scores for the forward direction
            motif2 = allMotifs[j]
            alignment_fwd, alignScore_fwd = local_align_motifs(motif1, motif2)
            dist_fwd = calcDistance(alignment_fwd[0], alignment_fwd[1])
            r_fwd = calcCorrelation(alignment_fwd[0], alignment_fwd[1])
            # calc scores for one motif reversed
            alignment_rev, alignScore_rev = local_align_motifs(motif1, revCompMotif(motif2))
            dist_rev = calcDistance(alignment_rev[0], alignment_rev[1])
            r_rev = calcCorrelation(alignment_rev[0], alignment_rev[1])

            dist = np.max([dist_fwd, dist_rev])
            r = np.max([r_fwd, r_rev])
            alignScore = np.max([alignScore_fwd, alignScore_rev])

            # add better score to matrix
            alignScores[i][j] = alignScore
            alignScores[j][i] = alignScore
            distances[i][j] = dist
            distances[j][i] = dist
            rs[i][j] = r
            rs[j][i] = r
    
    # fill in diagonal scores
    for i in range(len(allMotifs)):
        rs[i][i] = 1.0
        distances[i][i] = 1.0

    rs = np.array(rs)
    distances = np.array(distances)
    alignScores = np.array(alignScores)

    rsArray = []
    distancesArray = []
    alignScoresArray = []
    print('Creating visualizations...')
    # plot the score distribution
    for i in range(rs.shape[0] - 1):
        for j in range(i+1, rs.shape[0]):
            rsArray.append(rs[i][j])
            distancesArray.append(distances[i][j])
            alignScoresArray.append(alignScores[i][j])
    plt.hist(rsArray, bins = 20)
    plt.title("Rs Value Distribution")
    plt.savefig(outputPath + "/correlation_distribution.png")
    plt.close()
            
    plt.hist(alignScoresArray, bins = 20)
    plt.title("AlignScores Value Distribution")
    plt.savefig(outputPath + "/alignScores_distribution.png")
    plt.close()

    plt.hist(distancesArray, bins = 20)
    plt.title("Distances Value Distribution")
    plt.savefig(outputPath + "/distances_distribution.png")
    plt.close()

    # plot estimates for the number of motifs that would be merged according to threshold
    numMerged = []
    thresholds = []
    for t in np.arange(0,1,0.01):
        merged = set()
        thresholds.append(t)
        for i in range(rs.shape[0] - 1):
            for j in range(i+1, rs.shape[0]):
                if rs[i][j] > t:
                    merged.add(i)
                    merged.add(j)
        numMerged.append(len(merged))

    plt.plot(thresholds, numMerged)
    plt.xlim([min(thresholds),1.0])
    plt.xlabel("correlation threshold")
    plt.ylabel("number of motifs merged")
    plt.savefig(outputPath + "/correlationThreshold.png")
    plt.close()

    std = np.std(rsArray)
    mean = np.mean(rsArray)
    significance = np.arange(0.95,1,0.0005)
    numMerged = []
    thresholds = []

    for sig in significance:
        merged = set()
        zThreshold = st.norm.ppf(sig)
        threshold = min([1.0, zThreshold * std + mean]) # pearson can't exceed 1.0
        thresholds.append(threshold)
        for i in range(rs.shape[0] - 1):
            for j in range(i+1, rs.shape[0]):
                if rs[i][j] > threshold:
                    merged.add(i)
                    merged.add(j)
        numMerged.append(len(merged))


    plt.plot(significance, numMerged)
    plt.xlabel("significance threshold")
    plt.ylabel("number of motifs merged")
    plt.savefig(outputPath + "/significanceThreshold.png")
    plt.close()

    # save matrix
    print('Serializing scores... \ncorrelation should be used for clustering.')
    # pearson correlation for each pair of motifs
    np.savez_compressed(outputPath+"/correlation", rs, motifNames)
    # euclidian distance between each pair of motifs
