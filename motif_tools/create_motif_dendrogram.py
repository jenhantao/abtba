#!/usr/bin/env python
"""
"""

### imports ###
import argparse
import numpy as np
import os
import time
import scipy 
import matplotlib
from itertools import combinations
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
from scipy import ndimage
import pandas as pd
import Bio
from Bio import motifs
from motif_utilities import *
matplotlib.rcParams['savefig.dpi'] = 400

import warnings
warnings.filterwarnings("ignore")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Given a list of motifs and \
                                      an matrix of similarity scores, creates \
                                      a dendrogram of the motifs' )
    parser.add_argument("similarity_scores",
        help="path to a pickled numpy array containing motif similarity scores",
        type = str)
    parser.add_argument("output_path",
        help="directory where output file should be written",
        default="./", type=str)
    parser.add_argument("motif_files",
        help="list of motif files",
        type=str,
        nargs="+")
    parser.add_argument("-logos", action='store_true', 
        help="generate logos for dendrogram",
        default=False)

    # parse arguments
    args = parser.parse_args()

    similarity_scores_path= args.similarity_scores
    output_path = args.output_path
    motif_files = args.motif_files
    plot_logos = args.logos

    if not os.path.isdir(output_path):
        os.mkdir(output_path)

    score_frame = pd.read_csv(similarity_scores_path, sep='\t',index_col=0)
    correlations = score_frame.values
    motif_names = score_frame.columns.values

    # create directory to hold logos
    if plot_logos:
        if not os.path.isdir(output_path + '/logos'):
            os.mkdir(output_path + '/logos')
        # create logos
        for f in motif_files:
            logo_path = output_path + '/logos/' + '.'.join(f.split('/')[-1].split('.')[:-1])+'.png'
            with open(f) as mf:
                m = Bio.motifs.read(mf, 'jaspar')
            create_logo(m, logo_path, fmt='PNG')

    motif_path = output_path + '/logos/'
    sns.set_style('white')
    with sns.axes_style('white'):
        width = len(motif_names)/8
        height = width/2
        width = np.max([width, 4])
        height = np.max([height, 4])
        plt.figure(figsize=(width, height))
        correlations = np.clip(correlations, -1, 1)
        dissimilarity = 1-np.abs(correlations)
        coords = list(combinations(range(len(motif_names)),2))
        dissimilarity_as_pdist = [dissimilarity[x[0]][x[1]] for x in coords]

        Z=scipy.cluster.hierarchy.linkage(dissimilarity_as_pdist,
            method='centroid')
        
        if plot_logos:
            gs = matplotlib.gridspec.GridSpec(2, len(motif_names), wspace=0.0, hspace=0.0)
        else:
            gs = matplotlib.gridspec.GridSpec(1, len(motif_names), wspace=0.0, hspace=0.0)
        dendrogram_axis = plt.subplot(gs[0,:len(motif_names)])
        sns.despine()
        
        scipy.cluster.hierarchy.dendrogram(Z,
                                           color_threshold=0.2,
                                           ax=dendrogram_axis,
                                           labels=motif_names,)
        plt.axhline(0.2, linestyle='--', color='grey')
        plt.ylabel('Correlation Difference')
        
        sorted_motif_names = [x.get_text() for x in  dendrogram_axis.get_xticklabels()]
        if plot_logos:
            dendrogram_axis.set_xticklabels([])
        else:
            plt.xticks(rotation=90)

        if plot_logos:
            for i in range(len(motif_names)):
                current_axis = plt.subplot(gs[1, i])
                mn = sorted_motif_names[i]
                img = plt.imread(motif_path + '/' + mn + '.png')
                rotated_img = scipy.ndimage.rotate(img, 90)

                current_axis.imshow(rotated_img, origin = 'upper', extent=[0.0, 1.0, 0.0, 4.0])
                current_axis.set_xticks([])
                current_axis.set_yticks([])
                current_axis.set_xlabel(mn, rotation=90)
    plt.tight_layout()
    plt.savefig(output_path + '/motif_clustering_dendrogram.pdf', bbox_inches='tight')
