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
matplotlib.rcParams['savefig.dpi'] = 400



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

    # create directory to hold logos
    if not os.path.isdir(output_path + '/logos'):
        os.mkdir(output_path + '/logos')
    # create logos
    for f in motif_files:
        logo_path = output_path + '/logos/' + '.'.join(f.split('/')[-1].split('.')[:-1])
        if not os.path.isfile(logo_path + '.png'):
            os.system('motif2Logo.pl ' + f + ' -o ' + logo_path)

    correlation_data = np.load(similarity_scores_path)
    correlations = correlation_data['arr_0']
    motif_names = correlation_data['arr_1']

    motif_path = output_path + '/logos/'
    sns.set_style('white')
    with sns.axes_style('white'):
        width = len(motif_names)/8
        height = width/2
        plt.figure(figsize=(width, height))
        dissimilarity = 1-correlations
        coords = list(combinations(range(len(motif_names)),2))
        dissimilarity_as_pdist = [dissimilarity[x[0]][x[1]] for x in coords]

        Z=scipy.cluster.hierarchy.linkage(dissimilarity_as_pdist, 
                                         )
        gs = matplotlib.gridspec.GridSpec(2, len(motif_names), wspace=0.0, hspace=0.0)
        dendrogram_axis = plt.subplot(gs[0,:len(motif_names)])
        sns.despine()
        scipy.cluster.hierarchy.dendrogram(Z, 
                                           color_threshold=0.1,
                                           ax=dendrogram_axis,
                                           labels=motif_names)
        plt.axhline(0.1, linestyle='--', color='grey')
        plt.ylabel('Correlation Difference')
        
        sorted_motif_names = [x.get_text() for x in  dendrogram_axis.get_xticklabels()]
        for i in range(len(motif_names)):
            current_axis = plt.subplot(gs[1, i])
            mn = sorted_motif_names[i]
            img = plt.imread(motif_path + '/' + mn + '.png')
            rotated_img = scipy.ndimage.rotate(img, 90)
            current_axis.imshow(rotated_img)
            current_axis.set_xticks([])
            current_axis.set_yticks([])
            current_axis.axis('off')
    plt.savefig(output_path + '/motif_clustering_dendrogram.pdf', bbox_inches='tight')
