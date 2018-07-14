#!/usr/bin/env python
"""
Given standardized motif features, calculates the variance inflation factor
for each feature
"""

### imports ###
import argparse
import numpy as np
import os
import sys
import time
import pandas as pd
import sklearn
from sklearn import linear_model

### functions ###
def calculate_vif(features):
    '''
    calculates the VIF for each feature
    inputs: features, n X m (numSamples x numFeatures) vector of features
    output: VIFS, list of m VIFS
    '''
    vifs = []
    all_motifs = features.columns.values
    toolbar_width = 50

    # setup toolbar
    sys.stdout.write("[%s]" % (" " * toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line, after '['

    counter = 0
    num_motifs = len(all_motifs)
    for motif in all_motifs:
        current_motif_scores = features[[motif]]
        other_motif_scores = features[[x for x in all_motifs if not x == motif]]
        lr = sklearn.linear_model.LinearRegression(n_jobs=-1)
        lr.fit(other_motif_scores, current_motif_scores)
        
        # calculate the coefficient of determination
        coeff_det = lr.score(other_motif_scores, current_motif_scores)
        # calculate VIF
        if coeff_det == 1:
            vif = 10
        else:
            vif = 1/(1-coeff_det)
        vifs.append(vif)
        counter += 1
        sys.stdout.write('\r')
        sys.stdout.write("[%-50s] %d%%" % ('='*int(counter/num_motifs*toolbar_width), counter/num_motifs*100))
        sys.stdout.flush()
    

    sys.stdout.write("\n")
    to_return = pd.DataFrame({'Motif':all_motifs, 'VIF':vifs})
    return to_return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Given standardized motif features, \
        calculates the variance inflation factor')
    parser.add_argument("feature_path",
        help="path to a standardized feature set created by create_features.py",
        type = str)
    parser.add_argument("output_path",
        help="directory where output file should be written",
        default="./", type=str)

    # parse arguments
    args = parser.parse_args()

    feature_path = args.feature_path
    output_path = args.output_path

    # read in features
    feature_frame = pd.read_csv(feature_path, sep='\t', index_col=0)
    
    vifs = calculate_vif(feature_frame)

    vifs.to_csv(output_path, index = False, sep='\t')
