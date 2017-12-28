#!/usr/bin/env python
"""
Given standardized motif features and matching labels trains a classifier and 
returns performance metrics and model coefficients
"""

### imports ###
import argparse
import numpy as np
import os
import time
import pandas as pd
from sklearn import preprocessing
import sklearn
from sklearn import linear_model
from sklearn import cross_validation
import scipy

### functions ###
# split data into GC content matched training and test data
def get_split(features, labels, test_size):
    '''
    feature: 2D array (samples x features)
    labels: 1D boolean array (samples x)
    test_size: fraction of data to test on 
    '''
    
    ### match GC content of samples labelled True with those labelled False by thowing out False samples 
    # retrieve sequences using index of labels
    index_label_tuples = tuple(zip(labels.index.values, labels.values))
    
    true_ids = [x[0] for x in index_label_tuples if x[1]]
    
    false_ids = [x[0] for x in index_label_tuples if not x[1]]
       
    filtered_ids = true_ids + false_ids
    filtered_features = features[features.index.isin(filtered_ids)]
    filtered_labels = labels[labels.index.isin(filtered_ids)]

    if test_size <= 0.5: 
        training_indices, test_indices = next(iter(
                sklearn.cross_validation.StratifiedKFold(filtered_labels, int(1/test_size), shuffle=True)))
    else:   
        test_indices, training_indices = next( 
            iter(sklearn.cross_validation.StratifiedKFold(filtered_labels, int(1/(1-test_size)), shuffle=True)))
    training_ids = [filtered_ids[i] for i in training_indices]
    test_ids = [filtered_ids[i] for i in test_indices]
    
    training_features = filtered_features[filtered_features.index.isin(training_ids)]
    test_features = filtered_features[filtered_features.index.isin(test_ids)]
    training_labels = filtered_labels[filtered_labels.index.isin(training_ids)]
    test_labels = filtered_labels[filtered_labels.index.isin(test_ids)]
    
    return training_features, test_features, training_labels, test_labels

def calc_model_log_likelihood(probas, labels):
    log_likelihood = 0
    Y = labels.astype(float)
    for i in range(len(Y)):
        p_true = probas[i][1]
        p_false = probas[i][0]
        y = Y[i]
        prod = ((p_true)**y) * ((p_false) ** (1-y))
        log_prod = np.log(prod)
        log_likelihood += log_prod
    return log_likelihood

def calc_feature_pvals(features,
                       labels,
                       test_size=0.2,
                       num_iterations=5):
    pvals = []
    num_motifs = features.shape[1]
    # split data into training and test sets
    scaler = preprocessing.StandardScaler()

    # standardize features
    standardized_features = pd.DataFrame(scaler.fit_transform(features))
    standardized_features.columns = features.columns.values
    standardized_features.index = features.index.values

    for i in range(num_iterations):
        training_features, test_features, training_labels, test_labels = get_split(
            features, labels, test_size = test_size)
        
        # standardize training features
        standardized_training_features = pd.DataFrame(scaler.fit_transform(training_features))
        standardized_training_features.columns = training_features.columns.values
        standardized_training_features.index = training_features.index.values

        # standardize test features
        standardized_test_features = pd.DataFrame(scaler.fit_transform(test_features))
        standardized_test_features.columns = test_features.columns.values
        standardized_test_features.index = test_features.index.values
            
        #  Train affinity classifier
        classifier = sklearn.linear_model.LogisticRegression(penalty='l1', 
            solver='liblinear') 

        classifier.fit(standardized_training_features, training_labels)
        # score predictions

        probas = classifier.predict_proba(standardized_features) # [[p_false, p_true]...] 
        overall_log_likelihood = calc_model_log_likelihood(probas, labels)
        iter_pvals = []

        current_classifier = sklearn.linear_model.LogisticRegression(penalty='l1', 
            solver='liblinear')
        for motif in features.columns.values:
            print('testing', motif)
            current_features = standardized_features.drop(motif, axis=1, inplace=False)
            current_training_features = standardized_training_features.drop(motif, axis=1, inplace = False)
            current_classifier.fit(current_training_features, training_labels)
            current_probas = current_classifier.predict_proba(current_features)
            current_log_likelihood = calc_model_log_likelihood(current_probas, labels)

            stat = -2*(current_log_likelihood - overall_log_likelihood)
            p = scipy.stats.chisqprob(stat, df=1)

            iter_pvals.append(p)

        pvals.append(iter_pvals)
    return pvals
            

    
def read_labels(label_path):
    '''
    reads label files created by create_features.py and returns a pandas Series representation
    '''
    indices = []
    vals = []
    with open(label_path) as f:
        data = f.readlines()
    for line in data:
        tokens = line.strip().split()
        indices.append(tokens[0])
        if tokens[1] == '1':
            vals.append(True)
        else:
            vals.append(False)
    to_return = pd.Series(vals, index=indices)
    return to_return

def write_test_results(features, pvals, output_path):
    pval_dict = dict(zip(range(len(pvals)), pvals))
    pval_frame = pd.DataFrame(data = pval_dict, index = features.columns.values)
    pval_frame.to_csv(output_path + '/significance.tsv', sep = '\t')
    
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='given a fasta file and a list \
                                     and a list of motif files, calculates the \
                                     the best match to each motif for each \
                                     sequence' )
    parser.add_argument("feature_path",
        help="path to a standardized feature set created by create_features.py",
        type = str)
    parser.add_argument("label_path",
        help="path to a fasta_file containing negative sequences to score",
        type = str)
    parser.add_argument("output_path",
        help="directory where output file should be written",
        default="./", type=str)
    parser.add_argument("-num_iterations", 
        help="number of iterations to train classifier",
        type=int,
        default=5)
    parser.add_argument("-test_fraction", 
        help="fraction of data to use for testing classifier",
        type=float,
        default=0.5)

    # parse arguments
    args = parser.parse_args()

    feature_path = args.feature_path
    label_path = args.label_path
    output_path = args.output_path
    num_iterations = args.num_iterations
    test_fraction = args.test_fraction

    if not os.path.isdir(output_path):
        os.mkdir(output_path)

    # read in features
    print('reading features and labels')
    feature_frame = pd.read_csv(feature_path, sep='\t', index_col=0)

    # read in labels
    labels = read_labels(label_path)

    print('testing feature significance for', feature_path)
    pvals = calc_feature_pvals(feature_frame,
                            labels,
                            num_iterations=num_iterations,
                            test_size=0.2)
    print('writing results')
    write_test_results(feature_frame, pvals, output_path) 
