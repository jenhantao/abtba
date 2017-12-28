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
    
def train_classifier(features,
                     labels,
                     numIterations=5,
                     test_size=0.5):
    all_rocs = []
    all_precisions = []
    all_coefficients = []
    all_scores = []
    all_testLabels = []
    for i in range(numIterations):  
        # split data into training and test sets
        training_features, test_features, training_labels, test_labels = get_split(
            features, labels, test_size = test_size)

        #  Train affinity classifier
        classifier = sklearn.linear_model.LogisticRegression(penalty='l1', n_jobs=-1, tol=1e-8)
        classifier.fit(training_features, training_labels)

        # score predictions
        probas = classifier.predict_proba(test_features)
        current_roc = sklearn.metrics.roc_auc_score(test_labels, 
                                                              probas[:, 1], 
                                                              average = None)
        current_precision = sklearn.metrics.average_precision_score(test_labels,
                                                                probas[:, 1], 
                                                                average = None)

         # retrieve coefficients
        current_coefficients = classifier.coef_.flatten()
        
        all_rocs.append(current_roc)
        all_precisions.append(current_precision)
        all_coefficients.append(current_coefficients)
        all_scores.append(probas)
        all_testLabels.append(test_labels)
        
    # convert coefficients into data frame
    all_coefficients = pd.DataFrame(np.array(all_coefficients)).T
    all_coefficients.index = features.columns.values
    
    results = (all_rocs, 
               all_precisions, 
               all_coefficients,
               all_scores,
               all_testLabels)
    return results
    
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

def write_classifier_results(results, output_path):
    '''
    writes results of train_clasifier as tsv files
    '''
    all_rocs = results[0]
    all_precisions = results[1]
    all_coefficients = results[2]
    all_scores = results[3]
    all_test_labels = results[4]
    
    performance_frame = pd.DataFrame({'ROC':all_rocs, 
                                      'Precision Score': all_precisions})
    performance_frame.to_csv(output_path + '/performance.tsv', sep='\t', index=False)
    all_coefficients.to_csv(output_path + '/coefficients.tsv', sep='\t')

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

    # train classifier on each feature
    print('training individual classifiers')
    feature_performance_tuples = []
    for feature in feature_frame.columns.values:
        start = time.time()
        sub_features = feature_frame[[feature]]
        sub_results = train_classifier(sub_features,
                         labels,
                         numIterations=num_iterations,
                         test_size=test_fraction)
        sub_roc = np.mean(sub_results[0])
        sub_precision = np.mean(sub_results[1])
        feature_performance_tuples.append((feature, sub_roc, sub_precision))
        end = time.time()
        print('classifier using feature', feature, 'ROC:', sub_roc, 'Precision:', sub_precision, 'Training Time:', end-start)
    # use performance to sort individual features
    feature_performance_tuples.sort(key=lambda x:x[1], reverse=True)
    sorted_features = [x[0] for x in feature_performance_tuples]

    # write sorted performance of individual classifiers
    ind_performance_file = open(output_path + '/individual_performance.tsv', 'w')
    ind_performance_file.write('\t'.join(['Feature Name','ROC','Precision'])+'\n')
    for t in feature_performance_tuples:
        ind_performance_file.write('\t'.join([str(x) for x in t]) + '\n')
    ind_performance_file.close()
    
    print('Performing feed forward selection')
    # train classifier on increasingly large subsets of features
    numFeatures_performance_tuples = []
    for i in range(len(sorted_features)):
        num_features = i + 1
        sub_features = feature_frame[sorted_features[:num_features]]

        sub_results = train_classifier(sub_features,
                         labels,
                         numIterations=num_iterations,
                         test_size=test_fraction)

        sub_roc = np.mean(sub_results[0])
        sub_precision = np.mean(sub_results[1])
        numFeatures_performance_tuples.append((num_features, sub_roc, sub_precision))

    # write performance of classifiers
    ff_performance_file = open(output_path + '/feed_forward_performance.tsv', 'w')
    ff_performance_file.write('\t'.join(['Feature Name','ROC','Precision'])+'\n')
    for t in numFeatures_performance_tuples:
        ff_performance_file.write('\t'.join([str(x) for x in t]) + '\n')
    ff_performance_file.close()
